import os
import glob
import click
import pickle
import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN
from biotools import mash
from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def closerelatives_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        # First line of description should be number of close relatives desired.
        try:
            num_close_relatives = int(description[0])
        except ValueError:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Error! The first line of the description must be the number'
                                                ' of strains you want to find. The first line of your '
                                                'description was: {}'.format(description[0]),
                                          status_id=4)
            return

        # Second line of description should be the SEQID of what you want to find a close reference for.
        seqid = description[1]

        # Try to extract FASTA files for the specified SEQID.
        retrieve_nas_files(seqids=[seqid],
                           outdir=os.path.join(work_dir, 'fasta'),
                           filetype='fasta',
                           copyflag=False)
        if len(glob.glob(os.path.join(work_dir, 'fasta', '*.fasta'))) != 1:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Error! Could not find FASTA file for the specified SEQID. The SEQID'
                                                ' that you specified was: {}'.format(seqid),
                                          status_id=4)
            return

        # Run mash dist with the FASTQ file specified against the sketch of all our stuff.
        query_fasta = glob.glob(os.path.join(work_dir, 'fasta', '*.fasta'))[0]
        #mash.dist(query_fasta, '/mnt/nas2/databases/mash_all_updated2023_database/part_10_886.msh', ## Mathu changed the sketch files here with the updated one#####
        #mash.dist(query_fasta, '/mnt/nas2/redmine/bio_requests/14674/all_sequences.msh', #changed on 2024-08-20 AC
        mash.dist(query_fasta, '/mnt/nas2/databases/mash_nas_assemblies/olc_assemblies_mashsketch.msh',
                  threads=8, output_file=os.path.join(work_dir, 'distances.tab'))
        mash_results = mash.read_mash_output(os.path.join(work_dir, 'distances.tab'))
        result_dict = dict()
        # Put all the results into a dictionary, where the key is the sequence file and the value is mash distance
        # between query fastq and reference fastq.
        for item in mash_results:
            seq_name = os.path.split(item.query)[-1].split('_')[0]
            result_dict[seq_name] = item.distance

        # Sort the results, store the sorted dictionary keys in a list.
        sorted_distance_results = sorted(result_dict, key=result_dict.get)

        # Prepare a string that lists the top hit SEQIDs to be posted to redmine.
        upload_string = ''
        for i in range(num_close_relatives):
            upload_string = upload_string + sorted_distance_results[i].replace('.fasta', '') + ' (' + str(result_dict[sorted_distance_results[i]]) + ')\n'
            upload_string = upload_string + sorted_distance_results[i].replace('.fna.gz', '') + ' (' + str(result_dict[sorted_distance_results[i]]) + ')\n' #added by Ashley

        # Also make a CSV file of all results, in case someone wants to take a closer look.
        with open(os.path.join(work_dir, 'close_relatives_results.csv'), 'w') as f:
            f.write('Strain,MashDistance\n')
            for seq in sorted_distance_results:
                f.write('{},{}\n'.format(seq.replace('.fasta', ''), result_dict[seq]))
                f.write('{},{}\n'.format(seq.replace('.fna.gz', ''), result_dict[seq]))

        output_list = [
            {
                'path': os.path.join(work_dir, 'close_relatives_results.csv'),
                'filename': 'close_relatives_results.csv'
            }
        ]
        # Post the list of closely related SEQIDs to redmine, as well as the CSV result file.
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Process complete! Here is the list of the {num_relatives} closest strains '
                                            'to {query_strain} (mash distance between query and result in brackets):'
                                            '\n{upload_string}'.format(num_relatives=str(num_close_relatives),
                                                                       query_strain=seqid,
                                                                       upload_string=upload_string),
                                      status_id=4,
                                      uploads=output_list)

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')

if __name__ == '__main__':
    closerelatives_redmine()
