from accessoryFunctions.accessoryFunctions import make_path
from nastools.nastools import retrieve_nas_files
from biotools import mash
import pickle
import shutil
import click
import glob
import os

import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def staramr_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))
        # Parse description to get list of SeqIDs
        seqids = list()
        for i in range(0, len(description)):
            item = description[i]
            item = item.upper()
            # Minimal check to make sure IDs provided somewhat resemble a valid sample ID
            if item.isalpha():
                pass
            else:
                seqids.append(item)

        # Run Mash
        with open(os.path.join(work_dir, 'seqid.txt'), 'w') as f:
            for seqid in seqids:
                f.write(seqid + '\n')
        # Drop FASTA files into workdir
        retrieve_nas_files(seqids=seqids,
                           outdir=work_dir,
                           filetype='fasta',
                           copyflag=False)
        # Create output directory
        output_dir = os.path.join(work_dir, 'output')
        make_path(output_dir)
        # Get all of the FASTA files
        fasta_list = sorted(glob.glob(os.path.join(work_dir, '*.fasta')))
        # Set the folder to store all the PointFinder outputs
        staramr_output_dir = os.path.join(work_dir, 'staramr_outputs')
        # Initialise a dictionaries to store the mash-calculated, and pointfinder-formatted genus outputs for each strain
        genus_dict = dict()
        organism_dict = dict()
        # Create lists to store missing and unprocessed seqids
        unprocessed_seqs = list()
        missing_seqs = list()
        mash_fails = list()
        # Dictionary to convert the mash-calculated genus to the pointfinder format
        pointfinder_org_dict = {'Campylobacter': 'campylobacter',
                                'Escherichia': 'e.coli',
                                'Shigella': 'e.coli',
                                '‎Mycobacterium': 'tuberculosis',
                                'Neisseria': 'gonorrhoeae',
                                'Salmonella': 'salmonella'}
        # Reverse look-up dictionary
        rev_org_dict = {'campylobacter': 'Campylobacter',
                        'e.coli': 'Escherichia',
                        'tuberculosis': 'Mycobacterium',
                        'gonorrhoeae': 'Neisseria',
                        'salmonella': 'Salmonella'}

        # Run mash screen on each of the assemblies
        for item in fasta_list:
            seqid = os.path.splitext(os.path.basename(item))[0]
            screen_file = os.path.join(output_dir, '{seqid}_screen.tab'.format(seqid=seqid))
            mash.screen('/mnt/nas2/databases/confindr/databases/refseq.msh',
                        item,
                        threads=8,
                        w='',
                        i='0.95',
                        output_file=screen_file,
                        returncmd=True)
            screen_output = mash.read_mash_screen(screen_file)
            # Determine the genus from the screen output file
            for screen in screen_output:
                # Extract the genus from the mash results
                mash_organism = screen.query_id.split('/')[-3]
                # Use the organism as a key in the pointfinder database name conversion dictionary
                try:
                    mash_genus = pointfinder_org_dict[mash_organism]
                except KeyError:
                    mash_genus = 'NA'
                # Populate the dictionaries with the seqid, and the calculated genus/pointfinder name
                genus_dict[seqid] = mash_genus
                organism_dict[seqid] = mash_organism
        # Delete all of the FASTA files
        for fasta in fasta_list:
            os.remove(fasta)
        # # Delete the output folder
        # shutil.rmtree(output_dir)

        # Pointfinder
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/staramr'
        staramr_py = '/mnt/nas2/virtual_environments/staramr/bin/staramr'
        # List of organisms in the pointfinder database
        staramr_list = ['campylobacter', 'salmonella']
        os.makedirs(stararm_output_dir, exist_ok=True)

        genus_seqid_dict = dict()
        for seqid in sorted(seqids):
            try:
                seqid_genus = genus_dict[seqid]
                if seqid_genus not in genus_seqid_dict:
                    genus_seqid_dict[seqid_genus] = [seqid]
                else:
                    genus_seqid_dict[seqid_genus].append(seqid)
            except KeyError:  # Mash sometimes doesn't find a genus!
                mash_fails.append(seqid)

        for genus in genus_seqid_dict:
            if genus in staramr_list:
                assembly_folder = os.path.join(work_dir, genus)
                make_path(assembly_folder)
                retrieve_nas_files(seqids=genus_seqid_dict[genus],
                                   outdir=assembly_folder,
                                   filetype='fasta',
                                   copyflag=False)
                fastas = sorted(glob.glob(os.path.join(assembly_folder, '*.fasta')))
                outdir = os.path.join(staramr_output_dir, genus)
                cmd = '{py} search --pointfinder-organism {orgn} -o {output} ' \
                    .format(py=staramr_py,
                            orgn=genus,
                            output=outdir,
                            )
                for fasta in fastas:
                    cmd += fasta + ' '
                # Create another shell script to execute within the PlasmidExtractor conda environment
                template = "#!/bin/bash\n{} && {}".format(activate, cmd)
                pointfinder_script = os.path.join(work_dir, 'run_staramr.sh')
                with open(pointfinder_script, 'w+') as f:
                    f.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(pointfinder_script)
                # Run shell script
                os.system(pointfinder_script)
            else:
                for seqid in genus_seqid_dict[genus]:
                    unprocessed_seqs.append(seqid)

        # Zip output
        output_filename = 'staramr_output'
        zip_filepath = zip_folder(results_path=staramr_output_dir,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'
        # Prepare upload
        output_list = [
            {
                'filename': os.path.basename(zip_filepath),
                'path': zip_filepath
            }
        ]

        # Create a note to add to the updated Redmine issue
        notes = 'StarAMR process complete!'
        # If there are missing, or unprocessed sequences, add details to the note
        if unprocessed_seqs:
            seq_list = list()
            for sequence in unprocessed_seqs:
                seq_list.append('{seqid} ({organism})'.format(seqid=sequence,
                                                              organism=organism_dict[sequence]))
            if len(unprocessed_seqs) > 1:
                notes += '\n The following sequences were not processed, as they were determined to be genera not ' \
                         'present in the StarAMR database: {seqs}'.format(seqs=', '.join(seq_list))
            else:
                notes += '\n The following sequence was not processed, as it was determined to be a genus not ' \
                         'present in the StarAMR database: {seqs}'.format(seqs=', '.join(seq_list))
        if missing_seqs:
            if len(missing_seqs) > 1:
                notes += '\n The following sequences were not processed, as they could not be located in the strain ' \
                         'database: {seqs}'.format(seqs=', '.join(missing_seqs))
            else:
                notes += '\n The following sequence was not processed, as it could not be located in the strain database:' \
                         ' {seqs}'.format(seqs=', '.join(missing_seqs))
        if mash_fails:
            if len(mash_fails) > 1:
                notes += '\n The following sequences could not be processed by MASH screen: {seqs}'\
                    .format(seqs=', '.join(mash_fails))
            else:
                notes += '\n The following sequence could not be processed by MASH screen: {seqs}'\
                    .format(seqs=', '.join(mash_fails))
        # Create a list of all the folders - will be used to clean up the working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)
        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes=notes)
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')

def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)


def zip_folder(results_path, output_dir, output_filename):
    """
    Compress a folder
    :param results_path: The path of the folder to be compressed
    :param output_dir: The output directory
    :param output_filename: The output file name
    :return:
    """
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    staramr_redmine()
