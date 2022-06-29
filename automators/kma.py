import os
import glob
import click
import pickle
import shutil
from biotools import mash
from amrsummary import before_send
import sentry_sdk
from automator_settings import SENTRY_DSN
from nastools.nastools import retrieve_nas_files
from externalretrieve import upload_to_ftp
import fileinput
from pathlib import Path
import csv



@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def kma_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Current list of analysis types that KMA can use
    analyses = [
        'custom', 'amr', 'biocide', 'metal'
    ]
    # Current list of sequence types that KMA can analyse
    seqtypes = [
        'fasta', 'fastq'
    ]

    # Variable to hold supplied arguments
    argument_dict = {
        'analysis': str(),
        'seqtype': 'fasta',
        'nanopore': False,
    }

    # Set the database path for the analyses
    dbpath = '/mnt/nas2/databases/kma_v_1.4.2_db/'
    database_path = {
        'custom': os.path.join(work_dir, 'targets'),
        'amr': os.path.join(dbpath, 'NCBI-AMR-v3.10'),
        'biocide': os.path.join(dbpath, 'NCBI-biocide-v3.10'),
        'metal': os.path.join(dbpath, 'NCBI-metal-v3.10'),
    }
    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'ANALYSIS' in item:
                argument_dict['analysis'] = item.split('=')[1].lower()
                continue
            if 'SEQTYPE' in item:
                argument_dict['seqtype'] = item.split('=')[1].lower()
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        #store the custom target file and index it
        if argument_dict['analysis'] == 'custom':
            # Set and create the directory to store the custom targets
            target_dir = os.path.join(work_dir, 'targets')
            try:
                os.mkdir(target_dir)
            except FileExistsError:
                pass
            # Download the attached FASTA file.
            # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
            # out a better way to do it.
            attachment = redmine_instance.issue.get(issue.id, include='attachments')
            attachment_id = 0
            for item in attachment.attachments:
                attachment_id = item.id
            # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
            if attachment_id != 0:
                attachment = redmine_instance.attachment.get(attachment_id)
                attachment.download(savepath=target_dir, filename='targets.fasta')
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Analysis type custom requires an attached FASTA file of '
                                                    'targets. The automator could not find any attached files. '
                                                    'Please create a new issue with the FASTA file attached and try '
                                                    'again.',
                                              status_id=4)

        #TODO: if they upload a custom database, index it using kma_index
  

        # Ensure that the analysis type is provided
        if not argument_dict['analysis']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not identify an analysis type. '
                                                'Please ensure that the first line of the issue contains one'
                                                ' of the following keywords: {ats}'.format(ats=', '.join(analyses)),
                                          status_id=4)
            return
        elif argument_dict['analysis'] not in analyses:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied analysis type {at} currently not in the supported '
                                                'list of analyses: {ats}'.format(at=argument_dict['analysis'],
                                                                                 ats=', '.join(analyses)),
                                          status_id=4)
            return

        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No SEQIDs provided!',
                                          status_id=4)
            return

        #create a folder to hold the sequences
        seq_dir = os.path.join(work_dir, 'sequences')
        os.mkdir(seq_dir)

        #pull assemblies from nas
        if argument_dict['seqtype'] == 'fasta':
        # Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
            retrieve_nas_files(seqids=seqids,
                               outdir=seq_dir,
                               filetype='fasta',
                               copyflag=False)
            missing_fastas = verify_fasta_files_present(seqids, seq_dir)
            # Update the Redmine issue if one or more of the requested SEQIDs could not be located
            if missing_fastas:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastas))
        #    return
        
        #TODO: add function for fastq files
        #pull raw reads from nas
        #if argument_dict['seqtype'] == 'fastq':
        # Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
        #    retrieve_nas_files(seqids=seqids,
        #                       outdir=seq_dir,
        #                       filetype='fastq',
        #                       copyflag=False)
        #    missing_fastas = verify_fasta_files_present(seqids, seq_dir)
            # Update the Redmine issue if one or more of the requested SEQIDs could not be located
        #    if missing_fastas:
        #        redmine_instance.issue.update(resource_id=issue.id,
        #                                      notes='WARNING: Could not find the following requested SEQIDs on'
        #                                            ' the OLC NAS: {}'.format(missing_fastas))
        #    return


        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/card-rgi'
        kma_py = '/mnt/nas2/virtual_environments/card-rgi/bin/kma'

        # Run kma with the necessary arguments
        #assembly files
        for assembly in glob.glob(os.path.join(seq_dir, '*.fasta')):
            seqid = os.path.split(assembly)[1].split('.')[0]
            #prepare command for kma
            kma_cmd = 'kma -i {assembly} -o {seqid} -t_db {dbpath} -t 7 -nf'\
                .format(assembly=assembly,
                        seqid=seqid,
                        dbpath=database_path[argument_dict['analysis']])

            # Append the align and/or the unique flags are required
            kma_cmd += ' -bcNano' if argument_dict['nanopore'] else ''

            # Create another shell script to execute within the PlasmidExtractor conda environment
            template = "#!/bin/bash\n{} && cd {} && {}".format(activate, seq_dir, kma_cmd)
            kma_script = os.path.join(work_dir, 'run_kma.sh')
            with open(kma_script, 'w+') as file:
                file.write(template)
            # Modify the permissions of the script to allow it to be run on the node
            make_executable(kma_script)
            # Run shell script
            os.system(kma_script)

        #add the filename (which is the seqid) to the first column in all of the res files, then concatenate into a single output csv file
        outputfile = os.path.join(seq_dir, 'kma_output.csv')
        with open(outputfile, 'w', newline='') as file_output:
            csv_output = csv.writer(file_output)
            #for fname in glob.glob(os.path.basename(seq_dir, '*.res')):
            for fname in glob.glob(os.path.join(seq_dir, '*.res')):
                fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
                seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
                with open(fname, newline='') as f_input:
                    csv_input = csv.reader(f_input)

                    #Header processing
                    header = csv_input.__next__()
                    header.insert(0,"SeqID")

                    csv_output.writerow(header)
                    for row in csv_input:
                        print(row)
                        row.insert(0,seqname) #this adds the seqid to the file before concatenating
                        #row.insert(0,fname)
                        csv_output.writerow(row)

        #now create an output folder so we can delete all of the extra files we didn't need
        out_dir = os.path.join(work_dir, 'output')
        os.mkdir(out_dir)

        #now copy the output csv file to the new output directory... leaving as a directory so we can zip it if we add other functions and outputs later
        shutil.copyfile(os.path.join(seq_dir, 'kma_output.csv'), os.path.join(out_dir, 'kma_output.csv'))

           
        # Zip output
        output_filename = 'kma_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=out_dir,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'
        
        upload_successful = upload_to_ftp(local_file=zip_filepath)
        # Prepare upload
        if upload_successful:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='KMA analysis complete!\n\n'
                                                'Results are available at the following FTP address:\n'
                                                'ftp://ftp.agr.gc.ca/outgoing/cfia-ac/{}'
                                          .format(os.path.split(zip_filepath)[1]))
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity '
                                                'issues. '
                                                'Please try again later.')


        # Remove all the folders
        shutil.rmtree(seq_dir)
        # Remove the zip file
        os.remove(zip_filepath)
        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes='{at} analysis with KMA complete!'
                                      .format(at=argument_dict['analysis'].lower()))
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Please contact a bioinformatician '
                                            'to investigate.')




def verify_fasta_files_present(seqid_list, fasta_dir):
    """
    Makes sure that FASTQ files specified in seqid_list have been successfully copied/linked to directory specified
    by fastq_dir.
    :param seqid_list: List with SEQIDs.
    :param fasta_dir: Directory to which FASTA files should have been linked
    :return: List of SEQIDs that did not have files associated with them.
    """
    missing_fastas = list()
    for seqid in seqid_list:
        # Check forward.
        if len(glob.glob(os.path.join(fasta_dir, '{seqid}*fasta*'.format(seqid=seqid)))) == 0:
            missing_fastas.append(seqid)
    return missing_fastas


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
    kma_redmine()
