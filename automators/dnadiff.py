import os
import glob
import click
import pickle
import shutil
from biotools import mash
from amrsummary import before_send
import sentry_sdk
from automator_settings import FTP_FOLDER, SENTRY_DSN
from nastools.nastools import retrieve_nas_files
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
import fileinput
from pathlib import Path
import csv
import pandas




@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def dnadiff_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    # Current list of sequence types that KMA can analyse
    seqtypes = [
        'fasta', 'fastq', 'minionfastq'
    ]

    # Variable to hold supplied arguments
    argument_dict = {
        'reference': str(),
        'query': str(),
    }


    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.rstrip()
            if 'reference' in item:
                argument_dict['reference'] = item.split('=')[1].upper()
                seqids.append(argument_dict['reference'])
                continue
            if 'query' in item:
                argument_dict['query'] = item.split('=')[1].upper()
                seqids.append(argument_dict['query'])
                continue
            # Otherwise the item should be a SEQID
            #seqids.append(item)

        # Ensure that a reference sequence is provided
        if not argument_dict['reference']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No reference provided. ',
                                          status_id=4)
            return

        # Ensure that a query sequence is provided
        if not argument_dict['query']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No query provided',
                                          status_id=4)
            return


        #create a folder to hold the sequences
        seq_dir = os.path.join(work_dir, 'output')
        os.makedirs(seq_dir, exist_ok=True)

        #pull assemblies from nas
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

        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dRep'

        # Run dnadiff with the necessary arguments
        dnadiff_cmd = 'dnadiff -p {rid}_dnadiff {reference}.fasta {query}.fasta'.format(reference=argument_dict['reference'],
                                                                                        rid=issue.id,
                                                                                        query=argument_dict['query'])
        #need to use sudo?
        #pcmd = os.sytem('echo %s|sudo -S %s' % (CLUSTER_PASSWORD, dnadiff_cmd))

        # Create another shell script to execute within the KMA conda environment
        template = "#!/bin/bash\n{} && cd {} && {}".format(activate, seq_dir, dnadiff_cmd)
        dnadiff_script = os.path.join(work_dir, 'run_dnadiff.sh')
        with open(dnadiff_script, 'w+') as file:
            file.write(template)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(dnadiff_script)
        # Run shell script
        os.system(dnadiff_script)

        # Remove the fasta files from the seq_dir, which contains our outputs
        fastafiles = os.listdir(seq_dir)
        for file in fastafiles:
            if file.endswith(".fasta"):
                os.remove(file)
          
        # Zip output
        dnadiffout_filename = 'dnadiff_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=seq_dir,
                                  output_dir=work_dir,
                                  output_filename=dnadiffout_filename)
        zip_filepath += '.zip'
        
        # Upload the zip file to Dropbox
        download_link = upload_to_dropbox(
            access_token=DROPBOX_ACCESS_TOKEN,
            refresh_token=DROPBOX_REFRESH_TOKEN,
            app_key=DROPBOX_APP_KEY,
            app_secret=DROPBOX_APP_SECRET,
            local_file_path=zip_filepath
        )

        if download_link:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='dnadiff analysis complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=download_link)
            )
        else:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='Upload of results was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Please contact a bioinformatician '
                                            'to investigate: {e}.'.format(e=e))




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

def check_fastqs_present(fastq_list, fastq_dir):
    missing_fastqs = list()
    for seqid in fastq_list:
        if len(glob.glob(os.path.join(fastq_dir, seqid + '*.fastq.gz'))) < 2:
            # JAS adding an extra if statement here to allow for Nanopore (SE) reads
            if len(glob.glob(os.path.join(fastq_dir, seqid + ".fastq.gz"))) == 0:
                missing_fastqs.append(seqid)
    return missing_fastqs

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
    dnadiff_redmine()
