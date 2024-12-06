import os
from glob import glob
import click
import pickle
import shutil
from amrsummary import before_send
import sentry_sdk
from automator_settings import FTP_FOLDER, SENTRY_DSN
from nastools.nastools import retrieve_nas_files
import fileinput
from pathlib import Path
import csv
from upload_to_dropbox import upload_to_dropbox
from tokens import DROPBOX_ACCESS_TOKEN, DROPBOX_REFRESH_TOKEN, DROPBOX_APP_KEY, DROPBOX_APP_SECRET



@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def fastqc_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = []
        for item in description:
            item = item.upper().rstrip()
            seqids.append(item)

        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='WARNING: No SEQIDs provided!',
                status_id=4
            )
            return

        # make folder in the redmine biorequest working directory
        output = os.path.join(work_dir, 'fastqc')
        os.makedirs(output)
        raw_dir = os.path.join(work_dir, 'seqs')
        os.makedirs(raw_dir)

        #retrieve the fastq files
        retrieve_nas_files(
            seqids=seqids,
            outdir=raw_dir,
            filetype='fastq',
            copyflag=False
        )
        missing_fastqs = check_fastqs_present(seqids, raw_dir)
        if len(missing_fastqs) > 0:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='WARNING: Could not find the following requested FASTQ '
                'SEQIDs on the OLC NAS: {}'.format(missing_fastqs)
            )

        #now run fastqc on all of the files
        activate = 'source /home/ubuntu/miniconda3/bin/activate ' \
            '/mnt/nas2/virtual_environments/fastqc'
        fastqc_cmd = 'fastqc -o {out} {raw}/*.fastq.gz'.format(
            out=output,
            raw=raw_dir
        )

        #create another shell script to execute within the fastqc environment
        template = "#!/bin/bash\n{} && {}".format(activate, fastqc_cmd)
        fastqc_script = os.path.join(work_dir, 'run_fastqc.sh')

        with open(fastqc_script, 'w+') as file:
            file.write(template)
        
        # Modify the permissions of the script to allow it to be run on
        # the node
        make_executable(fastqc_script)
        # Run shell script on cluster
        os.system(fastqc_script)

        #now run multiqc on the fastqc outputs
        multiqc_cmd = 'multiqc {out} --title {req} --outdir {out}'.format(
            req=issue.id,
            out=output
        )
        #create another shell script to execute within the fastqc environment
        template2 = "#!/bin/bash\n{} && {}".format(activate, multiqc_cmd)
        multiqc_script = os.path.join(work_dir, 'run_multiqc.sh')

        with open(multiqc_script, 'w+') as file:
            file.write(template2)
        
        # Modify the permissions of the script to allow it to be run on
        # the node
        make_executable(multiqc_script)
        # Run shell script on cluster
        os.system(multiqc_script)

        # Now make a zip folder that we'll upload to Dropbox and to Redmine
        # (Dropbox as backup in case file is too big)
        output_filename = 'fastqc_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=output,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'

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
                status_id=2,
                notes='Fastqc & Multiqc process complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=download_link)
            )
        else:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=2,
                notes='Upload of result files was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )

        #upload some files to the redmine request
        output_list = []

        #multiqc report
        output_dict = {}
        report_name = '{}_multiqc_report.html'.format(issue.id)
        output_dict['path'] = os.path.join(work_dir, 'fastqc', report_name)
        output_dict['filename'] = report_name
        output_list.append(output_dict)

        #the zipped fastqc report folder
        output_dict = {}
        output_dict['path'] = zip_filepath
        output_dict['filename'] = '{}.zip'.format(output_filename)
        output_list.append(output_dict)

        #upload the files to redmine
        redmine_instance.issue.update(
            resource_id=issue.id,
            uploads=output_list,
            status_id=4,notes='Fastqc and Multiqc analysis complete!'
        )

        # And finally, do some file cleanup.
        shutil.rmtree(output)
        shutil.rmtree(raw_dir)
        os.remove(zip_filepath)

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes='Something went wrong! We log this automatically and will '
            'look into the problem and get back to you with a fix soon.'
        )



def verify_fasta_files_present(seqid_list, fasta_dir):
    """
    Makes sure that FASTQ files specified in seqid_list have been successfully
    copied/linked to directory specified by fastq_dir.
    :param seqid_list: List with SEQIDs.
    :param fasta_dir: Directory to which FASTA files should have been linked
    :return: List of SEQIDs that did not have files associated with them.
    """
    missing_fastas = []
    for seqid in seqid_list:
        # Check forward.
        forward = glob(
            os.path.join(
                fasta_dir,
                '{seqid}*fasta*'.format(seqid=seqid)
            )
        )
        if len(forward) == 0:
            missing_fastas.append(seqid)
    return missing_fastas

def check_fastqs_present(fastq_list, fastq_dir):
    missing_fastqs = []
    for seqid in fastq_list:
        number_fastqs = len(
            glob(
                os.path.join(
                    fastq_dir, seqid + '*.fastq.gz')
            )
        )
        if number_fastqs < 2:
            # JAS adding an extra if statement here to allow for Nanopore
            # (SE) reads
            if len(glob(os.path.join(fastq_dir, seqid + ".fastq.gz"))) == 0:
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
    fastqc_redmine()
