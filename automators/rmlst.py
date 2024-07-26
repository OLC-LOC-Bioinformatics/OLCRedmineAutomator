"""
Redmine automator for performing rMLST analyses. Includes downloading and
updating databases from PubMLST, as well as running the rMLST analysis.
"""

# Standard imports
import datetime
import glob
import os
import pickle
import shutil
import subprocess
from typing import Optional

# Third party imports
import click
from nastools.nastools import retrieve_nas_files
import requests
import sentry_sdk

# Local imports
from methods import (
    construct_file_paths,
    construct_log_file_paths,
    create_shell_script,
    extract_job_id,
    make_executable,
    process_allele_data,
)


@click.command()
@click.option(
    '--redmine_instance', help='Path to pickled Redmine API instance'
)
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def rmlst_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    
    # Set the database path for the analyses
    database_path = '/mnt/nas2/databases/rMLST'
    
    # Set the database update boolean
    update = False
    
    
    try:
        # Set the name and path of the log files
        stdout_log_path, stderr_log_path = construct_log_file_paths(work_dir)
        
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            # Check to see if the database should be updated
            if 'update' in item:
                update = True
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)
        
        # Find the most recent version of the database
        database_version_path = find_most_recent_database_version(
            genus_database_path=database_path
        )
        
        # If the path doesn't exist or update is True, download the database
        if not os.path.exists(database_path) or update or \
                not database_version_path:
            # Create the path
            os.makedirs(database_path, exist_ok=True)
            
            # Get the current date
            current_date = datetime.datetime.now()
            
            # Set the version of the database (YYMMDD)
            database_version_path = os.path.join(
                database_path,
                current_date.strftime('%y%m%d')
            )
            
            # Create the path
            os.makedirs(database_version_path, exist_ok=True)
            
            # Create a Python script to run the download
            rest_call = """
\"\"\"
Run the REST class in rest_auth_class to download rMLST profile and alleles.
Combine alleles with combinealleles in GET class

\"\"\"

# Standard imports
from argparse import ArgumentParser
from glob import glob
import logging
import os

# Local imports
from olctools.databasesetup.get_rmlst import Get
from olctools.databasesetup.rest_auth_class import REST

# Define the path to the credentials
credentials= '/mnt/nas2/databases/keys/rmlst/'

# Create arguments to feed into the rest_auth_class script
args = ArgumentParser
args.secret_file = os.path.join(credentials, 'secret.txt')
args.file_path = os.path.join(credentials)
args.output_path = '{database_version_path}'
args.logging = logging

# Create the REST object
rmlst = REST(args)

# Download the profile and alleles
rmlst.main()

# Get the new alleles into a list
alleles = glob(os.path.join('{database_version_path}', '*.tfa'))

# Create the combinedAlleles file
Get.combinealleles('{database_version_path}', alleles)
            """.format(database_version_path=database_version_path)
            
            # Set the name and path of the Python script
            python_download_script = os.path.join(
                work_dir, 'rmlst_download.py'
            )
            
            # Write the script to a file
            with open(python_download_script, 'w', encoding='utf-8') as file:
                file.write(rest_call)

            # Make the script executable
            make_executable(path=python_download_script)
            
            # Set up the download script
            activate = (
                'source /home/ubuntu/miniconda3/bin/activate '
                '/mnt/nas2/virtual_environments/plasmid_borne'
            )
            
            # Create the command to run the Python rMLST download
            rmlst_download_cmd = \
                'python {script}'.format(script=python_download_script)

            # Set the name and path of the download script
            download_script = os.path.join(
                work_dir, 'run_rmlst_download.sh'
            )

            # Create the script file and make it executable
            template = create_shell_script(
                activate=activate,
                cmd=rmlst_download_cmd,
                script_path=download_script
            )

            # Update the Redmine issue after attempting to download
            # definitions
            notes = (
                'Starting rMLST download with command:\n {template}'.format(
                    template=template
                )
            )
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes=notes,
                status_id=2
            )

            # Open the log files
            with open(stdout_log_path, 'a') as stdout_file, \
                open(stderr_log_path, 'a') as stderr_file:
                # Run shell script using subprocess.run and redirect output
                subprocess.run(
                    ['bash', download_script],
                    stdout=stdout_file, stderr=stderr_file,
                    check=True
                )
            
            # # Reformat the database to work with GeneSeekr
            # process_allele_data(db_dir=database_version_path)
        
        # Normalize the path to ensure it ends with a slash
        normalized_path = os.path.join(database_version_path, '')
        
        # Get the parent directory of the normalized path
        parent_dir = os.path.dirname(normalized_path)
        
        # Extract the version
        version = os.path.basename(parent_dir)
        
        # Ensure that the database is present
        if not version:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='ERROR: Could not create the rMLST database',
                status_id=4
            )
            return
        
        # Update the Redmine issue with the version
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes='Using database version: {version}'.format(version=version)
        )
        
        # Set the name and path of the query directory
        query_dir = os.path.join(
            work_dir, 'query'
        )
        
        # Create the query directory
        os.makedirs(query_dir, exist_ok=True)
        
        # Check to see if there were attached files
        
        # Retrieve the issue again, but include any attachments
        issue_with_attachments = redmine_instance.issue.get(
            issue.id, include='attachments'
        )

        # Check if the issue has attachments
        if issue_with_attachments.attachments:
            # Download the attached FASTA files
            for attachment in issue_with_attachments.attachments:
                # Check if the file has one of the desired extensions
                if attachment.filename.endswith(
                        ('.fasta', '.fa', '.txt', '.fas')):
                    # Attempt to get the attachment by ID
                    file = redmine_instance.attachment.get(attachment.id)
                    # Attempt to download the file
                    file.download(savepath=query_dir, filename=file.filename)
        
        # Run file linker to link local FASTA files to the query directory
        retrieve_nas_files(
            seqids=seqids,
            outdir=query_dir,
            filetype='fasta',
            copyflag=False
        )
        
        # Make sure that all FASTA files requested are present
        missing_fastas = verify_fasta_files_present(seqids, query_dir)
        
        # Update the Redmine issue if one or more of the requested SEQIDs
        # could not be located
        if missing_fastas:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='WARNING: Could not find the following requested SEQIDs '
                'on the OLC NAS: {missing}'.format(missing=missing_fastas)
        )

        # These unfortunate hard coded paths appear to be necessary
        activate = (
            'source /home/ubuntu/miniconda3/bin/activate '
            '/mnt/nas2/virtual_environments/geneseekr'
        )
        seekr_py = '/mnt/nas2/virtual_environments/geneseekr/bin/GeneSeekr'
        
        # Run GeneSeekr in rMLST mode
        seekr_cmd = (
            'python {py} blastn -s {queries} -r {reports} -t {db_path} -R'
            .format(
                py=seekr_py,
                queries=query_dir,
                reports=os.path.join(work_dir, 'reports'),
                db_path=database_version_path,
            )
        )

        # Set the name and path of the GeneSeekr script
        rmlst_script = os.path.join(work_dir, 'run_rmlst_geneseekr.sh')
        
        # Create the script file and make it executable
        template = create_shell_script(
            activate=activate,
            cmd=seekr_cmd,
            script_path=rmlst_script
        )

        # Update the Redmine issue after attempting to download sequences
        notes = (
            'Starting rMLST analyses with command:\n {template}'.format(
                template=template
            )
        )
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes=notes,
            status_id=2
        )

        # Open the log files
        with open(stdout_log_path, 'a') as stdout_file, \
                open(stderr_log_path, 'a') as stderr_file:
            # Run shell script using subprocess.run and redirect output
            subprocess.run(
                ['bash', rmlst_script],
                stdout=stdout_file, stderr=stderr_file,
                check=True
            )
        
        # Zip output
        output_filename = 'rmlst_output_{}'.format(issue.id)
        zip_filepath = zip_folder(
            results_path=os.path.join(query_dir, 'reports'),
            output_dir=work_dir,
            output_filename=output_filename
        )
        zip_filepath += '.zip'
        
        # Prepare upload
        output_list = [
            {
                'filename': os.path.basename(zip_filepath),
                'path': zip_filepath
            }
        ]

        # Create a list of all the folders - will be used to clean up the
        # working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)
        
        # Wrap up issue
        redmine_instance.issue.update(
            resource_id=issue.id,
            uploads=output_list,
            status_id=4,
            notes='rMLST analysis with GeneSeekr complete!'
        )
    except Exception as exc:
        sentry_sdk.capture_exception(exc)
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes='Something went wrong! Please see the logs in the issue '
            'directory ({work_dir}) for more info. Error: {exc}'.format(
                work_dir=work_dir,
                exc=exc
            ),
        )


def verify_fasta_files_present(seqid_list, fasta_dir):
    """
    Makes sure that FASTQ files specified in seqid_list have been successfully
    copied/linked to directory specified by fastq_dir.
    :param seqid_list: List with SEQIDs.
    :param fasta_dir: Directory to which FASTA files should have been linked
    :return: List of SEQIDs that did not have files associated with them.
    """
    missing_fastas = list()
    for seqid in seqid_list:
        # Check forward.
        if len(
            glob.glob(
                os.path.join(
                    fasta_dir, '{seqid}*fasta*'.format(seqid=seqid)
                    )
                )
            ) == 0:
            missing_fastas.append(seqid)
    return missing_fastas


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


def find_most_recent_database_version(
        genus_database_path: str) -> Optional[str]:
    """
    Finds the most recent database version folder within a specified path.

    This function searches for subdirectories within the given path and
    returns the path to the most recent one based on the naming convention
    that assumes folders are sorted in ascending order.

    Parameters:
    genus_database_path (str): The path to the directory containing database
    version folders.

    Returns:
    Optional[str]: The path to the most recent database version folder, or
    None if no folders are found.
    """
    # Create a sorted list of all folders in the database
    folders = sorted(glob.glob(os.path.join(genus_database_path, '*/')))
    
    # Check to make sure that at least one folder is present
    if folders:
        # Find the most recent version of the database
        database_version_path = folders[-1]
        return database_version_path
    else:
        return None

if __name__ == '__main__':
    rmlst_redmine()
