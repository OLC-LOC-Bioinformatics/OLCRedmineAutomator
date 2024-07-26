"""
Redmine automator for performing plasmid MLST analyses. Includes downloading
and updating databases from PubMLST, as well as running the plasmid MLST
analysis.
"""

# Standard imports
import datetime
import glob
import os
import pickle
import shutil
import subprocess
import traceback
from typing import Optional

# Third party imports
from Bio import SeqIO
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
def plasmid_mlst_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    
    # Set the database path for the analyses
    database_path = '/mnt/nas2/databases/pMLST'
    
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
        
        # Create a list to store all paths to the different schemes
        scheme_paths = []
        
        # Create a dictionary to store the URLs for the profiles and
        # alleles
        urls = {
            'IncI1': {
                'scheme': '1',
                'alleles': [
                    'repI1',
                    'ardA',
                    'trbA',
                    'sogS',
                    'pilL',
                ]
            },
            'IncHI2': {
                'scheme': '2',
                'alleles': [
                    'smr0018',
                    'smr0199',
                ]
            },
            'IncF': {
                'scheme': '3',
                'alleles': [
                    'FII',
                    'FIC',
                    'FIIK',
                    'FIIS',
                    'FIIY',
                    'FIA',
                    'FIB'
                ]
            },
            'IncN': {
                'scheme': '4',
                'alleles': [
                    'repN',
                    'traJ',
                    'korA',
                ]
            },
            'IncHI1': {
                'scheme': '5',
                'alleles': [
                    'HCM1_043',
                    'HCM1_064',
                    'HCM1_099',
                    'HCM1_116',
                    'HCM1_178ac',
                    'HCM1_259',
                ]
            },
            'IncAC': {
                'scheme': '6',
                'alleles': [
                    'repA',
                    'parA',
                    'parB',
                    'A053',
                ]
            },
            'IncAC_cg': {
                'scheme': '7',
                'alleles': [
                    'repA',
                    'parA',
                    'parB',
                    'A053',
                    'A002',
                    'A004',
                    'A005',
                    'A007',
                    'A008',
                    'A009',
                    'A010',
                    'A011',
                    'A015',
                    'A017',
                    'A050',
                    'A054',
                    'A055',
                    'A056',
                    'A058',
                    'A059',
                    'A157',
                    'A160',
                    'A165',
                    'A167',
                    'A169',
                    'A172',
                    'A173',
                    'A175'
                ]
            },
            'pBSSB1-family': {
                'scheme': '8',
                'alleles': [
                    'higB',
                    'mqsA',
                    'soj',
                ]
            },
            'Shigella_flexneri': {
                'scheme': '9',
                'alleles': [
                    'vapBC',
                    'parAB',
                    'yacAB',
                ]
            },
        }
        
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
            
            # Download the profile and alleles manually
            base_url = 'https://rest.pubmlst.org/db/pubmlst_plasmid_seqdef/'
        
            # Iterate over each scheme
            for url, scheme_info_dict in urls.items():
            
                # Set the profile URL
                profile_url = base_url + 'schemes/{url}/profiles_csv'.format(
                    url=scheme_info_dict['scheme']
                )
                    
                print(
                    'Downloading profile file: {url}'.format(url=profile_url)
                )
                    
                # Set the scheme path
                scheme_path = os.path.join(
                    database_version_path, scheme_info_dict['scheme']
                )
                
                # Create the scheme path if it doesn't exist
                os.makedirs(scheme_path, exist_ok=True)
                
                # Append the scheme path to the list of scheme paths
                scheme_paths.append(scheme_path)
                
                # Set the path to the local profile file
                profile_path = os.path.join(
                    scheme_path, 'profiles.txt'
                )
                
                # Create a list to store the allele file paths
                allele_paths = []
                
                # Download the profile CSV
                response = requests.get(profile_url)

                # Replace underscores with nothing in the content
                modified_content = \
                    response.content.decode('utf-8').replace('_', '')

                # Write the modified content to the file
                with open(profile_path, 'w', encoding='utf-8') as file:
                    file.write(modified_content)
                    
                # Create the allele URLs
                for allele in scheme_info_dict['alleles']:
                    # Set the allele URL
                    allele_url = base_url + \
                        'loci/{allele}/alleles_fasta'.format(
                            allele=allele)
                        
                    print('Downloading file: {url}'.format(url=allele_url))
                    
                    # Set the path to the local allele file
                    allele_path = os.path.join(
                        scheme_path, '{allele}.tfa'.format(
                            allele=allele.replace('_', '')
                        )
                    )
                    
                    # Add the path to the list of allele paths
                    allele_paths.append(allele_path)
                    
                    # Use requests to download the allele
                    response = requests.get(allele_url)
                    with open(allele_path, 'wb') as file:
                        file.write(response.content)
                            
                # Create the combinedtargets.fasta file
                combined_path = os.path.join(
                    scheme_path, 'combinedtargets.fasta'
                )
                
                # Combine the allele files into a single file
                with open(combined_path, 'w') as combined_file:
                    for allele_path in allele_paths:
                        for record in SeqIO.parse(allele_path, 'fasta'):
                            # Split the record ID by underscores
                            parts = record.id.rsplit('_', 1)
                            # Replace underscores with nothing in all but
                            # the last part
                            parts[0] = parts[0].replace('_', '')
                            # Join the parts back together
                            record.id = '_'.join(parts)
                            # Clear the description
                            record.description = ''
                            # Write the modified record to the combined file
                            SeqIO.write(record, combined_file, 'fasta')
                
        else:
            # Create a list to store the paths to the different schemes
            scheme_paths = glob.glob(os.path.join(database_version_path, '*/'))
        
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
                notes='ERROR: Could not create the database',
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
        
        # Run GeneSeekr on each scheme
        for scheme_path in sorted(scheme_paths):
       
            # Create a variable to store the scheme name
            scheme = str()
            
            # Extract the scheme name from the urls dictionary
            for scheme_name, info_dict in urls.items():
                last_directory = os.path.basename(
                    os.path.normpath(scheme_path)
                )
                if last_directory == info_dict['scheme']:
                    scheme = scheme_name
                    break
                    
            print('Scheme: {scheme}'.format(scheme=scheme))
            
            # I don't have a profile file for IncF, so run a simple BLAST
            # search instead
            if scheme == 'IncF':
                seekr_cmd = (
                    'python {py} blastn -s {queries} -r {reports} '
                    '-t {db_path} -c 99 -u '
                    .format(
                        py=seekr_py,
                        queries=query_dir,
                        reports=os.path.join(work_dir, 'reports'),
                        db_path=scheme_path,
                    )
                )
            
            else:
                
                # Run GeneSeekr in MLST mode
                seekr_cmd = (
                    'python {py} blastn -s {queries} -r {reports} -t {db_path} -M'
                    .format(
                        py=seekr_py,
                        queries=query_dir,
                        reports=os.path.join(work_dir, 'reports'),
                        db_path=scheme_path,
                    )
                )

            # Set the name and path of the download script
            mlst_script = os.path.join(
                work_dir, 'run_mlst_geneseekr_{scheme}.sh'.format(
                    scheme=scheme
                )
            )
            
            # Create the script file and make it executable
            template = create_shell_script(
                activate=activate,
                cmd=seekr_cmd,
                script_path=mlst_script
            )

            # # Update the Redmine issue after attempting to download sequences
            # notes = (
            #     'Starting MLST analyses with command:\n {template}'.format(
            #         template=template
            #     )
            # )
            # redmine_instance.issue.update(
            #     resource_id=issue.id,
            #     notes=notes,
            #     status_id=2
            # )

            # Open the log files
            with open(stdout_log_path, 'a') as stdout_file, \
                    open(stderr_log_path, 'a') as stderr_file:
                # Run shell script using subprocess.run and redirect output
                subprocess.run(
                    ['bash', mlst_script],
                    stdout=stdout_file, stderr=stderr_file,
                    check=True
                )
            
            # Define the results folder to store all the outputs
            results_path = os.path.join(work_dir, 'results')
            
            # Create results_path if required
            os.makedirs(results_path, exist_ok=True)
            
            if scheme == 'IncF':
                query_reports = os.path.join(work_dir, 'reports')
            else:
                # Define the query/reports folder
                query_reports = os.path.join(query_dir, 'reports')
            
            # Ensure that the query/reports folder exists
            if os.path.exists(query_reports):
                
                # Copy the query_reports folder to results folder. Rename it
                # to scheme_reports
                shutil.copytree(
                    query_reports,
                    os.path.join(
                        results_path,
                        '{scheme}_reports'.format(scheme=scheme)
                    )
                )
                
            # Delete the reports and the query/reports folders
            folders = glob.glob(os.path.join(query_dir, '*/'))
            for folder in folders:
                shutil.rmtree(folder)
            
            shutil.rmtree(os.path.join(work_dir, 'reports'))
                
         # Zip output
        output_filename = 'plasmid_mlst_output_{issue}'.format(
            issue=issue.id
        )
        
        zip_filepath = zip_folder(
            results_path=results_path,
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
        
        # Wrap up issue
        redmine_instance.issue.update(
            resource_id=issue.id,
            uploads=output_list,
            status_id=4,
            notes='Plasmid MLST analysis with GeneSeekr complete!'
        )
        # Create a list of all the folders - will be used to clean up the
        # working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)
    except Exception as exc:
        sentry_sdk.capture_exception(exc)
        # Get the full traceback
        exc_traceback = traceback.format_exc()
        print(exc_traceback)
        # Update the Redmine issue with the full traceback
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes='Something went wrong! Please see the logs in the issue '
                'directory ({work_dir}) for more info. Error: {exc_traceback}'.format(
                    work_dir=work_dir,
                    exc_traceback=exc_traceback
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
    plasmid_mlst_redmine()
