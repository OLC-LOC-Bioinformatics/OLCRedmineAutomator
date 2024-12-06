#!/usr/env/python3

"""
Automator for downloading assembled sequencing runs from FoodPort
"""

# Local imports
import csv
from glob import glob
import os
import pickle
import re
import shutil
import traceback
from typing import (
    Dict,
    List,
    Tuple
)

# Third-party imports
from biotools import mash
import click
from nastools.nastools import retrieve_nas_files
import sentry_sdk

# Local imports
from amrsummary import before_send
from automator_settings import SENTRY_DSN


@click.command()
@click.option(
    '--redmine_instance', help='Path to pickled Redmine API instance'
)
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def azure_download_redmine(redmine_instance, issue, work_dir, description):
    """
    The automator for downloading sequencing runs from an Azure container
    """
    # Initialise sentry
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)

    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    # Current list of supported sequencing machines
    machines = [
        'miseq', 'nextseq', 'other'
    ]

    # List of locations in which data can be saved
    locations = [
        'atcc', 'enterobase', 'ncbi', 'merged', 'other', 'refseq'
    ]

    # Variables to hold supplied arguments
    run_name = str()
    machine = str()
    location = False
    output_list = []

    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = []
        for item in description:
            # Remove any trailing whitespace
            item = item.rstrip()

            # Extract the supplied arguments
            if 'sequence_folder' in item:
                run_name = item.split('=')[1].lower()
                continue
            if 'machine' in item:
                machine = item.split('=')[1].lower()
                continue
            if 'location' in item:
                location = item.split('=')[1].lower()
                continue

            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that the sequence folder was provided
        if not run_name:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='WARNING: Could not identify an sequence folder/run '
                'name. Please ensure that the first line of the issue contains'
                ' sequence_folder=YYMMDD-m05722 (or -LAB)',
                status_id=4
            )
            return

        # Ensure that the sequencing machine was provided
        if not machine:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='WARNING: Could not identify a sequencer/machine type. '
                'Please ensure that the first line of the issue contains '
                ' machine= and one of the following keywords: {ats}'.format(
                    ats=', '.join(machines)
                ),
                status_id=4
            )
            return
        
        # Ensure that the supplied machine is supported
        elif machine not in machines:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='WARNING: supplied machine type {at} currently not in '
                'the supported list of machines: {ats}'.format(
                    at=machine,
                    ats=', '.join(machines)
                ),
                status_id=4
            )
            return
        
        # Ensure that if "other" is chosen, user also supplies the location
        if machine == 'other':

            # Ensure that the location argument was provided
            if not location:
                redmine_instance.issue.update(
                    resource_id=issue.id,
                    notes='WARNING: Could not identify a location to upload '
                    'data to. Please ensure that the issue contains location= '
                    'and one of the following keywords: {ats}'.format(
                        ats=', '.join(locations)
                    ),
                    status_id=4
                )
                return
            # Ensure that a valid location argument was provided
            elif location not in locations:
                redmine_instance.issue.update(
                    resource_id=issue.id,
                    notes='WARNING: supplied location {at} currently not in '
                    'the supported list of locations: {ats}'.format(
                        at=location,
                        ats=', '.join(locations)
                    ),
                    status_id=4
                )
                return

        # Create the local folder(s) that we'll need.
        local_folder, local_raw_folder = _set_local_folders(
            location=location,
            machine=machine,
            run_name=run_name
        )

        # Get the activation command
        activate_command = _create_activate_command()
        
        # Handle MiSeq or other downloads
        if machine in ['miseq', 'other']:
            # Create the local_raw_folder_directory
            os.makedirs(local_raw_folder, exist_ok=True)

            # Create the download command
            azure_download_command = _create_download_command(
                sequence_folder=run_name,
                work_dir=work_dir
            )
            
            # Create and run the download script
            _create_and_run_azure_script(
                activate_command=activate_command,
                azure_command=azure_download_command,
                command_type='Download',
                sequence_folder=run_name,
                work_dir=work_dir
            )

            # Move the download files to the appropriate folders
            output_list = _process_downloaded_data(
                local_folder=local_folder,
                local_raw_folder=local_raw_folder,
                machine=machine,
                run_name=run_name,
                work_dir=work_dir
            )
        
        # Handle NextSeq downloads
        if machine == 'nextseq':
            # Create a command to list the containers that match the supplied
            # sequence_folder argument
            azure_list_command, output_file = _list_azure_containers(
                sequence_folder=run_name,
                work_dir=work_dir
            )

            # Create and run the listing script
            _create_and_run_azure_script(
                activate_command=activate_command,
                azure_command=azure_list_command,
                command_type='List',
                sequence_folder=run_name,
                work_dir=work_dir
            )

            # Parse the outputs from AzureList
            nextseq_containers = _parse_azure_list_output(
                output_file=output_file
            )
            
            # Download the assemblies for each sub-run
            for iterator, container in enumerate(nextseq_containers):

                # Create the local folder(s) that we'll need.
                local_folder, local_raw_folder = _set_local_folders(
                    location=location,
                    machine=machine,
                    run_name=container
                )

                # Create directories
                os.makedirs(local_folder, exist_ok=True)
                os.makedirs(local_raw_folder, exist_ok=True)
                
                # Create the download command
                azure_download_command = _create_download_command(
                    sequence_folder=container,
                    work_dir=work_dir
                )
                
                # Create and run the download script
                _create_and_run_azure_script(
                    activate_command=activate_command,
                    azure_command=azure_download_command,
                    command_type='Download-{iterator}'.format(
                        iterator=str(iterator)
                    ),
                    sequence_folder=container,
                    work_dir=work_dir
                )

                # Process the downloaded data for this sub-run
                _process_downloaded_data(
                    local_folder=local_folder,
                    local_raw_folder=local_raw_folder,
                    machine=machine,
                    run_name=container,
                    work_dir=work_dir
                )

        # After processing all sub-runs, combine the SampleSheets and
        # metadata files
        output_list = _combine_nextseq_files(
            work_dir=work_dir,
            local_raw_folder=local_raw_folder
        )

        # Wrap up issue
        redmine_instance.issue.update(
            resource_id=issue.id,
            uploads=output_list,
            status_id=4,
            assigned_to_id=787,  # Assign to Bridgette Kelly
            notes='AzureDownload complete!\n\n'
            'Legacy combined metadata reports are attached. '
            'Other files can be downloaded using the report retrieve function.'
        )

    except Exception as exc:
        sentry_sdk.capture_exception(exc)

        # Capture the full traceback as a formatted string
        exc_info = traceback.format_exc()

        # Update the Redmine issue with the error message
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes='Something went wrong! Please contact a bioinformatician '
            'to investigate: {exc_info}'.format(exc_info=exc_info)
        )


def make_executable(
    *,
    path: str
) -> None:
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)


def _set_local_folders(
    *,
    location: str,
    machine: str,
    run_name: str
) -> Tuple[str, str]:
    """
    Create the local directories for the supplied sequence folder

    Args:
        location (str): Location of the data (ATCC, Enterobase, NCBI, etc.).
        Used when 'other' is selected as the machine type.
        machine (str): Type of sequencing machine used to generate data
        run_name (str): Name of the sequencing run

    Returns:
        local_folder (str): Path to the local folder for storing the assembled
        data
        local_raw_folder (str): Path to the local folder for storing the raw
        sequence data
    """
    
    # Base directories
    base_processed = '/mnt/nas2/processed_sequence_data'
    base_raw = '/mnt/nas2/raw_sequence_data'

    # Mapping for processed directories based on machine or location
    processed_map = {
        'miseq': 'miseq_assemblies',
        'nextseq': 'nextseq_assemblies',
        'atcc': 'atcc',
        'enterobase': 'enterobase_assemblies',
        'ncbi': 'ncbi',
        'merged': 'merged_assemblies',
        'other': 'other_assemblies',
        'refseq': 'refseq_assemblies'
    }

    # Mapping for raw directories based on location
    raw_map = {
        'miseq': 'miseq',
        'nextseq': 'nextseq',
        'atcc': 'atcc',
        'enterobase': 'enterobase',
        'ncbi': 'ncbi',
        'merged': 'merged_sequences',
        'other': 'other',
        'refseq': 'refseq'
    }

    # Determine processed_subdir and raw_subdir based on machine type
    if machine in ['miseq', 'nextseq']:
        location = machine
        processed_subdir = processed_map[machine]
    elif machine == 'other':
        location = location
        processed_subdir = processed_map[location]

    # Set the raw sequence data directory
    raw_subdir = raw_map[location]

    # Construct full paths
    local_folder = os.path.join(
        base_processed, processed_subdir
    )
    local_raw_folder = os.path.join(
        base_raw,
        raw_subdir,
        run_name
    )

    # Return the local_folder and local_raw_folder variables
    return local_folder, local_raw_folder


def _create_and_run_azure_script(
    *,
    activate_command: str,
    azure_command: str,
    command_type: str,
    sequence_folder: str,
    work_dir: str
):
    """
    Create and execute a shell script to run AzureStorage commands.

    Args:
        activate_command (str): Command to activate the Conda environment.
        azure_command (str): Azure command
        command_type (str): Type of Azure command (Download, List, etc.) and
        any other identifiers (e.g. Download-1, etc.)
        sequence_folder (str): The sequence folder name.
        work_dir (str): The working directory where the script will be created.

    Returns:
        None
    """

    # Create the shell script content
    template = (
        "#!/bin/bash\n{activate_command} && cd {work_dir} && "
        "{azure_command}".format(
            activate_command=activate_command,
            azure_command=azure_command,
            work_dir=work_dir
        )
    )

    # Set the path of the script 
    azure_script = os.path.join(
        work_dir,
        'run_Azure{command_type}.sh'.format(command_type=command_type)
    )

    # Write the shell script
    with open(azure_script, 'w+', encoding='utf-8') as file:
        file.write(template)

    # Make the script executable
    make_executable(path=azure_script)

    # Execute the script
    os.system(azure_script)


def _create_activate_command() -> str:
    """
    Creates the AzureStorage environment activation command

    Args:
        None

    Returns:
        str: string of the AzureStorage environment activation command
    """
    # Activation command
    activate_command = ('source /home/ubuntu/miniconda3/bin/activate '
        '/mnt/nas2/virtual_environments/azurestorage'
    )
   
    return activate_command


def _create_download_command(
    *,
    sequence_folder: str,
    work_dir: str
 ) -> str:
    """
    Create the command to download the sequencing run from Azure storage

    Args:
        sequence_folder (str): The name of the container in Azure storage
        to download
        work_dir (str): The path to the working directory for this issue

    Returns:
        str, str: Sting of the azure download command
    """
    
    # AzureDownload command
    azure_download_command = (
        'AzureDownload container -a carlingst01 -c {run} -o {output}'.format(
            output=work_dir,
            run=sequence_folder
        )
    )

    return azure_download_command


def _list_azure_containers(
    *,
    sequence_folder: str,
    work_dir: str
) -> Tuple[str, str]:
    """
    Create the command to list the containers in the Azure storage account that
    match the sequence_folder argument

    Args:
        sequence_folder (str): The name of the sequencing run
        work_dir (str): The path to the working directory for this issue

    Returns:
        Tuple[str, str]: string of the AzureList command, string of the path
        to the output file
    """
    # Set the name of the output file
    output_file = os.path.join(
        work_dir,
        '{sequence_folder}_containers.txt'.format(
            sequence_folder=sequence_folder
        )
    )

    # Delete the output file if it exists (this is only necessary when
    # testing, as the file should never be created more than once otherwise)
    if os.path.exists(output_file):
        os.remove(output_file)

    # AzureList command
    azure_list_command = (
        'AzureList container -a carlingst01 -o {output_file} '
        '{sequence_folder}\*'.format(
            output_file=output_file,
            sequence_folder=sequence_folder
        )
    )

    return azure_list_command, output_file


def _parse_azure_list_output(
    output_file: str
) -> List:
    """
    Extract the containers that match the NextSeq run name. Return all the
    sub-runs (e.g. 241112-nb501933-1, 241112-nb501933-2, etc. but not 
    241112-nb501933)

    Args:
        output_file (str): name and absolute path to the output file from
        AzureList with all the matching container names

    Returns:
        List: list of the extracted matching containers
    """
    # Regular expression pattern to match container names ending with -<number>
    pattern = re.compile(r'.+-\d+$')

    # Initialize a list to store the filtered sub-run names
    containers = []

    # Open the file and read its contents
    with open(output_file, 'r', encoding='utf-8') as file:
        raw_containers = file.readlines()

    # Iterate through each line/container name
    for container in raw_containers:
        # Strip whitespace and newline characters
        clean_container = container.strip()

        # Check if the container name matches the pattern
        if pattern.match(clean_container):
            containers.append(clean_container)

    return containers


def _process_downloaded_data(
    *,
    local_folder: str,
    local_raw_folder: str,
    run_name: str,
    work_dir: str,
    machine: str
) -> List[Dict[str, str]]:
    """
    Process the downloaded sequencing data by organizing files into their
    appropriate directories.

    Args:
        local_folder (str): Path to the processed sequence data directory.
        local_raw_folder (str): Path to the raw sequence data directory.
        run_name (str): Name of the sequencing run.
        work_dir (str): The working directory where files are downloaded.
        machine (str): Type of sequencing machine ('miseq', 'nextseq', or
        'other').

    Returns:
        List[Dict[str, str]]: A list containing information about the legacy
            report file copied to the Redmine request.
    """
    # Determine the downloaded folder path
    downloaded_folder = os.path.join(work_dir, run_name)

    # List files in the downloaded folder
    output_files = os.listdir(downloaded_folder)

    # Move .fastq.gz files to the raw data folder
    for file in output_files:
        if file.endswith(".fastq.gz"):
            shutil.move(
                os.path.join(downloaded_folder, file),
                os.path.join(local_raw_folder, file)
            )

    # Copy the legacy report file to the working directory
    legacy_report = 'legacy_combinedMetadata'
    legacy_report_src = os.path.join(
        downloaded_folder,
        'reports',
        legacy_report
    )
    legacy_report_dest = os.path.join(
        work_dir,
        '{legacy_report}_{run_name}.csv'.format(run_name=run_name)
    )

    # Check if the legacy report file exists before copying
    if os.path.exists(legacy_report_src):
        shutil.copyfile(legacy_report_src, legacy_report_dest)

    # Prepare the output list with legacy report information
    output_list = []
    if os.path.exists(legacy_report_dest):
        output_dict = {
            'path': legacy_report_dest,
            'filename': '{legacy_report}_{run_name}.csv'.format(
                legacy_report=legacy_report,
                run_name=run_name
            )
        }
        output_list.append(output_dict)

    # Move remaining files to the processed sequence data folder
    shutil.move(downloaded_folder, local_folder)

    return output_list


def _combine_nextseq_files(
    *,
    local_raw_folder: str,
    work_dir: str
) -> List[Dict[str, str]]:
    """
    Combine legacy_combinedMetadata.csv files from NextSeq sub-runs.

    Args:
        local_raw_folder (str): Path to the raw sequence data directory.
        work_dir (str): The working directory where files are downloaded.

    Returns:
        List[Dict[str, str]]: A list containing information about the combined
            legacy report file copied to the Redmine request.
    """

    # Combine the legacy_combinedMetadata.csv files
    metadata_files = glob(
        os.path.join(
            work_dir,
            'legacy_combinedMetadata_*.csv'
        )
    )
    combined_metadata = os.path.join(
        work_dir,
        'legacy_combinedMetadata.csv'
    )

    # Write all the information from the individual combined metadata files
    with open(combined_metadata, 'w', encoding='utf-8') as out_file:
        for iterator, metadata_file in enumerate(metadata_files):
            with open(metadata_file) as current_metadata:
                content = current_metadata.read()
                # Skip headers for subsequent files
                if iterator > 0:
                    content = '\n'.join(content.strip().split('\n')[1:])
                out_file.write(content)
                out_file.write('\n')

    # Prepare the output list with combined legacy report information
    output_list = []
    if os.path.exists(combined_metadata):
        output_dict = {
            'path': combined_metadata,
            'filename': 'legacy_combinedMetadata.csv'
        }
        output_list.append(output_dict)

    return output_list


if __name__ == '__main__':
    azure_download_redmine()
