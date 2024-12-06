"""
Redmine automator for downloading SRA files from the NCBI SRA database to the
NAS, uploading them to Azure, and assembling them on FoodPort
"""

# Standard imports
import os
import pickle
import re
import subprocess

# Third party imports
import click
import sentry_sdk

# Local imports
from automator_settings import (
    before_send,
    SENTRY_DSN
)
from methods import (
    api_request_populate,
    api_request_python,
    api_request_script,
    azure_automate,
    azure_upload,
    construct_file_paths,
    construct_log_file_paths,
    create_shell_script,
    extract_job_id,
    make_executable,
)

@click.command()
@click.option(
    '--redmine_instance', help='Path to pickled Redmine API instance'
)
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def sra_download_redmine(redmine_instance, issue, work_dir, description):
    """
    Download a list of SRA files from the NCBI SRA database to the NAS,
    upload them to Azure, and assembles them on FoodPort
    """
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)

    # Unpickle Redmine objects
    with open(redmine_instance, 'rb') as file:
        redmine_instance = pickle.load(file)

    with open(issue, 'rb') as file:
        issue = pickle.load(file)

    with open(description, 'rb') as file:
        description = pickle.load(file)

    # Variable to store sequence folder/container name
    container = str()
    
    # Variable to store email address
    email = str()
    
    # basic_assembly flag
    basic_assembly = False
    
    # preprocess flag
    preprocess = False

    try:
        # Set the name and path of the log files
        stdout_log_path, stderr_log_path = construct_log_file_paths(work_dir)
        
        # Initialise lists to store SEQID information
        seqids = []
        already_present_seqids = []
        download_seqids = []
        # Parse description to figure out what SRA files we need to run on.
        for item in description:

            # Remove whitespace from the items
            item = item.rstrip()
            if 'run_name' in item:
                # Extract the container name from the description. Set it
                # to lower case and replace underscores with hyphens
                container = item.split('=')[1].lower().replace('_', '-')
                continue
            if 'email' in item:
                # Extract the email address from the description
                email = item.split('=')[1]
                continue
            if 'basic_assembly' in item:
                # Set the basic_assembly flag to True
                basic_assembly = True
                continue
            if 'preprocess' in item:
                # Set the preprocess flag to True
                preprocess = True
                continue

            # Otherwise the item should be a SRA ID
            seqids.append(item)

        # Ensure that the analysis type is provided
        if not container:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='WARNING: no sequence folder/run name provided. '
                      'Please ensure that the first line of the issue contains'
                      ' run_name=YYMMDD-sra',
                status_id=3
            )
            return
        
        # Define the pattern for YYMMDD-sra naming scheme
        pattern = r'\d{6}-sra'

        # Use re to determine whether the pattern matches the container value
        if not re.match(pattern, container):
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='Invalid run name format. Please ensure the run name '
                      'follows the YYMMDD-sra naming scheme.',
                status_id=3
            )
            return

        # Ensure that both basic_assembly and preprocess are not True
        if basic_assembly and preprocess:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='Both basic_assembly and preprocess arguments were '
                      'provided. Please ensure you select a single option.',
                status_id=3
            )
            return

        # Download the files from SRA to the NAS

        # Set the absolute path to the download folder
        download_path = '/mnt/nas2/raw_sequence_data/ncbi/sra'

        # Check each seqid for existing files
        for seqid in seqids:
            # Construct file paths for R1 and R2
            r1_path, r2_path = construct_file_paths(
                download_path=download_path,
                seqid=seqid
            )

            # Check if both files exist
            if os.path.exists(r1_path) and os.path.exists(r2_path):
                # Log as already present
                already_present_seqids.append(seqid)
            else:
                # Add to download list
                download_seqids.append(seqid)

        # Update issue with already present seqids, if any
        if already_present_seqids:
            notes = (
                'The following seqids are already present in the '
                'SRA folder:  ' + ', '.join(already_present_seqids)
            )
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes=notes,
                status_id=2
            )

        # Only run the download script if there are SRA files to download
        if not download_seqids:
            # Update the Redmine issue if no new seqids need to be downloaded
            if already_present_seqids:
                # Construct the note about already present seqids
                notes = (
                    'All specified seqids are already present in the '
                    'SRA folder: ' + ', '.join(already_present_seqids)
                )
            else:
                # Construct a note indicating no seqids were specified
                # for download
                notes = 'No new seqids were specified for download.'

            # Update the Redmine issue
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes=notes,
                status_id=2
            )
        else:
            # Create the run info table for seqids to be downloaded
            run_info_table = os.path.join(
                download_path, '{issue}_run_info_table.csv'.format(
                    issue=issue.id
                )
            )

            # Write the SEQID information to the table
            with open(run_info_table, 'w', encoding='utf-8') as file:
                # Write header
                file.write('SampleName,Run\n')
                for seqid in download_seqids:
                    # Validate seqid before writing
                    if seqid:  # This checks if seqid is not empty or None
                        # Write each seqid to the run info table
                        file.write('{seqid},{seqid}\n'.format(seqid=seqid))

            # Set up the download script
            activate = (
                'source /home/ubuntu/miniconda3/bin/activate '
                '/mnt/nas2/virtual_environments/plasmid_borne'
            )

            # Create the download command
            sra_download_cmd = (
                'python -m genewrappers.biotools.sra_download '
                '-p {download_path} -r {run_info_table}'.format(
                    download_path=download_path,
                    run_info_table=run_info_table
                )
            )

            # Set the name and path of the download script
            download_script = os.path.join(work_dir, 'run_sra_download.sh')

            # Create the script file and make it executable
            template = create_shell_script(
                activate=activate,
                cmd=sra_download_cmd,
                script_path=download_script
            )

            # Update the Redmine issue after attempting to download sequences
            notes = (
                'Starting SRA download with command:\n {template}'.format(
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

        # Upload SRA files to Azure

        # Create a list of all the SRA files to be uploaded
        files_to_upload = []
        
        # Create a list to store any files not found
        missing_files = []

        for seqid in seqids:
            # Construct file paths for R1 and R2
            r1_path, r2_path = construct_file_paths(
                download_path=download_path,
                seqid=seqid
            )
            # Check if both files exist
            if os.path.exists(r1_path) and os.path.exists(r2_path):
                # Add both paths to the files_to_upload list
                files_to_upload.extend([r1_path, r2_path])
            else:
                # Add seqid to missing files list
                missing_files.append(seqid)

        # If there are missing files, update the Redmine issue
        if missing_files:
            missing_files_str = ', '.join(missing_files)
            notes = (
                "Missing files for the following SEQIDs: {missing_files_str}"
                .format(missing_files_str=missing_files_str)
            )
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes=notes,
                status_id=2
            )

        # Create a .tsv file for AzureAutomate
        azure_tsv = azure_automate(
            container=container,
            files_to_upload=files_to_upload,
            work_dir=work_dir
        )

        # Create the upload script
        azure_upload_script, template = azure_upload(
            azure_tsv=azure_tsv,
            work_dir=work_dir
        )

        # Update the Redmine issue after attempting to upload sequences
        notes = (
            'Uploading SRA files to container {container} with command:\n '
            '{template}'.format(
                container=container,
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
            result = subprocess.run(
                ['bash', azure_upload_script],
                stdout=stdout_file, stderr=stderr_file,
            )
        
        # Make an API request to FoodPort to run a research assembly on the
        # SRA files
        api_log_path, script_content, script_content_redacted = \
            api_request_populate(
            basic_assembly=basic_assembly,
            container=container,
            email=email,
            preprocess=preprocess,
            work_dir=work_dir
        )

        # Make a Python script to create the API call
        api_call_script = api_request_python(
            script_content=script_content,
            work_dir=work_dir
        )
        
        # Create a shell script to run the Python script
        foodport_api_script = api_request_script(
            api_call_script=api_call_script,
            work_dir=work_dir
        )

        # Open the log files
        with open(stdout_log_path, 'a') as stdout_file, \
             open(stderr_log_path, 'a') as stderr_file:
            # Run shell script using subprocess.run and redirect output
            subprocess.run(
                ['bash', foodport_api_script],
                stdout=stdout_file, stderr=stderr_file,
                check=True
            )
            
        # Open and read the .err log file
        with open(api_log_path, 'r') as file:
            log_content = file.read()

        # Define the error message to search for
        error_message = "Exception Value: duplicate key value violates " \
            "unique constraint"
        
        # Check if the specific error message is in the log file
        if error_message in log_content:
            # Define the notes to update the Redmine issue with
            notes = (
                "The supplied run_name, {container} already exists in "
                "FoodPort. Please check the details and submit a unique "
                "name.".format(container=container)
            )
            
            # Update the Redmine issue
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes=notes,
                status_id=3
            )
            return
            
        # Create notes for updating the Redmine issue after submitting the
        # API request
        notes = (
            'Running COWBAT on SRA files in container {container} with '
            'API call:\n {script_content_redacted}\nPlease follow the '
            'progress on FoodPort'.format(
                container=container,
                script_content_redacted=script_content_redacted
            )
        )
        
        # Add an additional note if an email address was provided
        if email:
            notes += (
                '\nAn email will be sent to {email} when the assembly is '
                'complete.'.format(email=email)
            )
        
        # Update the Redmine issue
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes=notes,
            status_id=4
        )
        
    except Exception as exc:
        sentry_sdk.capture_exception(exc)
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes=(
                'Something went wrong! Please contact a bioinformatician '
                'to investigate: {exc}'.format(
                    exc=str(exc)
                )
            ),
            status_id=3
        )


if __name__ == '__main__':
    sra_download_redmine()