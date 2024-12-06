"""
Collection of common methods for automators
"""

# Standard imports
from glob import glob
import os
import shutil
from typing import Tuple

# Third party imports
from Bio import SeqIO

# Local imports
from automator_settings import (
    POSTGRES_PASSWORD,
)

def construct_file_paths(download_path: str, seqid: str) -> Tuple[str, str]:
    """
    Construct the file paths for R1 and R2 fastq.gz files for a given seqid.

    This function takes a base download path and a sequence identifier (seqid)
    as inputs and returns the paths for the corresponding R1 and R2 fastq.gz
    files. These paths are constructed by appending the seqid and the
    appropriate suffix ('_R1.fastq.gz' for R1, '_R2.fastq.gz' for R2) to the
    download path.

    Parameters:
        download_path (str): The base path where files are downloaded.
        seqid (str): The sequence identifier for which to construct file paths.

    Returns:
        Tuple[str, str]: A tuple containing the file paths for the R1 and R2
                         files, in that order.

    Example:
        >>> construct_file_paths('/path/to/download', 'sample123')
        ('/path/to/download/sample123_R1.fastq.gz',
        '/path/to/download/sample123_R2.fastq.gz')
    """
    # Construct the file path for the R1 file using the seqid
    r1_path = os.path.join(
        download_path,
        '{seqid}_R1.fastq.gz'.format(
            seqid=seqid
        )
    )
    # Construct the file path for the R2 file using the seqid
    r2_path = os.path.join(
        download_path, 
        '{seqid}_R2.fastq.gz'.format(
            seqid=seqid
        )
    )

    return r1_path, r2_path


def create_shell_script(
        activate: str,
        cmd: str,
        script_path: str) -> str:
    """
    Creates a shell script with the specified activation and command,
    makes the script executable using the make_executable function, and returns
    the shell script contents.

    Args:
        activate (str): The command to activate the environment.
        cmd (str): The SRA download command to be executed.
        script_path (str): The path where the shell script will be saved.

    Returns:
        str: The contents of the shell script.

    Example:
        >>> activate = 'source /path/to/activate environment'
        >>> cmd = 'python -m module -arg'
        >>> script_path = '/path/to/script.sh'
        >>> script_content = create_shell_script(activate, cmd,
        script_path)
        >>> print(script_content)
    """
    # Create the shell script content
    template = '#!/bin/bash\n{activate} && {cmd}'.format(
        activate=activate,
        cmd=cmd
    )

    # Write the shell script to the specified path
    with open(script_path, 'w', encoding='utf-8') as file:
        file.write(template)

    # Make the shell script executable using the provided function
    make_executable(script_path)

    # Return the shell script contents
    return template


def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode

    # This copies read bits to execute bits
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)


def extract_job_id(work_dir):
    """
    Extracts the SLURM job ID from the working directory path.
    
    Assumes the last part of the working directory path is the SLURM job ID.
    
    :param work_dir: Path to the working directory
    :return: The extracted SLURM job ID as a string
    """
    # Extract the last component of the work_dir path
    job_id = os.path.basename(work_dir)
    return job_id


def construct_log_file_paths(work_dir):
    """
    Constructs the file paths for the SLURM job's standard output and error
    logs based on the working directory.
    
    :param work_dir: Path to the working directory
    :return: A tuple containing the paths to the standard output and error
             log files
    """
    job_id = extract_job_id(work_dir)
    # Construct the paths for stdout and stderr log files
    stdout_log_path = os.path.join(
        work_dir, "job_{job_id}.out".format(job_id=job_id)
    )
    stderr_log_path = os.path.join(
        work_dir, "job_{job_id}.err".format(job_id=job_id)
    )
    return stdout_log_path, stderr_log_path

def process_allele_data(db_dir):
    """
    Processes allele data by creating a dictionary of alleles for each gene
    and generating .tfa files for each gene in the specified directory.

    Parameters:
    db_dir (str): Directory where the profile copy and .tfa files will be
        saved.
    """
    
    # Set the name of the profile file
    profile_file = os.path.join(db_dir, 'profiles_csv')
    
    # Read gene names from the first line of the profile file, excluding 'ST'
    # and names containing 'clonal'
    with open(profile_file, 'r', encoding='utf-8') as profile:
        # Read the first line and split the headers by tab
        headers = profile.readline().strip().split('\t')
        
        # Create a list of gene names excluding 'ST' and 'clonal' names
        gene_names = [
            gene for gene in headers if gene != 'ST' and 'clonal' not in gene
        ]

    # Copy the profile file to db_dir named profile.txt
    shutil.copy(profile_file, os.path.join(db_dir, 'profile.txt'))

    # Set the glob pattern to match all .fasta files and then filter out
    # 'combinedtargets.fasta'
    allele_files = glob(os.path.join(db_dir, '*.fasta'))
    allele_file = next(
        (f for f in allele_files if 'combinedtargets.fasta' not in f), None
    )
    
    # Parse allele_file to create a dictionary of alleles for each gene
    alleles_dict = {gene: [] for gene in gene_names}
    for record in SeqIO.parse(allele_file, 'fasta'):
        # Allele IDs are formatted as gene_allele
        gene_id = record.id.split('_')[0]
        if gene_id in alleles_dict:
            alleles_dict[gene_id].append(record)

    # Create a .tfa file for each gene with its alleles in db_dir
    for gene, alleles in alleles_dict.items():
        
        # Set the path for the .tfa file for the current gene
        tfa_path = os.path.join(db_dir, '{}.tfa'.format(gene))
        with open(tfa_path, 'w', encoding='utf-8') as tfa_file:
            
            # Iterate over the alleles for the current gene and write them to
            # the .tfa file
            for allele in alleles:
                SeqIO.write(allele, tfa_file, 'fasta')
                
    # Copy the allele_file to 'combinedtargets.fasta' in the db_dir
    shutil.copy(allele_file, os.path.join(db_dir, 'combinedtargets.fasta'))


def azure_automate(
    container: str,
    files_to_upload: list,
    work_dir: str
):
    """
    Creates a .tsv file for uploading files to Azure Blob Storage.
    
    Parameters:
        container (str): The name of the Azure Blob Storage container to upload
            the files to.
        files_to_upload (list): A list of file paths to upload to Azure Blob
            Storage.
        work_dir (str): The working directory where the .tsv file will be
            created.
    
    Returns:
        str: The path to the .tsv file containing the information for uploading
            files to Azure Blob
    """
    
    # Set the path for the .tsv file
    azure_tsv = os.path.join(work_dir, 'azure_upload.tsv')

    # Write the necessary information to the .tsv file
    with open(azure_tsv, 'w', encoding='utf-8') as file:
        for upload_file in files_to_upload:
            file.write(
                '{container}\t{upload_file}\t""\n'.format(
                    container=container,
                    upload_file=upload_file
                )
            )
    
    return azure_tsv


def azure_upload(
    azure_tsv: str,
    work_dir: str
):
    """
    Creates a shell script for uploading files to Azure Blob Storage.
    
    Parameters:
        azure_tsv (str): The path to the .tsv file containing the information
            for uploading files to Azure Blob Storage.
        work_dir (str): The working directory where the shell script will be
            created.
            
    Returns:
        Tuple[str, str]: A tuple containing the path to the shell script and
        the contents of the shell script.
    """
    
    # Create the Python interpreter source command
    activate = (
        'source /home/ubuntu/miniconda3/bin/activate '
        '/mnt/nas2/virtual_environments/azurestorage'
    )

    # Create the Azure upload command
    azure_upload_cmd = (
        'AzureAutomate upload file -a carlingst01 -b {azure_tsv}'.format(
            azure_tsv=azure_tsv
        )
    )

    # Set the name and path of the upload script
    azure_upload_script = os.path.join(work_dir, 'sra_azure_upload.sh')

    # Create the script file and make it executable
    template = create_shell_script(
        activate=activate,
        cmd=azure_upload_cmd,
        script_path=azure_upload_script
    )
    
    return azure_upload_script, template


def api_request_populate(
    basic_assembly: str,
    container: str,
    email: str,
    preprocess: str,
    work_dir: str,
):
    """
    Populate the script content template with the necessary information for
    making an API request to the Research Assembly API.
    
    Parameters:
        basic_assembly (str): The name of the basic assembly.
        container (str): The name of the container.
        email (str): The email address to send the API response to.
        preprocess (str): The name of the preprocess.
        work_dir (str): The working directory.
        
    Returns:
        Tuple[str, str]: A tuple containing the script content and the
        redacted script content.
            
    """
    # Base script content template
    script_content_template = (
        "import json\n"
        "import requests\n\n"
        "API_ENDPOINT = "
        "\"http://10.148.57.4/en-ca/api/research_assembly/\"\n"
        "data = {\n"
        "    'run_name': '%s',\n"
        "    'basic_assembly': '%s',\n"
        "    'preprocess': '%s',\n"
        "    'nextseq': False,\n"
    )

    # Check if the email variable exists and is not empty
    if email:
        # Append the email line to the template
        script_content_template += "    'email': '{}',\n".format(email)
    
    # Define the API log path variable outside the template for clarity
    api_log_path = os.path.join(work_dir, 'api_response_log.txt')
    
    # Finish the template with the request, print statements, and log file
    # writing
    script_content_template += (
        "}\n"
        "headers = {'Content-Type': 'application/json'}\n"
        "response = requests.post(\n"
        "    API_ENDPOINT, data=json.dumps(data), headers=headers,\n"
        "    auth=('%s', '%s'),\n"
        "    verify=False, allow_redirects=True)\n\n"
        
        # Add print statements
        "print('Status Code:', response.status_code)\n"
        "print('Response Content:', response.content.decode('utf-8'))\n\n"
        
        # Open log file in append mode and write the response
        "with open('%s', 'a') as log_file:\n"
        "    log_file.write('Response Content: ' + \\\n"
        "    response.content.decode('utf-8') + '\\n\\n')\n"
    )

    # Format the script content with actual credentials for use
    script_content = script_content_template % (
        container,
        basic_assembly,
        preprocess,
        "olcbioinformatics@gmail.com",
        POSTGRES_PASSWORD,
        api_log_path
    )

    # Format the script content with 'REDACTED' for the Redmine issue
    # update
    script_content_redacted = script_content_template % (
        container,
        basic_assembly,
        preprocess,
        "REDACTED",
        "REDACTED",
        api_log_path
    )
    
    return api_log_path, script_content, script_content_redacted


def api_request_python(
    script_content: str,
    work_dir: str
):
    """
    Create a Python script for making an API request to the Research Assembly
    API.
    
    Parameters:
        script_content (str): The content of the Python script that makes the
            API request.
        work_dir (str): The working directory where the Python script will be
            created.

    Returns:
        str: The path to the Python script that makes the API request.
    """
    # Set the path for the Python script
    api_call_script = os.path.join(work_dir, 'foodport_api_call.py')

    # Write the script to a file
    with open(api_call_script, 'w', encoding='utf-8') as file:
        file.write(script_content)

    # Make the script executable
    make_executable(path=api_call_script)
    
    return api_call_script


def api_request_script(
    api_call_script: str,
    work_dir: str
):
    """
    Create a shell script for making an API request to the Research Assembly
    API.
    
    Parameters:
        api_call_script (str): The path to the Python script that makes the API
            request.
        work_dir (str): The working directory where the shell script will be
            created.
    
    Returns:
        Tuple[str, str]: A tuple containing the path to the shell script and
        the contents of the shell script.
    """
    
    # Create the Python interpreter source command
    activate = (
        'source /home/ubuntu/miniconda3/bin/activate '
        '/mnt/nas2/virtual_environments/azurestorage'
    )

    # Create the API request
    foodport_api_cmd = (
        'python {api_call_script}'.format(
            api_call_script=api_call_script
        )
    )

    # Set the name and path of the upload script
    foodport_api_script = os.path.join(work_dir, 'foodport_api.sh')

    # Create the script file and make it executable
    create_shell_script(
        activate=activate,
        cmd=foodport_api_cmd,
        script_path=foodport_api_script
    )
    
    return foodport_api_script