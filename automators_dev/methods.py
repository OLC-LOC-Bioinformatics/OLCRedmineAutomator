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
            and 'species' not in gene
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