"""
Redmine automator for performing MLST analyses. Includes downloading and
updating databases from PubMLST, as well as running the MLST analysis.
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
def mlst_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    
    # Set the database path for the analyses
    database_path = '/mnt/nas2/databases/MLST'
    
    # Set the organism to None
    organism = None
    
    # Set the database update boolean
    update = False
    
    species_list = [
        'Achromobacter spp.',
        'Acinetobacter baumannii#1',
        'Acinetobacter baumannii#2',
        'Aeromonas spp.',
        'Aggregatibacter actinomycetemcomitans',
        'Anaplasma phagocytophilum',
        'Arcobacter spp.',
        'Aspergillus fumigatus',
        'Bacillus cereus',
        'Bacillus licheniformis',
        'Bacillus subtilis',
        'Bacteroides fragilis',
        'Bartonella bacilliformis',
        'Bartonella henselae',
        'Bartonella washoensis',
        'Bordetella spp.',
        'Borrelia spp.',
        'Brachyspira hampsonii',
        'Brachyspira hyodysenteriae',
        'Brachyspira intermedia',
        'Brachyspira pilosicoli',
        'Brachyspira spp.',
        'Brucella spp.',
        'Burkholderia cepacia complex',
        'Burkholderia pseudomallei',
        'Campylobacter concisus/curvus',
        'Campylobacter fetus',
        'Campylobacter helveticus',
        'Campylobacter hyointestinalis',
        'Campylobacter insulaenigrae',
        'Campylobacter jejuni',
        'Campylobacter lanienae',
        'Campylobacter lari',
        'Campylobacter sputorum',
        'Campylobacter upsaliensis',
        'Candida albicans',
        'Candida glabrata',
        'Candida krusei',
        'Candida tropicalis',
        'Candidatus Liberibacter solanacearum',
        'Carnobacterium maltaromaticum',
        'Chlamydiales spp.',
        'Citrobacter freundii',
        'Clonorchis sinensis',
        'Clostridioides difficile',
        'Clostridium botulinum',
        'Clostridium perfringens',
        'Clostridium septicum',
        'Corynebacterium diphtheriae',
        'Cronobacter spp.',
        'Cutibacterium acnes',
        'Dichelobacter nodosus',
        'Edwardsiella spp.',
        'Enterobacter cloacae',
        'Enterococcus faecalis',
        'Enterococcus faecium',
        'Escherichia coli#1',
        'Escherichia coli#2',
        'Flavobacterium psychrophilum',
        'Gallibacterium anatis',
        'Geotrichum spp.',
        'Glaesserella parasuis',
        'Haemophilus influenzae',
        'Helicobacter cinaedi',
        'Helicobacter pylori',
        'Helicobacter suis',
        'Kingella kingae',
        'Klebsiella aerogenes',
        'Klebsiella oxytoca',
        'Klebsiella pneumoniae',
        'Kudoa septempunctata',
        'Lactobacillus salivarius',
        'Lactococcus lactis bacteriophage',
        'Leptospira spp.',
        'Leptospira spp.#2',
        'Leptospira spp.#3',
        'Listeria monocytogenes',
        'Macrococcus canis',
        'Macrococcus caseolyticus',
        'Mammaliicoccus sciuri',
        'Mannheimia haemolytica',
        'Melissococcus plutonius',
        'Moraxella catarrhalis',
        'Mycobacteria spp.',
        'Mycobacteroides abscessus',
        'Mycoplasma agalactiae',
        'Mycoplasma anserisalpingitidis',
        'Mycoplasma bovis',
        'Mycoplasma flocculare',
        'Mycoplasma gallisepticum#1',
        'Mycoplasma gallisepticum#2',
        'Mycoplasma hominis',
        'Mycoplasma hyopneumoniae',
        'Mycoplasma hyorhinis',
        'Mycoplasma iowae',
        'Mycoplasma pneumoniae',
        'Mycoplasma synoviae',
        'Neisseria spp.',
        'Orientia tsutsugamushi',
        'Ornithobacterium rhinotracheale',
        'Paenibacillus larvae',
        'Pasteurella multocida#1',
        'Pasteurella multocida#2',
        'Pediococcus pentosaceus',
        'Photobacterium damselae',
        'Piscirickettsia salmonis',
        'Porphyromonas gingivalis',
        'Proteus spp.',
        'Providencia spp.',
        'Pseudomonas aeruginosa',
        'Pseudomonas fluorescens',
        'Pseudomonas putida',
        'Rhodococcus spp.',
        'Riemerella anatipestifer',
        'Salmonella enterica',
        'Saprolegnia parasitica',
        'Serratia spp.',
        'Shewanella spp.',
        'Sinorhizobium spp.',
        'Staphylococcus aureus',
        'Staphylococcus chromogenes',
        'Staphylococcus epidermidis',
        'Staphylococcus haemolyticus',
        'Staphylococcus hominis',
        'Staphylococcus lugdunensis',
        'Staphylococcus pseudintermedius',
        'Stenotrophomonas maltophilia',
        'Streptococcus agalactiae',
        'Streptococcus bovis/equinus complex (SBSEC)',
        'Streptococcus canis',
        'Streptococcus dysgalactiae equisimilis',
        'Streptococcus gallolyticus',
        'Streptococcus oralis',
        'Streptococcus pneumoniae',
        'Streptococcus pyogenes',
        'Streptococcus suis',
        'Streptococcus thermophilus',
        'Streptococcus thermophilus#2',
        'Streptococcus uberis',
        'Streptococcus zooepidemicus',
        'Streptomyces spp',
        'Taylorella spp.',
        'Tenacibaculum spp.',
        'Treponema pallidum',
        'Trichomonas vaginalis',
        'Ureaplasma spp.',
        'Vibrio cholerae',
        'Vibrio cholerae#2',
        'Vibrio parahaemolyticus',
        'Vibrio spp.',
        'Vibrio tapetis',
        'Vibrio vulnificus',
        'Wolbachia',
        'Xylella fastidiosa',
        'Yersinia pseudotuberculosis',
        'Yersinia ruckeri',
    ]
    
    try:
        # Set the name and path of the log files
        stdout_log_path, stderr_log_path = construct_log_file_paths(work_dir)
        
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            # Determine the genus to use in the analyses
            if 'organism' in item:
                organism = item.split('=')[1].capitalize()
                continue
            # Check to see if the database should be updated
            if 'update' in item:
                update = True
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)
        
        # Ensure that the organism has been provided
        if not organism:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='ERROR: You must provide the genus to be used for the '
                      'analyses. Please create a new issue with '
                      'organism=organism e.g. organism=Escherichia coli#1 '
                      'included in the issue.',
                      status_id=4
                )
            return
        
        
        # Ensure that the organism database is present in PubMLST
        if organism not in species_list:
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes='ERROR: The organism {organism} is not present in '
                      'PubMLST. Please ensure that the organism is '
                      'spelled correctly and try again.'.format(
                          organism=organism
                      ),
                status_id=4
            )
            return
        
        # Replace any spaces with underscores in the organism name
        genus = organism.replace(' ', '_')
        
        # Set the genus database path
        genus_database_path = os.path.join(
            database_path, genus
        )
        
        # Find the most recent version of the database
        database_version_path = find_most_recent_database_version(
            genus_database_path=genus_database_path
        )
        
        # If the path doesn't exist or update is True, download the database
        if not os.path.exists(genus_database_path) or update or \
                not database_version_path:
            # Create the path
            os.makedirs(genus_database_path, exist_ok=True)
            
            # Get the current date
            current_date = datetime.datetime.now()
            
            # Set the version of the database (YYMMDD)
            database_version_path = os.path.join(
                genus_database_path,
                current_date.strftime('%y%m%d')
            )
            
            # Create the path
            os.makedirs(database_version_path, exist_ok=True)
            
            print(organism)
            
            # Treat Serratia, Proteus, and Providencia differently
            non_xml_organisms = [
                'Serratia spp.',
                'Proteus spp.',
                'Providencia spp.'
            ]
            
            if organism in non_xml_organisms:
                
                # Create a dictionary to store the URLs for the profiles and
                # alleles
                urls = {
                    'Serratia spp.': {
                        'identifier': 'serratia',
                        'alleles': [
                            'adk',
                            'fumC',
                            'gyrB',
                            'icd',
                            'mdh',
                            'recA'
                        ]
                    },
                    'Proteus spp.': {
                        'identifier': 'proteus',
                        'alleles': [
                            'atpD',
                            'dnaJ',
                            'mdh',
                            'pyrC',
                            'recA',
                            'rpoD'
                        ]
                    },
                    'Providencia spp.': {
                        'identifier': 'providencia',
                        'alleles': [
                            'fusA',
                            'gyrB',
                            'ileS',
                            'lepA',
                            'leuS'
                        ]
                    }
                }
                
                # Extract the identifier from the dictionary
                identifier = urls[organism]['identifier']
                
                # Download the profile and alleles manually
                base_url = \
                    'https://rest.pubmlst.org/db/' \
                    'pubmlst_{identifier}_seqdef/'.format(
                        identifier=identifier
                    )
                
                # Set the profile URL
                profile_url = base_url + 'schemes/1/profiles_csv'
                
                print(
                    'Downloading profile file: {url}'.format(url=profile_url)
                )
                
                # Set the path to the local profile file
                profile_path = os.path.join(
                    database_version_path, 'profiles.txt'
                )
                
                # Download the profile CSV
                response = requests.get(profile_url)
                with open(profile_path, 'wb') as file:
                    file.write(response.content)
                
                # Create a list to store the paths to the allele files
                allele_paths = []
                
                # Create the allele URLs
                for allele in urls[organism]['alleles']:
                    
                    # Set the allele URL
                    allele_url = base_url + \
                        'loci/{allele}/alleles_fasta'.format(
                            allele=allele)
                    
                    print('Downloading file: {url}'.format(url=allele_url))
                    
                    # Set the path to the local allele file
                    allele_path = os.path.join(
                        database_version_path, '{allele}.tfa'.format(
                            allele=allele
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
                    database_version_path, 'combinedtargets.fasta'
                )
                
                # Combine the allele files into a single file
                with open(combined_path, 'w') as combined_file:
                    for allele_path in allele_paths:
                        for record in SeqIO.parse(allele_path, 'fasta'):
                            SeqIO.write(record, combined_file, 'fasta')
            else:
                # Set up the download script
                activate = (
                    'source /home/ubuntu/miniconda3/bin/activate '
                    '/mnt/nas2/virtual_environments/plasmid_borne'
                )

                # Create the download command
                mlst_download_cmd = (
                    'python -m olctools.databasesetup.get_mlst '
                    '--genus "{organism}" '
                    '--repository_url http://pubmlst.org/data/dbases.xml '
                    '--path {db_path}'.format(
                        organism=organism,
                        db_path=database_version_path
                    )
                )

                # Set the name and path of the download script
                download_script = os.path.join(
                    work_dir, 'run_mlst_download.sh'
                )

                # Create the script file and make it executable
                template = create_shell_script(
                    activate=activate,
                    cmd=mlst_download_cmd,
                    script_path=download_script
                )

                # Update the Redmine issue after attempting to download
                # sequences
                notes = (
                    'Starting MLST download with command:\n {template}'.format(
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
                
                # Reformat the database to work with GeneSeekr
                process_allele_data(db_dir=database_version_path)
        
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
                notes='ERROR: Could not create the database for the organism '
                      '{organism}. Please ensure that supplied genus is '
                      'present in PubMLST and try again.'.format(
                          organism=organism
                    ),
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
        
        # Run GeneSeekr in MLST mode
        seekr_cmd = (
            'python {py} blastn -s {queries} -r {reports} -t {db_path} -M'
            .format(
                py=seekr_py,
                queries=query_dir,
                reports=os.path.join(work_dir, 'reports'),
                db_path=database_version_path,
            )
        )

        # Set the name and path of the download script
        mlst_script = os.path.join(work_dir, 'run_mlst_geneseekr.sh')
        
        # Create the script file and make it executable
        template = create_shell_script(
            activate=activate,
            cmd=seekr_cmd,
            script_path=mlst_script
        )

        # Update the Redmine issue after attempting to download sequences
        notes = (
            'Starting MLST analyses with command:\n {template}'.format(
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
                ['bash', mlst_script],
                stdout=stdout_file, stderr=stderr_file,
                check=True
            )
        
        # Zip output
        output_filename = 'mlst_output_{}'.format(issue.id)
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
            notes='MLST analysis with GeneSeekr complete!'
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
    mlst_redmine()
