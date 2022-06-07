import os
import glob
import click
import pickle
import shutil
import sentry_sdk
from automator_settings import COWBAT_DATABASES
from nastools.nastools import retrieve_nas_files


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def geneseekr_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Current list of analysis types that the GeneSeekr can perform
    analyses = [
        'custom', 'gdcs', 'genesippr', 'mlst', 'resfinder', 'rmlst', 'serosippr', 'sixteens', 'virulence'
    ]
    # Current BLAST analyses supported
    blasts = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
    # Variable to hold supplied arguments
    argument_dict = {
        'analysis': str(),
        'align': False,
        'blast': 'blastn',
        'cutoff': 70,
        'evalue': '1E-5',
        'unique': False,
        'organism': str(),
        'fasta': False,
    }
    # Dictionary of analysis types to argument flags to pass to the script
    argument_flags = {
        'custom': '',
        'gdcs': '-Q',
        'genesippr': '-G',
        'mlst': '-M',
        'resfinder': '-A',
        'rmlst': '-R',
        'serosippr': '-S',
        'sixteens': '-X',
        'virulence': '-V'
    }
    # Set the database path for the analyses
    dbpath = COWBAT_DATABASES
    database_path = {
        'custom': os.path.join(work_dir, 'targets'),
        'gdcs': os.path.join(dbpath, 'GDCS'),
        'genesippr': os.path.join(dbpath, 'genesippr'),
        'mlst': os.path.join(dbpath, 'MLST'),
        'resfinder': os.path.join(dbpath, 'resfinder'),
        'rmlst': os.path.join(dbpath, 'rMLST'),
        'serosippr': os.path.join(dbpath, 'serosippr', 'Escherichia'),
        'sixteens': os.path.join(dbpath, 'sixteens_full'),
        'virulence': os.path.join(dbpath, 'virulence'),
    }
    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'ALIGN' in item:
                argument_dict['align'] = True
                continue
            if 'BLAST' in item:
                argument_dict['blast'] = item.split('=')[1].lower()
                continue
            if 'CUTOFF' in item:
                argument_dict['cutoff'] = int(item.split('=')[1].lower())
                continue
            if 'EVALUE' in item:
                argument_dict['evalue'] = item.split('=')[1].lower()
                continue
            if 'UNIQUE' in item:
                argument_dict['unique'] = item.split('=')[1].lower()
                continue
            if 'ORGANISM' in item:
                argument_dict['organism'] = item.split('=')[1].capitalize()
                continue
            if 'ANALYSIS' in item:
                argument_dict['analysis'] = item.split('=')[1].lower()
                continue
            if 'FASTA' in item:
                argument_dict['fasta'] = True
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)
        # Use a cutoff value of 99% identity for allele-based analyses
        if 'mlst' in argument_dict['analysis']:
            argument_dict['cutoff'] = 99
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
                attachment.download(savepath=target_dir, filename='targets.tfa')
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Analysis type custom requires an attached FASTA file of '
                                                    'targets. The automator could not find any attached files. '
                                                    'Please create a new issue with the FASTA file attached and try '
                                                    'again.',
                                              status_id=4)
                return
        # Ensure that the organism has been provided for organism-specific analyses
        if argument_dict['analysis'] == 'gdcs' or argument_dict['analysis'] == 'mlst':
            if not argument_dict['organism']:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Analysis type {at} requires the genus to be used for the '
                                                    'analyses. Please create a new issue with organism=ORGANISM '
                                                    'included in the issue.'.format(at=argument_dict['analysis']),
                                              status_id=4)
                return
            else:
                database_path[argument_dict['analysis']] = os.path.join(database_path[argument_dict['analysis']],
                                                                        argument_dict['organism'])
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
                                          notes='WARNING: supplied analysis type {at} current not in the supported '
                                                'list of analyses: {ats}'.format(at=argument_dict['analysis'],
                                                                                 ats=', '.join(analyses)),
                                          status_id=4)
            return
        # Ensure that the requested BLAST analysis is valid
        if argument_dict['blast'] not in blasts:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: requested BLAST analysis, {bt}, is not one of the currently '
                                                'supported analyses: {blasts}'.format(bt=argument_dict['blast'],
                                                                                      blasts=', '.join(blasts)),
                                          status_id=4)
            return
        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No SEQIDs provided!',
                                          status_id=4)
            return
        # Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
        retrieve_nas_files(seqids=seqids,
                           outdir=work_dir,
                           filetype='fasta',
                           copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, work_dir)
        # Update the Redmine issue if one or more of the requested SEQIDs could not be located
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))

        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/geneseekr'
        seekr_py = '/mnt/nas2/virtual_environments/geneseekr/bin/GeneSeekr'
        # Run sipprverse with the necessary arguments
        seekr_cmd = 'python {seekr_py} {blast} -s {seqpath} -r {outpath} -t {dbpath} -c {cutoff} -e {evalue} {atf}'\
            .format(seekr_py=seekr_py,
                    blast=argument_dict['blast'],
                    seqpath=work_dir,
                    outpath=os.path.join(work_dir, 'reports'),
                    dbpath=database_path[argument_dict['analysis']],
                    cutoff=argument_dict['cutoff'],
                    evalue=argument_dict['evalue'],
                    atf=argument_flags[argument_dict['analysis']])
        # Append the align and/or the unique flags are required
        seekr_cmd += ' -a' if argument_dict['align'] else ''
        seekr_cmd += ' -u' if argument_dict['unique'] else ''
        seekr_cmd += ' -f' if argument_dict['fasta'] else ''
        # Update the issue with the GeneSeekr command
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='GeneSeekr command:\n {cmd}'.format(cmd=seekr_cmd))
        # Create another shell script to execute within the PlasmidExtractor conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, seekr_cmd)
        geneseekr_script = os.path.join(work_dir, 'run_geneseekr.sh')
        with open(geneseekr_script, 'w+') as file:
            file.write(template)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(geneseekr_script)
        # Run shell script
        os.system(geneseekr_script)

        # Zip output
        output_filename = 'geneseekr_output'
        zip_filepath = zip_folder(results_path=os.path.join(work_dir, 'reports'),
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'
        # Prepare upload
        output_list = [
            {
                'filename': os.path.basename(zip_filepath),
                'path': zip_filepath
            }
        ]

        # Create a list of all the folders - will be used to clean up the working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)
        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes='{at} analysis with GeneSeekr complete!'
                                      .format(at=argument_dict['analysis'].lower()))
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')


# def run_mash(seqids, output_dir):
#     """
#     Use MASH to determine the genus of strains when the requested analysis has a genus-specific database
#     :return: dictionary of MASH-calculated genera
#     """
#     # Dictionary to store the MASH results
#     genus_dict = dict()
#     # Run mash screen on each of the assemblies
#     for seqid in seqids:
#         screen_file = os.path.join(output_dir, '{seqid}_screen.tab'.format(seqid=seqid))
#         mash.screen('/mnt/nas2/databases/confindr/databases/refseq.msh',
#                     item,
#                     threads=8,
#                     w='',
#                     i='0.95',
#                     output_file=screen_file,
#                     returncmd=True)
#         screen_output = mash.read_mash_screen(screen_file)
#         # Determine the genus from the screen output file
#         for screen in screen_output:
#             # Extract the genus from the mash results
#             mash_organism = screen.query_id.split('/')[-3]
#             # Populate the dictionary with the seqid, and the calculated genus
#             genus_dict[seqid] = mash_organism
#     return genus_dict


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
    geneseekr_redmine()
