import os
import glob
import click
import pickle
import shutil
import ftplib
import sentry_sdk
from automator_settings import SENTRY_DSN
from amrsummary import before_send
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def mummer2circos_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    #commands supported by automator
    alignmethods = [
        'nucmer', 'promer', 'megablast'
    ]


    #variable to hold supplied arguments
    argument_dict = {
        'reference': str(),
        'analysis': '',
        'alignmethod': 'nucmer',
        'gaps': False,
    }

    argument_flags = {
        'custom': '',
        'gaps': '-g'


    #parse description to figure out analysis type, and find fastas
    try:
        #parse description
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'REFERENCE' in item:
               argument_dict['reference'] = item.split('=')[1]
               continue
            if 'ANALYSIS' in item:
                argument_dict['analysis'] = item.split('=')[1].lower()
                continue
            if 'ALIGNMETHOD' in item:
                argument_dict['alignmethod'] = item.split('=')[1].lower()
                continue
            if 'GAPS' in item:
                argument_dict['gaps'] = True
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that the reference is provided
        if not argument_dict['reference']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No reference sequence provided. '
                                                'Please ensure that the first line of the issue contains one',
                                          status_id=4)
            return

        # Ensure that the command type is acceptable
        if argument_dict['alignmethod'] not in alignmethods:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied alignment method {at} is not in the supported '
                                                'list of alignment methods: {ats}'.format(at=argument_dict['alignmethod'],
                                                                                   ats=', '.join(alignmethods)),
                                          status_id=4)
            return

        #if just a custom analysis, create the fastas folder
        if argument_dict['analysis'] == 'custom':
            # Make fasta file directory
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)

        #create the folder to hold the fasta files
        fasta_dir = os.path.join(work_dir, 'fastas')
        if not os.path.isdir(fasta_dir):
            os.makedirs(fasta_dir)

        #now retrieve the reference file
        ref_dir = os.path.join(work_dir, 'reference')
        if not os.path.isdir(ref_dir):
            os.makedirs(ref_dir)
        retrieve_nas_files(seqids=argument_dict['reference'],
                           outdir=ref_dir,
                           filetype='fasta',
                           copyflag=False)
        missing_fastas = verify_fasta_files_present(argument_dict['reference'], ref_dir)
        

        # Rest of Description should just be a list of SEQIDs. Get the fasta files associated with them extracted
	# Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
        retrieve_nas_files(seqids=seqids,
                           outdir=fasta_dir,
                           filetype='fasta',
                           copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, fasta_dir)
        
	# Verify that specified fasta files are actually there, warn user if they aren't.
        missing_fastas = verify_fasta_files_present(seqid_list=seqids,
                                                    fasta_dir=os.path.join(work_dir, 'fastas'))
        if len(missing_fastas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))
	# Make output dir
        output_dir = os.path.join(work_dir, 'mummer2circos_results')
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

	
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/mummer2circos'
        # Run mummer2circos with the necessary arguments
        cmd = 'mummer2circos -a {align} -l -r ref_dir/{ref} -q fasta_dir/*.fasta -o {ID}_circos '.format(align=argument_dict['alignmethod'],
                                                                                                         ref=argument_dict['reference'],
                                                                                                         ID=issue.id)
        cmd += ' -g' if argument_dict['gaps'] else ''


        #create another shell script to execute within the dRep conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, cmd)
        run_script = os.path.join(work_dir, 'run_mummer2circos.sh')
        with open(run_script, 'w+') as file:
            file.write(template)

        # Modify the permissions of the script to allow it to be run on the node
        make_executable(run_script)
        # Run shell script
        os.system(run_script)

	# Zip output folder and upload to the FTP, as sometimes files are >10Mb
        shutil.make_archive(root_dir=output_dir,
                            format='zip',
                            base_name=os.path.join(work_dir, str(issue.id)))

        # Upload the zip file to Dropbox
        download_link = upload_to_dropbox(
            access_token=DROPBOX_ACCESS_TOKEN,
            refresh_token=DROPBOX_REFRESH_TOKEN,
            app_key=DROPBOX_APP_KEY,
            app_secret=DROPBOX_APP_SECRET,
            local_file_path=os.path.join(work_dir, str(issue.id) + '.zip')
        )

        if download_link:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='dRep analysis complete!\n\n'
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

        # Create a list of all the folders - will be used to clean up the working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)

        os.remove(os.path.join(work_dir, str(issue.id) + '.zip')) #delete zip file
        # Wrap up issue
#        redmine_instance.issue.update(resource_id=issue.id,
#                                      uploads=output_list,
#                                      status_id=4,
#                                      notes='{at} analysis with dRep complete!'.format(at=argument_dict['analysis'].lower()))
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))

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
    mummer2circos_redmine()
