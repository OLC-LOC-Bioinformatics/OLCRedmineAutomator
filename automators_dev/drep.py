import os
import glob
import click
import pickle
import shutil
import ftplib
import sentry_sdk
from automator_settings import SENTRY_DSN
from amrsummary import before_send
from externalretrieve import upload_to_ftp
from automator_settings import FTP_USERNAME, FTP_PASSWORD
from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def drep_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # find fastas
    try:
	# Make fasta file directory
        fasta_dir = os.path.join(work_dir, 'fastas')
        if not os.path.isdir(fasta_dir):
            os.makedirs(fasta_dir)

        # Description should just be a list of SEQIDs. Get the fasta files associated with them extracted
	# Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
        retrieve_nas_files(seqids=description,
                           outdir=fasta_dir,
                           filetype='fasta',
                           copyflag=False)
        #missing_fastas = verify_fasta_files_present(seqids, work_dir)
        
	# Verify that specified fasta files are actually there, warn user if they aren't.
        missing_fastas = verify_fasta_files_present(seqid_list=description,
                                                    fasta_dir=os.path.join(work_dir, 'fastas'))
        if len(missing_fastas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))
	# Make output dir
        output_dir = os.path.join(work_dir, 'dRep_results')
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

	# These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/drep'
        dRep = '/mnt/nas2/virtual_environments/drep/bin/dRep'
        # Run dRep compare with the necessary arguments
        cmd = '{dRep} compare {outpath} -g {seqpath}/*.fasta'.format(dRep=dRep,
                                                                  seqpath=fasta_dir,
                                                                  outpath=output_dir)
        os.system(cmd)

	# Zip output folder and upload to the FTP
        shutil.make_archive(root_dir=output_dir,
                            format='zip',
                            base_name=os.path.join(work_dir, str(issue.id)))

        upload_successful = upload_to_ftp(local_file=os.path.join(work_dir, str(issue.id) + '.zip'))

        if upload_successful:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='dRep genome comparison complete!\n\n'
                                                'Results are available at the following FTP address:\n'
                                                'ftp://ftp.agr.gc.ca/outgoing/cfia-ac/{}'.format(str(issue.id) + '.zip'))
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity issues. '
                                                'Please try again later.')

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')

        #Upload to redmine
        #output_filename = 'dRep_output'
        #zip_filepath = zip_folder(results_path=os.path.join(work_dir, 'dRep_results'),
        #                          output_dir=work_dir,
        #                          output_filename=output_filename)
        #zip_filepath += '.zip'

	# Prepare upload
        #output_list = [
        #    {
        #        'filename': os.path.basename(zip_filepath),
        #        'path': zip_filepath
        #    }
        #]
	
	# Create a list of all the folders - will be used to clean up the working directory
        #folders = glob.glob(os.path.join(work_dir, 'seqs/'))
        # Remove the seqs folder
        #for folder in folders:
        #    if os.path.isdir(folder):
        #        shutil.rmtree(folder)
        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes='{at} analysis with dRep complete!'
                                      .format(at=argument_dict['analysis'].lower()))
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')

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
    drep_redmine()
