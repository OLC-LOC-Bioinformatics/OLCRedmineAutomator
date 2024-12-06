import os
import glob
import click
import ftplib
import pickle
import shutil
import socket
import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN
from nastools.nastools import retrieve_nas_files
#from nastoolsgfaretrieve import retrieve_nas_files
from automator_settings import FTP_USERNAME, FTP_PASSWORD


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def gfaretrieve_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    print('External retrieving!')
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        #os.makedirs(os.path.join(work_dir, str(issue.id)))
        # Parse description to figure out what SEQIDs we need to run on.
        gfa_list = list()
        for item in description:
            #item = item.upper().rstrip()
            item = item.rstrip()
            #if 'SRA' in item:
            #    sra = True
            #    continue
            gfa_list.append(item)
        
        #create an output directory
        gfa_files = os.path.join(work_dir, 'gfa_files')
        os.makedirs(gfa_files, exist_ok=True)

        
        #write our list to a file so that nastools will work
        listfile = os.path.join(work_dir, 'gfa_list.txt')
        with open(listfile, 'w') as outfile:
            for seqid in gfa_list:
                outfile.write(seqid)
                outfile.write("\n")

        # Use NAStools based gfa retrieve script (modified for use by AC 24-05-02) to put gfa files into our working dir.
        #retrieve_nas_files(seqids=gfa_list, outdir=gfa_files, filetype=gfa, copyflag=True)
        gfaretrpy = '/mnt/nas2/redmine/applications/OLCRedmineAutomator/automators/nastoolsgfaretrieve.py'
        retrieve_cmd = 'python {gfaretrpy} --file {seqids} --type gfa --outdir {out} --copy'\
            .format(gfaretrpy=gfaretrpy,seqids=listfile, out=gfa_files)
        retr_script = os.path.join(work_dir, 'retrieve_gfas.sh')
        with open(retr_script, 'w+') as file:
            file.write(retrieve_cmd)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(retr_script)
        # Run shell script
        os.system(retr_script)

        # Check that we got all the requested files.
        #missing_fastas = check_fastas_present(fasta_list, os.path.join(work_dir, str(issue.id)))
        missing_gfas = check_gfas_present(gfa_list, gfa_files)

        if len(missing_gfas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested gfa files on'
                                                ' the OLC NAS: {}. It is possible they do not exist. These files '
                                                ' were not created in older versions of the Hybrid Assembly Pipeline'.format(missing_gfas))

	# Zip output folder and upload to the FTP, as sometimes files are >10Mb
        zipfilename = 'gfa_files_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=gfa_files,
                                  output_dir=work_dir,
                                  output_filename=zipfilename)
        zip_filepath += '.zip'

        sas_url = upload_to_ftp(local_file=zip_filepath)


        # And finally, do some file cleanup.
        try:
            shutil.rmtree(gfa_files)
            os.remove(os.path.join(work_dir, zipfilename + '.zip')) #delete zip file
        except:
            pass

        if sas_url is False:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='There are connection issues. Unable to complete '
                      'external retrieve process. Please try again later.'
                )
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='External Retrieve process complete!\n\n'
                                                'Results are available at the following FTP address:\n'
						'ftp://ftp.agr.gc.ca/outgoing/cfia-ac/{}'
                                          .format(str(issue.id) + '.zip'))
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon: {}'.format(e))


def upload_to_ftp(local_file):
    """
    Since our FTP site has been misbehaving, we now get to have a special FTP uploader that tries to
    upload multiple times (up to 10).
    :param local_file: File that you want to upload to the FTP. Will be uploaded with the same name that
    the local file has.
    :return: True if upload ended up being successful, False if even after 10 tries the upload didn't work.
    """
    num_upload_attempts = 0
    upload_successful = False
    while num_upload_attempts < 10:
        # Try uploading - if timeout, check if the upload managed to complete but hang at the end, which happens
        # sometimes. If it did complete, we're good to go. Otherwise, try again.
        try:
            s = ftplib.FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD, timeout=30)
            s.cwd('outgoing/cfia-ac')
            f = open(local_file, 'rb')
            s.storbinary('STOR {}'.format(os.path.split(local_file)[1]), f)
            f.close()
            s.quit()
            upload_successful = True
            break
        except socket.timeout:
            s = ftplib.FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD, timeout=30)
            s.cwd('outgoing/cfia-ac')
            uploaded_file_size = s.size(os.path.split(local_file)[1])
            s.quit()
            if uploaded_file_size == os.path.getsize(local_file):
                upload_successful = True
                break
            num_upload_attempts += 1
    return upload_successful

def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)


def check_gfas_present(gfa_list, fasta_dir):
    missing_gfas = list()
    for seqid in gfa_list:
        if len(glob.glob(os.path.join(fasta_dir, seqid + '*.gfa'))) == 0:
            missing_gfas.append(seqid)
    return missing_gfas

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
    gfaretrieve_redmine()
