import os
import glob
import click
import ftplib
import pickle
import shutil
import socket

from iridaextractfiles import MassExtractor
from iridasequencefile import SequenceInfo

from wgsassembly import quit_ftp
from automator_settings import FTP_USERNAME, FTP_PASSWORD
SFTP_HOST = 'ftp.agr.gc.ca'
SFTP_PORT = 22

from automator_settings import SENTRY_DSN
import sentry_sdk
from amrsummary import before_send

# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)

permitted_users = [106, 429, 225, 529, 745] # Adam Julie Cathy Ray Brenna

# I'm going to make a Frankenstein of externalretrieve.py and Irida_Retrieve.py
# so I can do IRIDA retrieves from home more easily
@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def externalretrieve_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
 
    try:
        os.makedirs(os.path.join(work_dir, str(issue.id)))
        # Parse description
        sequences_info = list()
        for input_line in description:
            if input_line is not '':
                sequences_info.append(SequenceInfo(input_line))
            sequences_info = get_validated_seqids(sequences_info)
        
        # Use class from original IRIDA retrieve to extract fastq files
        extractor = MassExtractor(nas_mnt='/mnt/nas')
        missing_files, low_quality = extractor.move_files(sequences_info, os.path.join(work_dir, str(issue.id)))
        if len(low_quality) > 0:
            redmine_instance.issue.update(resource_id=issue.id, notes='WARNING: The following sample IDs had average read qualities below Q30 and were not retrieved: {}'.format(low_quality))
        redmine_instance.issue.update(resource_id=issue.id, notes='Files have been extracted and formatted for IRIDA! Files will now be zipped and uploaded to FTP.')

        # Now make a zip folder that we'll upload to the FTP.
        shutil.make_archive(root_dir=os.path.join(work_dir, str(issue.id)),
                            format='zip',
                            base_name=os.path.join(work_dir, str(issue.id)))

        # Now need to login to the FTP to upload the zipped folder.
        # Lots of FTP issues lately - in the event that upload does not work, a timeout will occur.
        # Allow for up to 10 attempts at uploading. If upload has completed and we stall at the end, allow.
        upload_successful = upload_to_ftp(local_file=os.path.join(work_dir, str(issue.id) + '.zip'))

                # Set the path of the zip file
        zip_filepath = os.path.join(work_dir, str(issue.id) + '.zip')
        
        # Upload the zip file to Dropbox
        download_link = upload_to_dropbox(
            access_token=DROPBOX_ACCESS_TOKEN,
            refresh_token=DROPBOX_REFRESH_TOKEN,
            app_key=DROPBOX_APP_KEY,
            app_secret=DROPBOX_APP_SECRET,
            local_file_path=zip_filepath
        )

        if download_link:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=2,
                notes='External IRIDA Retrieve process complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=download_link)
            )
        else:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=2,
                notes='Upload of files was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )


        # And finally, do some file cleanup.
        try:
            shutil.rmtree(os.path.join(work_dir, str(issue.id)))
            os.remove(os.path.join(work_dir, str(issue.id) + '.zip'))
        except:
            pass

    
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the'
                                            'problem and get back to you with a fix soon. Here is the traceback: {e}'.format(e=e))


def upload_to_sftp(local_file):
    """
    Uploads a file to the SFTP server with multiple attempts.
    
    :param local_file: File that you want to upload to the SFTP. Will be
    uploaded with the same name that the local file has.
    :return: True if upload ended up being successful, False if even after 10
    tries the upload didn't work.
    """
    num_upload_attempts = 0
    upload_successful = False
    while num_upload_attempts < 10:
        try:
            # Establish an SFTP connection
            transport = paramiko.Transport((SFTP_HOST, SFTP_PORT))
            transport.connect(username=FTP_USERNAME, password=FTP_PASSWORD)
            sftp = paramiko.SFTPClient.from_transport(transport)
            
            # Change to the target directory
            try:
                sftp.chdir('outgoing/cfia-ak')
            except IOError as exc:
                print("Failed to change directory: {exc}".format(exc=exc))
                raise
            
            # Upload the file
            try:
                sftp.put(local_file, os.path.split(local_file)[1])
            except IOError as exc:
                print("Failed to upload file: {exc}".format(exc=exc))
                raise
            
            # Verify the file size
            try:
                remote_file_size = sftp.stat(
                    os.path.split(local_file)[1]
                ).st_size
                local_file_size = os.path.getsize(local_file)
            except IOError as exc:
                print("Failed to stat file: {exc}".format(exc=exc))
                raise
            
            if remote_file_size == local_file_size:
                upload_successful = True
                sftp.close()
                transport.close()
                break
        except (socket.timeout, paramiko.SSHException, IOError) as exc:
            print(
                "Upload attempt {num_upload_attempts} failed: {exc}".format(
                    num_upload_attempts=num_upload_attempts + 1,
                    exc=exc
                )
            )
            num_upload_attempts += 1
        finally:
            try:
                sftp.close()
                transport.close()
            except:
                pass
    return upload_successful

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
            quit_ftp(s)
            # s.quit()
            upload_successful = True
            break
        except socket.timeout:
            s = ftplib.FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD, timeout=30)
            s.cwd('outgoing/cfia-ac')
            uploaded_file_size = s.size(os.path.split(local_file)[1])
            quit_ftp(s)
            # s.quit()
            if uploaded_file_size == os.path.getsize(local_file):
                upload_successful = True
                break
            num_upload_attempts += 1
    return upload_successful


def check_fastqs_present(fastq_list, fastq_dir):
    missing_fastqs = list()
    for seqid in fastq_list:
        if len(glob.glob(os.path.join(fastq_dir, seqid + '*.fastq.gz'))) < 2:
            # JAS adding an extra if statement here to allow for Nanopore (SE) reads
            if len(glob.glob(os.path.join(fastq_dir, seqid + ".fastq.gz"))) == 0:
                missing_fastqs.append(seqid)
    return missing_fastqs


def get_validated_seqids(sequences_list):
    """
    A inputted list is checked for Seq-ID format, each of the Elements that are validated are returned to the user
    sequences_list: list of Seq-IDs to be validated
    """

    validated_sequence_list = list()
    regex = r'^(2\d{3}-\w{2,10}-\d{3,4})$'
    import re
    for sequence in sequences_list:
        # if re.match(regex, sequence.sample_name):
        validated_sequence_list.append(sequence)
        # else:
        #     raise ValueError("Invalid seq-id \"%s\"" % sequence.sample_name)

    if len(validated_sequence_list) < 1:
        raise ValueError("Invalid format for redmine request. Couldn't find any fastas or fastqs to extract")

    return validated_sequence_list


if __name__ == '__main__':
    externalretrieve_redmine()
