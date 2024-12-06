import os
import glob
import click
import ftplib
import pickle
import shutil
import socket
import datetime
import paramiko
import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN
from nastools.nastools import retrieve_nas_files
from automator_settings import (
    AZURE_ACCOUNT_NAME,
    AZURE_ACCOUNT_KEY,
    AZURE_CONNECT_STRING,
    FTP_USERNAME,
    FTP_PASSWORD
)
SFTP_HOST = 'ftp.agr.gc.ca'
SFTP_PORT = 22

# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
# Azure
from azure.storage.blob import (
    BlobSasPermissions,
    BlobServiceClient,
    generate_blob_sas,
    RetentionPolicy
)


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def externalretrieve_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    print('External retrieving!')
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        os.makedirs(os.path.join(work_dir, str(issue.id)), exist_ok=True)
        # Parse description to figure out what SEQIDs we need to run on.
        fasta_list = list()
        fastq_list = list()
        fasta = False
        fastq = True
        sra = False
        for item in description:
            #item = item.upper().rstrip()
            item = item.rstrip()
            if 'SRA' in item:
                sra = True
                continue
            if 'fasta' in item:#added 2024-03 to allow lower case entries
                fasta = True
                fastq = False
                continue
            if 'fastq' in item:#added 2024-03 to allow lower case entries
                fastq = True
                fasta = False
                continue
            if 'FASTA' in item:
                fasta = True
                fastq = False
                continue
            if 'FASTQ' in item:
                fastq = True
                fasta = False
                continue
            if fasta:
                fasta_list.append(item)
            elif fastq:
                fastq_list.append(item)
                # Add a warning message to the issue if FASTA and SRA are selected
                if fasta and sra:
                    redmine_instance.issue.update(resource_id=issue.id,
                                                  notes='SRA option is not compatible with FASTA files')
                    # Set the files to be FASTQ
                    fastq_list = fasta_list
                    fasta_list = list()
        # Use NAStools to put FASTA and FASTQ files into our working dir.
        retrieve_nas_files(seqids=fasta_list,
                           outdir=os.path.join(work_dir, str(issue.id)),
                           filetype='fasta',
                           copyflag=True)

        retrieve_nas_files(seqids=fastq_list,
                           outdir=os.path.join(work_dir, str(issue.id)),
                           filetype='fastq',
                           copyflag=True)
        # Check that we got all the requested files.
        missing_fastas = check_fastas_present(fasta_list, os.path.join(work_dir, str(issue.id)))
        missing_fastqs = check_fastqs_present(fastq_list, os.path.join(work_dir, str(issue.id)))
        if len(missing_fastqs) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastqs))

        if len(missing_fastas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested FASTA SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))
        # Rename the files to be consistent with the desired SRA naming scheme
        # e.g. 2014-SEQ-0349_S11_L001_R1_001.fastq.gz renamed to 2014-SEQ-0349_R1.fastq.gz
        if sra:
            # Create the system call
            rename_cmd = "cd {wd} && rename 's/_S\d+_L001_R(\d)_001/_R$1/' *.gz"\
                .format(wd=os.path.join(work_dir, str(issue.id)))
            # Run the command
            os.system(rename_cmd)
        # Now make a zip folder that we'll upload to the FTP.
        shutil.make_archive(root_dir=os.path.join(work_dir, str(issue.id)),
                            format='zip',
                            base_name=os.path.join(work_dir, str(issue.id)))

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
                notes='External retrieve process complete!\n\n'
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
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes='Something went wrong! We log this automatically and will '
            'look into the problem and get back to you with a fix soon: {}.'
            .format(e)
        )


def upload_to_blob(local_file: str):
    """
    Upload files to cloud storage. Returns the URL of the uploaded file.
    """
    try:
        # Create a service client
        blob_service_client = BlobServiceClient.from_connection_string(
            AZURE_CONNECT_STRING
        )
        
        # Set the container name
        container_name = 'ftp'
        
        # Extract the file name from local file
        blob_file = os.path.basename(local_file)
        
        # Create a blob client for the current blob
        blob_client = blob_service_client.get_blob_client(
            container=container_name,
            blob=blob_file
        )
        
        # Read in the file data as binary
        with open(local_file, "rb") as data:
            # Upload the file data to the blob
            blob_client.upload_blob(data)
            # Set the storage tier
            blob_client.set_standard_blob_tier(
                standard_blob_tier='Hot')
        
        # Create the SAS token
        sas_token = generate_blob_sas(
            account_name=AZURE_ACCOUNT_NAME,
            container_name=container_name,
            blob_name=blob_file,
            account_key=AZURE_ACCOUNT_KEY,
            permission=BlobSasPermissions(read=True),
            start=datetime.datetime.now(datetime.timezone.utc) -
                datetime.timedelta(minutes=60),
            expiry=datetime.datetime.now(datetime.timezone.utc) +
                datetime.timedelta(days=8)
        )
        
        # Create a URL for the blob
        sas_url = 'https://{account_name}.blob.core.windows.net/' \
                  '{container_name}/{blob_name}?{sas_token}'.format(
                      account_name=AZURE_ACCOUNT_NAME,
                      container_name=container_name,
                      blob_name=blob_file,
                      sas_token=sas_token
                  )
        
        # Return the URL
        return sas_url
    except Exception as e:
        print(e)
        return False
    

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


def upload_to_sftp(local_file):
    """
    Uploads a file to the SFTP server with multiple attempts.
    
    :param local_file: File that you want to upload to the SFTP. Will be
    uploaded with the same name that the local file has.
    :return: True if upload ended up being successful, False if even after
    10 tries the upload didn't work.
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

def check_fastas_present(fasta_list, fasta_dir):
    missing_fastas = list()
    for seqid in fasta_list:
        if len(glob.glob(os.path.join(fasta_dir, seqid + '*.fasta'))) == 0:
            missing_fastas.append(seqid)
    return missing_fastas


def check_fastqs_present(fastq_list, fastq_dir):
    missing_fastqs = list()
    for seqid in fastq_list:
        if len(glob.glob(os.path.join(fastq_dir, seqid + '*.fastq.gz'))) < 2:
            # JAS adding an extra if statement here to allow for Nanopore (SE) reads
            if len(glob.glob(os.path.join(fastq_dir, seqid + ".fastq.gz"))) == 0:
                missing_fastqs.append(seqid)
    return missing_fastqs


if __name__ == '__main__':
    externalretrieve_redmine()
