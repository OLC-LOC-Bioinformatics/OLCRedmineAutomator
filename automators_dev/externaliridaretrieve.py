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

# I'm going to make a Frankenstein of externalretrieve.py and Irida_Retrieve.py
# so I can do IRIDA retrieves from home more easily
@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def externalretrieve_redmine(redmine_instance, issue, work_dir, description):
    print('External IRIDA retrieving!')
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

        # And finally, do some file cleanup.
        try:
            shutil.rmtree(os.path.join(work_dir, str(issue.id)))
            os.remove(os.path.join(work_dir, str(issue.id) + '.zip'))
        except:
            pass

        if upload_successful is False:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='There are connection issues with the FTP site. Unable to complete '
                                                'external retrieve process. Please try again later.')
        else:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='External IRIDA Retrieve process complete!\nNow Julie just needs to upload the folder to IRIDA\n\n'
                                                'Results are available at the following FTP address:\n'
                                                'ftp://ftp.agr.gc.ca/outgoing/cfia-ak/{}'.format(str(issue.id) + '.zip'))
    
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))


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
            s.cwd('outgoing/cfia-ak')
            f = open(local_file, 'rb')
            s.storbinary('STOR {}'.format(os.path.split(local_file)[1]), f)
            f.close()
            quit_ftp(s)
            # s.quit()
            upload_successful = True
            break
        except socket.timeout:
            s = ftplib.FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD, timeout=30)
            s.cwd('outgoing/cfia-ak')
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
