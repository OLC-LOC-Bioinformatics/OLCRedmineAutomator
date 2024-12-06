import os
import glob
import click
import pickle
import shutil
from nastools.nastools import retrieve_nas_files
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
import ftplib
from ftplib import FTP
from automator_settings import FTP_USERNAME, FTP_PASSWORD, FTP_FOLDER
import traceback
import pandas as pd

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def scoary_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        #analyses supported by the automator
        correctiontypes = [
            'I', 'B', 'BH', 'PW', 'EPW', 'P'
        ]

        # Variable to hold supplied arguments
        argument_dict = {
            'files_folder': str(),
            'p_value': 0.05,
            'correction': 'I'
        }

        # Parse description to figure out what SEQIDs we need to run on.
#        seqids = list()
        for item in description:
            item = item.rstrip()
            if 'files_folder' in item:
                argument_dict['files_folder'] = item.split('=')[1]
            if 'p_value' in item:
                argument_dict['p_value'] = item.split('=')[1]
                continue
            if 'correction' in item:
                argument_dict['correction'] = item.split('=')[1].upper
                continue
            return
#            # Otherwise the item should be a SEQID
#            seqids.append(item)


        # Sanity check for arguments - no folder on ftp
#        if not argument_dict['files_folder']:
#            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
#                                          notes='ERROR: Scoary requires a folder containing the traits.csv and '
#                                                'gene_presence_absence.csv files. The automator could not find '
#                                                'a folder on the ftp. '
#                                                'please upload your folder to '
#                                                'ftp://ftp.agr.gc.ca/incoming/cfia-ac '
#                                                'and include folder={foldername} as the first line of the '
#                                                'redmine request description. ')
#            return
        #incorrect correction type provided
        if argument_dict['correction'] not in correctiontypes:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied correction/filtration measure type {at} current not in the supported '
                                                'list of correction types: {ats}'.format(at=argument_dict['correction'],
                                                                                       ats=', '.join(correctiontypes)),
                                          status_id=4)
            return        

        # Create output folder
#        output_folder = os.path.join(work_dir, 'scoary_output')
#        os.makedirs(output_folder)

        #download the files

        #verify that the file_folder is there
        file_folder = argument_dict['file_folder']
#        local_folder = file_folder
        
        #do verification checks that files are present
        validation = verify_all_the_things(file_folder=argument_dict['file_folder'],
                                           issue=issue,
                                           work_dir=work_dir,
                                           redmine_instance=redmine_instance)
        if validation is False:
            return

        download_info_sheets(file_folder, work_dir)
        redmine_instance.issue.update(resource_id=issue.id, status_id=2,
                                      notes='Beginning download '
                                            'of files.')

        #create the local folder that we'll need.
        local_folder = os.path.join(work_dir, file_folder)

        if not os.path.isdir(local_folder):
            os.makedirs(local_folder)

        # download the folder from ftp recursively
        download_successful = download_dir(file_folder, local_folder)

        if download_successful is False:
            redmine_instance.issue.update(resource_id=issue.id,
                                          assigned_to_id=429,
                                          subject='Scoary: {}'.format(description[0]),
                                          notes='Download of files from FTP was not successful.')
            return
 
        # Once the folder has been downloaded, run scoary
        # Run scoary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/roary'
        scoary = '/mnt/nas2/virtual_environments/roary/bin/scoary'

        scoary_cmd = 'cd {output_dir} && scoary -g {gpa} -t {traits} -p {p_value} -c {correction}'.format(output_dir=local_folder, 
                                                                                                          gpa='gene_presence_absence.csv', 
                                                                                                          traits='traits.csv',
                                                                                                          p_value=argument_dict['p_value'],
                                                                                                          correction=argument_dict['correction'])
        redmine_instance.issue.update(resource_id=issue.id, notes='Going to run scoary with this command:{scoary}'.format(scoary=scoary_cmd))
        template = "#!/bin/bash\n{activate} && {scoary}".format(activate=activate, scoary=scoary_cmd)
        scoary_script = os.path.join(work_dir, 'run_scoary.sh')
        with open(scoary_script, 'w+') as file:
            file.write(template)
        make_executable(scoary_script)
        os.system(scoary_script)

        # Get scoary results uploaded
        # zip the scoary folder
        output_filename = 'scoary_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=local_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'
        # Prepare upload
        # This file can get too big to upload to Redmine, so we should put it on the FTP.
#        output_list = [
#            {
#                'filename': os.path.basename(zip_filepath),
#                'path': zip_filepath
#            }
#        ]
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
                status_id=4,
                notes='scoary analysis complete!\n\n'
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
        # Clean up files
#        shutil.rmtree(output_folder)
        os.remove(zip_filepath)
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))

def verify_folder_exists(file_folder):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd('incoming/cfia-ac')
    folders_present = ftp.nlst()
    if file_folder in folders_present:
        folder_exists = True
    else:
        folder_exists = False
    # ftp.quit()
    quit_ftp(ftp)
    return folder_exists

def download_info_sheets(file_folder, local_folder):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd(os.path.join('incoming/cfia-ac', file_folder))
    info_sheets = ['gene_presence_absence.csv', 'traits.csv']
    for sheet in info_sheets:
        try:
            f = open(os.path.join(local_folder, sheet), 'wb')
            ftp.retrbinary('RETR ' + sheet, f.write)
            f.close()
        except:
            pass
    # ftp.quit()
    quit_ftp(ftp)

def verify_all_the_things(file_folder, redmine_instance, issue, work_dir):
    validation = True
    # Verify that the folder specified does in fact exist. If it doesn't, give up.
    if verify_folder_exists(file_folder) is False:
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='ERROR: Could not find the folder ({}) specified in this issue on '
                                            'the FTP. Please ensure that it is uploaded correctly, create a new issue,'
                                            ' and try again.'.format(file_folder))
        validation = False
        return validation   # Can't check anything else if the folder doesn't exist, so stop here.
    # Next up, validate that gene_presence_absence.csv and traits.csv are present.
    missing_files = validate_files(file_folder)
    if len(missing_files) > 0:
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='ERROR: The following files were missing from the FTP'
                                            ' folder: {}\nPlease reupload to '
                                            'the FTP, including these files ('
                                            'spelling must be identical!) and'
                                            ' create a new issue.'.format(str(missing_files)))
        validation = False
    # Now, download the files to a temporary folder.
    if not os.path.isdir(os.path.join(work_dir, file_folder)):
        os.makedirs(os.path.join(work_dir, file_folder))
    download_info_sheets(file_folder, os.path.join(work_dir, file_folder))



def download_dir(ftp_dir, local_dir):
    all_downloads_successful = True
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd(os.path.join('incoming/cfia-ac', ftp_dir))
    present_in_folder = ftp.nlst()
    for item in present_in_folder:
        if check_if_file(item, ftp_dir):
            ftp_file = os.path.join(ftp_dir, item)
            download_successful = download_ftp_file(ftp_file=ftp_file, local_dir=local_dir)
            if download_successful is False:
                all_downloads_successful = False
        else:
            if not os.path.isdir(os.path.join(local_dir, item)):
                os.makedirs(os.path.join(local_dir, item))
            download_dir(os.path.join(ftp_dir, item), os.path.join(local_dir, item))
    quit_ftp(ftp)
    # ftp.quit()
    return all_downloads_successful

def check_if_file(file_name, ftp_dir):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd(os.path.join('incoming/cfia-ac', ftp_dir))
    is_file = True
    try:
        ftp.size(file_name)
    except:
        is_file = False
    # ftp.quit()
    quit_ftp(ftp)
    return is_file

def delete_ftp_dir(ftp_dir):
    """
    Cleaning up the FTP after things have finished is a good idea. This allows recursive deletion of an FTP dir.
    :param ftp_dir: Name of directory within incoming/cfia-ac you want deleted.
    """
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd(os.path.join('incoming/cfia-ac', ftp_dir))
    present_in_folder = ftp.nlst()
    for item in present_in_folder:
        if check_if_file(item, ftp_dir):
            ftp.delete(item)
        else:
            delete_ftp_dir(os.path.join(ftp_dir, item))
    ftp.cwd('..')
    ftp.rmd(os.path.split(ftp_dir)[-1])

def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)

def quit_ftp(ftp_object):
    """
    Apparently our connection to the FTP is so good that sometimes we can manage to get a timeout when
    trying to quit the FTP. This function will try to call quit(), will give up and not error out when a timeout happens
    :param ftp_object: An instantiated FTP object from python's ftplib
    """
    try:
        ftp_object.quit()
    except ftplib.error_temp:
        print('Timeout occurred when trying to close connection to the FTP. Ignoring the problem!')

def validate_files(file_name):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    missing_files = list()
    ftp.cwd(os.path.join('incoming/cfia-ac', file_name))
    files_present = ftp.nlst()
    if 'gene_presence_absence.csv' not in files_present:
        missing_files.append('gene_presence_absence.csv')
    if 'traits.csv' not in files_present:
        missing_files.append('traits.csv')
    # ftp.quit()
    quit_ftp(ftp)
    return missing_files

def verify_fasta_files_present(seqid_list, fasta_dir):
    missing_fastas = list()
    for seqid in seqid_list:
        if len(glob.glob(os.path.join(fasta_dir, seqid + '*.fasta'))) == 0:
            missing_fastas.append(seqid)
    return missing_fastas


def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    scoary_redmine()
