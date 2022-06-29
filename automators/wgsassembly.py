import os
import re
import csv
import glob
import click
import pickle
import shutil
import socket
import fnmatch
import xml.etree.ElementTree as et
from externalretrieve import upload_to_ftp
import ftplib
from ftplib import FTP
from automator_settings import FTP_USERNAME, FTP_PASSWORD, COWBAT_IMAGE, COWBAT_DATABASES
import traceback
import pandas as pd

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def wgsassembly_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        # Add Cathy as a watcher so that we can make sure things get done. Also add Ashley in case people
        # forget to assign the issue to me.
        issue.watcher.add(225)  # This is Cathy
        issue.watcher.add(222)  # This is Ashley
        issue.watcher.add(904)  # This is Liam
        #issue.watcher.add(429)  # This is Julie.
        # instead of folder on NAS.
        # Verify that sequence folder in description is named correctly.
        sequence_folder = description[0]

        # If the sequence folder looks like an absolute path on our NAS (implement a tougher check here),
        # we assume the user knows what they're doing and don't bother validating anything
        if len(sequence_folder.split('/')) > 1:
            redmine_instance.issue.update(resource_id=issue.id, status_id=2,
                                          notes='Attempting to use files already on NAS. Only use this option if you'
                                                ' really know what you\'re doing!')
            # os.path.split does not work as I thought it did. Apparently if a trailing slash is present, last element
            # it returns is '', not the final folder in the path. This fixes that.
            if sequence_folder.endswith('/'):
                sequence_folder = sequence_folder[:-1]
            local_folder = sequence_folder
            # samplesheet_seqids = get_seqids_from_samplesheet(os.path.join(sequence_folder, 'SampleSheet.csv'))
            # lab_id = samplesheet_seqids[0].split('-')[1]
            sequence_folder = os.path.split(local_folder)[1]

        # Otherwise, do all verification checks on the FTP and download files.
        else:
            validation = verify_all_the_things(sequence_folder=sequence_folder,
                                               issue=issue,
                                               work_dir=work_dir,
                                               redmine_instance=redmine_instance)

            # All checks that needed to be done should now be done. If any of them returned something bad,
            # we stop and boot the user. Otherwise, go ahead with downloading files.
            if validation is False:
                return

            download_info_sheets(sequence_folder, work_dir)
            redmine_instance.issue.update(resource_id=issue.id, status_id=2,
                                          notes='All validation checks passed - beginning download '
                                                'and assembly of sequence files.')

            # Create the local folder that we'll need.
            local_folder = os.path.join('/mnt/nas2/raw_sequence_data/miseq', sequence_folder)

            if not os.path.isdir(local_folder):
                os.makedirs(local_folder)

            # Download the folder, recursively!
            download_successful = download_dir(sequence_folder, local_folder)

            if download_successful is False:
                redmine_instance.issue.update(resource_id=issue.id,
                                              assigned_to_id=429,
                                              subject='WGS Assembly: {}'.format(description[0]),
                                              notes='Download of files from FTP was not successful.')
                return

        # Once the folder has been downloaded, copy it to the hdfs and start assembling using docker image.
        cmd = 'cp -r {local_folder} /hdfs'.format(local_folder=local_folder)
        os.system(cmd)
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dev/cowbat' #ashley added the /dev/
        cowbat_py = '/mnt/nas2/virtual_environments/dev/cowbat/bin/assembly_pipeline.py' #ashley added the /dev/
        # Run COWBAT with the necessary arguments
        cowbat_cmd = 'python {cowbat_py} -s {seqpath} -r {database}' \
            .format(cowbat_py=cowbat_py,
                    seqpath=os.path.join('/hdfs', sequence_folder),
                    database=COWBAT_DATABASES)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='COWBAT command:\n {cmd}'.format(cmd=cowbat_cmd))
        # Create another shell script to execute within the conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, cowbat_cmd)
        cowbat_script = os.path.join(work_dir, 'run_cowbat.sh')
        with open(cowbat_script, 'w+') as file:
            file.write(template)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(cowbat_script)
        # Run shell script
        os.system(cowbat_script)

        # Now need to move to an appropriate processed_sequence_data folder.
        local_wgs_spades_folder = os.path.join('/mnt/nas2/processed_sequence_data/miseq_assemblies', sequence_folder)
        cmd = 'mv {hdfs_folder} {wgsspades_folder}'.format(hdfs_folder=os.path.join('/hdfs', sequence_folder),
                                                           wgsspades_folder=local_wgs_spades_folder)
        print(cmd)
        os.system(cmd)

        # Remove the raw sequence files from processed_sequence_data, since we already have them in raw.
        cmd = 'rm {fastq_files}'.format(fastq_files=os.path.join(local_wgs_spades_folder, '*.fastq.gz'))
        print(cmd)
        os.system(cmd)

        # fix the legacy_combinedMetadata.csv file. Hopefully Adam will be able to remove this soon
        fix_legacy_metadata(os.path.join(local_wgs_spades_folder, 'reports/legacy_combinedMetadata.csv'))

        # Upload the results of the sequencing run to Redmine.
        cmd = 'cp {samplesheet} {reports_folder}'.format(samplesheet=os.path.join(local_wgs_spades_folder, 'SampleSheet.csv'),
                                                         reports_folder=os.path.join(local_wgs_spades_folder, 'reports'))
        os.system(cmd)
        shutil.make_archive(os.path.join(work_dir, sequence_folder), 'zip', os.path.join(local_wgs_spades_folder, 'reports'))
        output_list = list()
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, sequence_folder + '.zip')
        output_dict['filename'] = sequence_folder + '.zip'
        output_list.append(output_dict)

        # Apparently we're also supposed to be uploading assemblies - these will be too big to be Redmine attachments,
        # so we'll need to upload to the ftp.
        folder_to_upload = os.path.join(work_dir, 'reports_and_assemblies')
        os.makedirs(folder_to_upload)
        cmd = 'cp -r {best_assemblies} {upload_folder}'.format(best_assemblies=os.path.join(local_wgs_spades_folder, 'BestAssemblies'),
                                                               upload_folder=folder_to_upload)
        os.system(cmd)
        cmd = 'cp -r {reports} {upload_folder}'.format(reports=os.path.join(local_wgs_spades_folder, 'reports'),
                                                       upload_folder=folder_to_upload)
        os.system(cmd)
        shutil.make_archive(os.path.join(work_dir, str(issue.id)), 'zip', folder_to_upload)

        # At this point, zip folder has been created (hopefully) called issue_id.zip in biorequest dir. Upload that
        # to the FTP.
        upload_successful = upload_to_ftp(local_file=os.path.join(work_dir, str(issue.id) + '.zip'))

        # Make redmine tell Paul that a run has finished and that we should add things to our DB so things don't get missed
        # to be made. assinged_to_id to use is 226. Priority is 3 (High).
        if upload_successful:
            redmine_instance.issue.update(resource_id=issue.id,
                                          assigned_to_id=529, priority_id=3,
                                          subject='WGS Assembly: {}'.format(description[0]), # Add run name to subject
                                          notes='This run has finished assembly! Please add it to the OLC Database.\n'
                                                'Reports and assemblies uploaded to FTP at: ftp://ftp.agr.gc.ca/outgoing/'
                                                'cfia-ac/{}.zip'.format(issue.id))
        else:
            redmine_instance.issue.update(resource_id=issue.id,
                                          assigned_to_id=429,
                                          subject='WGS Assembly: {}'.format(description[0]),
                                          notes='Upload of result files was not successful. Upload them manually!')
        try:
            delete_ftp_dir(description[0])
        except:  # Hakuna matata if things don't get deleted. This is just a nice-to-have
            pass

    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))
        print(traceback.print_exc())


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


def download_ftp_file(ftp_file, local_dir):
    num_download_attempts = 0
    download_successful = False
    while num_download_attempts < 10:
        # Try downloading - if timeout, check if the download managed to complete but hang at the end, which happens
        # sometimes. If it did complete, we're good to go. Otherwise, try again.
        try:
            s = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD, timeout=30)
            s.cwd('incoming/cfia-ac')
            local_path = os.path.join(local_dir, os.path.split(ftp_file)[1])
            f = open(local_path, 'wb')
            s.retrbinary('RETR ' + ftp_file, f.write)
            f.close()
            # s.quit()
            quit_ftp(s)
            download_successful = True
            break
        except socket.timeout:
            local_path = os.path.join(local_dir, os.path.split(ftp_file)[1])
            s = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD, timeout=30)
            s.cwd('incoming/cfia-ac')
            ftp_file_size = s.size(ftp_file)
            # s.quit()
            quit_ftp(s)
            if os.path.isfile(local_path):
                if ftp_file_size == os.path.getsize(local_path):
                    download_successful = True
                    break
            num_download_attempts += 1
    return download_successful


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


def check_for_fastq_on_nas(samplesheet_seqids):
    # fastq_files_on_nas = glob.glob('/mnt/nas2/raw_sequence_data/miseq/*/*.fastq.gz')
    fastq_files_on_nas = glob.glob('/mnt/nas/MiSeq_Backup/*/*.fastq.gz')
    fastq_files_on_nas += glob.glob('/mnt/nas/External_MiSeq_Backup/*/*/*.fastq.gz')
    duplicate_samples = list()
    for seqid in samplesheet_seqids:
        forward_pattern = seqid + '*_R1*.gz'
        reverse_pattern = seqid + '*_R2*.gz'
        for item in fastq_files_on_nas:
            nas_fastq = os.path.split(item)[-1]
            if fnmatch.fnmatch(nas_fastq, forward_pattern) or fnmatch.fnmatch(nas_fastq, reverse_pattern) \
                    and seqid not in duplicate_samples:
                duplicate_samples.append(seqid)
    return duplicate_samples


def validate_fastq_run_stats(samplesheet_seqids, sequence_folder):
    seqids_in_xml = list()
    missing_seqids = list()
    tree = et.ElementTree(file=os.path.join(sequence_folder, 'GenerateFASTQRunStatistics.xml'))
    for element in tree.iter():
        if element.tag == 'SampleID':
            seqids_in_xml.append(element.text)
    for seqid in samplesheet_seqids:
        if seqid not in seqids_in_xml:
            missing_seqids.append(seqid)
    return missing_seqids


def ensure_samples_are_present(samplesheet_seqids, sequence_folder):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd(os.path.join('incoming/cfia-ac', sequence_folder))
    missing_samples = list()
    ftp_files = ftp.nlst()
    for seqid in samplesheet_seqids:
        forward_pattern = seqid + '*_R1*.gz'
        reverse_pattern = seqid + '*_R2*.gz'
        forward_found = False
        reverse_found = False
        for item in ftp_files:
            if fnmatch.fnmatch(item, forward_pattern):
                forward_found = True
            if fnmatch.fnmatch(item, reverse_pattern):
                reverse_found = True
        if forward_found is False or reverse_found is False:
            missing_samples.append(seqid)
    # ftp.quit()
    quit_ftp(ftp)
    return missing_samples



def get_seqids_from_samplesheet(samplesheet):
    regex = r'^(2\d{3}-\w{2,10}-\d{3,4})$'
    csv_names = list()
    found = False

    with open(samplesheet, 'r') as input_file:
        reader = csv.reader(input_file, delimiter=',')
        for row in reader:
            if row:
                # iterate through the document unto 'Sample_ID' is found in the first column
                if 'Sample_ID' in row[0]:
                    found = True
                    continue

                # once past the 'Sample_ID' add any SeqID to the list of files names
                if found:
                    if re.match(regex, row[0]):
                        csv_names.append(row[0].rstrip())
    return csv_names


# super legit fix for column mislabeling issue
def fix_legacy_metadata(infile):
    md = pd.read_csv(infile)
    colorder = md.columns
    md = md.rename(columns={"PipelineVersion": "assemblydate", "AssemblyDate": "confindrcontamsnvs", "ConfindrContamSNVs": "pipelineversion"})
    md = md.rename(columns={"pipelineversion": "PipelineVersion", "assemblydate": "AssemblyDate", "confindrcontamsnvs": "ConfindrContamSNVs"})
    md = md[colorder].set_index("SeqID")
    md.to_csv(infile)


def validate_files(file_name):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    missing_files = list()
    ftp.cwd(os.path.join('incoming/cfia-ac', file_name))
    files_present = ftp.nlst()
    if 'SampleSheet.csv' not in files_present:
        missing_files.append('SampleSheet.csv')
    if 'RunInfo.xml' not in files_present:
        missing_files.append('RunInfo.xml')
    if 'GenerateFASTQRunStatistics.xml' not in files_present:
        missing_files.append('GenerateFASTQRunStatistics.xml')
    # ftp.quit()
    quit_ftp(ftp)
    return missing_files


def download_info_sheets(sequence_folder, local_folder):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd(os.path.join('incoming/cfia-ac', sequence_folder))
    info_sheets = ['SampleSheet.csv', 'RunInfo.xml', 'GenerateFASTQRunStatistics.xml']
    for sheet in info_sheets:
        try:
            f = open(os.path.join(local_folder, sheet), 'wb')
            ftp.retrbinary('RETR ' + sheet, f.write)
            f.close()
        except:
            pass
    # ftp.quit()
    quit_ftp(ftp)


def verify_fastq_sizes(sequence_folder):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    tiny_fastqs = list()
    ftp.cwd(os.path.join('incoming/cfia-ac', sequence_folder))
    ftp_files = ftp.nlst()
    for item in ftp_files:
        if item.endswith('.gz') and 'Undetermined' not in item:
            file_size = ftp.size(item)
            if file_size < 1000:
                tiny_fastqs.append(item)
    # ftp.quit()
    quit_ftp(ftp)
    return tiny_fastqs


def verify_seqid_formatting(sequence_folder):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    badly_formatted_files = list()
    ftp.cwd(os.path.join('incoming/cfia-ac', sequence_folder))
    # Find all the files in the specified sequence folder.
    ftp_files = ftp.nlst()
    for item in ftp_files:
        # Anything ending in .gz is a FASTQ file, and needs to be checked for SEQID formatting. Ignore Undetermined,
        # it's special.
        if item.endswith('.gz') and 'Undetermined' not in item:
            wrong_formatting = False
            seqid = item.split('_')[0]  # SEQID should be everything before the first underscore.
            if len(seqid.split('-')) != 3:  # SEQID should be in format YYYY-LAB-####, so len should be three.
                wrong_formatting = True
            elif len(seqid.split('-')) == 3:
                # Get the year, lab, and sample number.
                year = seqid.split('-')[0]
                lab = seqid.split('-')[1]
                samplenum = seqid.split('-')[2]
                # Check that the year is a) four digits long, and b) a number.
                if len(year) != 4:
                    wrong_formatting = True
                try:
                    num = int(year)
                except ValueError:
                    wrong_formatting = True
                # We'll assume the lab is OK.
                # Check that the samplenum is four digits long and a number as well.
                if len(samplenum) != 4:
                    wrong_formatting = True
                try:
                    num = int(samplenum)
                except ValueError:
                    wrong_formatting = True
            if wrong_formatting:
                badly_formatted_files.append(item)
    quit_ftp(ftp)
    # ftp.quit()
    return badly_formatted_files


def verify_folder_exists(sequence_folder):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd('incoming/cfia-ac')
    folders_present = ftp.nlst()
    if sequence_folder in folders_present:
        folder_exists = True
    else:
        folder_exists = False
    # ftp.quit()
    quit_ftp(ftp)
    return folder_exists


def verify_folder_name(sequence_folder):
    # Check that the folder name given is in format YYMMDD_LAB
    # If anything else, boot the user.
    properly_formatted = True
    if len(sequence_folder.split('_')) != 2:
        properly_formatted = False
        return properly_formatted
    yymmdd = sequence_folder.split('_')[0]
    labid = sequence_folder.split('_')[1]
    if len(yymmdd) != 6:
        properly_formatted = False
    try:
        num = int(yymmdd)
    except ValueError:
        properly_formatted = False
    try:
        num = int(labid)
        properly_formatted = False
    except ValueError:
        pass
    return properly_formatted


def verify_all_the_things(sequence_folder, redmine_instance, issue, work_dir):
    validation = True
    if verify_folder_name(sequence_folder) is False:
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='ERROR: The folder name ({}) was not properly formatted. The correct'
                                            ' format is YYMMDD_LAB. Please create a new folder that is properly named'
                                            ' and create a new issue.'.format(sequence_folder))
        validation = False

    # Verify that the sequence folder specified does in fact exist. If it doesn't, give up.
    if verify_folder_exists(sequence_folder) is False:
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='ERROR: Could not find the folder ({}) specified in this issue on '
                                            'the FTP. Please ensure that it is uploaded correctly, create a new issue,'
                                            ' and try again.'.format(sequence_folder))
        validation = False
        return validation   # Can't check anything else if the folder doesn't exist, so stop here.

    # Check that SEQIDs are properly formatted.
    badly_formatted_fastqs = verify_seqid_formatting(sequence_folder)
    if len(badly_formatted_fastqs) > 0:
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='ERROR: The following FASTQ files did not have their SEQIDs formatted '
                                            'correctly: {}\n\nThe correct format is YYYY-LAB-####, where #### is the'
                                            ' 4-digit sample number. Please rename the files, reupload to the FTP, '
                                            'and try again.'.format(str(badly_formatted_fastqs)))
        validation = False

    # Verify that all the files uploaded that are .gz files are at least 100KB. Anything that is smaller than that
    # almost certainly didn't upload properly. Ignore undetermined.
    tiny_fastqs = verify_fastq_sizes(sequence_folder)
    if len(tiny_fastqs) > 0:
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='ERROR: The following FASTQ files had file sizes '
                                            'smaller than 1KB: {}\n\nThey likely did not upload '
                                            'properly. Please re-upload to the FTP and create a new '
                                            'issue.'.format(str(tiny_fastqs)))
        validation = False

    # Next up, validate that SampleSheet.csv, RunInfo, and GenerateFASTQRunStatistics are present.
    missing_files = validate_files(sequence_folder)
    if len(missing_files) > 0:
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='ERROR: The following files were missing from the FTP'
                                            ' folder: {}\nPlease reupload to '
                                            'the FTP, including these files ('
                                            'spelling must be identical!) and'
                                            ' create a new issue.'.format(str(missing_files)))
        validation = False

    # Now, download the info sheets (to a temporary folder) and make sure that SEQIDs that are present are good to go.
    if not os.path.isdir(os.path.join(work_dir, sequence_folder)):
        os.makedirs(os.path.join(work_dir, sequence_folder))
    download_info_sheets(sequence_folder, os.path.join(work_dir, sequence_folder))
    if 'SampleSheet.csv' in missing_files:
        return False
    else:
        samplesheet_seqids = get_seqids_from_samplesheet(os.path.join(work_dir, sequence_folder, 'SampleSheet.csv'))
        missing_seqids = ensure_samples_are_present(samplesheet_seqids, sequence_folder)
        if len(missing_seqids) > 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='ERROR: The following SEQIDs from SampleSheet.csv could not'
                                                ' be found in the folder you uploaded to the FTP: {}\nPlease re-upload'
                                                ' your files to the FTP and create a new issue.'.format(str(missing_seqids)))
            validation = False

    if 'GenerateFASTQRunStatistics.xml' not in missing_files:
        samplesheet_seqids = get_seqids_from_samplesheet(os.path.join(work_dir, sequence_folder, 'SampleSheet.csv'))
        missing_seqids = validate_fastq_run_stats(samplesheet_seqids, os.path.join(work_dir, sequence_folder))
        if len(missing_seqids) > 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='ERROR: The following SEQIDs from SampleSheet.csv could not'
                                                ' be found in GenerateFASTQRunStatistics.xml: {}\nPlease re-upload'
                                                ' your files to the FTP and create a new issue.'.format(str(missing_seqids)))
            validation = False

    # Now check that the SEQIDs from the SampleSheet don't already exist on the OLC NAS.
    samplesheet_seqids = get_seqids_from_samplesheet(os.path.join(work_dir, sequence_folder, 'SampleSheet.csv'))
    shutil.rmtree(os.path.join(work_dir, sequence_folder))
    duplicate_fastqs = check_for_fastq_on_nas(samplesheet_seqids)
    if len(duplicate_fastqs) > 0:
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='ERROR: The following SEQIDs already have FASTQ files on the OLC NAS: {}\n'
                                                'Please rename and reupload your files and create a '
                                                'new issue.'.format(str(duplicate_fastqs)))
        validation = False
    return validation


def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)


if __name__ == '__main__':
    wgsassembly_redmine()
