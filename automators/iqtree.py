#!/usr/bin/env python

import os
import glob
import click
import pickle
import shutil
import sentry_sdk
import subprocess
import shutil
import xml.etree.ElementTree as et
from externalretrieve import upload_to_ftp
from biotools import mash
from amrsummary import before_send
from strainchoosr import strainchoosr
import ftplib
from ftplib import FTP
from automator_settings import SENTRY_DSN, FTP_USERNAME, FTP_PASSWORD
import traceback
from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def iqtree_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))
        #analyses supported by the automator
        analyses = [
            'custom', 'enterobacterales', 'listeriaceae', 'enterobacter', 'attached', 'attached_ftp'
        ]

        # Variable to hold supplied arguments
        argument_dict = {
            'ftp_folder': False,
            'analysis': 'attached',
            'alignment': False,
            'genomesize': 500000,
        }

        # Parse description for SEQIDs, write list that file_extractor needs.
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'FTP_FOLDER' in item:
               argument_dict['ftp_folder'] = item.split('=')[1].lower()
               continue
            if 'ANALYSIS' in item:
               argument_dict['analysis'] = item.split('=')[1].lower()
               continue
            if 'ALIGNMENT' in item:
               argument_dict['alignment'] = item.split('=')[1].lower()
               continue
            if 'GENOMESIZE' in item:
               argument_dict['genomesize'] = item.split('=')[1]
               continue
            if 'MINDEPTH' in item:
               argument_dict['mindepth'] = item.split('=')[1]
               continue
            if 'SKETCH-SIZE' in item:
               argument_dict['sketch-size'] = item.split('=')[1]
               continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that the analysis type is provided
        if not argument_dict['analysis']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not identify an analysis type. '
                                                'Please ensure that the first line of the issue contains one'
                                                ' of the following keywords: {ats}'.format(ats=', '.join(analyses)),
                                          status_id=4)

        elif argument_dict['analysis'] not in analyses:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied analysis type {at} current not in the supported '
                                                'list of analyses: {ats}'.format(at=argument_dict['analysis'],
                                                                                 ats=', '.join(analyses)),
                                          status_id=4)


        #create needed folders
        input_folder = os.path.join(work_dir, 'inputs')
        os.mkdir(alignments_folder)

        output_folder = os.path.join(work_dir, 'output')
        os.mkdir(output_folder)

        # Create the local folder that we'll need to bring our ftp files into.
        local_folder = os.path.join(work_dir, 'holding_folder')

    # link to order enterobacterales fasta files if analysis = enterobacterales
        if argument_dict['analysis'] == 'enterobacterales':
            src = '/mnt/nas2/processed_sequence_data/ncbi/enterobacterales_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

    #if listeriaceae, link to listeria reference sequences
        if argument_dict['analysis'] == 'listeriaceae':
            src = '/mnt/nas2/processed_sequence_data/ncbi/listeriaceae_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

    #if enterobacter, link to enterobacter reference sequences
        if argument_dict['analysis'] == 'enterobacter':
            src = '/mnt/nas2/processed_sequence_data/ncbi/enterobacter_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)
            src2 = '/mnt/nas2/processed_sequence_data/ncbi/enterobacterales_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd2 = 'ln -s {src}ENTEROBACTER_*.fasta {dst}'.format(src=src2,
                                                                    dst=dst)
            os.system(lncmd2)

    #if just a custom analysis, create the fastas folder
        if argument_dict['analysis'] == 'custom':
            # Make fasta file directory
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
    

        # Drop FASTA files into workdir
        retrieve_nas_files(seqids=seqids,
                           outdir=os.path.join(work_dir, 'fastas'),
                           filetype='fasta',
                           copyflag=False)

        # Download the attached files.
        # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
        # out a better way to do it.
        attachment = redmine_instance.issue.get(issue.id, include='attachments')
        attachment_id = 0
        for item in attachment.attachments:
            attachment_id = item.id
        # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
        if attachment_id != 0:
            for item in attachment.attachments:
                attachment_id = item.id
                attachment = redmine_instance.attachment.get(attachment_id)
                attachment.download(savepath=work_dir)


        #see if there were fasta files attached and move them to the fasta files directory
        if os.path.isfile(os.path.join(work_dir,'*.fasta')):
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Attached fasta file found.',
                                          status_id=2)
            lncmd = 'ln -s {wrk}*.fasta {dst}'.format(wrk=work_dir,
                                                      dst=fasta_dir)
            os.system(lncmd)

        if argument_dict['analysis'] == 'attached':
            if not os.path.isfile(os.path.join(work_dir,'alignment*')):
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Requested phylogeny analysis requires an attached alignment '
                                                    'file. This file must start with "alignment" in the name. '
                                                    'The automator could not find an attached alignment file. '
                                                    'Please create a new issue with the alignment file attached and try '
                                                    'again.',
                                              status_id=4)
            else:
                if os.path.isfile(os.path.join(work_dir, 'alignment*')):
                    redmine_instance.issue.update(resource_id=issue.id,
                                                  notes='Attached alignment file(s) found.',
                                                  status_id=2)
                #move any downloaded alignment files to the specified folder
                all_files = os.listdir(work_dir)
                for file in all_files:
                    if file.startswith("alignment"):
                        shutil.move(os.path.join(work_dir, file), os.path.join(input_folder, file))

        if argument_dict['analysis'] == 'attached_ftp':
            if not argument_dict['ftp_folder']:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Requested phylogeny analysis requires a attached alignment '
                                                    'file(s). '
                                                    'The automator could not find a listed ftp-folder to download. '
                                                    'Please create a new issue with the ftp_folder indicated and try '
                                                    'again.',
                                              status_id=4)
                return
            else:
#                if not os.path.isdir(local_folder):
#                    os.makedirs(local_folder)
                ftp_fold=argument_dict['ftp_folder']
                download_dir(ftp_fold, local_folder)

#                ftp_fold=argument_dict['ftp_folder']
#                download_successful = download_dir(ftp_fold, local_folder)

#                if download_successful is False:
#                    redmine_instance.issue.update(resource_id=issue.id,
#                                                  notes='Download of files from FTP was not successful.',
#                                                  status_id=4)



            #move any downloaded alignment files to the specified folder
            ftpdownloaded_files = os.listdir(local_folder)
            for file in ftpdownloaded_files:
#                if file.startswith("alignment"):
                shutil.move(os.path.join(local_folder, file), os.path.join(input_folder, file))


	# Verify that specified fasta files are actually there, warn user if they aren't.
        missing_fastas = verify_fasta_files_present(seqid_list=seqids,
                                                    fasta_dir=os.path.join(work_dir, 'fastas'))
        if len(missing_fastas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))

        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/iqtree'
        # Run iqtree with the necessary arguments
        iqtree_cmd = 'iqtree -s {inp} --prefix iqtree-output -T 10 '.format(inp=input_folder)

        #create another shell script to execute within the dRep conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, iqtree_cmd)
        iqtree_script = os.path.join(work_dir, 'run_iqtree.sh')
        with open(iqtree_script, 'w+') as file:
            file.write(template)
        make_executable(iqtree_script)

        # Run shell script
        os.system(iqtree_script)

        #move the output files to the output folder because there doesnt seem to be an out option with iqtree
        work_files = os.listdir(work_dir)
        for file in work_files:
            if file.startswith("iqtree-output"):
                shutil.move(os.path.join(work_dir, file), os.path.join(output_folder, file))

        #zip up the output file(s)
        output_filename = 'iqtree_output'
        zip_filepath = zip_folder(results_path=os.path.join(work_dir, 'output'),
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

        try:
            shutil.rmtree(os.path.join(work_dir, 'fastas'))
        except FileNotFoundError:
            pass

        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes='{at} analysis with iqtree complete!'
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

if __name__ == '__main__':
    iqtree_redmine()
