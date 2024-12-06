import os
import re
import csv
import glob
import click
import pickle
import shutil
import fnmatch
from nastools.nastools import retrieve_nas_files

from ftplib import FTP
from automator_settings import FTP_USERNAME, FTP_PASSWORD
import traceback

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def hybridassembly_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        sequence_folder = description[0]

        hybrid_info = list()
        for i in range(1, len(description)):
            x = description[i].rstrip().split(',')
            minion_id = x[0]
            seqid = x[1]
            oln_id = x[2]
            hybrid_info.append([minion_id, seqid, oln_id])

        local_folder = os.path.join('/mnt/nas2/raw_sequence_data/nanopore', sequence_folder)

        download_dir(ftp_dir=sequence_folder,
                     local_dir=local_folder)

        seqids = list()
        for item in hybrid_info:
            seqids.append(item[1])

        # Link the FASTQs needed to run hybrid assemblies to our working dir.
        retrieve_nas_files(seqids=seqids,
                           outdir=os.path.join(work_dir, 'fastqs'),
                           filetype='fastq')

        # TODO: Figure out conda env activation, it's always a pain
        # Now need to run the hybrid_assembly.py script, which will run Unicycler.
        for item in hybrid_info:
            minion_reads = glob.glob(os.path.join(local_folder, item[0] + '*'))[0]
            illumina_forward = glob.glob(os.path.join(work_dir, 'fastqs', item[1] + '*_R1*'))[0]
            illumina_reverse = glob.glob(os.path.join(work_dir, 'fastqs', item[1] + '*_R1*'))[0]
            output_dir = os.path.join(work_dir, item[2])
            cmd = 'python /mnt/nas/Redmine/OLCRedmineAutomator/automators/hybrid_assembly.py -1 {forward} -2 {reverse} ' \
                  '-p {minion} -o {output} -n {name}'.format(forward=illumina_forward,
                                                             reverse=illumina_reverse,
                                                             minion=minion_reads,
                                                             output=output_dir,
                                                             name=item[2])
            os.system(cmd)

        # Now that this is done, need to make a report that looks at least a bit like our old combined metadata report


        # At this point, zip folder has been created (hopefully) called issue_id.zip in biorequest dir. Upload that
        # to the FTP.
        s = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
        s.cwd('outgoing/cfia-ac')
        f = open(os.path.join(work_dir, str(issue.id) + '.zip'), 'rb')
        s.storbinary('STOR {}.zip'.format(str(issue.id)), f)
        f.close()
        s.quit()

        # Make redmine tell Paul that a run has finished and that we should add things to our DB so things don't get missed
        # to be made. assinged_to_id to use is 226. Priority is 3 (High).
        redmine_instance.issue.update(resource_id=issue.id,
                                      assigned_to=226,
                                      notes='Hybrid assembly complete. Please add it to the OLC Database.\n')

    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))
        print(traceback.print_exc())


def download_dir(ftp_dir, local_dir):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd(os.path.join('incoming/cfia-ac', ftp_dir))
    present_in_folder = ftp.nlst()
    for item in present_in_folder:
        if check_if_file(item, ftp_dir):
            local_path = os.path.join(local_dir, item)
            print(local_path)
            f = open(local_path, 'wb')
            ftp.retrbinary('RETR ' + item, f.write)
            f.close()
        else:
            if not os.path.isdir(os.path.join(local_dir, item)):
                os.makedirs(os.path.join(local_dir, item))
            download_dir(os.path.join(ftp_dir, item), os.path.join(local_dir, item))
    ftp.quit()


def check_if_file(file_name, ftp_dir):
    ftp = FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
    ftp.cwd(os.path.join('incoming/cfia-ac', ftp_dir))
    is_file = True
    try:
        ftp.size(file_name)
    except:
        is_file = False
    ftp.quit()
    return is_file


if __name__ == '__main__':
    hybridassembly_redmine()

