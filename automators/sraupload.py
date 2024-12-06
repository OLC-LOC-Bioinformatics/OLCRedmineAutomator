import os
import glob
import click
import ftplib
import pickle
import shutil
from nastools.nastools import retrieve_nas_files


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def externalretrieve_redmine(redmine_instance, issue, work_dir, description):
    print('SRA upload')
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        # Parse description to figure out what SEQIDs we need to run on.
        ftp_user = description.pop(0).rstrip()
        ftp_pass = description.pop(0).rstrip()
        ftp_folder = description.pop(0).rstrip()
        fastq_list = list()
        for item in description:
            #item = item.upper()
            item = item.rstrip()
            fastq_list.append(item)

        # Use NAStools to put FASTQ files into our working dir.
        retrieve_nas_files(seqids=fastq_list,
                           outdir=work_dir,
                           filetype='fastq',
                           copyflag=True)

        # Check that we got all the requested files.
        missing_fastqs = check_fastqs_present(fastq_list, work_dir)
        if len(missing_fastqs) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastqs))


        # Rename files to _R1 _R2, without anything else
        renamed_fastqs = rename_files(fastq_dir=work_dir)

        # Now need to login to the FTP to upload the zipped folder.
        try:
            s = ftplib.FTP('ftp-private.ncbi.nlm.nih.gov', user=ftp_user, passwd=ftp_pass)
        except ftplib.error_perm:
            redmine_instance.issue.update(resource_id=issue.id,
                                          status_id=4,
                                          notes='ERROR: Could not connect to the NCBI FTP site with the username '
                                                'and password provided. Make sure the username provided to you is the '
                                                'first line of the description, and the password is the second.')
            return

        try:
            s.cwd(ftp_folder)
        except ftplib.error_perm:
            redmine_instance.issue.update(resource_id=issue.id,
                                          status_id=4,
                                          notes='ERROR: Could not find the directory you specified on the NCBI FTP '
                                                'site. Make sure the directory is specified on the third line of the '
                                                'description.')
            return
        s.mkd(str(issue.id))
        s.cwd(str(issue.id))
        for fastq in renamed_fastqs:
            f = open(fastq, 'rb')
            s.storbinary('STOR {}'.format(os.path.split(fastq)[1]), f)
            f.close()
        s.quit()

        # Finally, do some file cleanup.
        try:
            os.system('rm {}/*.fastq.gz'.format(work_dir))
        except:
            pass

        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='SRA Retrieve process complete! You should now be able to select the FTP '
                                            'folder called {} when prompted for your SRA submission.\n\n'.format(issue.id))
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))


def rename_files(fastq_dir):
    renamed_fastqs = list()
    fastq_files = glob.glob(os.path.join(fastq_dir, '*.fastq.gz'))
    for fastq in fastq_files:
        fastq_base = os.path.split(fastq)[1]
        if '_R1' in fastq:
            newname = fastq_base.split('_')[0] + '_R1.fastq.gz'
        elif '_R2' in fastq:
            newname = fastq_base.split('_')[0] + '_R2.fastq.gz'
        cmd = 'mv ' + fastq + ' ' + os.path.join(fastq_dir, newname)
        print(cmd)
        os.system(cmd)
        renamed_fastqs.append(os.path.join(fastq_dir, newname))
    return renamed_fastqs


def check_fastqs_present(fastq_list, fastq_dir):
    missing_fastqs = list()
    for seqid in fastq_list:
        if len(glob.glob(os.path.join(fastq_dir, seqid + '*.fastq.gz'))) < 2:
            missing_fastqs.append(seqid)
    return missing_fastqs


if __name__ == '__main__':
    externalretrieve_redmine()
