import os
import glob
import click
import pickle
import shutil
import pandas as pd
from pandas import ExcelWriter
from nastools.nastools import retrieve_nas_files


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def merge_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        # Download the attached excel file.
        # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
        # out a better way to do it.
        attachment = redmine_instance.issue.get(issue.id, include='attachments')
        attachment_id = 0
        for item in attachment.attachments:
            attachment_id = item.id

        # Now download, if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
        if attachment_id != 0:
            attachment = redmine_instance.attachment.get(attachment_id)
            attachment.download(savepath=work_dir, filename='merge.xlsx')
        else:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='ERROR: Did not find any attached files. Please create a new issue with '
                                                'the merge excel file attached and try again.',
                                          status_id=4)
            return

        # Now use convert_excel_file to make compatible with merger.py
        convert_excel_file(os.path.join(work_dir, 'merge.xlsx'), os.path.join(work_dir, 'Merge.xlsx'))

        # Make a SEQID list of files we'll need to extract.
        seqid_list = generate_seqid_list(os.path.join(work_dir, 'Merge.xlsx'))

        # Create links of all seqids needed in our working dir
        retrieve_nas_files(seqids=seqid_list,
                           outdir=work_dir,
                           filetype='fastq',
                           copyflag=False)

        # Run the merger script.
        cmd = 'python /mnt/nas2/redmine/applications/OLCRedmineAutomator/automators_dev/merger.py -f {} -d ";" {}'\
            .format(os.path.join(work_dir, 'Merge.xlsx'), work_dir)
        os.system(cmd)

        # Make a folder to put all the merged FASTQs in biorequest folder. and put the merged FASTQs there.
        os.makedirs(os.path.join(work_dir, 'merged_' + str(issue.id)))
        cmd = 'mv {merged_files} {merged_folder}'.format(merged_files=os.path.join(work_dir, '*MER*/*.fastq.gz'),
                                                         merged_folder=os.path.join(work_dir, 'merged_' + str(issue.id)))
        os.system(cmd)

        if len(glob.glob(os.path.join(work_dir, 'merged_' + str(issue.id), '*fastq.gz'))) == 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='ERROR: Something went wrong, no merged FASTQ files were created.',
                                          status_id=4)
            return
        # Now copy those merged FASTQS to merge backup and the hdfs folder so they can be assembled.
        cmd = 'cp {merged_files} /mnt/nas/merge_Backup'.format(merged_files=os.path.join(work_dir, 'merged_' + str(issue.id),
                                                                                         '*.fastq.gz'))
        os.system(cmd)

        cmd = 'cp -r {merged_folder} /hdfs'.format(merged_folder=os.path.join(work_dir, 'merged_' + str(issue.id)))
        os.system(cmd)

        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Merged FASTQ files created, beginning assembly of merged files.')
        # With files copied over to the HDFS, start the assembly process.
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dev/cowbat'
        assembly_pipeline = '/mnt/nas2/virtual_environments/dev/cowbat/bin/assembly_pipeline.py'
        # Run sipprverse with the necessary arguments
        run_cmd = 'python {assembly_pipeline} -s {seqpath} -r {dbpath}' \
            .format(assembly_pipeline=assembly_pipeline,
                    seqpath=os.path.join('/hdfs', 'merged_' + str(issue.id)),
                    dbpath=COWBAT_DATABASES)

        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='COWBAT command:\n {cmd}'.format(cmd=run_cmd))
        # Create another shell script to execute within the conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, run_cmd)
        script = os.path.join(work_dir, 'run_pipeline.sh')
        with open(script, 'w+') as file:
            file.write(template)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(script)
        # Run shell script
        os.system(script)

        # Move results to merge_WGSspades, and upload the results folder to redmine.
        cmd = 'mv {hdfs_folder} {merge_WGSspades}'.format(hdfs_folder=os.path.join('/hdfs', 'merged_' + str(issue.id)),
                                                          merge_WGSspades=os.path.join('/mnt/nas/merge_WGSspades',
                                                          'merged_' + str(issue.id) + '_Assembled'))
        os.system(cmd)
        shutil.make_archive(os.path.join(work_dir, 'reports'), 'zip', os.path.join('/mnt/nas/merge_WGSspades', 'merged_' + str(issue.id) + '_Assembled', 'reports'))
        output_list = list()
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'reports.zip')
        output_dict['filename'] = 'merged_' + str(issue.id) + '_reports.zip'
        output_list.append(output_dict)
        redmine_instance.issue.update(resource_id=issue.id, uploads=output_list, status_id=4,
                                      notes='Merge Process Complete! Reports attached.')
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))


def convert_excel_file(infile, outfile):
    df = pd.read_excel(infile)
    to_keep = ['SEQID', 'OtherName']
    for column in df:
        if column not in to_keep:
            df = df.drop(column, axis=1)
    df = df.rename(columns={'SEQID': 'Name', 'OtherName': 'Merge'})
    writer = ExcelWriter(outfile)
    df.to_excel(writer, 'Sheet1', index=False)
    writer.save()


def generate_seqid_list(mergefile):
    df = pd.read_excel(mergefile)
    seqid_list = list()
    seqids = list(df['Merge'])
    for row in seqids:
        for item in row.split(';'):
            seqid_list.append(item)
    return seqid_list

def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)


if __name__ == '__main__':
    merge_redmine()