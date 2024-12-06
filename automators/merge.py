# Standard imports
from datetime import datetime
import glob
import os
import pickle
import shutil
import subprocess

# Third party imports
import click
from nastools.nastools import retrieve_nas_files
import pandas as pd
from pandas import ExcelWriter
import sentry_sdk


# Local imports
from amrsummary import before_send
from automator_settings import (
    COWBAT_DATABASES,
    COWBAT_IMAGE,
    SENTRY_DSN
)

# from externalretrieve import upload_to_ftp

from methods import (
    api_request_populate,
    api_request_python,
    api_request_script,
    azure_automate,
    azure_upload,
    construct_file_paths,
    construct_log_file_paths,
    create_shell_script,
    extract_job_id,
    make_executable,
)


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def merge_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    # Get the current date
    current_date = datetime.now()

    # Format the date as YYMMDD
    formatted_date = current_date.strftime('%y%m%d')

    # Set the container name using .format
    container = "{fd}-merge".format(fd=formatted_date)
    
    # Initialize variables to store user-defined arguments
    email = None
    basic_assembly = False
    preprocess = False
    
    try:
        
        # Set the name and path of the log files
        stdout_log_path, stderr_log_path = construct_log_file_paths(work_dir)
        
        # Parse description to determine requested arguments
        for item in description:

            # Remove whitespace from the items
            item = item.rstrip()
            if 'run_name' in item:
                # Extract the container name from the description. Set it
                # to lower case and replace underscores with hyphens
                container = item.split('=')[1].lower().replace('_', '-')
                continue
            if 'email' in item:
                # Extract the email address from the description
                email = item.split('=')[1]
                continue
            if 'basic_assembly' in item:
                # Set the basic_assembly flag to True
                basic_assembly = True
                continue
            if 'preprocess' in item:
                # Set the preprocess flag to True
                preprocess = True
                continue
        
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
                                          notes='ERROR: Did not find any attached files. Please create a new issue '
                                                'with the merge excel file attached and try again.',
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

        files_to_upload = merge_files(
            mergefile=os.path.join(work_dir, 'Merge.xlsx'),
            work_dir=work_dir
        )
        
        # Upload the merged files to cloud storage
        # Create a .tsv file for AzureAutomate
        azure_tsv = azure_automate(
            container=container,
            files_to_upload=files_to_upload,
            work_dir=work_dir
        )

        # Create the upload script
        azure_upload_script, template = azure_upload(
            azure_tsv=azure_tsv,
            work_dir=work_dir
        )

        # Update the Redmine issue after attempting to upload sequences
        notes = (
            'Uploading merged files to container {container} with command:\n '
            '{template}'.format(
                container=container,
                template=template
            )
        )
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes=notes,
            status_id=2
        )

        # Open the log files
        with open(stdout_log_path, 'a') as stdout_file, \
             open(stderr_log_path, 'a') as stderr_file:
            # Run shell script using subprocess.run and redirect output
            result = subprocess.run(
                ['bash', azure_upload_script],
                stdout=stdout_file, stderr=stderr_file,
            )
        
        # Run the merger script.
        # cmd = 'python /mnt/nas/Redmine/OLCRedmineAutomator/automators/merger.py -f {} -d ";" {}'.format(
        #     os.path.join(work_dir, 'Merge.xlsx'), work_dir)
        # os.system(cmd)

        # issue.watcher.add(226)  # Add Paul so he can put results into DB.
        # Make a folder to put all the merged FASTQs in biorequest folder. and put the merged FASTQs there.
        os.makedirs(os.path.join(work_dir, 'merged_' + str(issue.id)))
        cmd = 'mv {merged_files} {merged_folder}'\
            .format(merged_files=os.path.join(work_dir, '*MER*.fastq.gz'),
                    merged_folder=os.path.join(work_dir, 'merged_' + str(issue.id)))
        os.system(cmd)

        if len(glob.glob(os.path.join(work_dir, 'merged_' + str(issue.id), '*fastq.gz'))) == 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='ERROR: Something went wrong, no merged FASTQ files were created.',
                                          status_id=4)
            return
        # Now copy those merged FASTQS to merge backup
        cmd = 'cp {merged_files} /mnt/nas2/raw_sequence_data/merged_sequences'\
            .format(merged_files=os.path.join(work_dir, 'merged_' + str(issue.id), '*.fastq.gz'))
        os.system(cmd)

        # cmd = 'cp -r {merged_folder} /hdfs'.format(merged_folder=os.path.join(work_dir, 'merged_' + str(issue.id)))
        # os.system(cmd)

        redmine_instance.issue.update(resource_id=issue.id,
                                      #notes='Merged FASTQ files created, beginning assembly of merged files.')
                                      notes='Merged FASTQ files created. Uploading to FoodPort.')
        
        
        # Make an API request to FoodPort to run a research assembly on the
        # merged files
        api_log_path, script_content, script_content_redacted = \
            api_request_populate(
            basic_assembly=basic_assembly,
            container=container,
            email=email,
            preprocess=preprocess,
            work_dir=work_dir
        )

        # Make a Python script to create the API call
        api_call_script = api_request_python(
            script_content=script_content,
            work_dir=work_dir
        )
        
        # Create a shell script to run the Python script
        foodport_api_script = api_request_script(
            api_call_script=api_call_script,
            work_dir=work_dir
        )

        # Open the log files
        with open(stdout_log_path, 'a') as stdout_file, \
             open(stderr_log_path, 'a') as stderr_file:
            # Run shell script using subprocess.run and redirect output
            subprocess.run(
                ['bash', foodport_api_script],
                stdout=stdout_file, stderr=stderr_file,
                check=True
            )
            
        # Open and read the .err log file
        with open(api_log_path, 'r') as file:
            log_content = file.read()

        # Define the error message to search for
        error_message = "Exception Value: duplicate key value violates " \
            "unique constraint"
        
        # Check if the specific error message is in the log file
        if error_message in log_content:
            # Define the notes to update the Redmine issue with
            notes = (
                "The supplied run_name, {container} already exists in "
                "FoodPort. Please check the details and submit a unique "
                "name.".format(container=container)
            )
            
            # Update the Redmine issue
            redmine_instance.issue.update(
                resource_id=issue.id,
                notes=notes,
                status_id=3
            )
            return
            
        # Create notes for updating the Redmine issue after submitting the
        # API request
        notes = (
            'Running COWBAT on merged files in container {container} with '
            'API call:\n {script_content_redacted}\nPlease follow the '
            'progress on FoodPort'.format(
                container=container,
                script_content_redacted=script_content_redacted
            )
        )
        
        # Add an additional note if an email address was provided
        if email:
            notes += (
                '\nAn email will be sent to {email} when the assembly is '
                'complete.'.format(email=email)
            )
        
        # Update the Redmine issue
        redmine_instance.issue.update(
            resource_id=issue.id,
            notes=notes,
            status_id=4
        )
        
        #MODIFIED ON 2024-05-10. Will now upload the fastq files to the ftp, so that the requestor can upload to foodport. (the pipeline on the nas is old and broken... for now.)

        #commented out below lines 2024-05-10 AC
        # With files copied over to the HDFS, start the assembly process (Now using new pipeline!)
        # These unfortunate hard coded paths appear to be necessary
        # activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dev/cowbat'
        # assembly_pipeline = '/mnt/nas2/virtual_environments/dev/cowbat/bin/assembly_pipeline.py'
        # Run sipprverse with the necessary arguments
        # run_cmd = 'python {assembly_pipeline} -s {seqpath} -r {dbpath}' \
        #     .format(assembly_pipeline=assembly_pipeline,
        #             seqpath=os.path.join('/hdfs', 'merged_' + str(issue.id)),
        #             dbpath=COWBAT_DATABASES)
        #
        # redmine_instance.issue.update(resource_id=issue.id,
        #                               notes='COWBAT command:\n {cmd}'.format(cmd=run_cmd))
        # # Create another shell script to execute within the conda environment
        # template = "#!/bin/bash\n{} && {}".format(activate, run_cmd)
        # script = os.path.join(work_dir, 'run_pipeline.sh')
        # with open(script, 'w+') as file:
        #     file.write(template)
        # # Modify the permissions of the script to allow it to be run on the node
        # make_executable(script)
        # # Run shell script
        # os.system(script)
        #
        # # Move results to merge_WGSspades, and upload the results folder to redmine.
        # cmd = 'mv {hdfs_folder} {merge_WGSspades}'\
        #     .format(hdfs_folder=os.path.join('/hdfs', 'merged_' + str(issue.id)),
        #             merge_WGSspades=os.path.join('/mnt/nas2/processed_sequence_data/merged_assemblies',
        #                                          'merged_' + str(issue.id) + '_Assembled'))
        # os.system(cmd)
        # shutil.make_archive(os.path.join(work_dir, 'reports'), 'zip',
        #                     os.path.join('/mnt/nas2/processed_sequence_data/merged_assemblies',
        #                                  'merged_' + str(issue.id) + '_Assembled', 'reports'))
        # output_list = list()
        # output_dict = dict()
        # output_dict['path'] = os.path.join(work_dir, 'reports.zip')
        # output_dict['filename'] = 'merged_' + str(issue.id) + '_reports.zip'
        # output_list.append(output_dict)
        # redmine_instance.issue.update(resource_id=issue.id,
        #                               uploads=output_list,
        #                               status_id=4,
        #                               notes='Merge Process Complete! Reports attached.')

        #instead, we will zip the output sequences file and upload to the ftp
        # Zip output
        # out_filename = 'Merged_sequences_{}'.format(issue.id)
        # merged_files=os.path.join(work_dir, 'merged_' + str(issue.id))
        # zip_filepath = zip_folder(results_path=merged_files,
        #                           output_dir=work_dir,
        #                           output_filename=out_filename)
        # zip_filepath += '.zip'
        
        # upload_successful = upload_to_sftp(local_file=zip_filepath)

        # # Prepare upload
        # if upload_successful:
        #     redmine_instance.issue.update(resource_id=issue.id, uploads=output_list, status_id=4,
        #                                   notes='Sequence merge complete!\n\n'
        #                                         'The merged sequences can be downloaded from:\n'
        #                                         'ftp://ftp.agr.gc.ca/outgoing/cfia-ac/{} \n'
        #                                         'Please download the raw sequences, then upload them to foodport as a Research Assembly'
        #                                   .format(os.path.split(zip_filepath)[1]))
        # else:
        #     redmine_instance.issue.update(resource_id=issue.id, status_id=4,
        #                                   notes='Upload of result files was unsuccessful due to FTP connectivity '
        #                                         'issues. '
        #                                         'Please try again later.')

        # Remove all the folders
        # shutil.rmtree(merged_files)
        # Remove the zip file
        # os.remove(zip_filepath)


    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong: \n {error}.\n '
                                            'We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.'.format(error=e))


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


def merge_files(mergefile, work_dir):
    df = pd.read_excel(mergefile)
    
    # Define a list to store the files to upload
    files_to_upload = []
    for i in range(len(df['Name'])):
        seqids = df['Merge'][i].split(';')
        merge_name = df['Name'][i]
        merge_forward_reads = os.path.join(work_dir, merge_name + '_S1_L001_R1_001.fastq.gz')
        merge_reverse_reads = os.path.join(work_dir, merge_name + '_S1_L001_R2_001.fastq.gz')
        
        # Add the forward and reverse files to the list of files to upload
        files_to_upload.extend([merge_forward_reads, merge_reverse_reads])
        for seqid in seqids:
            forward = glob.glob(os.path.join(work_dir, seqid + '*_R1*.fastq.gz'))[0]
            reverse = glob.glob(os.path.join(work_dir, seqid + '*_R2*.fastq.gz'))[0]
            
            cmd = 'cat {forward} >> {merge_forward_reads}'.format(forward=forward,
                                                                  merge_forward_reads=merge_forward_reads)
            os.system(cmd)
            cmd = 'cat {reverse} >> {merge_reverse_reads}'.format(reverse=reverse,
                                                                  merge_reverse_reads=merge_reverse_reads)
            os.system(cmd)
    
    return files_to_upload


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

def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    merge_redmine()
