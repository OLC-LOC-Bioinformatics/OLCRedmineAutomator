import os
import glob
import click
import pickle
import shutil
from biotools import mash
from amrsummary import before_send
import sentry_sdk
from automator_settings import SENTRY_DSN
from nastools.nastools import retrieve_nas_files
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
import fileinput
from pathlib import Path
import csv
import fnmatch



@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def subsampler_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    #list of listtypes (haha) users can select
    listtypes = [
        'attached', 'below'
    ]

    # Variable to hold supplied arguments
    argument_dict = {
        'list': 'attached',
        'coverage': 'default',
    }

    try:
        #create an empty list for the seqids
        seqids = list()
        # Parse description to figure out whether the user has attached a list, or put the items in the description
        for item in description:
            item = item.rstrip()
            if 'list' in item:
                argument_dict['list'] = item.split('=')[1]
                continue
            if 'coverage' in item:
                argument_dict['coverage'] = item.split('=')[1]
                continue
            # Otherwise the item should be a SEQID
            #seqids.append(item)

        # Ensure that SEQIDs were included
        #if not seqids:
        #    redmine_instance.issue.update(resource_id=issue.id,
        #                                  notes='WARNING: No SEQIDs provided!',
        #                                  status_id=4)
        #    return


        # If user specified attachment as the reference file, download it to our working directory.
        if argument_dict['list'] == 'attached':
            # Get the attachment ID, and download if it isn't equal to zero (meaning no attachment, so boot user with
            # appropriate error message)
            attachment = redmine_instance.issue.get(issue.id, include='attachments')
            attachment_id = 0
            for item in attachment.attachments:
                attachment_id = item.id

            # Download if we found an attachment, and use as our reference. Otherwise, exit and tell user to try again
            if attachment_id != 0:
                attachment = redmine_instance.attachment.get(attachment_id)
                #attachment.download(savepath=ref_dir, filename='reference.fasta')
                attachment.download(savepath=work_dir)
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: You specified that the list would be an attached file,'
                                                    ' but no attached file was found. Please create a new issue and '
                                                    'try again.',
                                              status_id=4)
                return

        #now we need to parse either the description or the attached file in order to get the seqids, etc, for the analysis
        #first get a list of attached files..
        attachedfiles = os.listdir(work_dir)
        #now parse through attached csv file to get a list of the seqids we need
        for file in attachedfiles:
            if file.endswith('.csv'):
                with open(os.path.join(work_dir,file), 'r') as csv_file:
                    reader = csv.reader(csv_file)
                    #look through each row, and pull out the second item (SEQID) and add it to our seqid list
                    for row in reader:
                        seqids.append(row[1])
        

        #create a folder to hold the sequences
        seq_dir = os.path.join(work_dir, 'sequences')
        os.makedirs(seq_dir, exist_ok=True)

        #now use the nastools to retrieve all of the sequences in our list
        retrieve_nas_files(seqids=seqids, outdir=seq_dir, filetype='fastq', copyflag=False)
        missing_fastqs = check_fastqs_present(seqids, seq_dir)
        if len(missing_fastqs) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastqs))


        #create an output folder
        out_dir = os.path.join(work_dir, 'output')
        os.makedirs(out_dir, exist_ok=True)

        #list of the coverage levels we are looking for
        if argument_dict['coverage'] != 'default':
            #have to change the users' list into a list of integers... see below
            covgs = argument_dict['coverage']
            covglist = covgs.split(",")
            li = []
            for i in covglist:
                #li.append(int(i))
                li.append(float(i))
                print(float(i))
            coverage_levels = li
        else:
            coverage_levels = [1, 2.5, 5, 7.5, 10, 20]


        #need to create a script that we will run on the cluster... first need the environment we will activate
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/hybrid_assembly'
        runbbmap = os.path.join(work_dir, 'runbbmap_reformat.sh')

        #there might be a better way to do this.... but will try reading through each row and doing manually for now..
        for file in attachedfiles:
            if file.endswith(".csv"):
                with open(os.path.join(work_dir,file), 'r') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        strain_name = row[0]
                        seqid = row[1]
                        total_assembly_length = int(row[2])
                        coverage_bases = [int(total_assembly_length * coverage) for coverage in coverage_levels]
                        #create output directory for that particlar strain
                        output_dir = os.path.join(out_dir, '{strain_name}_coverage_fastq'.format(strain_name=strain_name))
                        os.makedirs(output_dir, exist_ok=True)
                        #generate 10 unique .fastq files for each coverage level using BBMap 'reformat.sh'
                        for i, bases in enumerate(coverage_bases):
                            for j in range(1, 11):
                                seqname = "{seqid}.fastq.gz".format(seqid=seqid)
                                input_fastq = os.path.join(seq_dir,seqname)
                                output_fastq = os.path.join(output_dir, '{}_{}x_{}.fastq.gz'.format(strain_name,coverage_levels[i],j))
                                cmd = 'reformat.sh in={fastq} out={outfq} samplebasestarget={bases}'.format(fastq=input_fastq, outfq=output_fastq, bases=bases)
                                template = "#!/bin/bash\n {} && {}".format(activate, cmd)
                                #now write the above into the runbbmap file we created
                                with open(runbbmap, 'w+') as g:
                                    g.write(template)
                                make_executable(runbbmap)
                                #run the script on the cluster
                                os.system(runbbmap)

        #the error file contains all of the terminal outputs, so we will copy that to the output folder just in case someone needs it
        workdirfiles = os.listdir(work_dir)
        for file in workdirfiles:
            if file.endswith(".err"):
                shutil.copy(os.path.join(work_dir, file), os.path.join(out_dir, file))
        
                        
        # Zip output
        subsampler_filename = 'subsampler_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=out_dir,
                                  output_dir=work_dir,
                                  output_filename=subsampler_filename)
        zip_filepath += '.zip'
        
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
                notes='Subsampling complete!\n\n'
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

        # Remove all the folders
        shutil.rmtree(seq_dir)
        shutil.rmtree(out_dir)

        # Remove the zip file
        #os.remove(zip_filepath)

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='Something went wrong! Please contact a bioinformatician '
                                            'to investigate. Here is the traceback: {e}'.format(e=e))




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

def check_fastqs_present(fastq_list, fastq_dir):
    missing_fastqs = list()
    for seqid in fastq_list:
        if len(glob.glob(os.path.join(fastq_dir, seqid + '*.fastq.gz'))) < 2:
            # JAS adding an extra if statement here to allow for Nanopore (SE) reads
            if len(glob.glob(os.path.join(fastq_dir, seqid + ".fastq.gz"))) == 0:
                missing_fastqs.append(seqid)
    return missing_fastqs

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


if __name__ == '__main__':
    subsampler_redmine()
