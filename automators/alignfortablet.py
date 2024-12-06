import os
import glob
import click
import pickle
import shutil
import ftplib
import sentry_sdk
from automator_settings import SENTRY_DSN
from amrsummary import before_send
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def alignfortablet_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    # Variable to hold supplied arguments
    argument_dict = {
        'fragmentlength': 500,
    }


    #parse description to figure out analysis type, and find fastas
    try:
        query_list = list()
        reference = list()
        compare = False
        # Go through description to figure out what our query is and what the reference is.
        for item in description:
            item = item.upper()
            if item == '':
                continue
            if 'FRAGMENTLENGTH' in item:
                argument_dict['fragmentlength'] = item.split('=')[1]
                continue
            if 'COMPARE' in item:
                compare = True
                continue
            if compare:
                query_list.append(item)
            else:
                if 'REFERENCE' not in item:
                    reference.append(item)
        
        # Create fastas folder
        fasta_folder = os.path.join(work_dir, 'output')
        os.makedirs(fasta_folder)

        # Retrieve our reference file. Error user if they selected anything but one reference and don't continue.
        if len(reference) != 1:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='ERROR: You must specify one reference strain, and you '
                                                'specified {} reference strains. Please create a new'
                                                ' issue and try again.'.format(len(reference)), status_id=4)
            return

        if reference[0].upper() != 'ATTACHED':
            # Extract our reference file to our working directory.
            retrieve_nas_files(seqids=reference,
                               outdir=fasta_folder,
                               filetype='fasta',
                               copyflag=True)
            # Check that the file was successfully extracted. If it wasn't boot the user.
            if len(glob.glob(os.path.join(fasta_folder, '*fasta'))) == 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Could not find the specified reference file.'
                                                    ' Please verify it is a correct SEQID, create a new '
                                                    'issue, and try again.', status_id=4)
                return

        # If user specified attachment as the reference file, download it to our working directory.
        else:
            # Get the attachment ID, and download if it isn't equal to zero (meaning no attachment, so boot user with
            # appropriate error message)
            attachment = redmine_instance.issue.get(issue.id, include='attachments')
            attachment_id = 0
            for item in attachment.attachments:
                attachment_id = item.id

            # Download if we found an attachment, and use as our reference. Otherwise, exit and tell user to try again
            if attachment_id != 0:
                attachment = redmine_instance.attachment.get(attachment_id)
                attachment.download(savepath=work_dir, filename='reference.fasta')
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: You specified that the reference would be in attached file,'
                                                    ' but no attached file was found. Please create a new issue and '
                                                    'try again.',
                                              status_id=4)
                return


        # Create fastq folder
        fastq_folder = os.path.join(work_dir, 'fastqs')
        os.makedirs(fastq_folder)
   	#now extract our query file
        # Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
        retrieve_nas_files(seqids=query_list,
                           outdir=fastq_folder,
                           filetype='fastq',
                           copyflag=False)
        # With our query files extracted, verify that all the SEQIDs the user specified were able to be found.
        missing_fastqs = verify_fastqs_present(query_list, os.path.join(work_dir, 'fastqs'))
        if len(missing_fastqs) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Warning! Could not find the following requested query SEQIDs: '
                                                '{}. \nYou may want to verify the SEQIDs, create a new issue, and try'
                                                ' again.'.format(str(missing_fastqs)))

        #first create the bowtie2 index using our reference
        for reffile in glob.glob(os.path.join(fasta_folder, '*.fasta')):
            refid = os.path.split(reffile)[1].split('.')[0]
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dev/cowbat'
        for rawread in glob.glob(os.path.join(fastq_folder, '*_R1_001.fastq.gz')):
            seqid1 = os.path.split(rawread)[1].split('.')[0]
            seqid = os.path.split(seqid1)[1].split('_R')[0]
            #seqids
            forward = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
            reverse = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)

            # Run bowtie2 index with the necessary arguments
            bowtieindex_cmd = 'bowtie2-build {ref}.fasta indexedreference'.format(ref=refid)
            bowtie2_cmd = 'bowtie2 -X {frag} -x indexedreference -1 {read1} -2 {read2} -S paired_end_alignment.sam'\
                .format(frag=argument_dict['fragmentlength'],
                        read1=os.path.join(fastq_folder, forward), read2=os.path.join(fastq_folder, reverse))

            #create another shell script to execute within the dRep conda environment
            template = "#!/bin/bash\n{} && cd {} && {} && {}".format(activate, fasta_folder, bowtieindex_cmd, bowtie2_cmd)
            bowtie_script = os.path.join(work_dir, 'run_bowtie.sh')
            with open(bowtie_script, 'w+') as file:
                file.write(template)

            # Modify the permissions of the script to allow it to be run on the node
            make_executable(bowtie_script)
            # Run shell script
            os.system(bowtie_script)

        #now create a bam file of alignment
        samview_cmd = 'samtools view -b -S paired_end_alignment.sam -t {ref}.fasta > your_alignment.bam'.format(ref=reference)
        samsort_cmd = 'samtools sort your_alignment.bam -o your_alignment.sorted.bam'
        samindex_cmd = 'samtools index your_alignment.sorted.bam'

        #create another shell script to execute within the dRep conda environment
        template2 = "#!/bin/bash\n{} && cd {} && {} && {} && {}".format(activate, fasta_folder,samview_cmd, samsort_cmd, samindex_cmd)
        samtools_script = os.path.join(work_dir, 'run_samtools.sh')
        with open(samtools_script, 'w+') as file:
            file.write(template2)

        # Modify the permissions of the script to allow it to be run on the node
        make_executable(samtools_script)
        # Run shell script
        os.system(samtools_script)

        # Zip output
        alginout_filename = 'bowtiesamtools_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=fasta_folder,
                                  output_dir=work_dir,
                                  output_filename=alginout_filename)
        zip_filepath += '.zip'
        
                # Upload the zip file to Dropbox
        download_link = upload_to_dropbox(
            access_token=DROPBOX_ACCESS_TOKEN,
            refresh_token=DROPBOX_REFRESH_TOKEN,
            app_key=DROPBOX_APP_KEY,
            app_secret=DROPBOX_APP_SECRET,
            local_file_path=zip_filepath
        )
        # Prepare upload
        if download_link:
            redmine_instance.issue.update(
                resource_id=issue.id, status_id=4,
                notes='Bowtie2 and samtools analysis complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=download_link)
        else:
            redmine_instance.issue.update(
                resource_id=issue.id, status_id=4,
                notes='Upload of result files was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )



        # Create a list of all the folders - will be used to clean up the working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)

        #os.remove(os.path.join(work_dir, str(issue.id) + '.zip')) #delete zip file
        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes='{at} alignment with bowtie2 and samtools complete!'
                                      .format(at=argument_dict['analysis'].lower()))
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))

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


def verify_fastqs_present(query_list, fastq_folder):
    missing_fastqs = list()
    for query in query_list:
        # Check that forward reads are present.
        if len(glob.glob(os.path.join(fastq_folder, query + '*R1*fastq*'))) == 0:
            missing_fastqs.append(query)
        # Check that reverse reads are present, and add to list if forward reads weren't missing
        if len(glob.glob(os.path.join(fastq_folder, query + '*R2*fastq*'))) == 0 and query not in missing_fastqs:
            missing_fastqs.append(query)
    # Returns list of SEQIDs for which we couldn't find forward and/or reverse reads
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
    alignfortablet_redmine()
