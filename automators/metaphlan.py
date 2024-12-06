import os
import glob
import click
import pickle
import shutil
import subprocess
from biotools import mash
from nastools.nastools import retrieve_nas_files
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
import csv


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def metaphlan_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        #sequence types supported by the automator
        seqtypes = [
            'fastq', 'fasta', 'nanoporefastq',
        ]

        #analyses supported by the automator
        analyses = [
            'rel_ab', 'rel_ab_w_read_stats', 'reads_map', 'clade_profiles','marker_ab_table', 'marker_counts','marker_pres_table',
        ]

        # Variable to hold supplied arguments
        argument_dict = {
            'seqtype': 'fastq',
            'analysis': 'rel_ab_w_read_stats',
        }


        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        analysistype = "" # Just for the output at the end
        for item in description:
            item = item.rstrip()
            if 'seqtype' in item:
                argument_dict['seqtype'] = item.split('=')[1].lower()
                continue
            if 'analysis' in item:
                argument_dict['analysis'] = item.split('=')[1].lower()
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        
        # Ensure that SEQIDs were included
        if len(seqids) == 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You did not include any SEQIDs.')
            return

        # Create folder to drop sequence files
        sequences_folder = os.path.join(work_dir, 'sequences')
        os.makedirs(sequences_folder, exist_ok=True)

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
#            attachment.download(savepath=target_dir, filename='traits.tsv')
#                attachment.download(savepath=work_dir, filename=attachment_id)
                attachment.download(savepath=sequences_folder)


        # Create metaphlan4 output folder
        metaphlan4_folder = os.path.join(work_dir, 'metaphlan4_output')
        os.makedirs(metaphlan4_folder)


        #IF fastqs requested, check that they are all found on the nas
        if argument_dict['seqtype'] == 'fastq':
            #pull out the fastq files
            # Extract FASTQ files.
            retrieve_nas_files(seqids=seqids, outdir=sequences_folder, filetype='fastq', copyflag=False)
            missing_fastqs = check_fastqs_present(seqids, sequences_folder)
            if len(missing_fastqs) > 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastqs))

        if argument_dict['seqtype'] == 'nanoporefastq':
            #pull out the fastq files
            # Extract FASTQ files.
            retrieve_nas_files(seqids=seqids, outdir=sequences_folder, filetype='fastq', copyflag=False)
            missing_fastqs = check_fastqs_present(seqids, sequences_folder)
            if len(missing_fastqs) > 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastqs))
        if argument_dict['seqtype'] == 'fasta':
            # Extract FASTA files.
            retrieve_nas_files(seqids=seqids, outdir=sequences_folder, filetype='fasta', copyflag=False)
            missing_fastas = verify_fasta_files_present(seqids, seq_dir)
            # Update the Redmine issue if one or more of the requested SEQIDs could not be located
            if missing_fastas:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastas))

        #set the database path for the ChocoPhlan database (used by metaphlan)
        dbpath = '/mnt/nas2/databases/metaphlan4/'

        #Commands required for metaphlan
        activateemetaphlan = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/metaphlan4'

        #if fastq is chosen, run analysis using the metaphlan4 command
        if argument_dict['seqtype'] == 'fastq':
            #now run the analysis
            #do not want to repeat twice, so added in the _R1
            for metagenome in glob.glob(os.path.join(sequences_folder, '*_R1_001.fastq.gz')):
                seqid1 = os.path.split(metagenome)[1].split('.')[0]
                seqid = os.path.split(seqid1)[1].split('_R')[0]
                forwardseq = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
                reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
                bt2out = '{seqid}.bowtie2.bz2'.format(seqid=seqid)
                outputfile = '{seqid}_metaphlan4_report'.format(seqid=seqid)
                #prepare command
                metaphlan4_cmd = 'metaphlan --input_type fastq --bowtie2db {db} --nproc 24 --bowtie2out {bt2o} -o {out}'\
                                 ' -s {seq} -t {analysis} --sample_id {redmid}'\
                                 ' {forwardseq},{reverseseq}'.format(db=dbpath,
                                                                     bt2o=bt2out,out=outputfile,
                                                                     seq=seqid,analysis=argument_dict['analysis'], redmid=issue.id,
                                                                     forwardseq=forwardseq, reverseseq=reverseseq)

                # Create another shell script to execute within the pyseer conda environment
                template = "#!/bin/bash\n{ntasks}\n{mem}\n{activate} && cd {seqdir} && {mpa}"\
                    .format(ntasks="#SBATCH --ntasks 30",mem="#SBATCH --mem=190000",
                            activate=activateemetaphlan,seqdir=sequences_folder,
                            mpa=metaphlan4_cmd)
                mpa4_script = os.path.join(work_dir, 'run_metaphlan4.sh')
                with open(mpa4_script, 'w+') as file:
                    file.write(template)
                make_executable(mpa4_script)

                # Run shell script
                os.system(mpa4_script)

            #now move the files
            output_files = os.listdir(sequences_folder)
            for file in output_files:
                if file.endswith("_report"):
                    shutil.move(os.path.join(sequences_folder,file), os.path.join(metaphlan4_folder, file))

        # Zip metaphlan output
        output_filename = 'metaphlan4_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=metaphlan4_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
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
                status_id=2,
                notes='MetaPhlan4 analysis complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=download_link)
            )
        else:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=2,
                notes='Upload of results was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )
        # Remove the zip file
        #os.remove(zip_filepath)
        #remove the sequences folder
        #shutil.rmtree(sequences_folder)
        #remove the output folder
        #shutil.rmtree(kraken2_folder)

    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))


def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)

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


def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    metaphlan_redmine()
