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
import pandas




@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def seqsero2_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Current list of workflow types that SeqSero2 can use
    #analyses = [
    #    'allele',
    #]
    # Current list of sequence types that SeqSero2 can analyse
    seqtypes = [
        'assembly', 'rawpaired', 'minionfastq'
    ]

    # Variable to hold supplied arguments
    argument_dict = {
        #'analysis': str(),
        'seqtype': 'rawpaired',
        'nanopore': False,
    }


    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            #if 'ANALYSIS' in item:
            #    argument_dict['analysis'] = item.split('=')[1].lower()
            #    continue
            if 'SEQTYPE' in item:
                argument_dict['seqtype'] = item.split('=')[1].lower()
                continue
            if 'NANOPORE' in item:
                argument_dict['nanopore'] = True
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that the analysis type is provided
        #if not argument_dict['analysis']:
        #    redmine_instance.issue.update(resource_id=issue.id,
        #                                  notes='WARNING: Could not identify an analysis type. '
        #                                        'Please ensure that the first line of the issue contains one'
        #                                        ' of the following keywords: {ats}'.format(ats=', '.join(analyses)),
        #                                  status_id=4)
        #    return
        #elif argument_dict['analysis'] not in analyses:
        #    redmine_instance.issue.update(resource_id=issue.id,
        #                                  notes='WARNING: supplied analysis type {at} currently not in the supported '
        #                                        'list of analyses: {ats}'.format(at=argument_dict['analysis'],
        #                                                                         ats=', '.join(analyses)),
        #                                  status_id=4)
        #    return

        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No SEQIDs provided!',
                                          status_id=4)
            return

        # Set the database path for the analyses

        #create a folder to hold the sequences
        seq_dir = os.path.join(work_dir, 'sequences')
        os.makedirs(seq_dir, exist_ok=True)

        #now create an output folder so we can store all of the outputs
        out_dir = os.path.join(work_dir, 'output')
        os.makedirs(out_dir, exist_ok=True)

        #pull assemblies from nas
        if argument_dict['seqtype'] == 'assembly':
        # Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
            retrieve_nas_files(seqids=seqids,
                               outdir=seq_dir,
                               filetype='fasta',
                               copyflag=False)
            missing_fastas = verify_fasta_files_present(seqids, seq_dir)
            # Update the Redmine issue if one or more of the requested SEQIDs could not be located
            if missing_fastas:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastas))

            # These unfortunate hard coded paths appear to be necessary
            activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/seqsero2'

            # Run seqsero with the necessary arguments
            #assembly files
            for assembly in glob.glob(os.path.join(seq_dir, '*.fasta')):
                seqid = os.path.split(assembly)[1].split('.')[0]
                seqseroout = 'SeqSero_result_{}'.format(seqid)
                #prepare command for seqsero
                seqsero_cmd = 'SeqSero2_package.py -t 4 -i {assembly} -d {out} -m k -p 8'\
                    .format(assembly=assembly,
                            out=os.path.join(out_dir,seqseroout))

                # Create another shell script to execute within the KMA conda environment
                template = "#!/bin/bash\n{} && cd {} && {}".format(activate, seq_dir, seqsero_cmd)
                seqsero_script = os.path.join(work_dir, 'run_seqsero.sh')
                with open(seqsero_script, 'w+') as file:
                    file.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(seqsero_script)
                # Run shell script
                os.system(seqsero_script)

        #pull paired-end raw reads from nas
        if argument_dict['seqtype'] == 'rawpaired':
            # Extract FASTQ files.
            retrieve_nas_files(seqids=seqids, outdir=seq_dir, filetype='fastq', copyflag=False)
            missing_fastqs = check_fastqs_present(seqids, seq_dir)
            if len(missing_fastqs) > 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastqs))

            # These unfortunate hard coded paths appear to be necessary
            activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/seqsero2'

            #run analysis for each pair of fastq files
            for rawread in glob.glob(os.path.join(seq_dir, '*_R1_001.fastq.gz')):
                seqid1 = os.path.split(rawread)[1].split('.')[0]
                seqid = os.path.split(seqid1)[1].split('_S')[0]
                #seqidforfile = os.path.split(seqid1)[1].split('_S')[0]
                #forwardseq = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
                #reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
                #prepare command for SeqSero2
                seqseroout = 'SeqSero_result_{}'.format(seqid)
                #prepare command for seqsero
                seqsero_cmd = 'SeqSero2_package.py -t 2 -i {seqid}* -d {out} -m a -p 12'\
                    .format(seqid=seqid,
                            out=os.path.join(out_dir,seqseroout))

                # Create another shell script to execute within the KMA conda environment
                template = "#!/bin/bash\n{} && cd {} && {}".format(activate, seq_dir, seqsero_cmd)
                seqsero_script = os.path.join(work_dir, 'run_seqsero.sh')
                with open(seqsero_script, 'w+') as file:
                    file.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(seqsero_script)
                # Run shell script
                os.system(seqsero_script)


        #now concatenate the result files, which are named SeqSero_result.txt in the subfolders, to one file in the out_dir directory...
        resultsfiles = glob.glob("{}/*/SeqSero_result.txt".format(out_dir))
        #open a new outfile that will be the concatenation of the results files
        ssoutputfile = os.path.join(out_dir, 'SeqSero_results_all.txt')
        with open(ssoutputfile, 'w') as outfile:
            for f in resultsfiles:
                with open(f) as infile:
                    outfile.write(infile.read())
                outfile.write("\n")
        df = pandas.read_fwf(ssoutputfile)
        resultfilename = 'SeqSero_results_{}.csv'.format(issue.id)
        resultfileloc = os.path.join(out_dir,resultfilename)
        df.to_csv(resultfileloc, sep="\t", index=False)

        #upload the SeqSero results file to the redmine request
        output_list = list()
        output_dict = dict()
        output_dict['path'] = os.path.join(out_dir,resultfilename)
        output_dict['filename'] = resultfilename
        output_list.append(output_dict)
           
        # Zip output
        seqseroout_filename = 'SeqSero2_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=out_dir,
                                  output_dir=work_dir,
                                  output_filename=seqseroout_filename)
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
                uploads-output_list,
                notes='SeqSero analysis complete!\n\n'
                      'Results are attached, or full results can be found\n'
                      'at the following URL:\n'
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
        os.remove(zip_filepath)

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Please contact a bioinformatician '
                                            'to investigate: {}'.format(e))




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
    seqsero2_redmine()
