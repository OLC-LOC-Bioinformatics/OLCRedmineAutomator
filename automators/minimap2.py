import os
import glob
import click
import pickle
import shutil
import pandas
from nastools.nastools import retrieve_nas_files
from externalretrieve import upload_to_ftp


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def minimap2_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        # Variable to hold supplied arguments
        argument_dict = {
            'k-mer': 15
        }

        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        analysistype = "" # Just for the output at the end
        for item in description:
            item = item.rstrip()
            if 'k-mer' in item:
                argument_dict['k-mer'] = int(item.split('=')[1].lower())
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # No SEQIDs
        if len(seqids) == 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You did not include any SEQIDs.')
            return

        # Create folder to drop FASTQ files
        seq_folder = os.path.join(work_dir, 'sequences')
        os.mkdir(seq_folder)

        # Create output folder which will be zipped later
        output_folder = os.path.join(work_dir, 'output')
        os.makedirs(output_folder)

        # Extract FASTQ files.
        retrieve_nas_files(seqids=seqids, outdir=seq_folder, filetype='fastq', copyflag=False)
        missing_fastqs = check_fastqs_present(seqids, seq_folder)
        if missing_fastqs:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on '
                                                'the OLC NAS: {}'.format(missing_fastqs))

        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/minimap2_sarah'
        targets_file = '/mnt/nas2/virtual_environments/minimap2_sarah/serotype_combined_targets.fasta'

        #now for each fastq file in the list, going to run the minimap2 command using the chosen targets file
        for fastq in glob.glob(os.path.join(seq_folder, '*.fastq.gz')):
            seqid = os.path.split(fastq)[1].split('.')[0]
            # Prepare command
            cmd = 'minimap2 -ax map-ont {target} {fastq} -o {output_file}'\
                .format(output_file='{seqid}-aln.sam'.format(seqid=seqid),
                        fastq=fastq,
                        seqid=seqid,
                        target=targets_file)

            # Create another shell script to execute within the roary conda environment
            template = "#!/bin/bash\n{} && cd {} && {}".format(activate, seq_folder, cmd) #this line takes the activate  command (above) and the command for minimap2 and writes it into bash file
            minimap2_script = os.path.join(work_dir, 'run_minimap2.sh')
            with open(minimap2_script, 'w+') as file:
                file.write(template)
            make_executable(minimap2_script)

            # Run shell script
            os.system(minimap2_script)

            #move the alignment files to the output folder
            sam_file = '{seqid}-aln.sam'.format(seqid=seqid)
            output_sam = os.path.join(output_folder, sam_file)
            shutil.copyfile(os.path.join(seq_folder,sam_file), output_sam)

        # Zip output
        output_filename = 'minimap2_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=output_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'

        upload_successful = upload_to_ftp(local_file=zip_filepath)
        # Prepare upload
        if upload_successful:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Minimap2 process complete!\n\n'
                                                'Results are available at the following FTP address:\n'
                                                'ftp://ftp.agr.gc.ca/outgoing/cfia-ac/{}'
                                          .format(os.path.split(zip_filepath)[1]))
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity '
                                                'issues. '
                                                'Please try again later.')
        # Remove the zip file
#        os.remove(zip_filepath)

        # Clean up files
        shutil.rmtree(output_folder)
#        os.remove(zip_filepath)
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
    missing_fastas = list()
    for seqid in seqid_list:
        if len(glob.glob(os.path.join(fasta_dir, seqid + '*.fasta'))) == 0:
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
    minimap2_redmine()
