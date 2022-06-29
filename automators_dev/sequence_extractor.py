import os
import glob
import click
import pickle
import shutil
from biotools import mash
from automator_settings import COWBAT_DATABASES
from nastools.nastools import retrieve_nas_files


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def geneseekr_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Set and create the directory to store the custom targets
    details_dir = os.path.join(work_dir, 'fasta_details')
    details_file = os.path.join(details_dir, 'fasta_details.txt')
    sequence_dir = os.path.join(work_dir, 'sequences')
    try:
        # Download the attached file.
        # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
        # out a better way to do it.
        attachment = redmine_instance.issue.get(issue.id, include='attachments')
        attachment_id = 0
        for item in attachment.attachments:
            attachment_id = item.id
        seqids = list()
        # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
        if attachment_id != 0:
            attachment = redmine_instance.attachment.get(attachment_id)
            attachment.download(savepath=details_dir, filename='fasta_details.txt')
            # Parse the details file
            with open(details_file, 'r') as details:
                for line in details:
                    seqid, contig, start, stop = line.split(';')
                    seqids.append(seqid)
                    redmine_instance.issue.update(resource_id=issue.id,
                                                  notes='Extracted the following SEQIDs from the supplied list:\n'
                                                        '{seqids}'.format(seqids='\n'.join(seqids)))
        else:
            for item in description:
                seqid, contig, start, stop = item.split(';')
                seqid = seqid.upper()
                with open(details_file, 'w') as details:
                    details.write(f'{seqid};{contig};{start}{stop}')
                seqids.append(seqid)
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Created sequence extraction details file from provided details')

        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No SEQIDs provided!',
                                          status_id=4)
            return
        # Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
        retrieve_nas_files(seqids=seqids,
                           outdir=sequence_dir,
                           filetype='fasta',
                           copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, sequence_dir)
        # Update the Redmine issue if one or more of the requested SEQIDs could not be located
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))

        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/sequence_extractor'
        py = '/mnt/nas2/virtual_environments/sequence_extractor/lib/python3.9/site-packages/genemethods/SequenceExtractor/src/sequenceExtractor.py'
        # Run the command the necessary arguments
        cmd = 'python {py} -s {seqpath} -f {fasta_details}'\
            .format(py=py,
                    fasta_details=details_file,
                    seqpath=sequence_dir)
        # Update the issue with the command
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='SequenceExtractor command:\n {cmd}'.format(cmd=cmd))
        # Create another shell script to execute within the PlasmidExtractor conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, cmd)
        script = os.path.join(work_dir, 'run_sequence_extractor.sh')
        with open(script, 'w+') as file:
            file.write(template)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(script)
        # Run shell script
        os.system(script)

        # Zip output
        output_filename = 'sequence_extractor_output'
        zip_filepath = zip_folder(results_path=os.path.join(sequence_dir, 'output'),
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

        # Create a list of all the folders - will be used to clean up the working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)
        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes='Sequence extraction complete!')
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
    geneseekr_redmine()
