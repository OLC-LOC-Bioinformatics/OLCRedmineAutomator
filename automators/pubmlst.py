import os
import glob
import click
import pickle
import shutil
import subprocess
from biotools import mash
from nastools.nastools import retrieve_nas_files
from externalretrieve import upload_to_ftp
import csv


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def pubmlst_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.rstrip()
            # Otherwise the item should be a SEQID
            seqids.append(item)

        
        # Ensure that SEQIDs were included
        if len(seqids) == 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You did not include any SEQIDs.')
            return

        # Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
        retrieve_nas_files(seqids=seqids,
                           outdir=work_dir,
                           filetype='fasta',
                           copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, work_dir)
        # Update the Redmine issue if one or more of the requested SEQIDs could not be located
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))

        #activate environment
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/pubmlst'
        pubmlst_py = '/mnt/nas2/virtual_environments/pubmlst/species_identification.py'

        assembly = glob.glob(os.path.join(work_dir, '*.fasta'))[0]
        seqid = os.path.split(assembly)[1].split('.')[0]
        pubmlst_cmd = 'python {py} -f {assembly} > {work_dir}/output.txt'.format(py=pubmlst_py, assembly=assembly, work_dir=work_dir)

        #create another shell script to execute within the dRep conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, pubmlst_cmd)
        pubmlst_script = os.path.join(work_dir, 'run_pubmlst.sh')
        with open(pubmlst_script, 'w+') as file:
            file.write(template)

        # Modify the permissions of the script to allow it to be run on the node
        make_executable(pubmlst_script)
        #run command
        os.system(pubmlst_script)

        #results
        rank_line = [0]
        taxon_line = [1]
        support_line = [2]
        taxonomy_line = [3]
        o_file = open(os.path.join(work_dir, 'output.txt'))
        for index, line in enumerate(o_file):
            if index in rank_line:
                Rank = line.split(':')[1]
            if index in taxon_line:
                Taxon = line.split(':')[1]
            if index in support_line:
                Support = line.split(':')[1]
            if index in taxonomy_line:
                Taxonomy = line.split('axonomy')[1]

        # Update the issue with the pubmlst output
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Process complete! Here is the closest match from pubMLST:\n \n'
                                            ' Rank: {}\n'
                                            ' Taxon: {}\n'
                                            ' Support: {}\n'
                                            ' Taxonomy: {}'.format(Rank,Taxon,Support,Taxonomy),
                                      status_id=4)


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
    pubmlst_redmine()
