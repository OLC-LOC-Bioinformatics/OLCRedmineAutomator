#!/usr/bin/env python

import os
import glob
import click
import pickle
import sentry_sdk
import subprocess
import shutil
from biotools import mash
from amrsummary import before_send
from strainchoosr import strainchoosr
from automator_settings import SENTRY_DSN

from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def mashtree_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))
        #analyses supported by the automator
        analyses = [
            'custom', 'enterobacterales', 'listeriaceae', 'enterobacter'
        ]

        # Variable to hold supplied arguments
        argument_dict = {
            'analysis': str(),
            'genomesize': 500000,
            'mindepth': 5,
            'kmerlength': 21,
            'sketch-size': 10000,
        }

        # Parse description for SEQIDs, write list that file_extractor needs.
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'ANALYSIS' in item:
               argument_dict['analysis'] = item.split('=')[1].lower()
               continue
            if 'GENOMESIZE' in item:
               argument_dict['genomesize'] = item.split('=')[1]
               continue
            if 'MINDEPTH' in item:
               argument_dict['mindepth'] = item.split('=')[1]
               continue
            if 'SKETCH-SIZE' in item:
               argument_dict['sketch-size'] = item.split('=')[1]
               continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that the analysis type is provided
        if not argument_dict['analysis']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not identify an analysis type. '
                                                'Please ensure that the first line of the issue contains one'
                                                ' of the following keywords: {ats}'.format(ats=', '.join(analyses)),
                                          status_id=4)
            return
        elif argument_dict['analysis'] not in analyses:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied analysis type {at} current not in the supported '
                                                'list of analyses: {ats}'.format(at=argument_dict['analysis'],
                                                                                 ats=', '.join(analyses)),
                                          status_id=4)
            return
    # link to order enterobacterales fasta files if analysis = enterobacterales
        if argument_dict['analysis'] == 'enterobacterales':
            src = '/mnt/nas2/processed_sequence_data/ncbi/enterobacterales_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

    #if listeriaceae, link to listeria reference sequences
        if argument_dict['analysis'] == 'listeriaceae':
            src = '/mnt/nas2/processed_sequence_data/ncbi/listeriaceae_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

    #if enterobacter, link to enterobacter reference sequences
        if argument_dict['analysis'] == 'enterobacter':
            src = '/mnt/nas2/processed_sequence_data/ncbi/enterobacter_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)
            src2 = '/mnt/nas2/processed_sequence_data/ncbi/enterobacterales_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd2 = 'ln -s {src}ENTEROBACTER_*.fasta {dst}'.format(src=src2,
                                                                    dst=dst)
            os.system(lncmd2)

    #if just a custom analysis, create the fastas folder
        if argument_dict['analysis'] == 'custom':
            # Make fasta file directory
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
    

        # Drop FASTA files into workdir
        retrieve_nas_files(seqids=seqids,
                           outdir=os.path.join(work_dir, 'fastas'),
                           filetype='fasta',
                           copyflag=False)
	# Verify that specified fasta files are actually there, warn user if they aren't.
        missing_fastas = verify_fasta_files_present(seqid_list=seqids,
                                                    fasta_dir=os.path.join(work_dir, 'fastas'))
        if len(missing_fastas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))

        # Download the attached FASTA file.
        # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
        # out a better way to do it.
        attachment = redmine_instance.issue.get(issue.id, include='attachments')
        for item in attachment.attachments:
            # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
            if item.id != 0:
                attachment = redmine_instance.attachment.get(item.id)
                attachment.download(os.path.join(work_dir, 'fastas'), filename=item.filename)

        # Full paths needed here since SLURM doesn't give the $PATH of the host machine to the script for some reason
        if not os.path.isdir(os.path.join(work_dir, 'output')):
            os.makedirs(os.path.join(work_dir, 'output'))
        cmd = '/home/ubuntu/bin/mashtree --numcpus 24 --outtree {output_newick} --genomesize {g}'\
            ' --mindepth {depth} --kmerlength {k} --sketch-size {s} {input_fastas}'\
            .format(output_newick=os.path.join(work_dir, 'output', 'mashtree.tree'),
                    g=argument_dict['genomesize'],
                    depth=argument_dict['mindepth'],
                    k=argument_dict['kmerlength'],
                    s=argument_dict['sketch-size'],
                    input_fastas=os.path.join(work_dir, 'fastas', '*.fasta'))

        returncode = subprocess.call(cmd, shell=True, env={'PERL5LIB': '$PERL5LIB:/home/ubuntu/lib/perl5'})
        if returncode != 0:
            raise Exception('Tree creation command ({}) for {} had return code {}'.format(cmd, issue.id, returncode))

        #zip up the output file(s)
        output_filename = 'mashtree_output'
        zip_filepath = zip_folder(results_path=os.path.join(work_dir, 'output'),
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

        try:
            shutil.rmtree(os.path.join(work_dir, 'fastas'))
        except FileNotFoundError:
            pass

        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes='{at} analysis with mashtree complete!'
                                      .format(at=argument_dict['analysis'].lower()))
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')


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
    mashtree_redmine()
