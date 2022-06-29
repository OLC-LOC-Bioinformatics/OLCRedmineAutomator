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
def diversitree_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        if not os.path.isdir(os.path.join(work_dir, 'fastas')):
            os.makedirs(os.path.join(work_dir, 'fastas'))
        # Check that the first line of the request is a number. If it isn't, tell author they goofed and give up.
        try:
            desired_num_strains = int(description[0])
        except ValueError:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Error! The first line of your request must be the number of'
                                                ' strains you want picked from the tree.',
                                          status_id=4)
            return

        # Parse description for SEQIDs, write list that file_extractor needs.
        seqids = list()
        for i in range(1, len(description)):
            item = description[i].upper()
            seqids.append(item)

        if 'TREEPROGRAM' in seqids[-1]:
            treemaker = seqids[-1].split('=')[1].lower()
            seqids.pop()
            if treemaker not in ['parsnp', 'mashtree']:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='Error! Available tree creation programs are mashtree and parsnp. '
                                                    'Your choice was {}'.format(treemaker),
                                              status_id=4)
                return
        else:
            treemaker = 'parsnp'

        # Drop FASTA files into workdir
        retrieve_nas_files(seqids=seqids,
                           outdir=os.path.join(work_dir, 'fastas'),
                           filetype='fasta',
                           copyflag=False)
        # Download the attached FASTA file.
        # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
        # out a better way to do it.
        attachment = redmine_instance.issue.get(issue.id, include='attachments')
        for item in attachment.attachments:
            # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
            if item.id != 0:
                attachment = redmine_instance.attachment.get(item.id)
                attachment.download(os.path.join(work_dir, 'fastas'), filename=item.filename)
        # Run a mash to figure out if any strains are particularly far apart and likely to make PARSNP fail.
        reference_file = glob.glob(os.path.join(work_dir, 'fastas', '*.fasta'))[0]
        bad_fastas = check_distances(reference_file, os.path.join(work_dir, 'fastas'))
        if bad_fastas:
            outstr = ''
            for fasta in bad_fastas:
                fasta = os.path.split(fasta)[-1]
                outstr += fasta + '\n'
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Warning! MASH screening thinks that the following samples may be too'
                                                ' far from the reference: {samples}\nIn this case, the reference file'
                                                ' was {reference}. You may want to create a new issue and '
                                                'try again.'.format(samples=outstr,
                                                                    reference=os.path.split(reference_file)[-1]))

        # Remove distances.tab and sketch.msh from fastas folder, because sometimes they make
        # parsnp crash. Other times they don't. I have no idea why, so remove just to be safe.
        try:
            os.remove(os.path.join(work_dir, 'fastas', 'distances.tab'))
            os.remove(os.path.join(work_dir, 'fastas', 'sketch.msh'))
        except OSError:
            pass

        # Full paths needed here since SLURM doesn't give the $PATH of the host machine to the script for some reason
        if treemaker == 'parsnp':
            cmd = '/mnt/nas2/virtual_environments/mob_suite/bin/parsnp -r ! -d {input} -o {output} -p {threads}'\
                .format(input=os.path.join(work_dir, 'fastas'),
                        output=os.path.join(work_dir, 'output'),
                        threads=24)
        elif treemaker == 'mashtree':
            if not os.path.isdir(os.path.join(work_dir, 'output')):
                os.makedirs(os.path.join(work_dir, 'output'))
            cmd = '/home/ubuntu/bin/mashtree --numcpus 24 --outtree {output_newick} {input_fastas}'\
                .format(output_newick=os.path.join(work_dir, 'output', 'parsnp.tree'),
                        input_fastas=os.path.join(work_dir, 'fastas', '*.fasta'))
        returncode = subprocess.call(cmd, shell=True, env={'PERL5LIB': '$PERL5LIB:/home/ubuntu/lib/perl5'})
        if returncode != 0:
            raise Exception('Tree creation command ({}) for {} had return code {}'.format(cmd, issue.id, returncode))
        # Now use diversitree to pick the strains we actually want.
        # IMPORTANT NOTE TO ANYONE MAINTAINING THIS: Need to have xvfb installed on nodes in order to make this run.
        # StrainChoosr uses ete3 to draw trees, which uses PyQt, which needs some sort of display.
        cmd = 'xvfb-run strainchoosr --treefile {tree} --number {number} ' \
              '--output_name {output}'.format(tree=os.path.join(work_dir, 'output', 'parsnp.tree'),
                                              number=desired_num_strains,
                                              output=os.path.join(work_dir, 'diversitree_report'))
        returncode = subprocess.call(cmd, shell=True)
        if returncode != 0:
            raise Exception('StrainChoosr command ({}) for {} had return code {}'.format(cmd, issue.id, returncode))
        output_list = list()
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'output', 'parsnp.tree')
        output_dict['filename'] = 'tree.nwk'
        output_list.append(output_dict)
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'diversitree_report.html')
        output_dict['filename'] = 'diversitree_report.html'
        output_list.append(output_dict)
        redmine_instance.issue.update(resource_id=issue.id, uploads=output_list, status_id=4,
                                      notes='DiversiTree process complete!')
        try:
            shutil.rmtree(os.path.join(work_dir, 'fastas'))
        except FileNotFoundError:
            pass
    except Exception as e:
        sentry_sdk.capture_exception(e)
        try:
            shutil.rmtree(os.path.join(work_dir, 'fastas'))
        except FileNotFoundError:
            pass
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))


def check_distances(ref_fasta, fasta_folder):
    bad_fastqs = list()
    # fastqs = glob.glob(os.path.join(fastq_folder, '*R1*'))
    mash.sketch(os.path.join(fasta_folder, '*.fasta'), output_sketch=os.path.join(fasta_folder, 'sketch.msh'), threads=56)
    mash.dist(os.path.join(fasta_folder, 'sketch.msh'), ref_fasta, threads=56, output_file=os.path.join(fasta_folder, 'distances.tab'))
    mash_output = mash.read_mash_output(os.path.join(fasta_folder, 'distances.tab'))
    for item in mash_output:
        print(item.reference, item.query, str(item.distance))
        if item.distance > 0.06:  # May need to adjust this value.
            bad_fastqs.append(item.reference)
    return bad_fastqs


if __name__ == '__main__':
    diversitree_redmine()
