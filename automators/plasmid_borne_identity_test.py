import os
import glob
import click
import pickle
import shutil
import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN
from automator_settings import COWBAT_DATABASES
from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def plasmid_borne_identity(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Variable to hold supplied arguments
    argument_dict = {
        'analysis': 'custom',
        'blast': 'blastn',
        'cutoff': 1,
        'evalue': '1E-5',
    }
    # Set the database path for the analyses
    database_path = {
        'custom': os.path.join(work_dir, 'targets')
    }
    # Current BLAST analyses supported
    blasts = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']
    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'CUTOFF' in item:
                argument_dict['cutoff'] = int(item.split('=')[1].lower())
                continue
            if 'EVALUE' in item:
                argument_dict['evalue'] = item.split('=')[1].lower()
                continue
            if 'BLAST' in item:
                argument_dict['blast'] = item.split('=')[1].lower()
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        retrieve_nas_files(seqids=seqids,
                           outdir=work_dir,
                           filetype='fasta',
                           copyflag=False)

        missing_fastas = verify_fasta_files_present(seqids, work_dir)
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))

        # Set and create the directory to store the custom targets
        target_dir = os.path.join(work_dir, 'targets')
        os.makedirs(target_dir, exist_ok=True)
        # Download the attached FASTA file.
        # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
        # out a better way to do it.
        attachment = redmine_instance.issue.get(issue.id, include='attachments')
        attachment_id = 0
        for item in attachment.attachments:
            attachment_id = item.id
        # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
        if attachment_id != 0:
            attachment = redmine_instance.attachment.get(attachment_id)
            attachment.download(savepath=target_dir, filename='targets.tfa')
        else:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='ERROR: Analysis type custom requires an attached FASTA file of '
                                                'targets. The automator could not find any attached files. '
                                                'Please create a new issue with the FASTA file attached and try '
                                                'again.',
                                          status_id=4)
            return
        # Ensure that the requested BLAST analysis is valid
        if argument_dict['blast'] not in blasts:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: requested BLAST analysis, {bt}, is not one of the currently '
                                                'supported analyses: {blasts}'.format(bt=argument_dict['blast'],
                                                                                      blasts=', '.join(blasts)),
                                          status_id=4)
            return
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/cowbat'
        seekr_py = '/mnt/nas2/virtual_environments/geneseekr/bin/GeneSeekr'
        # Run sipprverse with the necessary arguments
        seekr_cmd = 'python {seekr_py} {blast} -s {seqpath} -r {outpath} -t {dbpath} -c {cutoff} -e {evalue} -u' \
            .format(seekr_py=seekr_py,
                    blast=argument_dict['blast'],
                    seqpath=work_dir,
                    outpath=os.path.join(work_dir, 'reports'),
                    dbpath=database_path[argument_dict['analysis']],
                    cutoff=argument_dict['cutoff'],
                    evalue=argument_dict['evalue'])
        # Update the issue with the GeneSeekr command
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='GeneSeekr command:\n {cmd}'.format(cmd=seekr_cmd))
        # Create another shell script to execute within the PlasmidExtractor conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, seekr_cmd)
        geneseekr_script = os.path.join(work_dir, 'run_geneseekr.sh')
        with open(geneseekr_script, 'w+') as file:
            file.write(template)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(geneseekr_script)
        # Run shell script
        os.system(geneseekr_script)
        # Run MOB Recon
        # Use the COWBAT_DATABASES variable as the database path
        db_path = COWBAT_DATABASES
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/cowbat'
        # Run sipprverse with the necessary arguments
        mob_cmd = 'python -m genemethods.assemblypipeline.mobrecon -s {seqfolder} -r {targetfolder} -a geneseekr -b {blast}' \
            .format(seqfolder=work_dir,
                    targetfolder=os.path.join(db_path, 'mobrecon'),
                    blast=argument_dict['blast'])
        # Update the issue with the MOB Recon command
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='MOB Recon command:\n {cmd}'.format(cmd=mob_cmd))
        # Create another shell script to execute within the PlasmidExtractor conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, mob_cmd)
        mob_script = os.path.join(work_dir, 'run_mob_recon.sh')
        with open(mob_script, 'w+') as file:
            file.write(template)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(mob_script)
        # Run shell script
        os.system(mob_script)
        # Get the output file uploaded.
        output_list = list()
        output_dict = dict()
        # Add the reports separately to the output list
        # GeneSeekr Excel-formatted report
        output_dict['path'] = os.path.join(work_dir, 'reports', 'geneseekr_{blast}.xlsx'
                                           .format(blast=argument_dict['blast']))
        output_dict['filename'] = 'geneseekr_{blast}.xlsx'.format(blast=argument_dict['blast'])
        output_list.append(output_dict)
        # Detailed GeneSeekr report
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'reports', 'geneseekr_{blast}_detailed.csv'
                                           .format(blast=argument_dict['blast']))
        output_dict['filename'] = 'geneseekr_{blast}_detailed.csv'.format(blast=argument_dict['blast'])
        output_list.append(output_dict)
        # Simple GeneSeekr report
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'reports', 'geneseekr_{blast}.csv'
                                           .format(blast=argument_dict['blast']))
        output_dict['filename'] = 'geneseekr_{blast}.csv'.format(blast=argument_dict['blast'])
        output_list.append(output_dict)
        # MOB Recon summary report
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'reports', 'mob_recon_summary.csv')
        output_dict['filename'] = 'mob_recon_summary.csv'
        output_list.append(output_dict)
        # Plasmid-borne summary report
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'reports', 'plasmid_borne_summary.csv')
        output_dict['filename'] = 'plasmid_borne_summary.csv'
        output_list.append(output_dict)
        redmine_instance.issue.update(resource_id=issue.id, uploads=output_list, status_id=4,
                                      notes='PlasmidBorne identity process complete!')

        # Clean up all FASTA/FASTQ files so we don't take up too much space on the NAS
        os.system('rm {workdir}/*fasta'.format(workdir=work_dir))
        try:
            # Remove all other folders
            for dirpath, dirnames, filenames in os.walk(work_dir):
                for dirname in dirnames:
                    shutil.rmtree(os.path.join(dirpath, dirname))
        except IOError:
            pass
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')


def verify_fasta_files_present(seqid_list, fasta_dir):
    missing_fastas = list()
    for seqid in seqid_list:
        if len(glob.glob(os.path.join(fasta_dir, seqid + '*.fasta'))) == 0:
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


if __name__ == '__main__':
    plasmid_borne_identity()
