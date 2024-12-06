import os
import re
import glob
import click
import pickle
import shutil
import sentry_sdk
from automator_settings import COWBAT_DATABASES, SENTRY_DSN
from nastools.nastools import retrieve_nas_files


def before_send(event, hint):
    """
    The cluster we send jobs to has its clock ahead by a bunch for reasons I don't remember.
    Because of this, whenver we retrieve any SeqIDs, Sentry freaks out because the timestamps attached
    are ahead of where they should be, and ends up sending error notifications each time we retrieve any sequence
    data. Adding this callback to the sentry init finds if the message is formatted like a SeqID, and if it is
    the error report won't get sent, as it isn't really an error.
    """
    message = event['logentry']['message']
    if re.match(r'\d{4}-[A-Z]+-\d{4}', message):
        return None
    else:
        return event


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def resfinder_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.upper()
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

        # Use the COWBAT_DATABASES variable as the database path
        db_path = COWBAT_DATABASES
        # Run ResFindr
        cmd = 'GeneSeekr blastn -s {seqfolder} -t {targetfolder} -r {reportdir} -A'\
            .format(seqfolder=work_dir,
                    targetfolder=os.path.join(db_path, 'resfinder'),
                    reportdir=os.path.join(work_dir, 'reports'))
        # Update the issue with the ResFinder command
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='ResFinder command:\n {cmd}'.format(cmd=cmd))
        os.system(cmd)
        # These unfortunate hard coded paths appear to be necessary
        #activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/cowbat'
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dev/cowbat'
        # Run sipprverse with the necessary arguments
        mob_cmd = 'python -m genemethods.assemblypipeline.mobrecon -s {seqfolder} -r {targetfolder}' \
            .format(seqfolder=work_dir,
                    targetfolder=os.path.join(db_path, 'mobrecon'))
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
        # Add the three reports separately to the output list
        output_dict['path'] = os.path.join(work_dir, 'reports', 'resfinder_blastn.xlsx')
        output_dict['filename'] = 'resfinder_blastn.xlsx'
        output_list.append(output_dict)
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'reports', 'mob_recon_summary.csv')
        output_dict['filename'] = 'mob_recon_summary.csv'
        output_list.append(output_dict)
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir, 'reports', 'amr_summary.csv')
        output_dict['filename'] = 'amr_summary.csv'
        output_list.append(output_dict)
        redmine_instance.issue.update(resource_id=issue.id, uploads=output_list, status_id=4,
                                      notes='resfinder process complete!')

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
    resfinder_redmine()
