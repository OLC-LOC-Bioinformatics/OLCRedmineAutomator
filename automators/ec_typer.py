import os
import glob
import click
import pickle
import shutil
import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN
from nastools.nastools import retrieve_nas_files


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def ec_typer_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        # Parse description to get list of SeqIDs
        seqids = []
        for i in range(0, len(description)):
            item = description[i]
            item = item.upper()

            # Minimal check to make sure IDs provided somewhat resemble a valid sample ID
            seqids.append(item)

        # Create folder to drop FASTQ files
        assemblies_folder = os.path.join(work_dir, 'assemblies')
        os.makedirs(assemblies_folder, exist_ok=True)

        # Create output folder
        output_folder = os.path.join(work_dir, 'output')
        os.makedirs(output_folder)

        # Extract FASTQ files.
        retrieve_nas_files(seqids=seqids, outdir=assemblies_folder, filetype='fasta', copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, assemblies_folder)
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on '
                                                'the OLC NAS: {}'.format(missing_fastas))

        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/cowbat'
        ectyper = '/mnt/nas2/virtual_environments/cowbat/bin/ectyper'

        # Prepare command
        cmd = '{ectyper} -i {input_folder} -o {output_folder}'.format(ectyper=ectyper,
                                                                      input_folder=assemblies_folder,
                                                                      output_folder=output_folder)

        # Create another shell script to execute within the PlasmidExtractor conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, cmd)
        ec_script = os.path.join(work_dir, 'run_ec_typer.sh')
        with open(ec_script, 'w+') as file:
            file.write(template)
        make_executable(ec_script)

        # Run shell script
        os.system(ec_script)

        # Get the output file uploaded.
        output_list = list()
        output_dict = dict()
        # Add the reports separately to the output list
        # GeneSeekr Excel-formatted report
        output_dict['path'] = os.path.join(output_folder, 'output.tsv')
        output_dict['filename'] = 'ec_typer_report.tsv'
        output_list.append(output_dict)
        redmine_instance.issue.update(resource_id=issue.id, uploads=output_list, status_id=4,
                                      notes='ECTyper analyses complete!')
        # Clean up files
        shutil.rmtree(output_folder)
        shutil.rmtree(assemblies_folder)
    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')


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


if __name__ == '__main__':
    ec_typer_redmine()
