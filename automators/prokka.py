import os
import glob
import click
import pickle
import shutil
import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN
from nastools.nastools import retrieve_nas_files
from externalretrieve import upload_to_ftp


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def prokka_redmine(redmine_instance, issue, work_dir, description):
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
        os.mkdir(assemblies_folder)

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
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/prokka'
        prokka = '/mnt/nas2/virtual_environments/prokka/bin/prokka'

        for assembly in glob.glob(os.path.join(assemblies_folder, '*.fasta')):
            seqid = os.path.split(assembly)[1].split('.')[0]
            # Prepare command
            cmd = '{prokka} --outdir {output_folder} --prefix {seqid} {assembly}'.format(prokka=prokka,
                                                                                         output_folder=os.path.join(output_folder, seqid),
                                                                                         seqid=seqid,
                                                                                         assembly=assembly)

            # Create another shell script to execute within the PlasmidExtractor conda environment
            template = "#!/bin/bash\n{} && {}".format(activate, cmd)
            prokka_script = os.path.join(work_dir, 'run_prokka.sh')
            with open(prokka_script, 'w+') as file:
                file.write(template)
            make_executable(prokka_script)

            # Run shell script
            os.system(prokka_script)

        # Zip output
        output_filename = 'prokka_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=output_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'

        upload_successful = upload_to_ftp(local_file=zip_filepath)
        # Prepare upload
        if upload_successful:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Prokka process complete!\n\n'
                                                'Results are available at the following FTP address:\n'
                                                'ftp://ftp.agr.gc.ca/outgoing/cfia-ac/{}'.format(os.path.split(zip_filepath)[1]))
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity issues. '
                                                'Please try again later.')
        # Clean up files
        shutil.rmtree(output_folder)
        os.remove(zip_filepath)
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


def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    prokka_redmine()
