import os
import click
import pickle
import shutil
from nastools.nastools import retrieve_nas_files


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def confindr_redmine(redmine_instance, issue, work_dir, description):
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
        if item.isalpha():
            pass
        else:
            seqids.append(item)

    # Create folder to drop FASTQ files
    raw_reads_folder = os.path.join(work_dir, 'raw_reads')
    os.makedirs(raw_reads_folder, exist_ok=True)

    # Create output folder
    output_folder = os.path.join(work_dir, 'output')
    os.makedirs(output_folder, exist_ok=True)

    # Extract FASTQ files.
    retrieve_nas_files(seqids=seqids, outdir=raw_reads_folder, filetype='fastq', copyflag=False)

    # These unfortunate hard coded paths appear to be necessary
    activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/confindr'
    confindr_py = '/mnt/nas2/virtual_environments/confindr/bin/confindr.py'

    # Database locations
    confindr_db = '/mnt/nas2/databases/confindr/databases/'

    # Prepare command
    cmd = '{confindr_py} --rmlst ' \
          '-i {raw_reads_folder} ' \
          '-o {output_folder} ' \
          '-d {confindr_db}'.format(confindr_py=confindr_py,
                                    raw_reads_folder=raw_reads_folder,
                                    output_folder=output_folder,
                                    confindr_db=confindr_db)

    # Create another shell script to execute within the confindr conda environment
    template = "#!/bin/bash\n{} && {}".format(activate, cmd)
    confindr_script = os.path.join(work_dir, 'run_confindr.sh')
    with open(confindr_script, 'w+') as file:
        file.write(template)
    make_executable(confindr_script)

    # Run shell script
    os.system(confindr_script)

    # Zip output
    output_filename = 'confindr_output'
    zip_filepath = zip_folder(results_path=output_folder,
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

    # Wrap up issue
    redmine_instance.issue.update(resource_id=issue.id, uploads=output_list, status_id=4,
                                  notes='ConFindr process complete!')


def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)


def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    confindr_redmine()
