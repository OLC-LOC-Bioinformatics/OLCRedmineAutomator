import os
import glob
import click
import pickle
import shutil
from nastools.nastools import retrieve_nas_files
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def bagel4_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        # Variable to hold supplied arguments
        argument_dict = {
            'analysis': False
        }

        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.rstrip()
            if 'analysis' in item:
                argument_dict['analysis'] = item.split('=')[1]
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)
 

        # No SEQIDs
        if len(seqids) == 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You did not include any SEQIDs.')


        # Create folder to drop FASTA files
        assemblies_folder = os.path.join(work_dir, 'assemblies')
        os.makedirs(assemblies_folder, exist_ok=True)

        # Create output folder
        output_folder = os.path.join(work_dir, 'output')
        os.makedirs(output_folder)



        # Extract FASTA files to assemblies folder.
        retrieve_nas_files(seqids=seqids, outdir=assemblies_folder, filetype='fasta', copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, assemblies_folder)
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on '
                                                'the OLC NAS: {}'.format(missing_fastas))

          
        # Run bcgtree on the .faa files
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/bagel4'
        bagel4 = '/mnt/nas2/virtual_environments/bagel4/bagel4_2021/bagel4_wrapper.pl'

        #this command moves to the folder containing the .faa files and the config.txt file, then runs bcgTree
        bagel4_cmd = 'cd {assemblies_dir} && {bagel4}'.format(assemblies_dir=assemblies_folder, 
                                                              bagel4=bagel4)


        # Create another shell script to execute within the bagel4 conda environment
        templateb = "#!/bin/bash\n{} && {}".format(activate, bagel4_cmd)

        bagel4_script = os.path.join(work_dir, 'run_bagel4.sh')
        with open(bagel4_script, 'w+') as file:
            file.write(templateb)
        make_executable(bagel4_script)

        # Run shell script
        os.system(bagel4_script)

        #move files to new subfolder
#        output_files = os.listdir(assemblies_dir)

#        for file in output_files:
#            if file.endswith(".aln"):
#                shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))


        # Get bagel4 results uploaded

        # zip the bagel4 folder
        output_filename = 'bagel4_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=assemblies_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'

                # Upload the zip file to Dropbox
        download_link = upload_to_dropbox(
            access_token=DROPBOX_ACCESS_TOKEN,
            refresh_token=DROPBOX_REFRESH_TOKEN,
            app_key=DROPBOX_APP_KEY,
            app_secret=DROPBOX_APP_SECRET,
            local_file_path=zip_filepath
        )

        if download_link:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=2,
                notes='Analysis with bagel4 complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=download_link)
            )
        else:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=2,
                notes='Upload of results was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )
        # Clean up files
#        shutil.rmtree(output_folder)
#        shutil.rmtree(assemblies_folder)

        #don't remove the bagel4 zip output.. if you want to remove, uncomment the line below
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


def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    bagel4_redmine()
