import os
import glob
import click
import pickle
import shutil
from nastools.nastools import retrieve_nas_files
from externalretrieve import upload_to_ftp


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def bcgtree_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))


        # Variable to hold supplied arguments
        argument_dict = {
            'bootstraps': 100,
            'min_proteomes': 2,
            'aa_substitution_model': 'AUTO'
        }

        #amino acid substitution models supported by automator
        aa_substitution_models = [
            'AUTO', 'DAYHOFF', 'DCMUT', 'JTT', 'MTREV', 'WAG','RTREV', 'CPREV', 'VT', 'BLOSUM62', 'MTMAM', 'LG', \
            'MTART', 'MTZOA', 'PMB', 'HIVB','HIVW', 'JTTDCMUT', 'FLU', 'STMTREV', 'DUMMY', 'DUMMY2', 'AUTO', \
            'LG4M', 'LG4X', 'PROT_FILE', 'GTR_UNLINKED', 'GTR'
        ]

        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.rstrip()
            if 'bootstraps' in item:
                argument_dict['bootstraps'] = item.split('=')[1]
                continue
            if 'min_proteomes' in item:
                argument_dict['aa_substitution_model'] = item.split('=')[1]
                continue
            if 'aa_substitution_model' in item:
                argument_dict['aa_substitution_model'] = item.split('=')[1].upper()
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)
 

        # Ensure that the amino acid substitution model is acceptable
        if argument_dict['aa_substitution_model'] not in aa_substitution_models:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied amino acid substitution model {at} current not in the supported '
                                                'list of models: {ats}'.format(at=argument_dict['aa_substitution_model'],
                                                                               ats=', '.join(aa_substitution_models)),
                                          status_id=4)
            return        
        # No SEQIDs
        if len(seqids) == 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You did not include any SEQIDs.')


        # Create folder to drop FASTA files
        assemblies_folder = os.path.join(work_dir, 'assemblies')
        os.mkdir(assemblies_folder)

        # Create output folder
        output_folder = os.path.join(work_dir, 'output')
        os.makedirs(output_folder)

        # Create the folder in which .fa files are to be stored
        bcgtree_folder = os.path.join(work_dir, 'bcgtree')
        os.makedirs(bcgtree_folder)
        # Extract FASTA files to assemblies folder.
        retrieve_nas_files(seqids=seqids, outdir=assemblies_folder, filetype='fasta', copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, assemblies_folder)
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on '
                                                'the OLC NAS: {}'.format(missing_fastas))

        # These unfortunate hard coded paths appear to be necessary
        activatep = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/roary'
        prokka = '/mnt/nas2/virtual_environments/roary/bin/prokka'

        for assembly in glob.glob(os.path.join(assemblies_folder, '*.fasta')):
            seqid = os.path.split(assembly)[1].split('.')[0]
            # Prepare command
            cmd = '{prokka} --outdir {output_folder} --prefix {seqid} {assembly}'\
                .format(prokka=prokka,
                        output_folder=os.path.join(output_folder, seqid),
                        seqid=seqid,
                        assembly=assembly)

            # Create another shell script to execute within the conda environment
            template = "#!/bin/bash\n{} && {}".format(activatep, cmd)
            prokka_script = os.path.join(work_dir, 'run_prokka.sh')
            with open(prokka_script, 'w+') as file:
                file.write(template)
            make_executable(prokka_script)

            # Run shell script
            os.system(prokka_script)
            fa_file = '{seqid}.fa'.format(seqid=seqid)
            output_fa = os.path.join(work_dir, 'bcgtree', fa_file)

            # Copy the .fa files to the output folder
            shutil.copyfile(os.path.join(output_folder, seqid, fa_file),
                            output_fa)

        #create config file to pass to bcgtree
        config_file = os.path.join(output_folder, 'config.txt')
        # with open(config_file, 'w+') as file:
        for proteome in glob.glob(os.path.join(output_folder, '*.fa')):
            proteomeid = os.path.split(proteome)[1].split('.')[0]
            # add a line to the config file for each proteome/seqid in the output folder            
            with open(config_file, 'a+') as file:
                file.write("--proteome {proteomeid}={proteome} ".format(proteomeid=proteomeid, proteome=proteome))
                return
            return


        # Run bcgtree on the .fa files
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/bcgTree'
        bcgtree = '/mnt/nas2/virtual_environments/bcgTree/bcgTree/bin/bcgTree.pl'

        #the command moves to the folder containing the .fa files and the config.txt file
        bcgtree_cmd = 'cd {bcgtree_dir} && {bcgtree} @{config}'.format(bcgtree_dir=output_folder, 
                                                                       bcgtree=bcgtree,
                                                                       config=config_file)

#        redmine_instance.issue.update(resource_id=issue.id,
#                                      notes='Prokka complete!\n'
#                                            'bcgtree command:{bcgtreecmd}'.format(bcgtreecmd=bcgtree_cmd)

        # Create another shell script to execute within the bcgtree conda environment
        templateb = "#!/bin/bash\n{} && {}".format(activate, bcgtree_cmd)
#        template = "#!/bin/bash\n{} && {}".format(activate, bcgtree_cmd)
        bcgtree_script = os.path.join(work_dir, 'run_bcgtree.sh')
        with open(bcgtree_script, 'w+') as file:
            file.write(templateb)
        make_executable(bcgtree_script)

        # Run shell script
        os.system(bcgtree_script)

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
                                                'ftp://ftp.agr.gc.ca/outgoing/cfia-ac/{}'
                                          .format(os.path.split(zip_filepath)[1]))
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity '
                                                'issues. '
                                                'Please try again later.')
        # Remove the zip file
#        os.remove(zip_filepath)

        # Get bcgtree results uploaded
#        for fa in glob.glob(os.path.join(bcgtree_folder, '*.fa')):
#            os.remove(fa)
        # zip the roary folder
        output_filename = 'bcgtree_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=bcgtree_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'
        # Prepare upload
        # This file can get too big to upload to Redmine, so we should put it on the FTP.
#        output_list = [
#            {
#                'filename': os.path.basename(zip_filepath),
#                'path': zip_filepath
#            }
#        ]
        upload_successful = upload_to_ftp(local_file=zip_filepath)
        # Prepare upload
        if upload_successful:
        # Wrap up issue
            redmine_instance.issue.update(resource_id=issue.id,
                                      status_id=4,
                                      notes='Analysis with bcgTree complete!\n\nResults are available at the following FTP address:\nftp://ftp.agr.gc.ca/outgoing/cfia-ac/{l}'.format(l=os.path.split(zip_filepath)[1]))
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity '
                                                'issues. '
                                                'Please try again later.')
        # Clean up files
#        shutil.rmtree(output_folder)
#        shutil.rmtree(roary_folder)
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
    bcgtree_redmine()
