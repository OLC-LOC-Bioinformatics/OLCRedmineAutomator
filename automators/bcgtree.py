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
            'analysis': False,
            'bootstraps': 100,
            'min_proteomes': 2,
            'aa_substitution_model': 'AUTO'
        }

        #amino acid substitution models supported by automator
        aa_substitution_models = [
            'AUTO', 'DAYHOFF', 'DCMUT', 'JTT', 'MTREV', 'WAG','RTREV', 'CPREV', 'VT', 'BLOSUM62', 'MTMAM', 'LG', \
            'MTART', 'MTZOA', 'PMB', 'HIVB','HIVW', 'JTTDCMUT', 'FLU', 'STMTREV', 'DUMMY', 'DUMMY2', \
            'LG4M', 'LG4X', 'PROT_FILE', 'GTR_UNLINKED', 'GTR'
        ]

        #analysis types supported by automator
        analyses = [
            'enterobacterales', 'listeriaceae', 'enterobacter'
        ]

        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.rstrip()
            if 'analysis' in item:
                argument_dict['analysis'] = item.split('=')[1]
                continue
            if 'bootstraps' in item:
                argument_dict['bootstraps'] = item.split('=')[1]
                continue
            if 'min_proteomes' in item:
                argument_dict['min_proteomes'] = item.split('=')[1]
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
            dst = assemblies_folder
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

    #if enterobacter, link to enterobacter reference sequences
        if argument_dict['analysis'] == 'enterobacter':
            src = '/mnt/nas2/processed_sequence_data/ncbi/enterobacter_references/BestAssemblies/'
            dst = assemblies_folder
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)
            src2 = '/mnt/nas2/processed_sequence_data/ncbi/enterobacterales_references/BestAssemblies/'
            dst = assemblies_folder
            lncmd2 = 'ln -s {src}ENTEROBACTER_*.fasta {dst}'.format(src=src2,
                                                                    dst=dst)
            os.system(lncmd2)

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
            faa_file = '{seqid}.faa'.format(seqid=seqid)
            output_faa = os.path.join(work_dir, 'bcgtree', faa_file)

            # Copy the .faa files to the output folder
            shutil.copyfile(os.path.join(output_folder, seqid, faa_file),
                            output_faa)

        #create config file to pass to bcgtree
        config_file = os.path.join(bcgtree_folder, 'config.txt')
        bootstraps = "--bootstraps={boots}".format(boots=argument_dict['bootstraps'])
        outdir = "--outdir={out}".format(out=bcgtree_folder)
        minproteomes = "--min-proteomes={minp}".format(minp=argument_dict['min_proteomes'])
        aasubmod = '--raxml-aa-substitution-model "{aamod}"'.format(aamod=argument_dict['aa_substitution_model'])
        configtemplate = "--threads=8\n{}\n{}\n{}\n{}".format(bootstraps, outdir, minproteomes, aasubmod)
        # with open(config_file, 'w+') as file:
        with open(config_file, 'w+') as file:
            file.write(configtemplate)
        for proteome in glob.glob(os.path.join(bcgtree_folder, '*.faa')):
            proteomeid = os.path.split(proteome)[1].split('.')[0]
            proteomefa = os.path.split(proteome)[1]
            # add a line to the config file for each proteome/seqid in the output folder            
            with open(config_file, 'a+') as file:
                file.write("\n--proteome {proteomeid}={proteome} ".format(proteomeid=proteomeid, proteome=proteomefa))
            


        # Run bcgtree on the .faa files
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/bcgTree'
        bcgtree = '/mnt/nas2/virtual_environments/bcgTree/bcgTree/bin/bcgTree.pl'

        #this command moves to the folder containing the .faa files and the config.txt file, then runs bcgTree
        bcgtree_cmd = 'cd {bcgtree_dir} && {bcgtree} @{config}'.format(bcgtree_dir=bcgtree_folder, 
                                                                       bcgtree=bcgtree,
                                                                       config=config_file)

#        redmine_instance.issue.update(resource_id=issue.id,
#                                      notes='Prokka complete!\n'
#                                            'Now running bcgtree on amino acid fasta files (.faa)')

        # Create another shell script to execute within the bcgtree conda environment
        templateb = "#!/bin/bash\n{} && {}".format(activate, bcgtree_cmd)
#        template = "#!/bin/bash\n{} && {}".format(activate, bcgtree_cmd)
        bcgtree_script = os.path.join(work_dir, 'run_bcgtree.sh')
        with open(bcgtree_script, 'w+') as file:
            file.write(templateb)
        make_executable(bcgtree_script)

        # Run shell script
        os.system(bcgtree_script)

        #move alignment files to new subfolder
        output_files = os.listdir(bcgtree_folder)
#        aln_files = [file for file in output_files if file.endswith(".aln*")]
        alignments_folder = os.path.join(bcgtree_folder, 'alignment_files')
        os.makedirs(alignments_folder)        
        for file in output_files:
            if file.endswith(".aln"):
                shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
        for file in output_files:
            if file.endswith(".aln-gb"):
                shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
        for file in output_files:
            if file.endswith(".aln-gb.comp"):
                shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
        for file in output_files:
            if file.endswith(".aln-gb.fa"):
                shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
        for file in output_files:
            if file.endswith(".aln-gb.htm"):
                shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
        #move id files to a new subfolder
        ids_folder = os.path.join(bcgtree_folder, 'ids_files')
        os.makedirs(ids_folder)        
        for file in output_files:
            if file.endswith(".ids"):
                shutil.move(os.path.join(bcgtree_folder,file), os.path.join(ids_folder,file))
        #move protein fasta files to a new subfolder
        fa_files = os.listdir(bcgtree_folder)
        fastas_folder = os.path.join(bcgtree_folder, 'fasta_files')
        os.makedirs(fastas_folder)        
        for file in fa_files:
            if file.endswith(".fa"):
                shutil.move(os.path.join(bcgtree_folder,file), os.path.join(fastas_folder,file))

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
        # Remove the prokka zip file
        os.remove(zip_filepath)

        # Get bcgtree results uploaded
#        for faa in glob.glob(os.path.join(bcgtree_folder, '*.faa')):
#            os.remove(faa)
        # zip the bcgtree folder
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
        shutil.rmtree(output_folder)
        shutil.rmtree(assemblies_folder)
#        shutil.rmtree(bcgtree_folder)
        #don't remove the bcgtree zip output.. if you want to remove, uncomment the line below
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
