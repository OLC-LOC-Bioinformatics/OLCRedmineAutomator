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
def roary_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))
        # Variable to hold supplied arguments
        argument_dict = {
            'union': False,
            'intersection': False,
            'complement': False,
            'difference': False,
            'set_one': list(),
            'set_two': list()
        }
        set_one = True
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        analysistype = "" # Just for the output at the end
        for item in description:
            item = item.upper().rstrip()
            if 'UNION' in item:
                argument_dict['union'] = True
                analysistype += "union "
                continue
            if 'INTERSECTION' in item:
                argument_dict['intersection'] = True
                analysistype += "intersection "
                continue
            if 'COMPLEMENT' in item:
                argument_dict['complement'] = True
                analysistype += "complement "
                continue
            if 'DIFFERENCE' in item:
                argument_dict['difference'] = True
                analysistype += "difference "
                continue
            if 'SET_TWO' in item:
                set_one = False
                continue
            if 'SET_ONE' in item:
                set_one = True
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)
            # Add the SEQID to the appropriate list
            if set_one:
                argument_dict['set_one'].append(item)
            else:
                argument_dict['set_two'].append(item)
        # Sanity check for arguments - too many analysis types
        if argument_dict['union'] + argument_dict['intersection'] + argument_dict['complement'] + \
                argument_dict['difference'] > 1:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You included too many analysis types.')
        # No analysis types
        elif argument_dict['union'] + argument_dict['intersection'] + argument_dict['complement'] + \
                argument_dict['difference'] == 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You did not include an analysis type.')
        # No SEQIDs
        if len(seqids) == 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You did not include any SEQIDs.')
        # Difference analysis type, but no input_set_two provided
        if argument_dict['difference'] and not argument_dict['set_two']:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You specified the "difference" analysis type, but did not include any '
                                                'SEQIDs for input_set_two.')
        # Create folder to drop FASTA files
        assemblies_folder = os.path.join(work_dir, 'assemblies')
        os.mkdir(assemblies_folder)

        # Create output folder
        output_folder = os.path.join(work_dir, 'output')
        os.makedirs(output_folder)
        # Create the folder in which .gff files are to be stored
        roary_folder = os.path.join(work_dir, 'roary')
        os.makedirs(roary_folder)
        # Extract FASTA files.
        retrieve_nas_files(seqids=seqids, outdir=assemblies_folder, filetype='fasta', copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, assemblies_folder)
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on '
                                                'the OLC NAS: {}'.format(missing_fastas))

        set_one_gff = list()
        set_two_gff = list()
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/roary'
        prokka = '/mnt/nas2/virtual_environments/roary/bin/prokka'

        for assembly in glob.glob(os.path.join(assemblies_folder, '*.fasta')):
            seqid = os.path.split(assembly)[1].split('.')[0]
            # Prepare command
            cmd = '{prokka} --outdir {output_folder} --prefix {seqid} {assembly}'\
                .format(prokka=prokka,
                        output_folder=os.path.join(output_folder, seqid),
                        seqid=seqid,
                        assembly=assembly)

            # Create another shell script to execute within the roary conda environment
            template = "#!/bin/bash\n{} && {}".format(activate, cmd)
            prokka_script = os.path.join(work_dir, 'run_prokka.sh')
            with open(prokka_script, 'w+') as file:
                file.write(template)
            make_executable(prokka_script)

            # Run shell script
            os.system(prokka_script)
            gff_file = '{seqid}.gff'.format(seqid=seqid)
            output_gff = os.path.join(work_dir, 'roary', gff_file)
            # Copy the .gff file to the roary folder
            shutil.copyfile(os.path.join(output_folder, seqid, gff_file),
                            output_gff)
            # Put the output files into appropriate lists
            if seqid in argument_dict['set_two']:
                set_two_gff.append(gff_file)
            else:
                set_one_gff.append(gff_file)
        # Run roary on the .gff files
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/roary'
        roary = '/mnt/nas2/virtual_environments/roary/bin/roary'
        query = '/mnt/nas2/virtual_environments/roary/bin/query_pan_genome'

        analysis_type = str()
        for key, value in argument_dict.items():
            if type(value) == bool:
                if value:
                    analysis_type = key
        roary_cmd = 'cd {roary_dir} && {roary} *.gff'.format(roary_dir=roary_folder,
                                                             roary=roary)
        if analysis_type != 'difference':
            query_cmd = 'cd {roary_dir} && {query} -a {analysis_type} -o {output_file} *.gff'\
                .format(query=query,
                        analysis_type=analysis_type,
                        output_file=os.path.join(roary_folder, 'pan_genome_results'),
                        roary_dir=roary_folder,
                        gff_files=','.join(set_one_gff))
        else:
            query_cmd = 'cd {roary_dir} && {query} -a {analysis_type} -o {output_file} -i {set_one_files} ' \
                        '-t {set_two_files}'\
                .format(roary_dir=roary_folder,
                        query=query,
                        analysis_type=analysis_type,
                        output_file=os.path.join(roary_folder, 'pan_genome_results'),
                        set_one_files=','.join(set_one_gff),
                        set_two_files=','.join(set_two_gff))
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Prokka complete!\n'
                                            'roary command:{roary}\n'
                                            'query command:{query}'.format(roary=roary_cmd,
                                                                           query=query_cmd))
        # Create another shell script to execute within the roary conda environment
        template = "#!/bin/bash\n{activate} && {roary} && {query}".format(activate=activate,
                                                                          roary=roary_cmd,
                                                                          query=query_cmd)
        roary_script = os.path.join(work_dir, 'run_roary.sh')
        with open(roary_script, 'w+') as file:
            file.write(template)
        make_executable(roary_script)

        # Run shell script
        os.system(roary_script)

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
                                                'ftp://ftp.agr.gc.ca/outgoing/cfia-ak/{}'
                                          .format(os.path.split(zip_filepath)[1]))
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity '
                                                'issues. '
                                                'Please try again later.')
        # Remove the zip file
        os.remove(zip_filepath)
        # Get roary results uploaded
        for gff in glob.glob(os.path.join(roary_folder, '*.gff')):
            os.remove(gff)
        # zip the roary folder
        output_filename = 'roary_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=roary_folder,
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
 
        # Wrap up issue
        upload_successful = upload_to_ftp(local_file=zip_filepath)
        if upload_successful:
            redmine_instance.issue.update(resource_id=issue.id,
                                      status_id=4,
                                      notes='{at} analysis with roary complete!\n\nResults are available at the following FTP address:\nftp://ftp.agr.gc.ca/outgoing/cfia-ak/{l}'.format(at=analysistype, l=os.path.split(zip_filepath)[1]))
        # Clean up files
        shutil.rmtree(output_folder)
        shutil.rmtree(roary_folder)
        os.remove(zip_filepath)
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
    roary_redmine()
