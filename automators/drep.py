import os
import glob
import click
import pickle
import shutil
import ftplib
import sentry_sdk
from automator_settings import SENTRY_DSN
from amrsummary import before_send
from externalretrieve import upload_to_ftp
from automator_settings import FTP_USERNAME, FTP_PASSWORD
from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def drep_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    #analyses supported by the automator
    analyses = [
        'custom', 'enterobacterales', 'listeriaceae', 'enterobacter'
    ]
    #commands supported by automator
    commands = [
        'compare', 'dereplicate'
    ]
    #clustering algorithms supported by the automator
    clusteralgorithms = [
        'median', 'weighted', 'single', 'complete', 'average', 'ward', 'centroid'
    ]

    #algorithm for secondary clustering comparisons:
    comparisonalgorithms = [
        'fastANI', 'ANImf', 'ANIn', 'gANI', 'goANI'
    ]

    #variable to hold supplied arguments
    argument_dict = {
        'analysis': str(),
        'command': 'compare',
        'debug': 'False',
        'comparisonalgorithm': 'ANImf',
        'clusteralgorithm': 'average',
        'P_ANI': 0.9,
        'S_ANI': 0.99,
        'coverage_threshold': 0.1,
    }


    #parse description to figure out analysis type, and find fastas
    try:
        #parse description
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'ANALYSIS' in item:
               argument_dict['analysis'] = item.split('=')[1].lower()
               continue
            if 'COMMAND' in item:
                argument_dict['command'] = item.split('=')[1].lower()
                continue
            if 'DEBUG' in item:
                argument_dict['debug'] = item.split('=')[1].capitalize()
                continue
            if 'CLUSTERALGORITHM' in item:
                argument_dict['clusteralgorithm'] = item.split('=')[1].lower()
                continue
            if 'COMPARISONALGORITHM' in item:
               argument_dict['comparisonalgorithm'] = item.split('=')[1].lower()
               continue
#            if 'MASH_ANI_THRESHOLD' in item:
#               argument_dict['P_ANI'] = int(float(item.split('=')[1].lower()))
#               continue
            if 'MASH_ANI_THRESHOLD' in item:
               argument_dict['P_ANI'] = str(float(item.split('=')[1].lower()))
               continue
            if 'S_ANI' in item:
               argument_dict['S_ANI'] = str(float(item.split('=')[1].lower()))
               continue
            if 'COVERAGE_THRESHOLD' in item:
               argument_dict['coverage_threshold'] = str(float(item.split('=')[1].lower()))
               continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that the analysis type is provided
        if not argument_dict['analysis']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not identify an analysis type. '
                                                'Please ensure that the first line of the issue contains one'
                                                ' of the following keywords: {ats}'.format(ats=', '.join(analyses)),
                                          status_id=4)
            return
        elif argument_dict['analysis'] not in analyses:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied analysis type {at} current not in the supported '
                                                'list of analyses: {ats}'.format(at=argument_dict['analysis'],
                                                                                 ats=', '.join(analyses)),
                                          status_id=4)
            return
        # Ensure that the command type is acceptable
        if argument_dict['command'] not in commands:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied command type {at} current not in the supported '
                                                'list of commands: {ats}'.format(at=argument_dict['command'],
                                                                                   ats=', '.join(commands)),
                                          status_id=4)
            return
        # Ensure that the comparison algorithm type is acceptable
        if argument_dict['comparisonalgorithm'] not in comparisonalgorithms:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied comparison algorithm type {at} current not in the supported '
                                                'list of algorithms: {ats}'.format(at=argument_dict['comparisonalgorithm'],
                                                                                 ats=', '.join(comparisonalgorithms)),
                                          status_id=4)
            return
        # Ensure that the cluster algorithm type is acceptable
        if argument_dict['clusteralgorithm'] not in clusteralgorithms:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied cluster algorithm type {at} current not in the supported '
                                                'list of algorithms: {ats}'.format(at=argument_dict['clusteralgorithm'],
                                                                                 ats=', '.join(clusteralgorithms)),
                                          status_id=4)
            return
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
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)
            #also do this for the atcc sequences
            src2 = '/mnt/nas2/processed_sequence_data/atcc/listeria/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst2 = fasta_dir
            lncmd2 = 'ln -s {src}*.fasta {dst}'.format(src=src2,
                                                      dst=dst2)
            os.system(lncmd2)

    #if enterobacter, link to enterobacter reference sequences
        if argument_dict['analysis'] == 'enterobacter':
            src = '/mnt/nas2/processed_sequence_data/ncbi/enterobacter_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

            src2 = '/mnt/nas2/processed_sequence_data/ncbi/enterobacterales_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd2 = 'ln -s {src}ENTEROBACTER_*.fasta {dst}'.format(src=src2,
                                                                    dst=dst)
            os.system(lncmd2)

    #if just a custom analysis, create the fastas folder
        if argument_dict['analysis'] == 'custom':
            # Make fasta file directory
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)

    # Rest of Description should just be a list of SEQIDs. Get the fasta files associated with them extracted
	# Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
        retrieve_nas_files(seqids=seqids,
                           outdir=fasta_dir,
                           filetype='fasta',
                           copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, fasta_dir)
        
	# Verify that specified fasta files are actually there, warn user if they aren't.
        missing_fastas = verify_fasta_files_present(seqid_list=seqids,
                                                    fasta_dir=os.path.join(work_dir, 'fastas'))
        if len(missing_fastas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))
	# Make output dir
        output_dir = os.path.join(work_dir, 'dRep_results')
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

	
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dRep'
        # Run dRep compare with the necessary arguments
        if argument_dict['command'] == 'compare':
            dRep_cmd = 'dRep {command} {outpath} -p 24 -g {seqpath}/*.fasta -pa {P_ANI} -sa {S_ANI} -nc {cov_thresh} ' \
                       '--clusterAlg {clusteralgorithm} --S_algorithm {comparisonalgorithm}'.format(command=argument_dict['command'],
                                                                                                    seqpath=fasta_dir,
                                                                                                    P_ANI=argument_dict['P_ANI'],
                                                                                                    S_ANI=argument_dict['S_ANI'],
                                                                                                    cov_thresh=argument_dict['coverage_threshold'],
                                                                                                    clusteralgorithm=argument_dict['clusteralgorithm'],
                                                                                                    comparisonalgorithm=argument_dict['comparisonalgorithm'],
                                                                                                    outpath=output_dir)


            #create another shell script to execute within the dRep conda environment
            template = "#!/bin/bash\n{} && {}".format(activate, dRep_cmd)
            dRep_script = os.path.join(work_dir, 'run_dRep.sh')
            with open(dRep_script, 'w+') as file:
                file.write(template)

        # Run dRep dreplicate with the necessary arguments
        if argument_dict['command'] == 'dereplicate':
            dRep_cmd = 'dRep {command} {outpath} -p 24 -g {seqpath}/*.fasta -pa {P_ANI} -sa {S_ANI} -nc {cov_thresh} ' \
                       '--clusterAlg {clusteralgorithm} --S_algorithm {comparisonalgorithm}'.format(command=argument_dict['command'],
                                                                                                    seqpath=fasta_dir,
                                                                                                    P_ANI=argument_dict['P_ANI'],
                                                                                                    S_ANI=argument_dict['S_ANI'],
                                                                                                    cov_thresh=argument_dict['coverage_threshold'],
                                                                                                    clusteralgorithm=argument_dict['clusteralgorithm'],
                                                                                                    comparisonalgorithm=argument_dict['comparisonalgorithm'],
                                                                                                    outpath=output_dir)


            #create another shell script to execute within the dRep conda environment
            template = "#!/bin/bash\n{} && {}".format(activate, dRep_cmd)
            dRep_script = os.path.join(work_dir, 'run_dRep.sh')
            with open(dRep_script, 'w+') as file:
                file.write(template)

        # Modify the permissions of the script to allow it to be run on the node
        make_executable(dRep_script)
        # Run shell script
        os.system(dRep_script)

	# Zip output folder and upload to the FTP, as sometimes files are >10Mb
        shutil.make_archive(root_dir=output_dir,
                            format='zip',
                            base_name=os.path.join(work_dir, str(issue.id)))

        sas_url = upload_to_ftp(
            local_file=os.path.join(work_dir, str(issue.id) + '.zip')
        )

        if sas_url:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='dRep complete!\n\n'
                'Results are available at the following URL:\n'
                '{url}'.format(url=sas_url)
            )
        else:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='Upload of result files was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )

        # Create a list of all the folders - will be used to clean up the working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)

        os.remove(os.path.join(work_dir, str(issue.id) + '.zip')) #delete zip file
        # Wrap up issue
#        redmine_instance.issue.update(resource_id=issue.id,
#                                      uploads=output_list,
#                                      status_id=4,
#                                      notes='{at} analysis with dRep complete!'.format(at=argument_dict['analysis'].lower()))
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))

def verify_fasta_files_present(seqid_list, fasta_dir):
    """
    Makes sure that FASTQ files specified in seqid_list have been successfully copied/linked to directory specified
    by fastq_dir.
    :param seqid_list: List with SEQIDs.
    :param fasta_dir: Directory to which FASTA files should have been linked
    :return: List of SEQIDs that did not have files associated with them.
    """
    missing_fastas = list()
    for seqid in seqid_list:
        # Check forward.
        if len(glob.glob(os.path.join(fasta_dir, '{seqid}*fasta*'.format(seqid=seqid)))) == 0:
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

def zip_folder(results_path, output_dir, output_filename):
    """
    Compress a folder
    :param results_path: The path of the folder to be compressed
    :param output_dir: The output directory
    :param output_filename: The output file name
    :return:
    """
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    drep_redmine()
