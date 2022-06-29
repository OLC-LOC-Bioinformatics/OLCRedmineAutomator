import os
import glob
import click
import pickle
import shutil
from nastools.nastools import retrieve_nas_files
from automator_settings import COWBAT_DATABASES


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def sipprverse_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Current list of analysis types that the sipprverse can perform
    analyses = [
        'custom', 'full', 'gdcs', 'genesippr', 'mash', 'mlst', 'pointfinder', 'resfinder', 'rmlst', 'serosippr',
        'sixteens', 'virulence'
    ]
    # Dictionary of analysis types to argument flags to pass to the script
    argument_flags = {
                         'custom': '-U {fasta}'.format(fasta=os.path.join(work_dir, 'targets.tfa')),
                         'gdcs': '-Q',
                         'genesippr': '-G',
                         'mash': '-C',
                         'mlst': '-M',
                         'pointfinder': '-P',
                         'resfinder': '-A',
                         'rmlst': '-R',
                         'serosippr': '-S',
                         'sixteens': '-X',
                         'virulence': '-V',
                         'full': '-F'
    }
    # Variable to hold supplied arguments
    argument_dict = {
        'analysis': str(),
        'averagedepth': 2,
        'kmersize': 19,
        'cutoff': 0.90,
        'allowsoftclips': False,
    }
    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'AVERAGEDEPTH' in item:
                argument_dict['averagedepth'] = int(item.split('=')[1].lower())
                continue
            if 'CUTOFF' in item:
                argument_dict['cutoff'] = float(item.split('=')[1].lower())
                continue
            if 'KMERSIZE' in item:
                argument_dict['kmersize'] = int(item.split('=')[1].lower())
                continue
            if 'ANALYSIS' in item:
                argument_dict['analysis'] = item.split('=')[1].lower()
                continue
            if 'ALLOWSOFTCLIPS' in item:
                argument_dict['allowsoftclips'] = True
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)
        if argument_dict['analysis'] == 'custom':
            # Download the attached FASTA file.
            # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
            # out a better way to do it.
            attachment = redmine_instance.issue.get(issue.id, include='attachments')
            attachment_id = 0
            for item in attachment.attachments:
                attachment_id = item.id
            # Set the name of and create the folder to store the targets
            dbpath = work_dir
            # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
            if attachment_id != 0:
                attachment = redmine_instance.attachment.get(attachment_id)
                attachment.download(savepath=work_dir, filename='targets.tfa')
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Analysis type custom requires an attached FASTA file of '
                                                    'targets. The automator could not find any attached files. '
                                                    'Please create a new issue with the FASTA file attached and try '
                                                    'again.',
                                              status_id=4)
                return
        else:
            dbpath = COWBAT_DATABASES
        # Ensure that the analysis type is provided
        if not argument_dict['analysis']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No analysis type provided. '
                                                'Please ensure that issue contains '
                                                '"analysistype=requested_analysis_type", where requested_analysis_type '
                                                ' is one of the following keywords: {ats}. See the the usage guide: '
                                                'https://olc-bioinformatics.github.io/redmine-docs/analysis/sipprverse/'
                                                ' for additional details'.format(ats=','.join(analyses)))
            return
        # Ensure that the supplied analysis type is valid
        if argument_dict['analysis'] not in analyses:
            redmine_instance.issue \
                .update(resource_id=issue.id,
                        notes='WARNING: Requested analysis type: {at} not in list of supported analyses: {ats}. Please '
                              'see https://olc-bioinformatics.github.io/redmine-docs/analysis/sipprverse/ '
                              'for additional details'.format(at=argument_dict['analysis'],
                                                              ats=','.join(analyses)))
            return
        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No SEQIDs provided!')
        # Run file linker and then make sure that all FASTQ files requested are present. Warn user if they
        # requested things that we don't have.
        retrieve_nas_files(seqids=seqids,
                           outdir=work_dir,
                           filetype='fastq',
                           copyflag=False)
        missing_fastqs = verify_fastq_files_present(seqids, work_dir)
        # Update the Redmine issue if one or more of the requested SEQIDs could not be located
        if missing_fastqs:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastqs))
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dev/sipprverse'
        sippr_py = '/mnt/nas2/virtual_environments/dev/sipprverse/bin/sippr.py'
        # Run sipprverse with the necessary arguments
        sippr_cmd = 'python {sippr_py} -s {seqpath} -o {outpath} -r {dbpath} -a {ad} -k {ks} -c {cut} {at}'\
            .format(sippr_py=sippr_py,
                    seqpath=work_dir,
                    outpath=work_dir,
                    dbpath=dbpath,
                    ad=argument_dict['averagedepth'],
                    ks=argument_dict['kmersize'],
                    cut=argument_dict['cutoff'],
                    at=argument_flags[argument_dict['analysis']])
        # Add the allow_soft_clips option if required
        sippr_cmd += ' -sc' if argument_dict['allowsoftclips'] else ''
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Sipprverse command:\n {cmd}'.format(cmd=sippr_cmd))
        # Create another shell script to execute within the conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, sippr_cmd)
        sipprverse_script = os.path.join(work_dir, 'run_sipprverse.sh')
        with open(sipprverse_script, 'w+') as file:
            file.write(template)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(sipprverse_script)
        # Run shell script
        os.system(sipprverse_script)

        # Zip output
        output_filename = 'sipprverse_output'
        zip_filepath = zip_folder(results_path=os.path.join(work_dir, 'reports'),
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

        # Create a list of all the folders - will be used to clean up the working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)
        # Wrap up issue
        redmine_instance.issue.update(resource_id=issue.id,
                                      uploads=output_list,
                                      status_id=4,
                                      notes='{at} analysis with sipprverse complete!'
                                      .format(at=argument_dict['analysis'].lower()))
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {error}\nAlternatively, please check out '
                                            'the usage guide: '
                                            'https://olc-bioinformatics.github.io/redmine-docs/analysis/sipprverse/'
                                      .format(error=e))


def verify_fastq_files_present(seqid_list, fastq_dir):
    """
    Makes sure that FASTQ files specified in seqid_list have been successfully copied/linked to directory specified
    by fastq_dir.
    :param seqid_list: List with SEQIDs.
    :param fastq_dir: Directory that FASTQ files (both forward and reverse reads should have been copied to
    :return: List of SEQIDs that did not have files associated with them.
    """
    missing_fastqs = list()
    for seqid in seqid_list:
        # Check forward.
        if len(glob.glob(fastq_dir + '/' + seqid + '*R1*fastq*')) == 0:
            missing_fastqs.append(seqid)
        # Check reverse. Only add to list of missing if forward wasn't present.
        if len(glob.glob(fastq_dir + '/' + seqid + '*R2*fastq*')) == 0 and seqid not in missing_fastqs:
            missing_fastqs.append(seqid)
    return missing_fastqs


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
    sipprverse_redmine()
