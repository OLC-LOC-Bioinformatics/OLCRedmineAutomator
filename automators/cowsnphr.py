import os
import glob
import click
import pickle
import shutil
import zipfile
import sentry_sdk
from automator_settings import SENTRY_DSN
from amrsummary import before_send
from biotools import mash
from externalretrieve import upload_to_ftp
from nastools.nastools import retrieve_nas_files


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def cowsnphr_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))
        #
        query_list = list()
        reference = list()
        compare = False
        # Go through description to figure out what our query is and what the reference is.
        for item in description:
            item = item.upper()
            if item == '':
                continue
            if 'COMPARE' in item:
                compare = True
                continue
            if compare:
                query_list.append(item)
            else:
                if 'REFERENCE' not in item:
                    reference.append(item)
        # Create output folder
        reference_folder = os.path.join(work_dir, 'ref')
        os.makedirs(reference_folder)
        # Retrieve our reference file. Error user if they selected anything but one reference and don't continue.
        if len(reference) != 1:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='ERROR: You must specify one reference strain, and you '
                                                'specified {} reference strains. Please create a new'
                                                ' issue and try again.'.format(len(reference)), status_id=4)
            return

        if reference[0].upper() != 'ATTACHED':
            # Extract our reference file to our working directory.
            retrieve_nas_files(seqids=reference,
                               outdir=reference_folder,
                               filetype='fasta',
                               copyflag=True)
            # Check that the file was successfully extracted. If it wasn't boot the user.
            if len(glob.glob(os.path.join(reference_folder, '*fasta'))) == 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Could not find the specified reference file.'
                                                    ' Please verify it is a correct SEQID, create a new '
                                                    'issue, and try again.', status_id=4)
                return

        # If user specified attachment as the reference file, download it to our working directory.
        else:
            # Get the attachment ID, and download if it isn't equal to zero (meaning no attachment, so boot user with
            # appropriate error message)
            attachment = redmine_instance.issue.get(issue.id, include='attachments')
            attachment_id = 0
            ref_name = 'reference.fasta'
            for item in attachment.attachments:
                attachment_id = item.id
                ref_name = item.filename
            # Download if we found an attachment, and use as our reference. Otherwise, exit and tell user to try again
            if attachment_id != 0:
                attachment = redmine_instance.attachment.get(attachment_id)
                attachment.download(savepath=reference_folder, filename=ref_name)
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: You specified that the reference would be in attached file,'
                                                    ' but no attached file was found. Please create a new issue and '
                                                    'try again.',
                                              status_id=4)
                return

        # PROKKA
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/prokka'
        prokka = '/mnt/nas2/virtual_environments/prokka/bin/prokka'

        # Prepare command
        ref_file = glob.glob(os.path.join(reference_folder, '*fasta'))[0]
        prefix = os.path.splitext(os.path.basename(ref_file))[0]
        cmd = '{prokka} --force --outdir {output_folder} --prefix {prefix} {ref_file}'\
            .format(prokka=prokka,
                    output_folder=reference_folder,
                    prefix=prefix,
                    ref_file=ref_file)
        # Update the issue with the GeneSeekr command
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Prokka command:\n {cmd}'.format(cmd=cmd))
        # Create another shell script to execute within the PlasmidExtractor conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, cmd)
        prokka_script = os.path.join(work_dir, 'run_prokka.sh')
        with open(prokka_script, 'w+') as file:
            file.write(template)
        make_executable(prokka_script)

        # Run shell script
        os.system(prokka_script)

        #
        seq_folder = os.path.join(work_dir, 'fastqs')
        # Now extract our query files.
        retrieve_nas_files(seqids=query_list,
                           outdir=seq_folder,
                           filetype='fastq',
                           copyflag=False)

        # With our query files extracted, verify that all the SEQIDs the user specified were able to be found.
        missing_fastqs = verify_fastqs_present(query_list, seq_folder)
        if len(missing_fastqs) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Warning! Could not find the following requested query SEQIDs: '
                                                '{}. \nYou may want to verify the SEQIDs, create a new issue, and try'
                                                ' again.'.format(str(missing_fastqs)))

        # Now check that the FASTQ files aren't too far away from the specified reference file.
        bad_fastqs = check_distances(ref_fasta=glob.glob(os.path.join(reference_folder, '*fasta'))[0],
                                     fastq_folder=seq_folder,
                                     work_dir=work_dir)
        if len(bad_fastqs) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Warning! The following SEQIDs were found to be fairly'
                                                ' divergent from the reference file specified:{} \nYou may'
                                                ' want to start a new COWSNPhR issue without them and try '
                                                'again.'.format(str(bad_fastqs)))

        # COWSNPhR
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/vsnp_dev'
        binary = '/mnt/nas2/virtual_environments/vSNP/cowsnphr_src/cowsnphr_src.py'

        # Prepare command
        cmd = '{bin} -s {seq_folder} -r {ref_folder} -w /mnt/nas2' \
            .format(bin=binary,
                    seq_folder=seq_folder,
                    ref_folder=reference_folder)
        # Update the issue with the GeneSeekr command
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='COWSNPhR command:\n {cmd}'.format(cmd=cmd))
        # Create another shell script to execute within the PlasmidExtractor conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, cmd)
        cowsnphr_script = os.path.join(work_dir, 'run_cowsnphr.sh')
        with open(cowsnphr_script, 'w+') as file:
            file.write(template)
        make_executable(cowsnphr_script)

        # Run shell script
        os.system(cowsnphr_script)

        # Zip output
        output_filename = 'cowsnphr_output_{}'.format(issue.id)
        zip_filepath = zip_folder(output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'
        #
        sas_url = upload_to_ftp(local_file=zip_filepath)
        # Prepare upload
        if sas_url:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='COWSNPhr process complete!\n\n'
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
        # Clean up files
        #shutil.rmtree(reference_folder)
        #shutil.rmtree(seq_folder)
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


def zip_folder(output_dir, output_filename):
    folder_base = ['alignments', 'summary_tables', 'tree_files', 'vcf_files']
    output_base = os.path.join(output_dir, 'fastqs', output_filename)
    with zipfile.ZipFile(output_base + '.zip', 'w', zipfile.ZIP_DEFLATED) as zippo:
        for folder in folder_base:
            to_zip = os.path.join(output_dir, 'fastqs', folder)
            for root, dirs, files in os.walk(to_zip):
                for filename in files:
                    dir_name = os.path.split(root)[-1]
                    zippo.write(os.path.join(root, filename), arcname=os.path.join(dir_name, filename))
    return output_base


def check_distances(ref_fasta, fastq_folder, work_dir):
    bad_fastqs = list()
    # fastqs = glob.glob(os.path.join(fastq_folder, '*R1*'))
    mash.sketch(os.path.join(fastq_folder, '*R1*'), output_sketch=os.path.join(work_dir, 'sketch.msh'), threads=5)
    mash.dist(os.path.join(work_dir, 'sketch.msh'), ref_fasta, threads=48,
              output_file=os.path.join(work_dir, 'distances.tab'))
    mash_output = mash.read_mash_output(os.path.join(work_dir, 'distances.tab'))
    for item in mash_output:
        if item.distance > 0.15:  # Moved value from 0.06 to 0.15 - was definitely too conservative before.
            bad_fastqs.append(item.reference)
    return bad_fastqs


def verify_fastqs_present(query_list, fastq_folder):
    missing_fastqs = list()
    for query in query_list:
        # Check that forward reads are present.
        if len(glob.glob(os.path.join(fastq_folder, query + '*R1*fastq*'))) == 0:
            missing_fastqs.append(query)
        # Check that reverse reads are present, and add to list if forward reads weren't missing
        if len(glob.glob(os.path.join(fastq_folder, query + '*R2*fastq*'))) == 0 and query not in missing_fastqs:
            missing_fastqs.append(query)
    # Returns list of SEQIDs for which we couldn't find forward and/or reverse reads
    return missing_fastqs


if __name__ == '__main__':
    cowsnphr_redmine()
