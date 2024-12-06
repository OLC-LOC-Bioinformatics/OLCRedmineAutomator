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
def mob_suite(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    """
    """
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        # Description should just be a list of SEQIDs. Get the fasta files associated with them extracted
        # to the bio_request dir
        retrieve_nas_files(seqids=description,
                           outdir=os.path.join(work_dir, 'fastas'),
                           filetype='fasta',
                           copyflag=True)  # Since we're docker-ing need to copy. Files get cleaned up at end of process
        # Now we need to run mob_recon (and typing!) on each of the fasta files requested. Put all results into one
        # folder (this will need to be uploaded to FTP - will overwhelm max (10MB) file size limit on Redmine

        fasta_files = glob.glob(os.path.join(work_dir, 'fastas', '*.fasta'))
        # Verify that specified fasta files are actually there, warn user if they aren't.
        missing_fastas = verify_fasta_files_present(seqid_list=description,
                                                    fasta_dir=os.path.join(work_dir, 'fastas'))
        if len(missing_fastas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))

        # Make output dir
        output_dir = os.path.join(work_dir, 'mob_suite_results')
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        for fasta in fasta_files:
            seqid = os.path.split(fasta)[-1].split('.')[0]
            # Run mobsuite via docker, since I can't seem to make it work with slurm any other way.
            # changed the virtual environment from /mob_suite to /dev/cowbat on April 30, 2021 as this contains mobsuite version 3.0.0
            cmd = 'docker run --rm -i -u $(id -u) -v /mnt/nas2:/mnt/nas2 mob_suite:latest /bin/bash -c "source activate ' \
                  #'/mnt/nas2/virtual_environments/dev/cowbat && mob_recon -i {input_fasta} -o {output_dir} ' \
                  '/mnt/nas2/virtual_environments/mobsuite3 && mob_recon -i {input_fasta} -s {seqid} -p {seqid} -o {output_dir} ' \
                  '--run_typer"'.format(input_fasta=fasta,seqid=seqid,
                                        output_dir=os.path.join(output_dir, seqid))
            os.system(cmd)

        # With mobsuite done, zip up the results folder and upload to the FTP.
        shutil.make_archive(root_dir=output_dir,
                            format='zip',
                            base_name=os.path.join(work_dir, str(issue.id)))

        sas_url = upload_to_ftp(
            local_file=os.path.join(work_dir, str(issue.id) + '.zip')
        )

        # And finally, do some file cleanup.
        try:
            shutil.rmtree(output_dir)
            shutil.rmtree(os.path.join(work_dir, 'fastas'))
            os.remove(os.path.join(work_dir, str(issue.id) + '.zip'))
        except:
            pass

        if sas_url:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='Mob-suite process complete!\n\n'
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

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')


def verify_fasta_files_present(seqid_list, fasta_dir):
    missing_fastas = list()
    for seqid in seqid_list:
        if len(glob.glob(os.path.join(fasta_dir, seqid + '*.fasta'))) == 0:
            missing_fastas.append(seqid)
    return missing_fastas


if __name__ == '__main__':
    mob_suite()
