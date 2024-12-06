import os
import glob
import click
import pickle
import shutil
import sentry_sdk
from automator_settings import SENTRY_DSN
from amrsummary import before_send
from nastools.nastools import retrieve_nas_files
from externalretrieve import check_fastas_present
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
def psortb_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))
        # Parse description to get list of SeqIDs
        seqids = list()
        for i in range(0, len(description)):
            item = description[i].rstrip()
            item = item.upper()
            seqids.append(item)

        assemblies_folder = os.path.join(work_dir, 'assemblies')
        retrieve_nas_files(seqids=seqids,
                           outdir=assemblies_folder,
                           filetype='fasta',
                           copyflag=False)
        missing_fastas = check_fastas_present(seqids, assemblies_folder)
        if len(missing_fastas) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested FASTA SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastas))
        for fasta in missing_fastas:
            seqids.remove(fasta)
        # Steps to follow here:
        # 1) Generate protein file for sequence(s) of interest - run prokka so proteins get named nicely.

        prokka_folder = os.path.join(work_dir, 'prokka')
        os.makedirs(prokka_folder)
        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/prokka'
        prokka = '/mnt/nas2/virtual_environments/prokka/bin/prokka'

        for assembly in glob.glob(os.path.join(assemblies_folder, '*.fasta')):
            seqid = os.path.split(assembly)[1].split('.')[0]
            # Prepare command
            cmd = '{prokka} --outdir {output_folder} --prefix {seqid} {assembly}'.format(prokka=prokka,
                                                                                         output_folder=os.path.join(prokka_folder, seqid),
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
        # 2) Figure out for each sequence if gram positive or negative
        # Psort says using Omp85 works pretty well for determining. Use Omp85 proteins from Neisseria, Thermosipho,
        # Synechoccus, and Thermus. If any hits with e value < 10^-3 gram positive, otherwise gram negative. There
        # are exceptions to this, but I don't think any of them are organisms that we work with/care about
        gram_pos_neg_dict = dict()
        protein_files = glob.glob(os.path.join(prokka_folder, '*', '*.faa'))
        for protein_file in protein_files:
            seqid = os.path.split(protein_file)[1].replace('.faa', '')
            # Make a blast DB from proteins.
            cmd = 'makeblastdb -in {} -dbtype prot'.format(protein_file)
            os.system(cmd)
            # Now BLAST our OMP85 proteins against the genome proteins.
            blast_result_file = protein_file.replace('.faa', '_blast.tsv')
            omp85_proteins = '/mnt/nas2/redmine/applications/OLCRedmineAutomator/data_and_stuff/omp85_proteins.fasta'
            cmd = 'blastp -db {} -query {} -out {} -outfmt "6 qseqid sseqid evalue"'.format(protein_file,
                                                                                            omp85_proteins,
                                                                                            blast_result_file)
            os.system(cmd)
            # Now parse through blast report to find if gram positive or negative
            has_omp_85 = False
            with open(blast_result_file) as f:
                for line in f:
                    evalue = float(line.rstrip().split()[-1])
                    if evalue < 0.0001:
                        has_omp_85 = True

            if has_omp_85 is True:
                gram_pos_neg_dict[seqid] = 'Negative'
            else:
                gram_pos_neg_dict[seqid] = 'Positive'

        """
        IMPORTANT NOTES ON GETTING PSORTB TO RUN:
        You'll need to 
        1) have pulled a docker image of PSORTB to each of the 
        nodes (docker pull brinkmanlab/psortb_commandline:1.0.2 should do the trick)
        2) Put a psortb executable into the data_and_stuff folder: 
        wget -O data_and_stuff/psortb https://raw.githubusercontent.com/brinkmanlab/psortb_commandline_docker/master/psortb
        3) chmod the psortb executable to actually make it executable.
        4) remove the 'sudo' from lines 35 and 61 of the psortb executable, otherwise nodes get unhappy, and make the -it
        in the commands into just a -i
        """

        # 3) Run PsortB!
        for seqid in seqids:
            protein_file = os.path.join(prokka_folder, seqid, seqid + '.faa')
            output_dir = os.path.join(prokka_folder, seqid)
            psortb_executable = '/mnt/nas2/redmine/applications/OLCRedmineAutomator/data_and_stuff/psortb'
            cmd = '{} -i {} -r {} '.format(psortb_executable, protein_file, output_dir)
            if gram_pos_neg_dict[seqid] == 'Negative':
                cmd += '--negative'
            else:
                cmd += '--positive'
            os.system(cmd)

        # Now need to: upload results, do file cleanup.
        report_dir = os.path.join(work_dir, 'psortb_reports_{}'.format(issue.id))
        os.makedirs(report_dir)
        raw_reports = glob.glob(os.path.join(prokka_folder, '*', '*psortb*.txt'))
        for raw_report in raw_reports:
            new_name = raw_report.split('/')[-2] + '_' + os.path.split(raw_report)[1]
            cmd = 'cp {} {}'.format(raw_report, os.path.join(report_dir, new_name))
            os.system(cmd)
        cmd = 'cp {} {}'.format(os.path.join(prokka_folder, '*', '*.faa'), report_dir)
        os.system(cmd)

        shutil.make_archive(report_dir, 'zip', report_dir)
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
                status_id=4,
                notes='PsortB analysis complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=download_link)
            )
        else:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='Upload of results was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )
        
        shutil.rmtree(assemblies_folder)
        shutil.rmtree(prokka_folder)
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


def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    psortb_redmine()
