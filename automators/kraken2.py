import os
import glob
import click
import pickle
import shutil
import subprocess
from biotools import mash
from nastools.nastools import retrieve_nas_files
from externalretrieve import upload_to_ftp
import csv


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def kraken2_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        #analyses supported by the automator
        seqtypes = [
            'paired', 'nanopore', 'assembly',
        ]

        #databases supported by the automator
        databases = [
            'kraken2', 'greengenes','plusPF','rdp','silva',
        ]

        #bracken filter levels supported by the automator
        brackenfilters = [
            'S', 'G','F','O','C','P','K',
        ]

        # Variable to hold supplied arguments
        argument_dict = {
            'seqtype': 'paired',
            'database': 'kraken2',
            'brackenfilter': 'S',
        }


        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        analysistype = "" # Just for the output at the end
        for item in description:
            item = item.rstrip()
            if 'seqtype' in item:
                argument_dict['seqtype'] = item.split('=')[1].lower()
                continue
            if 'database' in item:
                argument_dict['database'] = item.split('=')[1]
                continue
            if 'brackenfilter' in item:
                argument_dict['brackenfilter'] = item.split('=')[1].upper()
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        
        # Ensure that SEQIDs were included
        if len(seqids) == 0:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='You did not include any SEQIDs.')
            return

        # Download the attached files.
        # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
        # out a better way to do it.
        attachment = redmine_instance.issue.get(issue.id, include='attachments')
        attachment_id = 0
        for item in attachment.attachments:
            attachment_id = item.id
        # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
        if attachment_id != 0:
            for item in attachment.attachments:
                attachment_id = item.id
                attachment = redmine_instance.attachment.get(attachment_id)
#            attachment.download(savepath=target_dir, filename='traits.tsv')
#                attachment.download(savepath=work_dir, filename=attachment_id)
                attachment.download(savepath=work_dir)

        # Create folder to drop FASTQ files
        sequences_folder = os.path.join(work_dir, 'sequences')
        os.makedirs(sequences_folder, exist_ok=True)

        # Create Kraken2 output folder
        kraken2_folder = os.path.join(work_dir, 'kraken2_output')
        os.makedirs(kraken2_folder)

        #pull out the fastq files
        # Extract FASTQ files.
        retrieve_nas_files(seqids=seqids, outdir=sequences_folder, filetype='fastq', copyflag=False)
        missing_fastqs = check_fastqs_present(seqids, sequences_folder)
        if len(missing_fastqs) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastqs))
        #set the database path for the analyses
        dbpath = '/mnt/nas2/databases/kraken2/'
        database_path = {
            'kraken2': os.path.join(dbpath, 'k2_standard_20230605'),
            'greengenes': os.path.join(dbpath, 'greengenes'),
            'plusPF': os.path.join(dbpath, 'plusPF_20230605'),
            'rdp': os.path.join(dbpath, 'rdp'),
            'silva': os.path.join(dbpath, 'silva')
        }

        #Commands required for kraken2
        activateckraken2 = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kraken2'

        #if paired end is chosen, run analysis
        if argument_dict['seqtype'] == 'paired':
            # Rename the files to be consistent with the desired SRA naming scheme
            # e.g. 2014-SEQ-0349_S11_L001_R1_001.fastq.gz renamed to 2014-SEQ-0349_R1.fastq.gz
            #create the system call
            #rename_cmd = "cd {seqs} && rename 's/_S\d+_L001_R(\d)_001/_R$1/' *.gz".format(seqs=sequences_folder)
            #run the command
            #os.system(rename_cmd)
            #now run the analysis
            #do not want to repeat twice, so added in the _R1
            for metagenome in glob.glob(os.path.join(sequences_folder, '*_R1_001.fastq.gz')):
                seqid1 = os.path.split(metagenome)[1].split('.')[0]
                seqid = os.path.split(seqid1)[1].split('_R')[0]
                forwardseq = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
                reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
                #prepare command
                classified_out = '{seqid}#_classified'.format(seqid=seqid)
                unclassified_out = '{seqid}#_unclassified'.format(seqid=seqid)
                report_out = '{seqid}_{db}_report'.format(seqid=seqid, db=argument_dict['database'])
                kraken2_cmd = 'kraken2 --db {database} --threads 24 --paired --output {outdir} --gzip-compressed'\
                              ' --classified-out {outdir}/{outfilec} --unclassified-out {outdir}/{outfileu}'\
                              ' --report {outdir}/{reportout} --report-zero-counts'\
                              ' {forwardseq} {reverseseq}'.format(database=database_path[argument_dict['database']],
                                                                  outdir=kraken2_folder,
                                                                  outfilec=classified_out, outfileu=unclassified_out,
                                                                  reportout=report_out,
                                                                  forwardseq=forwardseq, reverseseq=reverseseq)

                # Create another shell script to execute within the pyseer conda environment
                template = "#!/bin/bash\n{ntasks}\n{mem}\n{activate} && cd {seqdir} && {kraken2}"\
                    .format(ntasks="#SBATCH --ntasks 30",mem="#SBATCH --mem=190000",
                            activate=activateckraken2,seqdir=sequences_folder,
                            kraken2=kraken2_cmd)
                kraken2_script = os.path.join(work_dir, 'run_kraken2.sh')
                with open(kraken2_script, 'w+') as file:
                    file.write(template)
                make_executable(kraken2_script)

                # Run shell script
                os.system(kraken2_script)

            for metagenome in glob.glob(os.path.join(sequences_folder, '*_R1.fastq.gz')):
                seqid1 = os.path.split(metagenome)[1].split('.')[0]
                seqid = os.path.split(seqid1)[1].split('_R')[0]
                forwardseq = '{seqid}_R1.fastq.gz'.format(seqid=seqid)
                reverseseq = '{seqid}_R2.fastq.gz'.format(seqid=seqid)
                #prepare command
                classified_out = '{seqid}#_classified'.format(seqid=seqid)
                unclassified_out = '{seqid}#_unclassified'.format(seqid=seqid)
                report_out = '{seqid}_{db}_report'.format(seqid=seqid, db=argument_dict['database'])
                kraken2_cmd = 'kraken2 --db {database} --threads 24 --paired --output {outdir} --gzip-compressed'\
                              ' --classified-out {outdir}/{outfilec} --unclassified-out {outdir}/{outfileu}'\
                              ' --report {outdir}/{reportout} --report-zero-counts'\
                              ' {forwardseq} {reverseseq}'.format(database=database_path[argument_dict['database']],
                                                                  outdir=kraken2_folder,
                                                                  outfilec=classified_out, outfileu=unclassified_out,
                                                                  reportout=report_out,
                                                                  forwardseq=forwardseq, reverseseq=reverseseq)

                # Create another shell script to execute within the pyseer conda environment
                template = "#!/bin/bash\n{ntasks}\n{mem}\n{activate} && cd {seqdir} && {kraken2}"\
                    .format(ntasks="#SBATCH --ntasks 30",mem="#SBATCH --mem=190000",
                            activate=activateckraken2,seqdir=sequences_folder,
                            kraken2=kraken2_cmd)
                kraken2_script = os.path.join(work_dir, 'run_kraken2.sh')
                with open(kraken2_script, 'w+') as file:
                    file.write(template)
                make_executable(kraken2_script)

                # Run shell script
                os.system(kraken2_script)

            #now move the files
            output_files = os.listdir(kraken2_folder)
            # Create reports, classified and unclassified folders in the output folder
            reports_folder = os.path.join(kraken2_folder, 'reports')
            os.makedirs(reports_folder)
            classified_folder = os.path.join(kraken2_folder, 'classified')
            os.makedirs(classified_folder)
            unclassified_folder = os.path.join(kraken2_folder, 'unclassified')
            os.makedirs(unclassified_folder)
            for file in output_files:
                if file.endswith("_report"):
                    shutil.move(os.path.join(kraken2_folder,file), os.path.join(reports_folder, file))
            for file in output_files:
                if file.endswith("_classified"):
                    shutil.move(os.path.join(kraken2_folder, file), os.path.join(classified_folder, file))
            for file in output_files:
                if file.endswith("_unclassified"):
                    shutil.move(os.path.join(kraken2_folder,file), os.path.join(unclassified_folder, file))

        #run bracken on reports
        for report in glob.glob(os.path.join(reports_folder, '*_report')):
            reportname = os.path.split(report)[1].split('_re')[0]
            outname = '{reportname}_bracken'.format(reportname=reportname)
            activatebracken = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/bracken'
            bracken_cmd = 'bracken -d {db} -i {report} -o {outname} -l {cat}'.format(db=database_path[argument_dict['database']],
                                                                                     report=report, outname=outname,
                                                                                     cat=argument_dict['brackenfilter'])
            # Create another shell script to execute within the pyseer conda environment
            templatebr = "#!/bin/bash\n{activate} && cd {rdir} && {bracken}"\
                .format(activate=activatebracken, rdir=reports_folder, bracken=bracken_cmd)
            bracken_script = os.path.join(work_dir, 'run_bracken.sh')
            with open(bracken_script, 'w+') as file:
                file.write(templatebr)
            make_executable(bracken_script)

            # Run shell script
            os.system(bracken_script)

        #now move the files
        output_files2 = os.listdir(reports_folder)
        # Create reports, classified and unclassified folders in the output folder
        bracken_folder = os.path.join(kraken2_folder, 'bracken_reports')
        os.makedirs(bracken_folder)
        bracken_fileinfo = os.path.join(kraken2_folder, 'bracken_fileinfo')
        os.makedirs(bracken_fileinfo)
        for file in output_files2:
            if file.endswith("_bracken"):
                shutil.move(os.path.join(reports_folder,file), os.path.join(bracken_fileinfo, file))
        for file in output_files2:
            if file.endswith("_bracken_*"):
                shutil.move(os.path.join(reports_folder,file), os.path.join(bracken_folder, file))

        # Zip kraken2 output
        output_filename = 'kraken2_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=kraken2_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'

        sas_url = upload_to_ftp(local_file=zip_filepath)
        # Prepare upload
        if sas_url:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='Kraken2 analysis complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=sas_url)
                )
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity '
                                                'issues. '
                                                'Please try again later.')
        # Remove the zip file
        #os.remove(zip_filepath)
        #remove the sequences folder
        #shutil.rmtree(sequences_folder)
        #remove the output folder
        #shutil.rmtree(kraken2_folder)

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

def check_fastqs_present(fastq_list, fastq_dir):
    missing_fastqs = list()
    for seqid in fastq_list:
        if len(glob.glob(os.path.join(fastq_dir, seqid + '*.fastq.gz'))) < 2:
            # JAS adding an extra if statement here to allow for Nanopore (SE) reads
            if len(glob.glob(os.path.join(fastq_dir, seqid + ".fastq.gz"))) == 0:
                missing_fastqs.append(seqid)
    return missing_fastqs


def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    kraken2_redmine()
