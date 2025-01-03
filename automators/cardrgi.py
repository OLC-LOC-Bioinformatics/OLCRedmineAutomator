import os
import glob
import click
import pickle
import shutil
import subprocess
from biotools import mash
from nastools.nastools import retrieve_nas_files
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
import csv
import pandas
from pathlib import Path
#from openpyxl import load_workbook


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def cardrgi_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        #analyses supported by the automator
        analysistypes = [
            'isolate', 'metagenome'
        ]


        # Variable to hold supplied arguments
        argument_dict = {
            'analysistype': 'isolate',
            'fastaattached': False,
            'loosehits': False,
            'partialgenes': False,
        }


        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        analysistype = "" # Just for the output at the end
        for item in description:
            item = item.rstrip()
            if 'analysistype' in item:
                argument_dict['analysistype'] = item.split('=')[1].lower()
                continue
            if 'loosehits' in item:
                argument_dict['loosehits'] = True
                continue
            if 'fastaattached' in item:
                argument_dict['fastaattached'] = True
                continue
            if 'partialgenes' in item:
                argument_dict['partialgenes'] = True
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        
        # Ensure that SEQIDs were included/... now that we are added the attachment option, need to remove this
        #if len(seqids) == 0:
        #    redmine_instance.issue.update(resource_id=issue.id, status_id=4,
        #                                  notes='You did not include any SEQIDs.')
        #    return


        # Create folder to drop FASTA files
        sequences_folder = os.path.join(work_dir, 'sequences')
        os.makedirs(sequences_folder, exist_ok=True)

        # Create cardrgi output folder
        cardrgi_folder = os.path.join(work_dir, 'card_output')
        os.makedirs(cardrgi_folder)

        #if using CARD for isolates, extract fasta files:
        if argument_dict['analysistype'] == 'isolate':
            # Extract FASTA files.
            retrieve_nas_files(seqids=seqids, outdir=sequences_folder, filetype='fasta', copyflag=False)
            missing_fastas = verify_fasta_files_present(seqids, sequences_folder)
            if missing_fastas:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested SEQIDs on '
                                                    'the OLC NAS: {}'.format(missing_fastas))

            # Download the attached files if included
            if argument_dict['fastaattached'] == True:
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
                        attachment.download(savepath=sequences_folder)
                #NCBI sequences are often .fna... lets rename the file(s)
                for f in glob.glob(os.path.join(sequences_folder, '*.fna')):
                    #Split the file into the filename and the extension saving as separate variables
                    filename, ext=os.path.splittext(f)
                    if "." in filename:
                        #If a '.' is in the name, rename, appending the suffix (ext) to the new file
                        new_name = filename.replace(".","_")
                        os.rename(
                            os.path.join(sequences_folder, f),
                            os.path.join(sequences_folder, new_name + ext))
                    #rename the suffix to .fasta
                    f.rename(f.with_suffix('.fasta'))

        #if metagenome is chosen, pull fastq files
        if argument_dict['analysistype'] == 'metagenome':
            # Extract FASTQ files.
            retrieve_nas_files(seqids=seqids, outdir=sequences_folder, filetype='fastq', copyflag=False)
            missing_fastqs = check_fastqs_present(seqids, sequences_folder)
            if len(missing_fastqs) > 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastqs))

        #activating the correct environment
        activatecardrgi = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/card-rgi-6'


        #if running for isolates
        if argument_dict['analysistype'] == 'isolate':
            for assembly in glob.glob(os.path.join(sequences_folder, '*.fasta')):
                seqid = os.path.split(assembly)[1].split('.')[0]
                #prepare command
                cardrgiload_cmd = 'rgi load -i /mnt/nas2/databases/card-rgi-6/card.json --card_annotation /mnt/nas2/databases/card-rgi-6/card_database_v3.2.5.fasta --local'
                outputfile = '{seqid}_rgiresults'.format(seqid=seqid)
                cardrgi_cmd = 'rgi main -i {assembly} -o {outfile} -n 16 --local'.format(assembly=assembly,
                                                                                         outfile=outputfile)
                #Append other flags if required
                cardrgi_cmd += ' --include-loose' if argument_dict['loosehits'] else ''
                cardrgi_cmd += ' --low_quality' if argument_dict['partialgenes'] else ''

                # Create another shell script to execute within the pyseer conda environment
                templatergi = "#!/bin/bash\n{activate} && cd {seqdir} && {cardload} && {cardrgi} --clean".format(activate=activatecardrgi,
                                                                                                                 seqdir=sequences_folder,
                                                                                                                 cardload=cardrgiload_cmd,
                                                                                                                 cardrgi=cardrgi_cmd)
                cardrgi_script = os.path.join(work_dir, 'run_card-rgi.sh')
                with open(cardrgi_script, 'w+') as file:
                    file.write(templatergi)
                make_executable(cardrgi_script)

                # Run shell script
                os.system(cardrgi_script)

            #move all card-rgi files to the output folder
            json_folder = os.path.join(cardrgi_folder, 'json_files')
            os.makedirs(json_folder)
            #now move the files
            output_files = os.listdir(sequences_folder)
            for file in output_files:
                if file.endswith(".txt"):
                    shutil.move(os.path.join(sequences_folder,file), os.path.join(cardrgi_folder, file))
            for file in output_files:
                if file.endswith(".json"):
                    shutil.move(os.path.join(sequences_folder,file), os.path.join(json_folder, file))

            #adding in the heatmap just in case users want it
            rgiheatmap_cmd = 'rgi heatmap -i {jsonfolder}'.format(jsonfolder=json_folder)
            # Create another shell script to execute within the pyseer conda environment
            templatergihmp = "#!/bin/bash\n{activate} && cd {outfolder} && {heatmap}".format(activate=activatecardrgi, outfolder=cardrgi_folder, heatmap=rgiheatmap_cmd)
                
            cardhmp_script = os.path.join(work_dir, 'run_card-rgi-heatmap.sh')
            with open(cardhmp_script, 'w+') as file:
                file.write(templatergihmp)
            make_executable(cardhmp_script)

            # Run shell script
            os.system(cardhmp_script)
            
            #add the filename (which is the seqid) to the first column in all of the card report files, then concatenate into a single output csv file
            outputfile = os.path.join(cardrgi_folder, 'CARDRGI_output.csv')
            with open(outputfile, 'w', newline='') as file_output:
                csv_output = csv.writer(file_output, delimiter='\t')
                #for fname in glob.glob(os.path.basename(seq_dir, '*.res')):
                for fname in glob.glob(os.path.join(cardrgi_folder, '*.txt')):
                    fbasename = os.path.basename(fname) #this is to just get the seqid.txt name of file
                    seqname1 = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
                    seqname = os.path.split(seqname1)[1].split('_r')[0] #this is to just pull out the seqid
                    with open(fname, newline='') as f_input:
                        csv_input = csv.reader(f_input, delimiter='\t')
                        #Header processing
                        header = csv_input.__next__()
                        header.insert(0,"SeqID")

                        csv_output.writerow(header)
                        for row in csv_input:
                            print(row)
                            row.insert(0,seqname) #this adds the seqid to the file before concatenating
                            #row.insert(0,fname)
                            csv_output.writerow(row)
            #write the csv file into an excel file so our window's users have easier viewing
            carddf = pandas.read_csv(os.path.join(cardrgi_folder, 'CARDRGI_output.csv'), sep='\t')
            #write it to output xlsx file
            cardexcelfile = os.path.join(cardrgi_folder, 'CARD_RGI_results_{}.xlsx'.format(issue.id))
            carddf.to_excel(cardexcelfile, sheet_name='Isolate_resistance', index=False)
            #remove the csv file
            os.remove(os.path.join(cardrgi_folder, 'CARDRGI_output.csv'))
            

        #if running for metagenomes
        if argument_dict['analysistype'] == 'metagenome':
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
                print(forwardseq)
                reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
                #prepare command
                cardrgiload_cmd = 'rgi load -i /mnt/nas2/databases/card-rgi-6/card.json --card_annotation /mnt/nas2/databases/card-rgi-6/card_database_v3.2.5.fasta --local'
                outputfile = '{seqid}_rgiresults'.format(seqid=seqid)
                cardrgi_cmd = 'rgi bwt -1 {forward} -2 {reverse} -a kma -o {outfile} -n 16 --local'.format(forward=forwardseq,reverse=reverseseq,
                                                                                                           outfile=outputfile)
                #Append other flags if required
                cardrgi_cmd += ' --include-loose' if argument_dict['loosehits'] else ''
                cardrgi_cmd += ' --low_quality' if argument_dict['partialgenes'] else ''

                # Create another shell script to execute within the pyseer conda environment
                templatergi = "#!/bin/bash\n{activate} && cd {seqdir} && {cardload} && {cardrgi} --clean".format(activate=activatecardrgi,
                                                                                                                 seqdir=sequences_folder,
                                                                                                                 cardload=cardrgiload_cmd,
                                                                                                                 cardrgi=cardrgi_cmd)
                cardrgi_script = os.path.join(work_dir, 'run_card-rgi.sh')
                with open(cardrgi_script, 'w+') as file:
                    file.write(templatergi)
                make_executable(cardrgi_script)

                # Run shell script
                os.system(cardrgi_script)

            #move all card-rgi files to the output folder
            json_folder = os.path.join(cardrgi_folder, 'json_files')
            os.makedirs(json_folder)
            #now move the files
            output_files = os.listdir(sequences_folder)
            for file in output_files:
                if file.endswith(".txt"):
                    shutil.move(os.path.join(sequences_folder,file), os.path.join(cardrgi_folder, file))
            for file in output_files:
                if file.endswith(".json"):
                    shutil.move(os.path.join(sequences_folder,file), os.path.join(json_folder, file))

            #add the filename (which is the seqid) to the first column in all of the report files, then concatenate into a single output csv file
            outputfile = os.path.join(cardrgi_folder, 'CARDRGI_gene_mapping_output.csv')
            with open(outputfile, 'w', newline='') as file_output:
                csv_output = csv.writer(file_output, delimiter='\t')
                #for fname in glob.glob(os.path.basename(seq_dir, '*.res')):
                for fname in glob.glob(os.path.join(cardrgi_folder, '*.gene_mapping_data.txt')):
                    fbasename = os.path.basename(fname) #this is to just get the seqid.txt name of file
                    seqname1 = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
                    seqname = os.path.split(seqname1)[1].split('_r')[0] #this is to just pull out the seqid
                    with open(fname, newline='') as f_input:
                        csv_input = csv.reader(f_input, delimiter='\t')
                        #Header processing
                        header = csv_input.__next__()
                        header.insert(0,"SeqID")

                        csv_output.writerow(header)
                        for row in csv_input:
                            print(row)
                            row.insert(0,seqname) #this adds the seqid to the file before concatenating
                            #row.insert(0,fname)
                            csv_output.writerow(row)

            #add the filename (which is the seqid) to the first column in all of the res files, then concatenate into a single output csv file
            outputfile2 = os.path.join(cardrgi_folder, 'CARDRGI_allele_mapping_output.csv')
            with open(outputfile2, 'w', newline='') as file_output:
                csv_output = csv.writer(file_output, delimiter='\t')
                #for fname in glob.glob(os.path.basename(seq_dir, '*.res')):
                for fname in glob.glob(os.path.join(cardrgi_folder, '*.allele_mapping_data.txt')):
                    fbasename = os.path.basename(fname) #this is to just get the seqid.txt name of file
                    seqname1 = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
                    seqname = os.path.split(seqname1)[1].split('_r')[0] #this is to just pull out the seqid
                    with open(fname, newline='') as f_input:
                        csv_input = csv.reader(f_input, delimiter='\t')
                        #Header processing
                        header = csv_input.__next__()
                        header.insert(0,"SeqID")

                        csv_output.writerow(header)
                        for row in csv_input:
                            print(row)
                            row.insert(0,seqname) #this adds the seqid to the file before concatenating
                            #row.insert(0,fname)
                            csv_output.writerow(row)

            #write the csv files into an excel file so our window's users have easier viewing
            cardexcelfile = os.path.join(cardrgi_folder, 'CARD_RGI_results_{}.xlsx'.format(issue.id))            
            carddf = pandas.read_csv(os.path.join(cardrgi_folder, 'CARDRGI_gene_mapping_output.csv'), sep='\t')
            cardalleledf = pandas.read_csv(os.path.join(cardrgi_folder, 'CARDRGI_allele_mapping_output.csv'), sep='\t')
            #write them to output xlsx file
            with pandas.ExcelWriter(cardexcelfile) as writer:
                carddf.to_excel(cardexcelfile, sheet_name='Gene_mapping_results', index=False)
                cardalleledf.to_excel(cardexcelfile, sheet_name='Allele_mapping_results', index=False)
            #remove the csv files
            os.remove(os.path.join(cardrgi_folder, 'CARDRGI_gene_mapping_output.csv'))
            os.remove(os.path.join(cardrgi_folder, 'CARDRGI_allele_mapping_output.csv'))

        # Zip card-rgi output
        output_filename = 'card-rgi_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=cardrgi_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'

        #Get the outputs uploaded
        output_list = list()
        output_dict = dict()
        #upload the CARDRGI excel file
        cardxlsxfilename = 'CARD_RGI_results_{}.xlsx'.format(issue.id)
        output_dict['path'] = os.path.join(cardrgi_folder,cardxlsxfilename)
        output_dict['filename'] = cardxlsxfilename
        output_list.append(output_dict)
        #and the zipfile
        output_dict = dict()
        output_dict['path'] = os.path.join(work_dir,output_filename)
        output_dict['filename'] = output_filename
        output_list.append(output_dict)


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
                notes='CARD-RGI analysis complete!\n\n'
                      'CARD-RGI version 6 with database v3.2.5\n\n'
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


        #Upload attachments to the redmine issue    
        redmine_instance.issue.update(resource_id=issue.id, 
                                      uploads=output_list, status_id=4,
                                      notes='CARD-RGI analysis complete!\n\n'
                                            'CARD-RGI version 6 with database v3.2.5\n'
                                            'Results are attached.'
                                     )
        #os.remove(zip_filepath)
        #remove the sequences folder
        shutil.rmtree(sequences_folder)
        #remove the output folder
        shutil.rmtree(cardrgi_folder)

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
    cardrgi_redmine()
