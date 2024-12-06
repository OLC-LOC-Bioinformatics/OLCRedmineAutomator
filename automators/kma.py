import os
import glob
import click
import pickle
import shutil
from biotools import mash
from amrsummary import before_send
import sentry_sdk
from automator_settings import SENTRY_DSN
from nastools.nastools import retrieve_nas_files
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
import fileinput
from pathlib import Path
import csv
import pandas




@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def kma_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Current list of analysis types that KMA can use
    analyses = [
        'custom', 'amr', 'biocide', 'metal', 'bacmet', 'all_ices', 'aice', 'cime', 'ime', 't4ss','verotoxin'
    ]
    # Current list of sequence types that KMA can analyse
    seqtypes = [
        'fasta', 'fastq', 'minionfastq', 'singlefastq'
    ]

    # Variable to hold supplied arguments
    argument_dict = {
        'analysis': str(),
        'targetsfile': False,
        'seqtype': 'fasta',
        'nanopore': False,
        'readcount': False,
        'pairedmethod': False,
        'vcf': False,
        'align': False,
        'hmm': False,
        'min_ID': False,
    }


    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            #item = item.upper().rstrip()
            item = item.rstrip()
            if 'analysis' in item:
                argument_dict['analysis'] = item.split('=')[1].lower()
                continue
            if 'targetsfile' in item:
                argument_dict['targetsfile'] = item.split('=')[1]
                continue
            if 'seqtype' in item:
                argument_dict['seqtype'] = item.split('=')[1].lower()
                continue
            if 'nanopore' in item:
                argument_dict['nanopore'] = True
                continue
            if 'readcount' in item:
                argument_dict['readcount'] = True
                continue
            if 'vcf' in item:
                argument_dict['vcf'] = True
                continue
            if 'align' in item:
                argument_dict['align'] = True
                continue
            if 'pairedmethod' in item:
                argument_dict['pairedmethod'] = True
                continue
            if 'hmm' in item:
                argument_dict['hmm'] = True
                continue
            if 'min_ID' in item:
                argument_dict['min_ID'] = int(item.split('=')[1].lower())
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
                                          notes='WARNING: supplied analysis type {at} currently not in the supported '
                                                'list of analyses: {ats}'.format(at=argument_dict['analysis'],
                                                                                 ats=', '.join(analyses)),
                                          status_id=4)
            return

        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No SEQIDs provided!',
                                          status_id=4)
            return

        # Set the database path for the analyses
        dbpath = '/mnt/nas2/databases/kma_v1.4.9/'
        #bacmetpath = '/mnt/nas2/databases/bacmet/'
        icebergpath = '/mnt/nas2/databases/ICEberg/'
        database_path = {
        #    'custom': os.path.join(target_dir, 'targets_KMA'),
            #'amr': os.path.join(dbpath, 'NCBI-AMR-v3.10'),
            #'biocide': os.path.join(dbpath, 'NCBI-biocide-v3.10'),
            #'metal': os.path.join(dbpath, 'NCBI-metal-v3.10'),
            'amr': os.path.join(dbpath, 'NCBI-AMR'),
            'biocide': os.path.join(dbpath, 'NCBI-BIOCIDE'),
            'metal': os.path.join(dbpath, 'NCBI-METAL'),
            'verotoxin': os.path.join(dbpath, 'VEROTOXIN'),
            'bacmet': os.path.join(dbpath, 'bacmet_2018-03-11-v2_renamed_kma'),
            'all_ices': os.path.join(icebergpath, 'ICEs_kmav1.4.9'),
            'aice': os.path.join(icebergpath, 'AICEs_kmav1.4.9'),
            'cime': os.path.join(icebergpath, 'CIMEs_kmav1.4.9'),
            'ime': os.path.join(icebergpath, 'IMEs_kmav1.4.9'),
            't4ss': os.path.join(icebergpath, 'T4SS_kmav1.4.9'),
        }

        #store the custom target file and index it
        if argument_dict['analysis'] == 'custom':
            #make sure the user provided a name for their target file
            if not argument_dict['targetsfile']:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Name of targets fasta file not provided. '
                                                    'If requesting a custom analysis, please ensure that the name of your targets.fasta '
                                                    ' file is provided (including the suffix .fasta). E.g. genes.fasta, telluriumtarget.fasta',
                                              status_id=4)
                return
            # Set and create the directory to store the custom targets
            target_dir = os.path.join(work_dir, 'targets')
            os.makedirs(target_dir, exist_ok=True)
            # Download the attached FASTA file.
            # First, get the attachment id - this seems like a kind of hacky way to do this, but I have yet to figure
            # out a better way to do it.
            attachment = redmine_instance.issue.get(issue.id, include='attachments')
            attachment_id = 0
            for item in attachment.attachments:
                attachment_id = item.id
            # Download if attachment id is not 0, which indicates that we didn't find anything attached to the issue.
            if attachment_id != 0:
                attachment = redmine_instance.attachment.get(attachment_id)
                attachment.download(savepath=target_dir, filename='{}'.format(argument_dict['targetsfile']))
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Analysis type custom requires an attached FASTA file of '
                                                    'targets. The automator could not find any attached files. '
                                                    'Please create a new issue with the .fasta file attached and try '
                                                    'again.',
                                              status_id=4)

            #If the user uploads a custom database, index it using kma_index
            activatei = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kma_v149'
            kmaindx_cmd = 'kma_index -i {targets} -o targets_KMA'.format(targets=argument_dict['targetsfile'])

            # Create another shell script to execute within the KMA conda environment
            templatei = "#!/bin/bash\n{} && cd {} && {}".format(activatei, target_dir, kmaindx_cmd)
            kmaindex_script = os.path.join(work_dir, 'run_kmaindex.sh')
            with open(kmaindex_script, 'w+') as file:
                file.write(templatei)
            # Modify the permissions of the script to allow it to be run on the node
            make_executable(kmaindex_script)
            # Run shell script
            os.system(kmaindex_script)

            # Set the database path for the analyses
            dbpath = '/mnt/nas2/databases/kma_v1.4.9/'
            #bacmetpath = '/mnt/nas2/databases/bacmet/'
            icebergpath = '/mnt/nas2/databases/ICEberg/'
            database_path = {
                'custom': os.path.join(target_dir, 'targets_KMA'),
                'amr': os.path.join(dbpath, 'NCBI-AMR'),
                'biocide': os.path.join(dbpath, 'NCBI-BIOCIDE'),
                'metal': os.path.join(dbpath, 'NCBI-METAL'),
                'verotoxin': os.path.join(dbpath, 'VEROTOXIN'),
                'bacmet': os.path.join(dbpath, 'bacmet_2018-03-11-v2_renamed_kma'),
                'all_ices': os.path.join(icebergpath, 'ICEs_kmav1.4.9'),
                'aice': os.path.join(icebergpath, 'AICEs_kmav1.4.9'),
                'cime': os.path.join(icebergpath, 'CIMEs_kmav1.4.9'),
                'ime': os.path.join(icebergpath, 'IMEs_kmav1.4.9'),
                't4ss': os.path.join(icebergpath, 'T4SS_kmav1.4.9'),
            }

        #create a folder to hold the sequences
        seq_dir = os.path.join(work_dir, 'sequences')
        os.makedirs(seq_dir, exist_ok=True)

        #now create an output folder so we can delete all of the extra files we didn't need
        out_dir = os.path.join(work_dir, 'output')
        os.makedirs(out_dir, exist_ok=True)

        #pull assemblies from nas
        if argument_dict['seqtype'] == 'fasta':
        # Run file linker and then make sure that all FASTA files requested are present. Warn user if they
        # requested things that we don't have.
            retrieve_nas_files(seqids=seqids,
                               outdir=seq_dir,
                               filetype='fasta',
                               copyflag=False)
            missing_fastas = verify_fasta_files_present(seqids, seq_dir)
            # Update the Redmine issue if one or more of the requested SEQIDs could not be located
            if missing_fastas:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastas))

            # These unfortunate hard coded paths appear to be necessary
            activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kma_v149'
            kma_py = '/mnt/nas2/virtual_environments/kma_v149/bin/kma'

            # Run kma with the necessary arguments
            #assembly files
            for assembly in glob.glob(os.path.join(seq_dir, '*.fasta')):
                seqid = os.path.split(assembly)[1].split('.')[0]
                #prepare command for kma
                kma_cmd = 'kma -i {assembly} -o {seqid} -t_db {dbpath} -t 7 -nf'\
                    .format(assembly=assembly,
                            seqid=seqid,
                            dbpath=database_path[argument_dict['analysis']])

                # Append the align and/or the unique flags are required
                kma_cmd += ' -bcNano' if argument_dict['nanopore'] else ''
                kma_cmd += ' -ef' if argument_dict['readcount'] else ''
                kma_cmd += ' -vcf 2' if argument_dict['vcf'] else ''
                kma_cmd += ' -hmm' if argument_dict['hmm'] else ''
                kma_cmd += ' -ID {}'.format(argument_dict['min_ID']) if argument_dict['min_ID'] else ''

                # Create another shell script to execute within the KMA conda environment
                template = "#!/bin/bash\n{} && cd {} && {}".format(activate, seq_dir, kma_cmd)
                kma_script = os.path.join(work_dir, 'run_kma.sh')
                with open(kma_script, 'w+') as file:
                    file.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(kma_script)
                # Run shell script
                os.system(kma_script)


        #pull raw reads from nas
        if argument_dict['seqtype'] == 'fastq':
            # Extract FASTQ files.
            retrieve_nas_files(seqids=seqids, outdir=seq_dir, filetype='fastq', copyflag=False)
            missing_fastqs = check_fastqs_present(seqids, seq_dir)
            if len(missing_fastqs) > 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastqs))

            #run analysis for each pair of fastq files
            if glob.glob(os.path.join(seq_dir, '*_R1_001.fastq.gz')):
                for rawread in glob.glob(os.path.join(seq_dir, '*_R1_001.fastq.gz')):
                    seqid1 = os.path.split(rawread)[1].split('.')[0]
                    seqid = os.path.split(seqid1)[1].split('_R')[0]
                    #seqidforfile = os.path.split(seqid1)[1].split('_S')[0]
                    forwardseq = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
                    reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
                    #prepare command for kma
                    activatepe = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kma_v149'
                    kma_py = '/mnt/nas2/virtual_environments/kma_v149/bin/kma'
                    kma_cmdpe = 'kma -ipe {forward} {reverse} -o {seqid} -t_db {db} -t 7 -nf'\
                        .format(forward=forwardseq,
                                reverse=reverseseq,
                                seqid=seqid,
                                db=database_path[argument_dict['analysis']])
                    # Append the align and/or the unique flags are required
                    kma_cmdpe += ' -ef' if argument_dict['readcount'] else ''
                    kma_cmdpe += ' -nc' if argument_dict['align'] == False else ''
                    kma_cmdpe += ' -vcf 2' if argument_dict['vcf'] else ''
                    kma_cmdpe += ' -hmm' if argument_dict['hmm'] else ''
                    kma_cmdpe += ' -ID {}'.format(argument_dict['min_ID']) if argument_dict['min_ID'] else ''

                    # Create another shell script to execute within the KMA conda environment
                    template_paired = "#!/bin/bash\n{} && cd {} && {}".format(activatepe, seq_dir, kma_cmdpe)
                    kma_script_paired = os.path.join(work_dir, 'run_kma.sh')
                    with open(kma_script_paired, 'w+') as file:
                        file.write(template_paired)
                    # Modify the permissions of the script to allow it to be run on the node
                    make_executable(kma_script_paired)
                    # Run shell script
                    os.system(kma_script_paired)

            #TODO: fix this... there is definitely a better way to account for files without the _001, but I dont want to look it up atm
            if glob.glob(os.path.join(seq_dir, '*_R1.fastq.gz')):
                for rawread in glob.glob(os.path.join(seq_dir, '*_R1.fastq.gz')):
                    seqid1 = os.path.split(rawread)[1].split('.')[0]
                    seqid = os.path.split(seqid1)[1].split('_R')[0]
                    #seqidforfile = os.path.split(seqid1)[1].split('_S')[0]
                    forwardseq = '{seqid}_R1.fastq.gz'.format(seqid=seqid)
                    reverseseq = '{seqid}_R2.fastq.gz'.format(seqid=seqid)
                    #prepare command for kma
                    activatepe = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kma_v149'
                    kma_py = '/mnt/nas2/virtual_environments/kma_v149/bin/kma'
                    kma_cmdpe = 'kma -ipe {forward} {reverse} -o {seqid} -t_db {db} -t 7 -nf'\
                        .format(forward=forwardseq,
                                reverse=reverseseq,
                                seqid=seqid,
                                db=database_path[argument_dict['analysis']])
                    # Append the align and/or the unique flags are required
                    kma_cmdpe += ' -ef' if argument_dict['readcount'] else ''
                    kma_cmdpe += ' -nc' if argument_dict['align'] == False else ''
                    kma_cmdpe += ' -vcf 2' if argument_dict['vcf'] else ''
                    kma_cmdpe += ' -hmm' if argument_dict['hmm'] else ''
                    kma_cmdpe += ' -ID {}'.format(argument_dict['min_ID']) if argument_dict['min_ID'] else ''

                    # Create another shell script to execute within the KMA conda environment
                    template_paired = "#!/bin/bash\n{} && cd {} && {}".format(activatepe, seq_dir, kma_cmdpe)
                    kma_script_paired = os.path.join(work_dir, 'run_kma.sh')
                    with open(kma_script_paired, 'w+') as file:
                        file.write(template_paired)
                    # Modify the permissions of the script to allow it to be run on the node
                    make_executable(kma_script_paired)
                    # Run shell script
                    os.system(kma_script_paired)

        #the minion raw data is single ended, not paired-end, so need a different command
        if argument_dict['seqtype'] == 'singlefastq':
            # Extract FASTQ files.
            retrieve_nas_files(seqids=seqids, outdir=seq_dir, filetype='fastq', copyflag=False)
            missing_fastqs = check_fastqs_present(seqids, seq_dir)
            if len(missing_fastqs) > 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastqs))
            # Run kma with the necessary arguments
            #fastq minion files
            for rawread in glob.glob(os.path.join(seq_dir, '*.fastq.gz')):
                seqid = os.path.split(rawread)[1].split('.')[0]
                #print(rawread)
                #print(seqid)
                #prepare command for kma
                activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kma_v149'
                kma_cmd = 'kma -i {rawread} -o {seqid} -t_db {dbpath} -t 7 -nf'\
                    .format(rawread=rawread,
                            seqid=seqid,
                            dbpath=database_path[argument_dict['analysis']])

                # Append the align and/or the unique flags are required
                #kma_cmd += ' -bcNano' if argument_dict['nanopore'] else ''
                kma_cmd += ' -ef' if argument_dict['readcount'] else ''
                kma_cmd += ' -vcf 2' if argument_dict['vcf'] else ''
                kma_cmd += ' -hmm' if argument_dict['hmm'] else ''
                kma_cmd += ' -ID {}'.format(argument_dict['min_ID']) if argument_dict['min_ID'] else ''

                # Create another shell script to execute within the KMA conda environment
                template = "#!/bin/bash\n{} && cd {} && {}".format(activate, seq_dir, kma_cmd)
                kma_script = os.path.join(work_dir, 'run_kma.sh')
                with open(kma_script, 'w+') as file:
                    file.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(kma_script)
                # Run shell script
                os.system(kma_script)

        #the minion raw data is single ended, not paired-end, so need a different command
        if argument_dict['seqtype'] == 'minionfastq':
            # Extract FASTQ files.
            retrieve_nas_files(seqids=seqids, outdir=seq_dir, filetype='fastq', copyflag=False)
            missing_fastqs = check_fastqs_present(seqids, seq_dir)
            if len(missing_fastqs) > 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastqs))
            # Run kma with the necessary arguments
            #fastq minion files
            for rawread in glob.glob(os.path.join(seq_dir, '*.fastq.gz')):
                seqid = os.path.split(rawread)[1].split('.')[0]
                #print(rawread)
                #print(seqid)
                #prepare command for kma
                activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kma_v149'
                kma_cmd = 'kma -i {rawread} -o {seqid} -t_db {dbpath} -t 7 -nf'\
                    .format(rawread=rawread,
                            seqid=seqid,
                            dbpath=database_path[argument_dict['analysis']])

                # Append the align and/or the unique flags are required
                kma_cmd += ' -bcNano' if argument_dict['nanopore'] else ''
                kma_cmd += ' -ef' if argument_dict['readcount'] else ''
                kma_cmd += ' -vcf 2' if argument_dict['vcf'] else ''
                kma_cmd += ' -hmm' if argument_dict['hmm'] else ''
                kma_cmd += ' -ID {}'.format(argument_dict['min_ID']) if argument_dict['min_ID'] else ''

                # Create another shell script to execute within the KMA conda environment
                template = "#!/bin/bash\n{} && cd {} && {}".format(activate, seq_dir, kma_cmd)
                kma_script = os.path.join(work_dir, 'run_kma.sh')
                with open(kma_script, 'w+') as file:
                    file.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(kma_script)
                # Run shell script
                os.system(kma_script)

        #add the filename (which is the seqid) to the first column in all of the res files, then concatenate into a single output csv file
        outputfile = os.path.join(seq_dir, 'kma_output.csv')
        with open(outputfile, 'w', newline='') as file_output:
            csv_output = csv.writer(file_output)
            #for fname in glob.glob(os.path.basename(seq_dir, '*.res')):
            if argument_dict['seqtype'] == 'fasta' or argument_dict['seqtype'] == 'minionfastq':
                for fname in glob.glob(os.path.join(seq_dir, '*.res')):
                    fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
                    seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
                    with open(fname, newline='') as f_input:
                        csv_input = csv.reader(f_input)
                        #Header processing
                        header = csv_input.__next__()
                        header.insert(0,"SeqID")

                        csv_output.writerow(header)
                        for row in csv_input:
                            print(row)
                            row.insert(0,seqname) #this adds the seqid to the file before concatenating
                            #row.insert(0,fname)
                            csv_output.writerow(row)
            elif argument_dict['seqtype'] == 'fastq':
                for fname in glob.glob(os.path.join(seq_dir, '*.res')):
                    fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
                    seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
                    seqname2 = os.path.split(seqname)[1].split('_')[0] #to remove the _SXX_L001 from seqid 
                    with open(fname, newline='') as f_input:
                        csv_input = csv.reader(f_input)
                        #Header processing
                        header = csv_input.__next__()
                        header.insert(0,"SeqID")

                        csv_output.writerow(header)
                        for row in csv_input:
                            print(row)
                            row.insert(0,seqname2) #this adds the seqid to the file before concatenating
                            #row.insert(0,fname)
                            csv_output.writerow(row)


        #KMA has finished, so we will modify the output
        #lets modify the output kma files for the databases Ashley curated so we end up with an excel file
        if argument_dict['analysis'] == 'amr' or argument_dict['analysis'] == 'metal' or argument_dict['analysis'] == 'biocide':
            kmadf = pandas.read_csv(os.path.join(seq_dir,'kma_output.csv'), sep=',|\t', engine='python')
            #rename the #Template column
            kmadf2 = kmadf.rename(columns={'#Template': 'Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass'})
            #split the newly renamed column
            kmadf2[['Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass']] = kmadf2['Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass'].str.split('|',expand=True)
            #drop the unwanted column
            kmadf2.drop('Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass', axis=1, inplace=True)
            #reorder the dataframe
            kmadf3 = kmadf2[['SeqID','Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value']]
            #write it to an excel file
            kmaexcelfile = os.path.join(seq_dir, 'KMA_results_redmine{}.xlsx'.format(issue.id))
            kmadf3.to_excel(kmaexcelfile,sheet_name='{}_Resistance'.format(argument_dict['analysis']), index=False)

        #modifications if bacmet is used
        elif argument_dict['analysis'] == 'bacmet':
            kmadf = pandas.read_csv(os.path.join(seq_dir,'kma_output.csv'), sep=',|\t', engine='python')
            #rename the #Template column
            kmadf2 = kmadf.rename(columns={'#Template': 'Gene|Bacmet_accession|Resistance_description'})
            #split the newly renamed column
            kmadf2[['Gene','Bacmet_accession','Resistance_description']] = kmadf2['Gene|Bacmet_accession|Resistance_description'].str.split('|',expand=True)
            #drop the unwanted column
            kmadf2.drop('Gene|Bacmet_accession|Resistance_description', axis=1, inplace=True)
            #reorder the dataframe
            kmadf3 = kmadf2[['SeqID','Gene','Bacmet_accession','Resistance_description','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value']]
            #write it to an excel file
            kmaexcelfile = os.path.join(seq_dir, 'KMA_results_redmine{}.xlsx'.format(issue.id))
            kmadf3.to_excel(kmaexcelfile,sheet_name='bacmet_Resistance', index=False)
        
        #modifications for everything else (custom, ICE, verotoxin)
        else:
            kmadf = pandas.read_csv(os.path.join(seq_dir,'kma_output.csv'), sep=',|\t', engine='python')
            #rename the #Template column
            kmadf2 = kmadf.rename(columns={'#Template': 'Query'})
            #write it to an excel file
            kmaexcelfile = os.path.join(seq_dir, 'KMA_results_redmine{}.xlsx'.format(issue.id))
            kmadf2.to_excel(kmaexcelfile,sheet_name='Targets_query', index=False)


        #now copy the output csv file to the new output directory... leaving as a directory so we can zip it if we add other functions and outputs later
        #shutil.copyfile(os.path.join(seq_dir, 'kma_output1.csv'), os.path.join(out_dir, 'kma_output1.csv'))
        shutil.copyfile(os.path.join(seq_dir, 'KMA_results_redmine{}.xlsx'.format(issue.id)), os.path.join(out_dir, 'KMA_results_redmine{}.xlsx'.format(issue.id)))

        #move all of the .mapstat files for readcounts to their own folder
        if argument_dict['readcount']:
            #make a combined csv of the mapstat results.... like we did for the res files above
            readcountsoutputfile = os.path.join(seq_dir, 'kma_readcounts_output.csv')
            with open(readcountsoutputfile, 'w', newline='') as readsfile_output:
                csv_output = csv.writer(readsfile_output)
                for rfname in glob.glob(os.path.join(seq_dir, '*.mapstat')):
                    rfbasename = os.path.basename(rfname) #this is to just get the seqid.res name of file
                    seqname = os.path.split(rfbasename)[1].split('.')[0] #this is to just pull out the seqid
                    seqname2 = os.path.split(seqname)[1].split('_')[0] #to remove the _SXX_L001 from seqid 
                    with open(rfname, newline='') as f_input:
                        csv_input = csv.reader(f_input)

                        #Header processing
                        header = csv_input.__next__()
                        header.insert(0,"SeqID")

                        csv_output.writerow(header)
                        for row in csv_input:
                            print(row)
                            row.insert(0,seqname2) #this adds the seqid to the file before concatenating
                            #row.insert(0,fname)
                            csv_output.writerow(row)

            #TODO: convert the file to a csv (comma-delimited)... currently not working I think because of ## in the file

            #now copy the output csv file to the new output directory... leaving as a directory so we can zip it if we add other functions and outputs later
            shutil.copyfile(os.path.join(seq_dir, 'kma_readcounts_output.csv'), os.path.join(out_dir, 'kma_readcounts_output.csv'))

            #now move all the individual readcount files to their own folder
            readcount_folder = os.path.join(out_dir,'readcount_files')
            os.makedirs(readcount_folder)
            output_files = os.listdir(seq_dir)
            for file in output_files:
                if file.endswith(".mapstat"):
                    shutil.move(os.path.join(seq_dir, file), os.path.join(readcount_folder, file))

        #move all of the .vcf files to their own folder
        if argument_dict['vcf']:
            vcf_folder = os.path.join(out_dir,'vcf_files')
            os.makedirs(vcf_folder)
            output_files = os.listdir(seq_dir)
            for file in output_files:
                if file.endswith(".vcf.gz"):
                    shutil.move(os.path.join(seq_dir, file), os.path.join(vcf_folder, file))


        #move all of the alignment files to their own folder
        if argument_dict['align']:
            align_folder = os.path.join(out_dir,'alignment_files')
            os.makedirs(align_folder)
            output_files = os.listdir(seq_dir)
            for file in output_files:
                if file.endswith(".aln"):
                    shutil.move(os.path.join(seq_dir, file), os.path.join(align_folder, file))
                if file.endswith(".fsa"):
                    shutil.move(os.path.join(seq_dir, file), os.path.join(align_folder, file))

        #upload the KMA excel file to the redmine request
        kmafilename = 'KMA_results_redmine{}.xlsx'.format(issue.id)
        output_list = list()
        output_dict = dict()
        output_dict['path'] = os.path.join(out_dir,kmafilename)
        output_dict['filename'] = kmafilename
        output_list.append(output_dict)

        redmine_instance.issue.update(resource_id=issue.id, uploads=output_list,
                                      notes='KMA analysis complete!\n\n'
                                            'Results are attached. \n'
                                            'The ftp is currently down, so you may only be able to access the attached report file.')

           
        # Zip output
        kmaout_filename = 'kma_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=out_dir,
                                  output_dir=work_dir,
                                  output_filename=kmaout_filename)
        zip_filepath += '.zip'
        
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
                notes='KMA analysis complete!\n\n'
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

        # Remove all the folders
        shutil.rmtree(seq_dir)
        # Remove the zip file
        os.remove(zip_filepath)

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Please contact a bioinformatician '
                                            'to investigate: {}'.format(e))




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

def check_fastqs_present(fastq_list, fastq_dir):
    missing_fastqs = list()
    for seqid in fastq_list:
        if len(glob.glob(os.path.join(fastq_dir, seqid + '*.fastq.gz'))) < 2:
            # JAS adding an extra if statement here to allow for Nanopore (SE) reads
            if len(glob.glob(os.path.join(fastq_dir, seqid + ".fastq.gz"))) == 0:
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
    kma_redmine()
