import os
import glob
import click
import pickle
import shutil
import ftplib
import sentry_sdk
import subprocess
import fnmatch
import pandas as pd
from pathlib import Path
from Bio import Phylo
from Bio import SeqIO
from biotools import mash
from strainchoosr import strainchoosr
from collections import OrderedDict
from automator_settings import SENTRY_DSN
from amrsummary import before_send
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)
from nastools.nastools import retrieve_nas_files

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def unknownisolate_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    #analyses supported by the automator
    analyses = [
        'custom','unknown', 'enterobacterales', 'listeriaceae', 'enterobacter', 'enterococcus', 'bacillus', 'mycobacterium','genus'
    ]

    #clustering algorithms supported by the automator
    clusteralgorithms = [
         'median', 'weighted', 'single', 'complete', 'average', 'ward', 'centroid'
    ]


    #variable to hold supplied arguments
    argument_dict = {
        'analysis': str(),
        'query': str(),
        'genus': False,
    }


    #parse description to figure out analysis type, and find fastas
    try:
        #parse description
        seqids = list()
        query_list = list()
        for item in description:
            item = item.rstrip()
            if 'analysis' in item:
               argument_dict['analysis'] = item.split('=')[1].lower()
               continue
            if 'query' in item:
               argument_dict['query'] = item.split('=')[1]
               query_list.append(argument_dict['query'])
               continue
            if 'genus' in item:
               argument_dict['genus'] = item.split('=')[1].upper()
               continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that the analysis type is provided
        if not argument_dict['analysis']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not identify an analysis type. '
                                                'Please ensure that the first line of the issue contains "analysis=" one'
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

        # Ensure that the user includes a query sequence
        if not argument_dict['query']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not identify query sequence. '
                                                'Please ensure that the request includes a line: "query=YYYY-SEQ-NNNN".'.format(ats=', '.join(analyses)),
                                          status_id=4)
            return

        #if they choose the genus option, make sure they actually include a genus in their request
        if argument_dict['analysis'] == 'genus':
            if not argument_dict['genus']:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Analysis type "genus" was selected. '
                                              'No genus argument was found. Please ensure that the request includes a line: "genus=GENUS".',
                                              status_id=4)

        #TODO: maybe a better way to do this is to create a sketch of all organisms, then compare the assembly to that sketch? should be faster
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

            #link to the atcc enterobacterales
            enterobacteraleses = ['ATLANTIBACTER','BUTTIAUXELLA','CITROBACTER',
                                'CRONOBACTER','EDWARDSIELLA','ENTEROBACTER',
                                'ESCHERICHIA','HAFNIA',
                                'KLEBSIELLA',
                                'LECLERCIA','LELLIOTTIA',
                                'PANTOEA','PASTEURELLA','PECTOBACTERIUM','PHYTOBACTER',
                                'PLESIOMONAS','PLURALIBACTER',
                                'RAOULTELLA','SALMONELLA','SERRATIA','TATUMELLA',
                                'SHIGELLA','SHIMWELLIA','YERSINIA',
                                'YOKENELLA']
            #list of enterobacterales genera that we don't currently have ATCC sequences for:
            #'PRAGIA','CEDECEA','DICKEYA','LEVINEA','BUDVICIA','KLUYVERA','RAHNELLA','AVERYELLA','EWINGELLA','KOSAKONIA','KOSERELLA',
            #'GIBBSIELLA','IZHAKIELLA','JEJUBACTER','EDAPHOVIRGA','LEMINORELLA','NISSABACTER','SICCIBACTER','LIMNOBACULUM','SCANDINAVIUM',
            #'TRABULSIELLA','ENTEROBACILLUS','FRANCONIBACTER','JINSHANIBACTER','MANGROVIBACTER','PHASEOLIBACTER','ROSENBERGIELLA',
            #'SACCHAROBACTER','ERYTHROBACILLUS',,'INSECTIHABITANS','INTESTINIRHABDUS','CAYMMATOBACTERIUM','PSEUDOENTEROBACTER','PSEUDOESCHERICHIA','SUPERFICIEIBACTER',
            for genus in enterobacteraleses:
                src3 = '/mnt/nas2/processed_sequence_data/atcc/all_genera/BestAssemblies/'
                dst = fasta_dir
                lncmd3 = 'ln -s {src}{genus}_*.fasta {dst}'.format(src=src3,genus=genus,
                                                                    dst=dst)
                os.system(lncmd3)

                #remove empty files.. because of my weird list above we end up with some *_.fasta.. these weren't empty, which is probably why this didnt work
                #fastafiles = os.listdir(fasta_dir)
                #for filename in glob.glob(os.path.join(work_dir, 'fastas',"*.fasta")):
                #for file in fastafiles:
                    #if os.stat(file).st_size == 0:
                    #    os.remove(file)
                    #filewithoutsuffix=os.path.splitext(file)[0]
                    #if filename=='{}_\*.fasta':
                    #    os.remove(filename)
    
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
            src2 = '/mnt/nas2/processed_sequence_data/atcc/all_genera/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst2 = fasta_dir
            lncmd2 = 'ln -s {src}LISTERIA_*.fasta {dst}'.format(src=src2,
                                                      dst=dst2)
            os.system(lncmd2)

            #remove the -BAA sequences because they make the automator angry...
            fastafiles = os.listdir(fasta_dir)
            for file in fastafiles:
                #if fnmatch.fnmatch(file, '*-BAA*'):
                file_without_suffix=os.path.splitext(file)[0]
                if "*-BAA*" in file_without_suffix:
                    os.remove(file)
                    print(file_without_suffix)

            #remove empty files.. because of my weird list above we end up with some *_.fasta
            #fastafiles = os.listdir(fasta_dir)
            #for file in fastafiles:
            #    if os.stat(file).st_size == 0:
            #        os.remove(file)

        #if enterococcus, link to enterococcus reference sequences
        if argument_dict['analysis'] == 'enterococcus':
            src = '/mnt/nas2/processed_sequence_data/ncbi/enterococcus_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

            #also do this for the atcc sequences
            src2 = '/mnt/nas2/processed_sequence_data/atcc/all_genera/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst2 = fasta_dir
            lncmd2 = 'ln -s {src}ENTEROCOCCUS_*.fasta {dst}'.format(src=src2,
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

            #the enterobacterales ncbi reference Enterobacters
            src2 = '/mnt/nas2/processed_sequence_data/ncbi/enterobacterales_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd2 = 'ln -s {src}ENTEROBACTER_*.fasta {dst}'.format(src=src2,
                                                                    dst=dst)
            os.system(lncmd2)

            #link to the atcc enterobacters
            src3 = '/mnt/nas2/processed_sequence_data/atcc/all_genera/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd3 = 'ln -s {src}ENTEROBACTER_*.fasta {dst}'.format(src=src3,
                                                                    dst=dst)
            os.system(lncmd3)

        #if bacillus, link to bacillus reference sequences
        if argument_dict['analysis'] == 'bacillus':
            src = '/mnt/nas2/processed_sequence_data/ncbi/bacillus_refseq/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

            #link to the atcc bacillus
            src3 = '/mnt/nas2/processed_sequence_data/atcc/all_genera/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd3 = 'ln -s {src}BACILLUS_*.fasta {dst}'.format(src=src3,
                                                                    dst=dst)
            os.system(lncmd3)

        #if mycobacterium, link to mycobacterium reference sequences
        if argument_dict['analysis'] == 'mycobacterium':
            src = '/mnt/nas2/processed_sequence_data/ncbi/mycobacterium/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

            #link to the atcc mycobacterium
            src3 = '/mnt/nas2/processed_sequence_data/atcc/all_genera/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd3 = 'ln -s {src}MYCOBACTERIUM_*.fasta {dst}'.format(src=src3,
                                                                    dst=dst)
            os.system(lncmd3)

        #if a custom analysis, create the fastas folder
        if argument_dict['analysis'] == 'custom':
            #also warn user that custom analysis should only be used when they have reason to believe their selected sequences are closely related
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: analysis type "custom" was selected. \n'
                                                'This analysis will use only the sequences listed in the description for comparison to your query. \n'
                                                'This analysis type should only be used if you would like to compare specific '
                                                'sequences, as by excluding reference sequences you are likely to miss '
                                                'something that may be closely related.',
                                          status_id=2)
            # Make fasta file directory
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            #now add all of the seqids from the list into the fasta_dir
            retrieve_nas_files(seqids=seqids,
                               outdir=fasta_dir,
                               filetype='fasta',
                               copyflag=False)
            missing_fastas = verify_fasta_files_present(seqids, fasta_dir)

        #if genus, take that genus name and pull from all reference folders (ATCC and NCBI) 
        if argument_dict['analysis'] == 'genus':
            genus_list = list()
            atccpath = '/mnt/nas2/processed_sequence_data/atcc/all_genera/BestAssemblies/'
            atcclist = os.listdir(atccpath)
            for file in atcclist:
                if file.startswith("{}".format(argument_dict['genus'])):
                    genus_list.append(os.path.splitext(file)[0])

            #ncbipath = '/mnt/nas2/processed_sequence_data/ncbi/'
            paths_list = glob.glob('/mnt/nas2/processed_sequence_data/ncbi/*/BestAssemblies/')
            dirs_list = [path for path in paths_list if os.path.isdir(path)]
            #search each directory
            for d in dirs_list:
                files_list = os.listdir(d)
                for file in files_list:
                    if file.startswith("{}".format(argument_dict['genus'])):
                        genus_list.append(os.path.splitext(file)[0])
            #now use nastools to link
            # Make fasta file directory
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            #now add all of the seqids
            retrieve_nas_files(seqids=genus_list,
                               outdir=fasta_dir,
                               filetype='fasta',
                               copyflag=False)

        #if a completely unknown isolate, pull all of the ATCC sequences for analysis
        if argument_dict['analysis'] == 'unknown':
            src = '/mnt/nas2/processed_sequence_data/atcc/all_genera/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd = 'ln -s {src}*.fasta {dst}'.format(src=src,
                                                      dst=dst)
            os.system(lncmd)

            #pull the ncbi enterobacterales and listeriaceae refseq sequences, as there are some weirdos not in ATCC
            src2 = '/mnt/nas2/processed_sequence_data/ncbi/enterobacterales_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd2 = 'ln -s {src}*.fasta {dst}'.format(src=src2,
                                                      dst=dst)

            os.system(lncmd2)

            src3 = '/mnt/nas2/processed_sequence_data/ncbi/listeriaceae_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd3 = 'ln -s {src}*.fasta {dst}'.format(src=src3,
                                                      dst=dst)

            os.system(lncmd3)
            #pull the enterococcus ncbi refseq sequences
            src4 = '/mnt/nas2/processed_sequence_data/ncbi/enterococcus_references/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd4 = 'ln -s {src}*.fasta {dst}'.format(src=src4,
                                                      dst=dst)

            os.system(lncmd4)
            #pull the bacillus ncbi refseq sequences
            src5 = '/mnt/nas2/processed_sequence_data/ncbi/bacillus_refseq/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd5 = 'ln -s {src}*.fasta {dst}'.format(src=src5,
                                                      dst=dst)
            os.system(lncmd5)

            #pull the mycobacterium ncbi refseq sequences
            src5 = '/mnt/nas2/processed_sequence_data/ncbi/mycobacterium/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd5 = 'ln -s {src}*.fasta {dst}'.format(src=src5,
                                                      dst=dst)
            os.system(lncmd5)

            #pull the other random ncbi refseq sequences
            src6 = '/mnt/nas2/processed_sequence_data/ncbi/other_refseq_for_unknownisolates/BestAssemblies/'
            fasta_dir = os.path.join(work_dir, 'fastas')
            if not os.path.isdir(fasta_dir):
                os.makedirs(fasta_dir)
            dst = fasta_dir
            lncmd6 = 'ln -s {src}*.fasta {dst}'.format(src=src6,
                                                      dst=dst)
            os.system(lncmd6)


            #also warn user that they should use the output from this to do a more targeted analysis
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: analysis type "unknown" was selected. \n'
                                                'This analysis will use >2600 ATCC reference genomes and hundreds of refseq sequences for comparison to your query. \n'
                                                'This will take a while.',
                                          status_id=2)


        # Only allowed to have one query file - boot the user if they tried to specify too many queries.
        if len(query_list) > 1:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='ERROR: You specified {query_list_len} query files ({query_list}). Only'
                                                ' one query file is supposed to be specified. '
                                                'Please try again.'.format(query_list_len=len(query_list),
                                                                           query_list=query_list))

        #allow users to upload a fasta file as the query...
        if argument_dict['query'] == 'attached': 
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
                    attachment.download(savepath=work_dir)
            #TODO: get the id of this file... may need to modify rest of script so we aren't removing and re-adding all the time? copyflag=True?
            for assembly in glob.glob(os.path.join(work_dir, '*.fasta')):
                seqid = os.path.split(assembly)[1].split('.')[0]
        else:
            # retrieve the query file if it is not an attachment
            retrieve_nas_files(seqids=query_list,
                               outdir=work_dir,
                               filetype='fasta',
                               copyflag=False)
            verify_fasta_files_present(query_list, work_dir)


        #run the sequence through pubmlst and print the output to redmine
        #activate environment
        pubmlstactivate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/pubmlst'
        pubmlst_py = '/mnt/nas2/virtual_environments/pubmlst/species_identification.py'

        assembly = glob.glob(os.path.join(work_dir, '*.fasta'))[0]
        queryseqid = os.path.split(assembly)[1].split('.')[0]
        pubmlst_cmd = 'python {py} -f {assembly} > {work_dir}/output.txt'.format(py=pubmlst_py, assembly=assembly,work_dir=work_dir)

        #create another shell script to execute within the dRep conda environment
        template = "#!/bin/bash\n{} && {}".format(pubmlstactivate, pubmlst_cmd)
        pubmlst_script = os.path.join(work_dir, 'run_pubmlst.sh')
        with open(pubmlst_script, 'w+') as file:
            file.write(template)

        # Modify the permissions of the script to allow it to be run on the node
        make_executable(pubmlst_script)
        #run command
        os.system(pubmlst_script)

        #check if the pubmlst file output is empty... sometimes the server disconnects us.
        if os.stat(os.path.join(work_dir, 'output.txt')).st_size == 0:
            #warn the user that it is empty, but continue the pipeline
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='PubMLST typing analysis (ribosomal MLST) failed. \n'
                                                'The connection to the server was likely severed. Please re-submit. \n'
                                                'If problem persists, contact a bioinformatician.')
        else:
            #parse the pubmlst results
            rank_line = [0]
            taxon_line = [1]
            support_line = [2]
            taxonomy_line = [3]
            o_file = open(os.path.join(work_dir, 'output.txt'))
            for index, line in enumerate(o_file):
                if index in rank_line:
                    Rank = line.split(':')[1]
                if index in taxon_line:
                    Taxon = line.split(':')[1]
                if index in support_line:
                    Support = line.split(':')[1]
                if index in taxonomy_line:
                    Taxonomy = line.split('axonomy')[1]

            # Update the issue with the pubmlst output
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='PubMLST typing analysis (ribosomal MLST) complete! Here is the closest match from pubMLST:\n \n'
                                                ' Rank: {}\n'
                                                ' Taxon: {}\n'
                                                ' Support: {}\n'
                                                ' Taxonomy: {}\n'
                                                'Now running MASH, and dRep compare'.format(Rank,Taxon,Support,Taxonomy))

        #TODO: run a 16S analysis of the sequence
        # sixteensdb = '/mnt/nas2/databases/assemblydatabases/0.5.0.18/sixteens_full'
        # if not os.path.isdir(os.path.join(work_dir, '16S_results')):
        #     os.makedirs(os.path.join(work_dir, '16S_results'))
        #
        # # These unfortunate hard coded paths appear to be necessary
        # activateseekr = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/geneseekr'
        # seekr_py = '/mnt/nas2/virtual_environments/geneseekr/bin/GeneSeekr'
        # # Run sipprverse with the necessary arguments
        # seekr_cmd = 'python {seekr_py} blastn -s {seqpath} -r {outpath} -t {dbpath} -c {cutoff} -e {evalue} -f'\
        #     .format(seekr_py=seekr_py,
        #             seqpath=work_dir,
        #             outpath=os.path.join(work_dir, '16S_results'),
        #             dbpath=sixteensdb,
        #             cutoff=70,
        #             evalue='1E-15')
        # # Update the issue with the GeneSeekr command
        # redmine_instance.issue.update(resource_id=issue.id,
        #                               notes='GeneSeekr command:\n {cmd}'.format(cmd=seekr_cmd))
        # # Create another shell script to execute within the PlasmidExtractor conda environment
        # templateseekr = "#!/bin/bash\n{} && {}".format(activateseekr, seekr_cmd)
        # geneseekr_script = os.path.join(work_dir, 'run_geneseekr.sh')
        # with open(geneseekr_script, 'w+') as file:
        #     file.write(templateseekr)
        # # Modify the permissions of the script to allow it to be run on the node
        # make_executable(geneseekr_script)
        # # Run shell script
        # os.system(geneseekr_script)


	    #Run a mash to figure out if any strains are particularly far apart and likely to make PARSNP fail.
	    # Make the reference file
        reference_file = glob.glob(os.path.join(work_dir, '*.fasta'))[0]
        make_ref(reference_file, os.path.join(work_dir, 'reference.fasta'))

        #Only report this if the analysis type is not unknown, or else we'll end up with thousands of lines in the redmine output
        if argument_dict['analysis'] != 'unknown':
            bad_fastas = check_distances(reference_file, os.path.join(work_dir, 'fastas')) #the float/string issue was because some sequences had spaces in their filenames
            print(bad_fastas)
            if bad_fastas:
                outstr = ''
                for fasta in bad_fastas:
                    fasta = os.path.split(fasta)[-1]
                    outstr += fasta + '\n'
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='Warning! MASH screening thinks that the following samples may be too'
                                                    ' far from the reference: {samples}\nIn this case, the reference file'
                                                    ' was {reference}.'.format(samples=outstr,
                                                                        reference=os.path.split(reference_file)[-1]))

        #report the closest matches from MASH to the user
        close_fastas = closest_distances(reference_file, os.path.join(work_dir, 'fastas'))
        close_list = list()
        if close_fastas:
            outstr = ''
            for fasta in close_fastas:
                fasta = os.path.split(fasta)[-1]
                fbasename = os.path.basename(fasta)
                seqname = os.path.split(fbasename)[1].split('.')[0]
                close_list.append(seqname) #add the seq name to a new list
                outstr += fasta + '\n'
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='MASH screening thinks that the following samples are'
                                                ' close-ish to the reference:\n {samples}\nIn this case, the reference file'
                                                ' was {reference}. The automator will run with these '
                                                ' sequences only.'.format(samples=outstr,
                                                                    reference=os.path.split(reference_file)[-1]))


        if not os.path.isdir(os.path.join(work_dir, 'mash_output')):
            os.makedirs(os.path.join(work_dir, 'mash_output'))
        # Move distances.tab and sketch.msh from fastas folder, because sometimes they make
        # parsnp crash. Other times they don't. I have no idea why, so remove just to be safe.
        #now going to move the mash files to the mash_output folder
        shutil.move(os.path.join(work_dir, 'fastas', 'distances.tab'), os.path.join(work_dir,'mash_output', 'all_sequences_distances.tab'))
        shutil.move(os.path.join(work_dir, 'fastas', 'sketch.msh'), os.path.join(work_dir,'mash_output', 'all_sequences_sketch.msh'))
        #try:
        #    os.remove(os.path.join(work_dir, 'fastas', 'distances.tab'))
        #    os.remove(os.path.join(work_dir, 'fastas', 'sketch.msh'))
        #except OSError:
        #    pass

        #now remove the reference.fasta file
        os.remove(os.path.join(work_dir, 'reference.fasta'))

        #add the query sequence to the fastas folder for use with mashtree
        if argument_dict['query'] == 'attached':
            #query_file = glob.glob(os.path.join(work_dir, '*.fasta'))[0]
            #fbasename = os.path.basename(query_file)
            #seqname = os.path.split(fbasename)[1].split('.')[0]
            #print(seqname)
            #query_file = glob.glob(os.path.join(work_dir, '*.fasta'))
            #query_id = os.path.split(query_file)[1].split('.')[0]
            #shutil.copy(os.path.join(work_dir, '{}.fasta'.format(seqname)), os.path.join(fasta_dir, '{}.fasta'.format(seqname))) ###here is problem
            #attempt 2
            #there HAS to be a better way to get the filename and copy the file? whyyy
            fileslist = os.listdir(work_dir)
            print(fileslist)
            for file in fileslist:
                if file.endswith(".fasta"):
                    shutil.copy(os.path.join(work_dir, file), os.path.join(fasta_dir,file))
        else:
            retrieve_nas_files(seqids=query_list,
                               outdir=fasta_dir,
                               filetype='fasta',
                               copyflag=False)

        #run mashtree
        # Full paths needed here since SLURM doesn't give the $PATH of the host machine to the script for some reason
        cmdm = '/home/ubuntu/bin/mashtree --numcpus 24 --outtree {output_newick} {input_fastas}'\
            .format(output_newick=os.path.join(work_dir, 'mash_output', 'all_sequences_mashtree.tree'),
                    input_fastas=os.path.join(work_dir, 'fastas', '*.fasta'))


        returncode = subprocess.call(cmdm, shell=True, env={'PERL5LIB': '$PERL5LIB:/home/ubuntu/lib/perl5'})
        if returncode != 0:
            raise Exception('Tree creation command ({}) for {} had return code {}'.format(cmdm, issue.id, returncode))



        #now that we have a mashtree, remove the fastas folder and re-build with only those in the "close" list
        shutil.rmtree(fasta_dir)

        #re-make fasta_dir
        fasta_dir = os.path.join(work_dir, 'fastas')
        if not os.path.isdir(fasta_dir):
            os.makedirs(fasta_dir)
        #now add all of the seqids from the close-ish matches list into the fasta_dir
        retrieve_nas_files(seqids=close_list,
                           outdir=fasta_dir,
                           filetype='fasta',
                           copyflag=False)


        #will try with mashtree instead, since parsnp won't play nice. The whole point of this is just to get a list of close isolates to run ANI on
        tree = Phylo.read(os.path.join(work_dir, 'mash_output', 'all_sequences_mashtree.tree'), 'newick')
        if argument_dict['query'] != 'attached':
            ref_clades = tree.find_clades('{}'.format(argument_dict['query']))
        else:
            for file in fileslist:
                if file.endswith(".fasta"):
                    fbasename = os.path.basename(file)
                    query_id = os.path.split(fbasename)[1].split('.')[0]
            print('query_id:',query_id)
            #ref_clades = tree.find_clades('{}'.format(query_id))
            ref_clades = tree.find_clades('{}'.format(seqname)) #we created the seqname above when we copied the file
        #now find the closest matching tree tips
        for clade in ref_clades:
            ref_clade = clade
        clades = tree.get_terminals()
        distance_dict = dict()
        for clade in clades:
            distance = tree.distance(clade, ref_clade)
            distance_dict[clade.name] = distance

        # Use some stackoverflow magic to sort dict https://stackoverflow.com/questions/613183/how-do-i-sort-a-dictionary-by-value
        sorted_dict = OrderedDict(sorted(distance_dict.items(), key=lambda x: x[1]))

        #changing this number if mycobacterium is used, or else it takes FOREVER for pyani to complete
        if argument_dict['analysis'] == 'mycobacterium':
            desired_num_strains = 15
        else:
            desired_num_strains = 40 #the neartree automator allows users to select number of strains... I will set it here to simplify the process


        i = 0
        outstr = ''
        nearlist = list()
        for key in sorted_dict:
            if 'reference' not in key and i < desired_num_strains: #might need to edit this to argument_dict['query'] not in key...
                seqname = os.path.split(key)[1].split('.')[0]
                nearlist.append(seqname) #add the seq name to a new list
                outstr += key.replace('.fasta', '') + '\n'
                i += 1

        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='MASHtree process complete! Closest strains are:\n {}\n'
                                            'Now running dRep compare, pyani (ANIb and ANIm) for ANI values'.format(outstr))

        #now that we have a mashtree, remove the fastas folder and re-build with only those in the "close" list
        shutil.rmtree(fasta_dir)

        #make a close matches fastas directory
        fasta_dir = os.path.join(work_dir, 'fastas')
        if not os.path.isdir(fasta_dir):
            os.makedirs(fasta_dir)

        #now add all of the seqids from the close-ish matches list into the fasta_dir
        retrieve_nas_files(seqids=nearlist,
                           outdir=fasta_dir,
                           filetype='fasta',
                           copyflag=False)
        #add the query sequence as well
        if argument_dict['query'] == 'attached':
            source = work_dir
            dst = fasta_dir
            lncmd = 'ln -s {src}/*.fasta {dst}'.format(src=source,dst=dst)
            os.system(lncmd)
            #assembly_file = os.path.join(work_dir, '*.fasta')
            #shutil.copy(os.path.join(work_dir, assembly_file), os.path.join(fasta_dir, assembly_file))
        else:
            retrieve_nas_files(seqids=query_list,
                               outdir=fasta_dir,
                               filetype='fasta',
                               copyflag=False)

        #run dRep compare using only the files from the mashtree closest matches!
        output_drep = os.path.join(work_dir, 'dRep_output')
        if not os.path.isdir(output_drep):
            os.makedirs(output_drep)

        # These unfortunate hard coded paths appear to be necessary
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/dRep'
        dRep_cmd = 'dRep compare {outpath} -p 24 -g {seqpath}/*.fasta -pa 0.8 -sa 0.8 -nc 0.1 ' \
                   '--clusterAlg ward'.format(seqpath=fasta_dir, outpath=output_drep)

        #create another shell script to execute within the dRep conda environment
        template = "#!/bin/bash\n{} && {}".format(activate, dRep_cmd)
        dRep_script = os.path.join(work_dir, 'run_dRep.sh')
        with open(dRep_script, 'w+') as file:
            file.write(template)

        # Modify the permissions of the script to allow it to be run on the node
        make_executable(dRep_script)
        # Run shell script
        os.system(dRep_script)

        #now run ANIm and ANIb (in the pyani conda env) to get ANI values for the sequences
        output_ani = os.path.join(work_dir, 'pyani_output')
        if not os.path.isdir(output_ani):
            os.makedirs(output_ani)

        #these hard coded paths are necessary
        activatepyani = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/pyani2'
        pyani_b_cmd = 'average_nucleotide_identity.py -i {seqpath} -m ANIb -o {outpath}/ANIb -g --gformat pdf,png'.\
            format(seqpath=fasta_dir, outpath=output_ani)
        pyani_m_cmd = 'average_nucleotide_identity.py -i {seqpath} -m ANIm -o {outpath}/ANIm -g --gformat pdf,png'.\
            format(seqpath=fasta_dir, outpath=output_ani)

        #create another shell script to execute within the dRep conda environment
        templatepyani = "#!/bin/bash\n{} && {} && {}".format(activatepyani, pyani_b_cmd, pyani_m_cmd)
        pyani_script = os.path.join(work_dir, 'run_pyani.sh')
        with open(pyani_script, 'w+') as file:
            file.write(templatepyani)

        # Modify the permissions of the script to allow it to be run on the node
        make_executable(pyani_script)
        # Run shell script
        os.system(pyani_script)

        #create a pdf file using the GROBI script in the automator folder on the nas
        activatefpdf = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/pythonreport'
        grobipy = '/mnt/nas2/redmine/applications/OLCRedmineAutomator/automators/grobi_pdf_report.py'
        grobireport = os.path.join(work_dir, 'GROBI_{}.pdf'.format(issue.id))
        grobi_cmd = 'python {grpy} -s {seq} -i_rmlst {pub} -i_anib {anib} -anib_mat {bmat} -i_anim {anim} ' \
                    '-anim_mat {mmat} -i_mash {mash} -drep_tree {tree} -o {out}'\
            .format(grpy=grobipy,seq=queryseqid,
                    pub=os.path.join(work_dir,'output.txt'),
                    anib=os.path.join(work_dir, 'pyani_output','ANIb','ANIb_percentage_identity.tab'),
                    bmat=os.path.join(work_dir, 'pyani_output','ANIb','ANIb_percentage_identity.png'),
                    anim=os.path.join(work_dir, 'pyani_output','ANIm','ANIm_percentage_identity.tab'),
                    mmat=os.path.join(work_dir, 'pyani_output','ANIm','ANIm_percentage_identity.png'),
                    mash=os.path.join(work_dir, 'mash_output', 'all_sequences_distances.tab'),
                    tree=os.path.join(work_dir,'dRep_output','figures','Primary_clustering_dendrogram.pdf'),
                    out=grobireport)


        #create another shell script to execute within the conda environment
        templategrobi = "#!/bin/bash\n{} && cd {} && {}".format(activatefpdf, work_dir, grobi_cmd)
        grobi_script = os.path.join(work_dir, 'create_GROBI.sh')
        with open(grobi_script, 'w+') as file:
            file.write(templategrobi)

        # Modify the permissions of the script to allow it to be run on the node
        make_executable(grobi_script)
        # Run shell script
        os.system(grobi_script)

        #make an output directory
        output_dir = os.path.join(work_dir, 'Results')
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        #move folders into the output directory
        folders_to_move = [os.path.join(work_dir, 'pyani_output'),
                           os.path.join(work_dir, 'dRep_output'),
                           os.path.join(work_dir, 'mash_output'),
                           #os.path.join(work_dir, 'parsnp_output')
                           ]
        for folder in folders_to_move:
            shutil.move(folder, output_dir)

        # move the pubmlst output file to its own folder
        rMLST_output = os.path.join(work_dir, 'Results','rMLST_output')
        os.makedirs(rMLST_output)
        # move the output file
        shutil.copy(os.path.join(work_dir, 'output.txt'), os.path.join(rMLST_output,'output.txt'))

        #copy the grobi report and table file to the results folder
        grobfile = 'GROBI_{}.pdf'.format(issue.id)
        shutil.copy(os.path.join(work_dir, grobfile), os.path.join(output_dir, grobfile))

        grobtable = 'GROBI_table.pdf'
        shutil.copy(os.path.join(work_dir, grobtable), os.path.join(output_dir, grobtable))

        #attach the GROBI report to the redmine issue
        output_list = list()
        output_dict = dict()
        grobifile = 'GROBI_{}.pdf'.format(issue.id)
        output_dict['path'] = os.path.join(work_dir, grobifile)
        output_dict['filename'] = grobifile
        output_list.append(output_dict)


	    # Zip output folder and upload as sometimes files are >10Mb
        zipfilename = 'unknown_isolate_results_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=output_dir,
                                  output_dir=work_dir,
                                  output_filename=zipfilename)
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
                uploads=output_list,
                notes='Analysis complete!\n\n'
                      'GROBI report is attached\n'
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

        # Create a list of all the folders - will be used to clean up the working directory
        folders = glob.glob(os.path.join(work_dir, '*/'))
        # Remove all the folders
        for folder in folders:
            if os.path.isdir(folder):
                shutil.rmtree(folder)

        #remove unnecessary files, but leave the report file just in case
        os.remove(os.path.join(work_dir, 'GROBI_ANI_matrices.pdf'))
        os.remove(os.path.join(work_dir, 'GROBI_table.pdf'))
        os.remove(os.path.join(work_dir, 'GROBI_table.png'))
        os.remove(os.path.join(work_dir, 'Primary_clustering_dendrogram0.png'))
        os.remove(os.path.join(work_dir, 'reportlab_coverpage.pdf'))
        os.remove(os.path.join(work_dir, '{}.fasta'.format(queryseqid)))
#        os.remove(os.path.join(work_dir, 'reference.fasta'))

        # os.remove(os.path.join(work_dir, zipfilename + '.zip')) #delete zip file
        # Wrap up issue
        #redmine_instance.issue.update(resource_id=issue.id,
        #                              #uploads=output_list,
        #                              status_id=4,
        #                              notes='Analysis complete!')
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))

def check_distances(ref_fasta, fasta_folder):
    bad_fastqs = list()
    # fastqs = glob.glob(os.path.join(fastq_folder, '*R1*'))
    mash.sketch(os.path.join(fasta_folder, '*.fasta'), output_sketch=os.path.join(fasta_folder, 'sketch.msh'), threads=46)
    mash.dist(os.path.join(fasta_folder, 'sketch.msh'), ref_fasta, threads=46, output_file=os.path.join(fasta_folder, 'distances.tab'))
    mash_output = mash.read_mash_output(os.path.join(fasta_folder, 'distances.tab'))
    for item in mash_output:
        print(item.reference, item.query) #trying to troubleshoot string to float issue... I have no idea why this fixed the problem for one isolate?.. is it the 1 values?
        print(item.reference, item.query, str(item.distance))
        #print(item.reference, item.query, float(item.distance))
        if item.distance > 0.20:  # May need to adjust this value. Ashley adjusted it from 0.06...
            bad_fastqs.append(item.reference)
    return bad_fastqs

def closest_distances(ref_fasta, fasta_folder):
    close_fastas = list()
    # fastqs = glob.glob(os.path.join(fastq_folder, '*R1*'))
    mash.sketch(os.path.join(fasta_folder, '*.fasta'), output_sketch=os.path.join(fasta_folder, 'sketch.msh'), threads=46)
    mash.dist(os.path.join(fasta_folder, 'sketch.msh'), ref_fasta, threads=46, output_file=os.path.join(fasta_folder, 'distances.tab'))
    mash_output = mash.read_mash_output(os.path.join(fasta_folder, 'distances.tab'))
    for item in mash_output:
        print(item.reference, item.query) #trying to troubleshoot string to float issue
        print(item.reference, item.query, str(item.distance))
        #print(item.reference, item.query, float(item.distance))
        if item.distance < 0.20:  # May need to adjust this value. Ashley adjusted it from 0.06... Adjusted to 0.20 after it failed for an enterococcus isolate
            close_fastas.append(item.reference)
    return close_fastas
    #in test runs, Salmonella distance to the refseq Salmonella was 0.009 to 0.094
    #the distance in these tests of Salmonella to Citrobacter (closely related), was 0.159... and some shigella and E. coli were 0.165 to 0.18... will set to 0.16

#def closest_distances_parsnpfailed(ref_fasta, fasta_folder):
#    close_fastas = list()
    # fastqs = glob.glob(os.path.join(fastq_folder, '*R1*'))
    #mash.sketch(os.path.join(fasta_folder, '*.fasta'), output_sketch=os.path.join(fasta_folder, 'sketch.msh'), threads=46)
    #mash.dist(os.path.join(fasta_folder, 'sketch.msh'), ref_fasta, threads=46, output_file=os.path.join(fasta_folder, 'distances.tab'))
    #mash_output = mash.read_mash_output(os.path.join(fasta_folder, 'distances.tab'))
#    mashdf = pd.read_csv(os.path.join(fasta_folder, 'distances.tab'), sep='\t', names=["Reference-ID","Query-ID","MASH-Distance","P-value","Matching-hashes"])
    #mashsorted = mashdf.set_index('MASH-Distance',ascending=True).groupby(
#    mashsorted = mashdf.set_index(['Reference-ID']).groupby('Query-ID')
#    for item in mashdf:
#        print(item.reference, item.query) #trying to troubleshoot string to float issue
#        print(item.reference, item.query, str(item.distance))
        #print(item.reference, item.query, float(item.distance))
#        if item.distance < 0.20:  # May need to adjust this value. Ashley adjusted it from 0.06... Adjusted to 0.20 after it failed for an enterococcus isolate
#            close_fastas.append(item.reference)
#    return close_fastas

def make_ref(input_file, output_file):
    """
    Parsnp doesn't like having multi-fastas as input files, so this method puts all contigs in a multi-fasta
    into a single contig.
    :param input_file: Path to your input multi-fasta
    :param output_file: Path to output fasta. Overwrites the file if it already exists.
    """
    contigs = SeqIO.parse(input_file, 'fasta')
    with open(output_file, 'w') as f:
        f.write('>reference\n')
        for s in contigs:
            f.write(str(s.seq) + '\n')

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
    unknownisolate_redmine()
