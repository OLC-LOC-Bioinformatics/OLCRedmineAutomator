import os
import glob
import click
import pickle
import shutil
from biotools import mash
from amrsummary import before_send
import sentry_sdk
from automator_settings import FTP_FOLDER, SENTRY_DSN
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
#import openpyxl
#from openpyxl import load_workbook




@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def metagenome_assembly_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Current list of analysis types that KMA can use
    #analyses = [
    #    'amr', 'biocide', 'metal', 'bacmet'
    #]
    # Current list of sequence types that this pipeline can analyse
    seqtypes = [
        'paired', 'minionfastq'
    ]

    # Variable to hold supplied arguments
    argument_dict = {
    #    'analysis': str(),
        'seqtype': 'paired',
        'kraken': False,
    }


    try:
        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        for item in description:
            item = item.upper().rstrip()
            if 'SEQTYPE' in item:
                argument_dict['seqtype'] = item.split('=')[1].lower()
                continue
            if 'KRAKEN' in item:
                argument_dict['kraken'] = True
                continue
            # if 'READCOUNT' in item:
            #     argument_dict['readcount'] = True
            #     continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No SEQIDs provided!',
                                          status_id=4)
            return

        # Set the database path for the analyses
        dbpath = '/mnt/nas2/databases/kma_v_1.4.2_db/'
        bacmetpath = '/mnt/nas2/databases/bacmet/'
        #icebergpath = '/mnt/nas2/databases/ICEberg/'
        database_path = {
        #    'custom': os.path.join(target_dir, 'targets_KMA'),
            'amr': os.path.join(dbpath, 'NCBI-AMR-v3.10'),
            'biocide': os.path.join(dbpath, 'NCBI-biocide-v3.10'),
            'metal': os.path.join(dbpath, 'NCBI-metal-v3.10'),
            'bacmet': os.path.join(bacmetpath, 'bacmet_2018-03-11-v2_renamed_kma'),
        }

        # create an output folder so we can delete all of the extra files we didn't need once analyses are finished
        out_dir = os.path.join(work_dir, 'output')
        os.makedirs(out_dir, exist_ok=True)

        # #create a folder to hold the sequences
        seq_dir = os.path.join(work_dir, 'sequences')
        os.makedirs(seq_dir, exist_ok=True)
        #retrieve the raw data files
        retrieve_nas_files(seqids=seqids, outdir=seq_dir, filetype='fastq', copyflag=False)
        missing_fastqs = check_fastqs_present(seqids, seq_dir)
        if len(missing_fastqs) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                ' the OLC NAS: {}'.format(missing_fastqs))

        # create a folder for assemblies
        assemblies_dir = os.path.join(work_dir, 'assemblies')
        os.makedirs(assemblies_dir, exist_ok=True)

        #create an output file to append assembly metrics to
        covgsummaryfile = os.path.join(work_dir,'Assembly_coverage_stats.csv')
        with open(covgsummaryfile, 'w') as file:
            csv_file = csv.writer(file, delimiter=',')
            with open(covgsummaryfile, 'a+') as file:
                file.write("SEQID,Average_coverage,SD,Mapped_reads,Mapped_bases,Ref_scaffolds,Ref_bases,Percent_mapped,"
                           "Percent_proper_pairs,Percent_scaffolds_with_any_coverage,Percent_of_reference_bases_covered\n")

        #assemble the metagenome using megahit
        activatemeg = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/megahit'
        activatebbmap = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/bbmap'
        bbwrap = '/mnt/nas2/virtual_environments/bbmap/bin/bbwrap.sh'
        bbmap = '/mnt/nas2/virtual_environments/bbmap/bin/bbwrap.sh'
        bbmappacbio = '/mnt/nas2/virtual_environments/bbmap/bin/mapPacBio.sh'
        bbpileup = '/mnt/nas2/virtual_environments/bbmap/bin/pileup.sh'

        # run analysis for each pair of fastq files
        if argument_dict['seqtype'] == 'paired':
            for rawread in glob.glob(os.path.join(seq_dir, '*_R1_001.fastq.gz')):
                seqid1 = os.path.split(rawread)[1].split('.')[0]
                seqid = os.path.split(seqid1)[1].split('_R')[0]
                # seqidforfile = os.path.split(seqid1)[1].split('_S')[0]
                forwardseq = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
                reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
                # prepare command for megahit
                megahit_cmdpe = 'megahit -1 {forward} -2 {reverse} -o {seqid} --out-prefix {seqid}'\
                    .format(forward=forwardseq,
                            reverse=reverseseq,
                            seqid=seqid)
                # Create another shell script to execute within the KMA conda environment
                template = "#!/bin/bash\n{} && cd {} && {}".format(activatemeg, seq_dir, megahit_cmdpe)
                megahit_script = os.path.join(work_dir, 'run_megahit.sh')
                with open(megahit_script, 'w+') as file:
                    file.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(megahit_script)
                # Run shell script
                os.system(megahit_script)

                #copy the assemblies to the assembly folder for use with KMA later
                assembly_file = '{seqid}.contigs.fa'.format(seqid=seqid)
                output_assembly = os.path.join(work_dir, 'assemblies', assembly_file)
                #copy the .fa file to the assemblies folder and the out_folder
                shutil.copyfile(os.path.join(seq_dir,seqid,assembly_file),
                                output_assembly)
                shutil.copyfile(os.path.join(seq_dir,seqid,assembly_file),
                                os.path.join(out_dir,assembly_file))

                #determine coverage of assembly using samtools, and print output to file
                #bbwrapcmd = '{bbwrap} ref={outassembly} in={fwr} in2={rev} out={seqid}.aln.sam.gz'\
                #    .format(bbwrap=bbwrap,outassembly=output_assembly,fwr=forwardseq,rev=reverseseq,seqid=seqid)
                bbmapcmd = '{bbwrap} ref={outassembly} in={fwr} in2={rev} out={seqid}.aln.sam.gz idtag ' \
                           'ordered slow scafstats={seqid}_scaffoldstats'\
                    .format(bbwrap=bbmap,outassembly=output_assembly,fwr=forwardseq,rev=reverseseq,seqid=seqid)
                pilecmd = '{pile} in={seqid}.aln.sam.gz out={seqid}_coverage.txt 2> {seqid}_covg_summary.txt'\
                    .format(pile=bbpileup,seqid=seqid)

                #create a shell script to execute within the environment
                templatecov = "#!/bin/bash\n{} && cd {} && {} && {}".format(activatebbmap, seq_dir, bbmapcmd, pilecmd)
                covg_script = os.path.join(work_dir, 'determine_coverage.sh')
                with open(covg_script, 'w+') as file:
                    file.write(templatecov)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(covg_script)
                # Run shell script
                os.system(covg_script)


                #parse the output coverage summary file, and append to the file created above
                #there is probably a better way to do this....
                covgfile=os.path.join(seq_dir,'{}_covg_summary.txt'.format(seqid))
                covgfilein=pandas.read_csv(covgfile,sep='\t', skiprows=3,header=None,engine='python')
                totalreads=covgfilein.loc[0][1]
                mappedreads=covgfilein.loc[1][1]
                mappedbases=covgfilein.loc[2][1]
                refbases=covgfilein.loc[4][1]
                refscaffolds=covgfilein.loc[3][1]
                percentmapped=covgfilein.loc[5][1]
                avgcovg=covgfilein.loc[7][1]
                avgcovgdels=covgfilein.loc[8][1]
                percproperpairs=covgfilein.loc[6][1]
                percscaffoldswithcoverage=covgfilein.loc[10][1]
                percrefbasescovered=covgfilein.loc[11][1]
                stdev=covgfilein.loc[9][1]
                with open(covgsummaryfile, 'a+') as file:
                    file.write('{},{},{},{},{},{},{},{},{},{},{}\n'.format(seqid,avgcovg,stdev,mappedreads,mappedbases,refscaffolds,refbases,percentmapped,percproperpairs,percscaffoldswithcoverage,percrefbasescovered))



        #run analysis if it is a single ended nanopore sequence
        if argument_dict['seqtype'] == 'minionfastq':
            #fastq minion files
            for rawread in glob.glob(os.path.join(seq_dir, '*.fastq.gz')):
                seqid = os.path.split(rawread)[1].split('.')[0]
                #prepare command for megahit
                megahit_cmd = 'megahit -r {rawread} -o {seqid} --out-prefix {seqid}'.format(rawread=rawread,seqid=seqid)

                # Create another shell script to execute within the megahit conda environment
                template = "#!/bin/bash\n{} && cd {} && {}".format(activatemeg, seq_dir, megahit_cmd)
                megahit_script = os.path.join(work_dir, 'run_megahit.sh')
                with open(megahit_script, 'w+') as file:
                    file.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(megahit_script)
                # Run shell script
                os.system(megahit_script)

                #copy the assembly files to the assembly folder
                assembly_file = '{seqid}.contigs.fa'.format(seqid=seqid)
                output_assembly = os.path.join(work_dir, 'assemblies', assembly_file)
                #copy the .fa file to the assemblies folder and outputfolder
                shutil.copyfile(os.path.join(seq_dir,seqid,assembly_file),
                                output_assembly)
                shutil.copyfile(os.path.join(seq_dir,seqid,assembly_file),
                                os.path.join(out_dir, assembly_file))

                #determine coverage of assembly using samtools, and print output to file
                #for the mappacbio function of bbmap, you can use a maxlength of up to 6000, but this site suggests
                #https://www.biostars.org/p/483592/
                #that using 1000 resulted in a higher mapping rate
                bbwrapcmd = '{bbwrap} ref={outassembly} maxlen=1000 minlen=200 in={raw} out={seqid}.aln.sam.gz ' \
                            'idtag ignorebadquality ordered slow scafstats={seqid}_scaffoldstats'\
                    .format(bbwrap=bbmappacbio,outassembly=output_assembly,raw=rawread,seqid=seqid)
                pilecmd = '{pile} in={seqid}.aln.sam.gz out={seqid}_coverage.txt 2> {seqid}_covg_summary.txt'\
                    .format(pile=bbpileup,seqid=seqid)

                #create a shell script to execute within the environment
                templatecov = "#!/bin/bash\n{} && cd {} && {} && {}".format(activatebbmap, seq_dir, bbwrapcmd, pilecmd)
                covg_script = os.path.join(work_dir, 'determine_coverage.sh')
                with open(covg_script, 'w+') as file:
                    file.write(templatecov)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(covg_script)
                # Run shell script
                os.system(covg_script)

                #TODO: FIX THIS
                #parse the output coverage summary file
                covgfile=os.path.join(seq_dir,'{}_covg_summary.txt'.format(seqid))
                covgfilein=pandas.read_csv(covgfile,sep='\t', skiprows=3,header=None,engine='python')
                totalreads=covgfilein.loc[0][1]
                mappedreads=covgfilein.loc[1][1]
                mappedbases=covgfilein.loc[2][1]
                refbases=covgfilein.loc[4][1]
                refscaffolds=covgfilein.loc[3][1]
                percentmapped=covgfilein.loc[5][1]
                avgcovg=covgfilein.loc[7][1]
                avgcovgdels=covgfilein.loc[8][1]
                percproperpairs=covgfilein.loc[6][1]
                percscaffoldswithcoverage=covgfilein.loc[10][1]
                percrefbasescovered=covgfilein.loc[11][1]
                stdev=covgfilein.loc[9][1]
                with open(covgsummaryfile, 'a+') as file:
                    file.write('{},{},{},{},{},{},{},{},{},{},{}\n'.format(seqid,avgcovg,stdev,mappedreads,mappedbases,refscaffolds,refbases,percentmapped,percproperpairs,percscaffoldswithcoverage,percrefbasescovered))

        #write the assembly statistics to an excel file
        covgdf = pandas.read_csv(covgsummaryfile, sep=',', engine='python')
        covgexcelfile = os.path.join(out_dir, 'Metagenome_Assembly_stats_redmine{}.xlsx'.format(issue.id))
        with pandas.ExcelWriter(covgexcelfile) as writer:
            covgdf.to_excel(writer,sheet_name='Assembly_metrics', index=False)


        #run KMA on the metagenome assemblies and raw data. put it into an excel file
        #I call a different kma python script, listed as kmapy below
        #TODO: add functionality for custom KMA analysis?
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kma_v149'
        kmapy = '/mnt/nas2/redmine/applications/OLCRedmineAutomator/automators/kma_metagenomes_analysis.py'
        kmacmd = 'python {kmapy} -a {assm} -w {workd} -s {seqdir} -type {t} -i {inp}'\
                  .format(kmapy=kmapy, assm=assemblies_dir, workd=work_dir, seqdir=seq_dir,
                          t=argument_dict['seqtype'], inp=issue.id)
        templatekma = "#!/bin/bash\n{} && {}".format(activate, kmacmd)
        kma_script = os.path.join(work_dir, 'kma_python.sh')
        with open(kma_script, 'w+') as file:
            file.write(templatekma)
        make_executable(kma_script)
        # Run shell script
        os.system(kma_script)
        

        #run metaphlan4
        #make a metaphlan folder
        metaphlan_fld = os.path.join(out_dir, 'Metaphlan4')
        os.makedirs(metaphlan_fld, exist_ok=True)
        metaphlan_asmbls = os.path.join(metaphlan_fld, 'Assembly_analysis')
        os.makedirs(metaphlan_asmbls, exist_ok=True)
        metaphlan_raw = os.path.join(metaphlan_fld, 'RawRead_analysis')
        os.makedirs(metaphlan_raw, exist_ok=True)
        #set the database path for the ChocoPhlan database (used by metaphlan)
        dbpath = '/mnt/nas2/databases/metaphlan4/'

        #Commands required for metaphlan
        activateemetaphlan = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/metaphlan4'

        for assembly in glob.glob(os.path.join(assemblies_dir, '*.contigs.fa')):
            seqid = os.path.split(assembly)[1].split('.')[0]
            bt2out = '{seqid}.bowtie2.bz2'.format(seqid=seqid)
            outputfile = '{seqid}_metaphlan4_assembly_report'.format(seqid=seqid)
            analysis='rel_ab_w_read_stats'
            #prepare command
            metaphlan4_cmd = 'metaphlan --input_type fasta --bowtie2db {db} --nproc 24 --bowtie2out {bt2o} -o {out}'\
                             ' -s {seq} -t {anls} --sample_id {seqid}'\
                             ' {assembly}'.format(db=dbpath,
                                                  bt2o=bt2out,out=outputfile,
                                                  seq=seqid,anls=analysis, seqid=seqid,
                                                  assembly=assembly)

            # Create another shell script to execute within the pyseer conda environment
            template = "#!/bin/bash\n{ntasks}\n{mem}\n{activate} && cd {seqdir} && {mpa}"\
                .format(ntasks="#SBATCH --ntasks 30",mem="#SBATCH --mem=190000",
                        activate=activateemetaphlan,seqdir=assemblies_dir,
                        mpa=metaphlan4_cmd)
            mpa4_script = os.path.join(work_dir, 'run_metaphlan4.sh')
            with open(mpa4_script, 'w+') as file:
                file.write(template)
            make_executable(mpa4_script)

            # Run shell script
            os.system(mpa4_script)

            #now move the assembly metaphlan files to the metaphlan assembly folder
            output_files = os.listdir(assemblies_dir)
            for file in output_files:
                if file.endswith("_assembly_report"):
                    shutil.move(os.path.join(assemblies_dir,file), os.path.join(metaphlan_asmbls, file))

        #merge metaphlan tables for assembled data
        #metaphlan_merge = '/mnt/nas2/virtual_environments/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/merge_metaphlan_tables.py'
        mergecmd = 'merge_metaphlan_tables.py *_assembly_report > merged_abundance_table_assembled.txt'
        template4 = "#!/bin/bash\n{activate} && cd {metaphlandir} && {mpa}"\
            .format(activate=activateemetaphlan,metaphlandir=metaphlan_asmbls,
                    mpa=mergecmd)
        mpa4m_script = os.path.join(work_dir, 'merge_metaphlan4.sh')
        with open(mpa4m_script, 'w+') as file:
            file.write(template4)
        make_executable(mpa4m_script)

        # Run shell script
        os.system(mpa4m_script)

        #run metaphlan on raw reads
        for metagenome in glob.glob(os.path.join(seq_dir, '*.fastq.gz')):
            seqid1 = os.path.split(metagenome)[1].split('.')[0]
            seqid = os.path.split(seqid1)[1].split('_R')[0]
            samoutfile = '{}_sam_metaphlan'.format(seqid)
            forwardseq = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
            reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
            bt2out = '{seqid}_metaphlan.bowtie2.bz2'.format(seqid=seqid)
            outputfile = '{seqid}_metaphlan4_rawread_report'.format(seqid=seqid)
            analysis='rel_ab_w_read_stats'
            #prepare command
            metaphlan4_cmd = 'metaphlan --input_type fastq --bowtie2db {db} --nproc 24 --bowtie2out {bt2o} -o {out}'\
                             ' -s {sam} -t {anls} --sample_id {rid}'\
                             ' {forwardseq},{reverseseq}'.format(db=dbpath,
                                                                 bt2o=bt2out,out=outputfile,
                                                                 sam=samoutfile,anls=analysis, rid=seqid,
                                                                 forwardseq=forwardseq, reverseseq=reverseseq)

            # Create another shell script to execute within the pyseer conda environment
            template = "#!/bin/bash\n{ntasks}\n{mem}\n{activate} && cd {seqdir} && {mpa}"\
                .format(ntasks="#SBATCH --ntasks 30",mem="#SBATCH --mem=190000",
                        activate=activateemetaphlan,seqdir=seq_dir,
                        mpa=metaphlan4_cmd)
            mpa4_script = os.path.join(work_dir, 'run_metaphlan4.sh')
            with open(mpa4_script, 'w+') as file:
                file.write(template)
            make_executable(mpa4_script)

            # Run shell script
            os.system(mpa4_script)

        #now move the rawread metaphlan files to the metaphlan rawread folder
        output_files = os.listdir(seq_dir)
        for file in output_files:
            if file.endswith("_rawread_report"):
                shutil.move(os.path.join(seq_dir,file), os.path.join(metaphlan_raw, file))

        #merge metaphlan tables for raw data
        #metaphlan_merge = '/mnt/nas2/virtual_environments/metaphlan4/lib/python3.7/site-packages/metaphlan/utils/merge_metaphlan_tables.py'
        mergecmd = 'merge_metaphlan_tables.py *_rawread_report > merged_abundance_table_raw.txt'

        #also merged the mapped reads data from the raw analyses
        mergereads = '/mnt/nas2/redmine/applications/OLCRedmineAutomator/automators/merge_metaphlanestimatedreadcounts_tables.py'
        mergereadscmd = 'python {} *_rawread_report -o merged_reads_table_raw.txt'.format(mergereads)

        template4 = "#!/bin/bash\n{activate} && cd {metaphlandir} && {mpa} && {mr}"\
            .format(activate=activateemetaphlan,metaphlandir=metaphlan_raw,
                    mpa=mergecmd, mr=mergereadscmd)
        mpa4m_script = os.path.join(work_dir, 'merge_metaphlan4.sh')
        with open(mpa4m_script, 'w+') as file:
            file.write(template4)
        make_executable(mpa4m_script)

        # Run shell script
        os.system(mpa4m_script)

        #write the metaphlan results to an excel file
        metarawdf = pandas.read_csv(os.path.join(metaphlan_raw,'merged_abundance_table_raw.txt'), sep=',|\t',skiprows=1, engine='python')
        metarawdfreads = pandas.read_csv(os.path.join(metaphlan_raw,'merged_reads_table_raw.txt'), sep=',|\t',engine='python')
        metaambldf = pandas.read_csv(os.path.join(metaphlan_asmbls,'merged_abundance_table_assembled.txt'), sep=',|\t',skiprows=1, engine='python')
        metaphlanexcelfile = os.path.join(metaphlan_fld, 'Metaphlan4_results_redmine{}.xlsx'.format(issue.id))
        with pandas.ExcelWriter(metaphlanexcelfile) as writer:
            metarawdf.to_excel(writer,sheet_name='Rawreads_abundance', index=False)
            metaambldf.to_excel(writer, sheet_name='Assembly_abundance', index=False)
            metarawdfreads.to_excel(writer,sheet_name='Rawreads_readcounts', index=False)

        #run kraken2 and bracken on the files if the user requests it
        if argument_dict['kraken']:
            activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kraken2'
            krknpy = '/mnt/nas2/redmine/applications/OLCRedmineAutomator/automators/kraken2_metagenomeassembly_pipeline.py'
            krkncmd = 'python {krknpy} -w {workd} -s {seqdir} -type {t} -i {iss}'\
                      .format(krknpy=krknpy, workd=work_dir, seqdir=seq_dir,
                              t=argument_dict['seqtype'], iss=issue.id)
            templatekrkn = "#!/bin/bash\n{} && {}".format(activate, krkncmd)
            krkn_script = os.path.join(work_dir, 'krkn_python.sh')
            with open(krkn_script, 'w+') as file:
                file.write(templatekrkn)
            make_executable(krkn_script)
            # Run shell script
            os.system(krkn_script)


        #now copy the output excel file to the new output directory... leaving as a directory so we can zip it if we add other functions and outputs later
        kmafilename = 'KMA_results_redmine{}.xlsx'.format(issue.id)
        kmareadsfilename = 'KMA_readcounts_redmine{}.xlsx'.format(issue.id)
        metafilename = 'Metaphlan4_results_redmine{}.xlsx'.format(issue.id)
        covgfilename = 'Metagenome_Assembly_stats_redmine{}.xlsx'.format(issue.id)
        shutil.copyfile(os.path.join(assemblies_dir, kmafilename), os.path.join(out_dir, kmafilename))
        shutil.copyfile(os.path.join(assemblies_dir, kmareadsfilename), os.path.join(out_dir, kmareadsfilename))
        shutil.copyfile(os.path.join(metaphlan_fld, metafilename), os.path.join(out_dir, metafilename))

        #upload some files to the redmine request (the results excel files, the rest will go to the ftp due to size)
        output_list = list()
        output_dict = dict()
        output_dict['path'] = os.path.join(out_dir,kmafilename)
        output_dict['filename'] = kmafilename
        output_list.append(output_dict)

        output_dict = dict()
        output_dict['path'] = os.path.join(out_dir,kmareadsfilename)
        output_dict['filename'] = kmareadsfilename
        output_list.append(output_dict)

        output_dict = dict()
        output_dict['path'] = os.path.join(out_dir,metafilename)
        output_dict['filename'] = metafilename
        output_list.append(output_dict)

        output_dict = dict()
        output_dict['path'] = os.path.join(out_dir,covgfilename)
        output_dict['filename'] = covgfilename
        output_list.append(output_dict)

        #if user asked for kraken2 and bracken analysis, attach these files as well
        if argument_dict['kraken']:
            krakenfilename='Kraken2andBracken_results_redmine{}.xlsx'.format(issue.id)
            output_dict = dict()
            output_dict['path'] = os.path.join(out_dir,'kraken2',krakenfilename)
            output_dict['filename'] = krakenfilename
            output_list.append(output_dict)


        #
        # Zip output
        out_filename = 'metagenome_assembly_pipeline_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=out_dir,
                                  output_dir=work_dir,
                                  output_filename=out_filename)
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
                notes='Metagenome assembly and analysis complete!\n\n'
                      'KMA and Metaphlan4 results have been attached to your '
                      'request.\nFull results are available at the following '
                      'URL:\n {url}'.format(url=download_link)
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
        shutil.rmtree(assemblies_dir)
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
    metagenome_assembly_redmine()
