import os
import glob
import shutil
import fileinput
from pathlib import Path
import csv
import pandas
import argparse
import argcomplete

#descriptions and argument assignment
parser = argparse.ArgumentParser(description='details',
         usage='use "%(prog)s --help" for more information',
         formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--info', default=None,
                    help='''
                    For use with redmine automator to create GROBI reports
                    ''')
#parser.add_argument('-a', dest='assemblies_dir', type=str, required=True, help="""Assemblies Directory""")
parser.add_argument('-w', dest='work_dir', type=str, required=True, help="""Working Directory""")
parser.add_argument('-s', dest='seq_dir', type=str, required=True, help="""Metagenome Sequences Directory""")
parser.add_argument('-type', dest='seqtype', type=str, required=True, help="""Sequence Type (fastq or minionfastq)""")
parser.add_argument('-i', dest='issueid', type=str, required=True, help="""Redmine issue id""")
argcomplete.autocomplete(parser)
args = parser.parse_args()

def make_executable(path):
    """
    Takes a shell script and makes it executable (chmod +x)
    :param path: path to shell script
    """
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2
    os.chmod(path, mode)

#set the database path for the analyses
dbpath = '/mnt/nas2/databases/kraken2/'
database_path = {
    'kraken2': os.path.join(dbpath, 'k2_standard_20230605'),
    'greengenes': os.path.join(dbpath, 'greengenes'),
    'plusPF': os.path.join(dbpath, 'plusPF_20230605'),
    'rdp': os.path.join(dbpath, 'rdp'),
    'silva': os.path.join(dbpath, 'silva')
    }

# These unfortunate hard coded paths appear to be necessary
activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kraken2'

kraken2_folder = os.path.join(args.work_dir, 'output','kraken2')
os.makedirs(kraken2_folder, exist_ok=True)

if args.seqtype == 'paired':
    for metagenome in glob.glob(os.path.join(args.seq_dir, '*_R1_001.fastq.gz')):
        seqid1 = os.path.split(metagenome)[1].split('.')[0]
        seqid = os.path.split(seqid1)[1].split('_R')[0]
        forwardseq = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
        reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
        #prepare command
        classified_out = '{seqid}#_classified'.format(seqid=seqid)
        unclassified_out = '{seqid}#_unclassified'.format(seqid=seqid)
        report_out = '{seqid}_report'.format(seqid=seqid)
        kraken2_cmd = 'kraken2 --db {database} --threads 24 --paired --output {outdir} --gzip-compressed'\
                      ' --classified-out {outdir}/{outfilec} --unclassified-out {outdir}/{outfileu}'\
                      ' --report {outdir}/{reportout} --report-zero-counts'\
                      ' {forwardseq} {reverseseq}'.format(database=database_path['plusPF'],
                                                          outdir=kraken2_folder,
                                                          outfilec=classified_out, outfileu=unclassified_out,
                                                          reportout=report_out,
                                                          forwardseq=forwardseq, reverseseq=reverseseq)

        # Create another shell script to execute within the pyseer conda environment
        template = "#!/bin/bash\n{ntasks}\n{mem}\n{activate} && cd {seqdir} && {kraken2}"\
            .format(ntasks="#SBATCH --ntasks 30",mem="#SBATCH --mem=190000",
                    activate=activate,seqdir=args.seq_dir,
                    kraken2=kraken2_cmd)
        kraken2_script = os.path.join(args.work_dir, 'run_kraken2.sh')
        with open(kraken2_script, 'w+') as file:
            file.write(template)
        make_executable(kraken2_script)

        # Run shell script
        os.system(kraken2_script)

    for metagenome in glob.glob(os.path.join(args.seq_dir, '*_R1.fastq.gz')):
        seqid1 = os.path.split(metagenome)[1].split('.')[0]
        seqid = os.path.split(seqid1)[1].split('_R')[0]
        forwardseq = '{seqid}_R1.fastq.gz'.format(seqid=seqid)
        reverseseq = '{seqid}_R2.fastq.gz'.format(seqid=seqid)
        #prepare command
        classified_out = '{seqid}#_classified'.format(seqid=seqid)
        unclassified_out = '{seqid}#_unclassified'.format(seqid=seqid)
        report_out = '{seqid}_report'.format(seqid=seqid)
        kraken2_cmd = 'kraken2 --db {database} --threads 24 --paired --output {outdir} --gzip-compressed'\
                      ' --classified-out {outdir}/{outfilec} --unclassified-out {outdir}/{outfileu}'\
                      ' --report {outdir}/{reportout} --report-zero-counts'\
                      ' {forwardseq} {reverseseq}'.format(database=database_path['plusPF'],
                                                          outdir=kraken2_folder,
                                                          outfilec=classified_out, outfileu=unclassified_out,
                                                          reportout=report_out,
                                                          forwardseq=forwardseq, reverseseq=reverseseq)

        # Create another shell script to execute within the pyseer conda environment
        template = "#!/bin/bash\n{ntasks}\n{mem}\n{activate} && cd {seqdir} && {kraken2}"\
            .format(ntasks="#SBATCH --ntasks 30",mem="#SBATCH --mem=190000",
                    activate=activate,seqdir=args.seq_dir,
                    kraken2=kraken2_cmd)
        kraken2_script = os.path.join(args.work_dir, 'run_kraken2.sh')
        with open(kraken2_script, 'w+') as file:
            file.write(template)
        make_executable(kraken2_script)

        # Run shell script
        os.system(kraken2_script)
#TODO: elif args.seqtype == 'minionfastq':

#run bracken on kraken2 reports
for report in glob.glob(os.path.join(kraken2_folder, '*_report')):
    reportname = os.path.split(report)[1].split('_re')[0]
    outname = '{}_bracken'.format(reportname)
    activatebracken = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/bracken'
    bracken_cmd = 'bracken -d {db} -i {report} -o {outname} -l S'.format(db=database_path['plusPF'],report=report,
                                                                         outname=outname)
    #create another shell script to run the bracken analysis of the reports
    templatebr = "#!/bin/bash\n{activate} && cd {rdir} && {bracken}" \
        .format(activate=activatebracken, rdir=kraken2_folder, bracken=bracken_cmd)
    bracken_script = os.path.join(args.work_dir, 'run_bracken.sh')
    with open(bracken_script, 'w+') as file:
        file.write(templatebr)
    make_executable(bracken_script)

    # Run shell script
    os.system(bracken_script)

#now convert the kraken and bracken reports to mpa format (like metaphlan) just because they are easier to merge and view for users
activatektools = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/krakentools'
kraken2mpa = '/mnt/nas2/virtual_environments/krakentools/bin/kreport2mpa.py'
mergekrakenmpas = '/mnt/nas2/virtual_environments/krakentools/bin/combine_mpa.py'

#convert the report files into mpa format
for krknreport in glob.glob(os.path.join(kraken2_folder, '*_report')):
    reportname = os.path.split(report)[1].split('_re')[0]
    outname = '{}_kraken_mpa'.format(reportname)
    outnameab = '{}_kraken_abnd_mpa'.format(reportname)
    convertkmpa = 'python {k2m} -r {report} -o {outname} --display-header --read_count'.format(k2m=kraken2mpa, report=krknreport,
                                                                              outname=outname)
    convertkmpa2 = 'python {k2m} -r {report} -o {outname} --display-header --percentages'.format(k2m=kraken2mpa, report=krknreport,
                                                                                outname=outnameab)
    template = "#!/bin/bash\n{activate} && cd {rdir} && {cmd} && {cmd2}" \
        .format(activate=activatektools, rdir=kraken2_folder, cmd=convertkmpa, cmd2=convertkmpa2)
    krkmpascript = os.path.join(args.work_dir, 'kraken2mpa.sh')
    with open(krkmpascript, 'w+') as file:
        file.write(template)
    make_executable(krkmpascript)
    # Run shell script
    os.system(krkmpascript)

for brknreport in glob.glob(os.path.join(kraken2_folder, '*_report_bracken_species')):
    reportname = os.path.split(report)[1].split('_br')[0]
    outname = '{}_bracken_mpa'.format(reportname)
    outnameab = '{}_bracken_abnd_mpa'.format(reportname)
    convertkmpa = 'python {k2m} -r {report} -o {outname}  --display-header --read_count'.format(k2m=kraken2mpa, report=brknreport, outname=outname)
    convertkmpa2 = 'python {k2m} -r {report} -o {outname}  --display-header --percentages'.format(k2m=kraken2mpa, report=brknreport, outname=outnameab)
    template = "#!/bin/bash\n{activate} && cd {rdir} && {cmd} && {cmd2}" \
        .format(activate=activatektools, rdir=kraken2_folder, cmd=convertkmpa, cmd2=convertkmpa2)
    krkmpascript = os.path.join(args.work_dir, 'kraken2mpa.sh')
    with open(krkmpascript, 'w+') as file:
        file.write(template)
    make_executable(krkmpascript)
    # Run shell script
    os.system(krkmpascript)

#merge the mpa reports together for kraken and bracken
mergekrk = 'python {} -i *_kraken_mpa -o merged_kraken_readcounts'.format(mergekrakenmpas)
mergeabkrk = 'python {} -i *_kraken_abnd_mpa -o merged_kraken_abundance'.format(mergekrakenmpas)
mergebrk = 'python {} -i *_bracken_mpa -o merged_braken_readcounts'.format(mergekrakenmpas)
mergeabbrk = 'python {} -i *_bracken_abnd_mpa -o merged_braken_abundance'.format(mergekrakenmpas)
#comand for cluster
template3 = "#!/bin/bash\n{activate} && cd {rdir} && {cmd} && {cmd2} && {cmd3} && {cmd4}".format(activate=activatektools,
                                                                                                 rdir=kraken2_folder,
                                                                                                 cmd=mergekrk,
                                                                                                 cmd2=mergeabkrk,
                                                                                                 cmd3=mergebrk,
                                                                                                 cmd4=mergeabbrk)
mergempascript = os.path.join(args.work_dir, 'mergekrakenmpa.sh')
with open(mergempascript, 'w+') as file:
    file.write(template3)
make_executable(mergempascript)
# Run shell script
os.system(mergempascript)

#convert the outputs to csv objects so we can add them to excel file
krkmpa = pandas.read_csv(os.path.join(kraken2_folder,'merged_kraken_readcounts'),sep=',|\t', engine='python')
krkabmpa = pandas.read_csv(os.path.join(kraken2_folder,'merged_kraken_abundance'),sep=',|\t', engine='python')
brkmpa = pandas.read_csv(os.path.join(kraken2_folder,'merged_braken_readcounts'),sep=',|\t', engine='python')
brkabmpa = pandas.read_csv(os.path.join(kraken2_folder,'merged_braken_abundance'),sep=',|\t', engine='python')

#write the kraken and bracken analysis results to a dataframe
krakenexcelfile = os.path.join(kraken2_folder, 'Kraken2andBracken_results_redmine{}.xlsx'.format(args.issueid))
with pandas.ExcelWriter(krakenexcelfile) as writer:
    krkmpa.to_excel(writer,sheet_name='Kraken2_readcounts', index=False)
    krkabmpa.to_excel(writer,sheet_name='Kraken2_abundance', index=False)
    brkmpa.to_excel(writer,sheet_name='Bracken_readcounts', index=False)
    brkabmpa.to_excel(writer,sheet_name='Bracken_abundance', index=False)

