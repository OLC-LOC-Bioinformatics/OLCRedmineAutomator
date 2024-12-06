import os
import glob
#import click
#import pickle
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
parser.add_argument('-a', dest='assemblies_dir', type=str, required=True, help="""Assemblies Directory""")
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

#database locations
#dbpath = '/mnt/nas2/databases/kma_v_1.4.2_db/'
dbpath = '/mnt/nas2/databases/kma_v1.4.9/'
bacmetpath = '/mnt/nas2/databases/bacmet/'
database_path = {
    #'custom': os.path.join(target_dir, 'targets_KMA'),
    'amr': os.path.join(dbpath, 'NCBI-AMR'),
    'biocide': os.path.join(dbpath, 'NCBI-BIOCIDE'),
    'metal': os.path.join(dbpath, 'NCBI-METAL'),
    'bacmet': os.path.join(bacmetpath, 'bacmet_2018-03-11-v2_renamed_kma'),
    }

# These unfortunate hard coded paths appear to be necessary
activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/kma_v149'
#kma_py = '/mnt/nas2/virtual_environments/card-rgi/bin/kma'

# # Run kma with the necessary arguments
# #assembly files AMR
for assembly in glob.glob(os.path.join(args.assemblies_dir, '*.contigs.fa')):
    seqid = os.path.split(assembly)[1].split('.')[0]
    #prepare command for kma
    kma_cmd = 'kma -i {assembly} -o {seqid}_amr -t_db {dbpath} -t 7 -nf'\
        .format(assembly=assembly,
                seqid=seqid,
                dbpath=database_path['amr'])

    # Append the align and/or the unique flags are required
    kma_cmd += ' -bcNano' if args.seqtype == 'minionfastq' else ''

    #prepare command for kma_metals
    kmametal_cmd = 'kma -i {assembly} -o {seqid}_metal -t_db {dbpath} -t 7 -nf'\
        .format(assembly=assembly,
                seqid=seqid,
                dbpath=database_path['metal'])

    # Append the align and/or the unique flags are required
    kmametal_cmd += ' -bcNano' if args.seqtype == 'minionfastq' else ''

    #prepare command for kma_biocides
    kmabiocide_cmd = 'kma -i {assembly} -o {seqid}_biocide -t_db {dbpath} -t 7 -nf'\
        .format(assembly=assembly,
                seqid=seqid,
                dbpath=database_path['biocide'])

    # Append the align and/or the unique flags are required
    kmabiocide_cmd += ' -bcNano' if args.seqtype == 'minionfastq' else ''

    # Create another shell script to execute within the KMA conda environment
    template = "#!/bin/bash\n{} && cd {} && {} && {} && {}".format(activate, args.assemblies_dir, kma_cmd, kmametal_cmd, kmabiocide_cmd)
    kma_script = os.path.join(args.work_dir, 'run_kma.sh')
    with open(kma_script, 'w+') as file:
        file.write(template)
    # Modify the permissions of the script to allow it to be run on the node
    make_executable(kma_script)
    # Run shell script
    os.system(kma_script)

#create a header list for when we combine the .res files into a single csv file
headerList = ['SeqID','Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value']

#add the filename (which is the seqid) to the first column in all of the res files, then concatenate into a single output csv file
outputfile = os.path.join(args.assemblies_dir, 'kma_amr_output.csv')
with open(outputfile, 'w', newline='') as file_output:
    csv_output = csv.writer(file_output)
    #gather info from each individual res file
    for fname in glob.glob(os.path.join(args.assemblies_dir, '*_amr.res')):
        fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
        seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
        seqname2 = os.path.split(seqname)[1].split('_')[0] #to remove the _SXX_L001 from seqid
        #first, find and replace the #Template column header
        with open(fname, 'r') as file:
            #read in the file
            data = file.read()
            #search and replace the text
            data = data.replace("#Template", "Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass")
        #open the file in write only mode to write the replaced content
        with open(fname, 'w') as file:
            file.write(data)
        with open(fname, newline='') as f_input:
            csv_input = csv.reader(f_input)
            #next(csv_input) #this line will allow us to skip the first line (header) of each res file, which starts with #Template
            for row in csv_input:
                print(row)
                row.insert(0,seqname2) #this adds the seqid to the file before concatenating
                csv_output.writerow(row)
df = pandas.read_csv(outputfile, sep=',|\t', header=None)
df.to_csv(outputfile, header=headerList, index=False)

#metal resistance file
outputfilemetal = os.path.join(args.assemblies_dir, 'kma_metal_output.csv')
with open(outputfilemetal, 'w', newline='') as file_output:
    csv_output = csv.writer(file_output)
    #gather info from each individual res file
    for fname in glob.glob(os.path.join(args.assemblies_dir, '*_metal.res')):
        fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
        seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
        seqname2 = os.path.split(seqname)[1].split('_')[0] #to remove the _SXX_L001 from seqid
        #first, find and replace the #Template column header
        with open(fname, 'r') as file:
            #read in the file
            data = file.read()
            #search and replace the text
            data = data.replace("#Template", "Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass")
        #open the file in write only mode to write the replaced content
        with open(fname, 'w') as file:
            file.write(data)
        with open(fname, newline='') as f_input:
            csv_input = csv.reader(f_input)
            #next(csv_input) #this line will allow us to skip the first line (header) of each res file, which starts with #Template
            for row in csv_input:
                print(row)
                row.insert(0,seqname2) #this adds the seqid to the file before concatenating
                #row.insert(0,fname)
                csv_output.writerow(row)
df = pandas.read_csv(outputfilemetal, sep=',|\t', header=None)
df.to_csv(outputfilemetal, header=headerList, index=False)

#biocide resistance file
outputfilebiocide = os.path.join(args.assemblies_dir, 'kma_biocide_output.csv')
with open(outputfilebiocide, 'w', newline='') as file_output:
    csv_output = csv.writer(file_output)
    #gather info from each individual res file
    for fname in glob.glob(os.path.join(args.assemblies_dir, '*_biocide.res')):
        fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
        seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
        seqname2 = os.path.split(seqname)[1].split('_')[0] #to remove the _SXX_L001 from seqid
        #first, find and replace the #Template column header
        with open(fname, 'r') as file:
            #read in the file
            data = file.read()
            #search and replace the text
            data = data.replace("#Template", "Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass")
        #open the file in write only mode to write the replaced content
        with open(fname, 'w') as file:
            file.write(data)
        with open(fname, newline='') as f_input:
            csv_input = csv.reader(f_input)
            #next(csv_input) #this line will allow us to skip the first line (header) of each res file, which starts with #Template
            for row in csv_input:
                print(row)
                row.insert(0,seqname2) #this adds the seqid to the file before concatenating
                #row.insert(0,fname)
                csv_output.writerow(row)
df = pandas.read_csv(outputfilebiocide, sep=',|\t', header=None)
df.to_csv(outputfilebiocide, header=headerList, index=False)

#lets modify the output kma files so we end up with an excel table
kmadf2 = pandas.read_csv(os.path.join(args.assemblies_dir,'kma_amr_output.csv'), sep=',|\t', engine='python')
#split the fasta header renamed column
kmadf2[['Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass']] = kmadf2['Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass'].str.split('|',expand=True)
#reorder the dataframe
kmadf3 = kmadf2[['SeqID','Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value','Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass']]
#drop the unwanted column
kmadf3.drop('Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass', axis=1, inplace=True)

#repeat for metal resistance files
kmadfm2 = pandas.read_csv(os.path.join(args.assemblies_dir,'kma_metal_output.csv'), sep=',|\t', engine='python')
#split the newly renamed column
kmadfm2[['Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass']] = kmadfm2['Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass'].str.split('|',expand=True)
#reorder the dataframe
kmadfm3 = kmadfm2[['SeqID','Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value','Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass']]
#drop the unwanted column
kmadfm3.drop('Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass', axis=1, inplace=True)

#and again for biocide resistance files
kmadfb2 = pandas.read_csv(os.path.join(args.assemblies_dir,'kma_biocide_output.csv'), sep=',|\t', engine='python')
#split the newly renamed column
kmadfb2[['Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass']] = kmadfb2['Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass'].str.split('|',expand=True)
#reorder the dataframe
kmadfb3 = kmadfb2[['SeqID','Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value','Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass']]
#drop the unwanted column
kmadfb3.drop('Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass', axis=1, inplace=True)


#repeat all of the above, but for raw data
if args.seqtype != 'minionfastq':
    for metagenome in glob.glob(os.path.join(args.seq_dir, '*.fastq.gz')):
        seqid1 = os.path.split(metagenome)[1].split('.')[0]
        seqid = os.path.split(seqid1)[1].split('_R')[0]
        forwardseq = '{seqid}_R1_001.fastq.gz'.format(seqid=seqid)
        reverseseq = '{seqid}_R2_001.fastq.gz'.format(seqid=seqid)
        kma_amrpe = 'kma -ipe {forward} {reverse} -o {seqid}_amr -t_db {db} -t 7 -nf -nc -ef'\
            .format(forward=forwardseq, reverse=reverseseq, seqid=seqid, db=database_path['amr'])
        kma_metrpe = 'kma -ipe {forward} {reverse} -o {seqid}_metal -t_db {db} -t 7 -nf -nc -ef'\
            .format(forward=forwardseq, reverse=reverseseq, seqid=seqid, db=database_path['metal'])
        kma_biocpe = 'kma -ipe {forward} {reverse} -o {seqid}_biocide -t_db {db} -t 7 -nf -nc -ef'\
            .format(forward=forwardseq, reverse=reverseseq, seqid=seqid, db=database_path['biocide'])
        # Create another shell script to execute within the KMA conda environment
        template_paired = "#!/bin/bash\n{} && cd {} && {} && {} && {}".format(activate, args.seq_dir, kma_amrpe, kma_metrpe, kma_biocpe)
        kma_script_paired = os.path.join(args.work_dir, 'run_kma_pe.sh')
        with open(kma_script_paired, 'w+') as file:
            file.write(template_paired)
        # Modify the permissions of the script to allow it to be run on the node
        make_executable(kma_script_paired)
        # Run shell script
        os.system(kma_script_paired)

#TODO: elif args.seqtype == 'minionfastq':
#    for metagenome in glob.glob(os.path.join(args.seq_dir, '*.fastq.gz')):

#now combine the results files
#add the filename (which is the seqid) to the first column in all of the res files, then concatenate into a single output csv file
outputfilepeamr = os.path.join(args.seq_dir, 'kmape_amr_output.csv')
outputreadsamr = os.path.join(args.seq_dir, 'kmape_amr_readcounts.csv')
with open(outputfilepeamr, 'w', newline='') as file_output:
    csv_output = csv.writer(file_output)
    #gather the info from all res files into one csv
    for fname in glob.glob(os.path.join(args.seq_dir, '*_amr.res')):
        fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
        seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
        seqname2 = os.path.split(seqname)[1].split('_')[0] #to remove the _SXX_L001 from seqid
        #first, find and replace the #Template column header
        with open(fname, 'r') as file:
            #read in the file
            data = file.read()
            #search and replace the text
            data = data.replace("#Template", "Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass")
        #open the file in write only mode to write the replaced content
        with open(fname, 'w') as file:
            file.write(data)
        with open(fname, newline='') as f_input:
            csv_input = csv.reader(f_input)
            #next(csv_input)
            for row in csv_input:
                print(row)
                row.insert(0,seqname2) #this adds the seqid to the file before concatenating
                csv_output.writerow(row)
df = pandas.read_csv(outputfilepeamr,sep=',|\t',header=None)
df.to_csv(outputfilepeamr, header=headerList,index=False)


outputfilepemetal = os.path.join(args.seq_dir, 'kmape_metal_output.csv')
outputreadsmetal = os.path.join(args.seq_dir, 'kmape_metal_readcounts.csv')
with open(outputfilepemetal, 'w', newline='') as file_output:
    csv_output = csv.writer(file_output)
    for fname in glob.glob(os.path.join(args.seq_dir, '*_metal.res')):
        fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
        seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
        seqname2 = os.path.split(seqname)[1].split('_')[0] #to remove the _SXX_L001 from seqid
        #first, find and replace the #Template column header
        with open(fname, 'r') as file:
            #read in the file
            data = file.read()
            #search and replace the text
            data = data.replace("#Template", "Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass")
        #open the file in write only mode to write the replaced content
        with open(fname, 'w') as file:
            file.write(data)
        with open(fname, newline='') as f_input:
            csv_input = csv.reader(f_input)
            #next(csv_input)
            for row in csv_input:
                print(row)
                row.insert(0,seqname2) #this adds the seqid to the file before concatenating
                #row.insert(0,fname)
                csv_output.writerow(row)
df = pandas.read_csv(outputfilepemetal,sep=',|\t',header=None)
df.to_csv(outputfilepemetal, header=headerList,index=False)

outputfilepebiocide = os.path.join(args.seq_dir, 'kmape_biocide_output.csv')
outputreadsbiocide = os.path.join(args.seq_dir, 'kmape_biocide_readcounts.csv')
with open(outputfilepebiocide, 'w', newline='') as file_output:
    csv_output = csv.writer(file_output)
    for fname in glob.glob(os.path.join(args.seq_dir, '*_biocide.res')):
        fbasename = os.path.basename(fname) #this is to just get the seqid.res name of file
        seqname = os.path.split(fbasename)[1].split('.')[0] #this is to just pull out the seqid
        seqname2 = os.path.split(seqname)[1].split('_')[0] #to remove the _SXX_L001 from seqid
        #first, find and replace the #Template column header
        with open(fname, 'r') as file:
            #read in the file
            data = file.read()
            #search and replace the text
            data = data.replace("#Template", "Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass")
        #open the file in write only mode to write the replaced content
        with open(fname, 'w') as file:
            file.write(data)
        with open(fname, newline='') as f_input:
            csv_input = csv.reader(f_input)
            #next(csv_input)
            for row in csv_input:
                print(row)
                row.insert(0,seqname2) #this adds the seqid to the file before concatenating
                #row.insert(0,fname)
                csv_output.writerow(row)
df = pandas.read_csv(outputfilepebiocide,sep=',|\t',header=None)
df.to_csv(outputfilepebiocide, header=headerList,index=False)


#now make dataframes of all the KMA hits to the raw reads
#lets modify the output kma files so we end up with an excel table
kmapedf2 = pandas.read_csv(os.path.join(args.seq_dir,'kmape_amr_output.csv'), sep=',|\t', engine='python')
#split the newly renamed column
kmapedf2[['Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass']] = kmapedf2['Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass'].str.split('|',expand=True)
#reorder the dataframe
kmapedf3 = kmapedf2[['SeqID','Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value','Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass']]
#drop the unwanted column
kmapedf3.drop('Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass', axis=1, inplace=True)

#metal files
kmapemetdf2 = pandas.read_csv(os.path.join(args.seq_dir,'kmape_metal_output.csv'), sep=',|\t', engine='python')
#split the newly renamed column
kmapemetdf2[['Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass']] = kmapemetdf2['Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass'].str.split('|',expand=True)
#reorder the dataframe
kmapemetdf3 = kmapemetdf2[['SeqID','Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value','Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass']]
#drop the unwanted column
kmapemetdf3.drop('Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass', axis=1, inplace=True)

#biocide files
kmapebiocidedf2 = pandas.read_csv(os.path.join(args.seq_dir,'kmape_biocide_output.csv'), sep=',|\t', engine='python')
#split the newly renamed column
kmapebiocidedf2[['Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass']] = kmapebiocidedf2['Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass'].str.split('|',expand=True)
#reorder the dataframe
kmapebiocidedf3 = kmapebiocidedf2[['SeqID','Allele','Gene_family','Description','ResistanceType','ResistanceClass','ResistanceSubClass','Score','Expected','Template_length','Template_Identity','Template_Coverage','Query_Identity','Query_Coverage','Depth','q_value','p_value','Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass']]
#drop the unwanted column
kmapebiocidedf3.drop('Allele|Gene_family|Description|ResistanceType|ResistanceClass|ResistanceSubClass', axis=1, inplace=True)

#write the KMA AMR, metal, biocide results to an excel file
kmaexcelfile = os.path.join(args.assemblies_dir, 'KMA_results_redmine{}.xlsx'.format(args.issueid))
with pandas.ExcelWriter(kmaexcelfile) as writer:
    kmadf3.to_excel(writer,sheet_name='AMR_assembly', index=False)
    kmadfm3.to_excel(writer, sheet_name='MetalR_assembly', index=False)
    kmadfb3.to_excel(writer, sheet_name='BiocideR_assembly', index=False)
    kmapedf3.to_excel(writer,sheet_name='AMR_rawreads', index=False)
    kmapemetdf3.to_excel(writer, sheet_name='MetalR_rawread', index=False)
    kmapebiocidedf3.to_excel(writer, sheet_name='BiocideR_rawreads', index=False)



#TODO: now combine the readcounts (the mapstat files) using the merge script Ashley created (based on metaphlan's merge script)
mergereadcounts = '/mnt/nas2/redmine/applications/OLCRedmineAutomator/automators/merge_kmamapstat_tables.py'
mergeamrcmd = 'python {} *_amr.mapstat -o {}'.format(mergereadcounts, outputreadsamr)
mergemetcmd = 'python {} *_metal.mapstat -o {}'.format(mergereadcounts, outputreadsmetal)
mergebioccmd = 'python {} *_metal.mapstat -o {}'.format(mergereadcounts, outputreadsbiocide)

# Create another shell script to execute within the KMA conda environment
template2 = "#!/bin/bash\n{} && cd {} && {} && {} && {}".format(activate, args.seq_dir, mergeamrcmd, mergemetcmd, mergebioccmd)
kmareads_script = os.path.join(args.work_dir, 'merge_kma_reads.sh')
with open(kmareads_script, 'w+') as file:
    file.write(template2)
# Modify the permissions of the script to allow it to be run on the node
make_executable(kmareads_script)
# Run shell script
os.system(kmareads_script)

#write the KMA AMR, metal, biocide read count results to an excel file
#first need to read in the files as data frames
outputreadsamr2 = pandas.read_csv(os.path.join(args.seq_dir, 'kmape_amr_readcounts.csv'), sep=',|\t', engine='python')
outputreadsmetal2 = pandas.read_csv(os.path.join(args.seq_dir, 'kmape_metal_readcounts.csv'), sep=',|\t', engine='python')
outputreadsbiocide2 = pandas.read_csv(os.path.join(args.seq_dir, 'kmape_biocide_readcounts.csv'), sep=',|\t', engine='python')

kmareadsfile = os.path.join(args.assemblies_dir, 'KMA_readcounts_redmine{}.xlsx'.format(args.issueid))
with pandas.ExcelWriter(kmareadsfile) as writer:
    outputreadsamr2.to_excel(writer,sheet_name='AMR_readcounts', index=False)
    outputreadsmetal2.to_excel(writer, sheet_name='MetalR_readcounts', index=False)
    outputreadsbiocide2.to_excel(writer, sheet_name='BiocideR_readcounts', index=False)



