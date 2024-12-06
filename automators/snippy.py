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
import fnmatch



@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def snippy_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    # Current list of sequence types that KMA can analyse
    querytypes = [
        'fasta', 'fastq'
    ]

    # Variable to hold supplied arguments
    argument_dict = {
        'reference': str(),
        'cleanup_&_tree': False,
        'querytype': 'fastq',
    }

    try:
        # Parse description to figure out what SEQIDs we need to run on.
        reference = list()
        seqids = list()
        for item in description:
            item = item.rstrip()
            if 'reference' in item:
                argument_dict['reference'] = item.split('=')[1]
                reference.append(argument_dict['reference'])
                continue
            if 'cleanup_&_tree' in item:
                argument_dict['cleanup_&_tree'] = True
                continue
            if 'querytype' in item:
                argument_dict['querytype'] = item.split('=')[1]
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)
        #print(reference)
        #print(seqids)
        # Ensure that a reference is provided
        if not argument_dict['reference']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No reference provided. ',
                                          status_id=4)
            return

        # Ensure that SEQIDs were included
        if not seqids:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: No SEQIDs provided!',
                                          status_id=4)
            return

        #create a folder to hold the sequences
        seq_dir = os.path.join(work_dir, 'sequences')
        os.makedirs(seq_dir, exist_ok=True)

        #create a folder to hold the reference
        ref_dir = os.path.join(work_dir, 'reference')
        os.makedirs(ref_dir, exist_ok=True)

        # Extract FASTA file to reference folder.
        if argument_dict['reference'] != 'attached':
            # Extract our reference file to our working directory.
            retrieve_nas_files(seqids=reference, outdir=ref_dir, filetype='fasta', copyflag=False)
            missing_fastas = verify_fasta_files_present(reference, ref_dir)
            if missing_fastas:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested SEQIDs on '
                                                    'the OLC NAS: {}'.format(missing_fastas))

        # If user specified attachment as the reference file, download it to our working directory.
        else:
            # Get the attachment ID, and download if it isn't equal to zero (meaning no attachment, so boot user with
            # appropriate error message)
            attachment = redmine_instance.issue.get(issue.id, include='attachments')
            attachment_id = 0
            for item in attachment.attachments:
                attachment_id = item.id

            # Download if we found an attachment, and use as our reference. Otherwise, exit and tell user to try again
            if attachment_id != 0:
                attachment = redmine_instance.attachment.get(attachment_id)
                #attachment.download(savepath=ref_dir, filename='reference.fasta')
                attachment.download(savepath=ref_dir)
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: You specified that the reference would be in attached file,'
                                                    ' but no attached file was found. Please create a new issue and '
                                                    'try again.',
                                              status_id=4)
                return
        #now create an output folder so we can delete all of the extra files we didn't need
        out_dir = os.path.join(work_dir, 'output')
        os.makedirs(out_dir, exist_ok=True)

        #if they specify fasta, then pull the fasta files into the query folder
        if argument_dict['querytype'] == 'fasta':
            # Extract our query fasta files to our working directory.
            retrieve_nas_files(seqids=seqids, outdir=seq_dir, filetype='fasta', copyflag=False)
            missing_fastas = verify_fasta_files_present(seqids, seq_dir)
            if missing_fastas:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested SEQIDs on '
                                                    'the OLC NAS: {}'.format(missing_fastas))
            #create a list of files in the query folder
            alllist = os.listdir(seq_dir)
            #get the first entry from the list for the snippy-core code
            first = alllist[0]
            firstloc = os.path.join(out_dir, first)

            # create the final line of our run  file so we can add each SEQID and then append to the end of our run file
            snippycore = os.path.join(work_dir, 'snippycore.txt')
            with open(snippycore, 'w+') as file:
                file.write("snippy-core --ref '{se}/ref.fa' ".format(se=firstloc))

            # now need to create the first sections of the snippy-multi run script, these will be done separately for single-end (se) and paired-end (pe) sequences
            refloc = os.path.join(ref_dir, argument_dict['reference'])
            runsnippy = os.path.join(out_dir, 'runmultisnippy.sh')

            for fname in glob.glob(os.path.join(seq_dir, '*.fasta')):
                fbasename = os.path.basename(fname)
                seqname = os.path.split(fbasename)[1].split('.')[0]
                outdir = os.path.join(out_dir, seqname)
                with open(runsnippy, 'a+') as file:
                    file.write("snippy --outdir '{ID}' --ctgs '{fnm}' --ref {ref}.fasta --cpus 20\n".format(ID=outdir,
                                                                                                            fnm=fname,
                                                                                                            ref=refloc))
                with open(snippycore, 'a+') as file:
                    file.write('{ID} '.format(ID=seqname))

            # now add the snippycore command to the runfile
            with open(runsnippy, 'a+') as file:
                snipend = open(snippycore, 'r')
                file.write(snipend.read())

        if argument_dict['querytype'] == 'fastq':
            # Extract the query FASTQ files to sequences folder.
            retrieve_nas_files(seqids=seqids, outdir=seq_dir, filetype='fastq', copyflag=False)
            missing_fastqs = check_fastqs_present(seqids, seq_dir)
            if len(missing_fastqs) > 0:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='WARNING: Could not find the following requested FASTQ SEQIDs on'
                                                    ' the OLC NAS: {}'.format(missing_fastqs))

            #instead of doing the below, will make the .sh files using our own script/coding. This allows us to include MIN single end sequences
            #create the tab-delimited file for snippy-multi
            #add the filename (which is the seqid) to the first column in all of the files, then concatenate into a single output csv file
            #tabfile = os.path.join(out_dir, 'multi_input.tsv')
            #with open(tabfile, 'w') as file:
            #    tsv_file = csv.writer(file, delimiter='\t')
            #    for fname in glob.glob(os.path.join(seq_dir, '*_R1_001.fastq.gz')):
            #        fbasename = os.path.basename(fname) #this is to just get the seqid.fastq.gz name of file
            #        seqname = os.path.split(fbasename)[1].split('_R')[0] #this is to just pull out the seqid
            #        forwrd = "{seq}_R1_001.fastq.gz".format(seq=seqname)
            #        forwardname = os.path.join(seq_dir, forwrd)
            #        reverse = "{seq}_R2_001.fastq.gz".format(seq=seqname)
            #        reversename = os.path.join(seq_dir, reverse)
            #        outdir = os.path.join(out_dir, seqname)
            #        with open(tabfile, 'a+') as file:
            #            file.write("{ID}\t{R1}\t{R2}\n".format(ID=outdir,R1=forwardname,R2=reversename))

            #now going to run snippy
            # These unfortunate hard coded paths appear to be necessary
            #activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/snippy'

            #cmd = 'snippy-multi multi_input.tsv --ref {seq}.fasta --cpus 20 > runmultisnippy.sh'.format(seq=os.path.join(ref_dir,argument_dict['reference']))

            # Create another shell script to execute within the conda environment
            #template = "#!/bin/bash\ncd {} && {} && {}".format(out_dir,activate, cmd)
            #snippymulti_script = os.path.join(work_dir, 'create_snippy_multi_tabfile_cmd.sh')
            #with open(snippymulti_script, 'w+') as file:
            #    file.write(template)
            #make_executable(snippymulti_script)

            # Run shell script
            #os.system(snippymulti_script)

            #will make lists of the raw sequence files... separate out the minion (single-end) and illumina (paired-end) into separate lists
            raw_files = os.listdir(seq_dir)
            alllist = list()
            selist = list()
            pelist = list()
            for raw in raw_files:
                if fnmatch.fnmatch(raw, '*_R1_*'):
                    seqname = os.path.split(raw)[1].split('_R')[0]
                    pelist.append(seqname)
                    alllist.append(seqname)
                elif fnmatch.fnmatch(raw, '*MIN*'):
                    seqname = os.path.split(raw)[1].split('.')[0]
                    selist.append(raw)
                    alllist.append(seqname)

            #get the first entry from the list for the snippy-core code
            first = alllist[0]
            firstloc = os.path.join(out_dir, first)

            #create the final line of our run  file so we can add each SEQID and then append to the end of our run file
            snippycore = os.path.join(work_dir, 'snippycore.txt')
            with open(snippycore, 'w+') as file:
                file.write("snippy-core --ref '{se}/ref.fa' ".format(se=firstloc))

            #now need to create the first sections of the snippy-multi run script, these will be done separately for single-end (se) and paired-end (pe) sequences
            refloc = os.path.join(ref_dir, argument_dict['reference'])
            runsnippy = os.path.join(out_dir, 'runmultisnippy.sh')
            if len(glob.glob(os.path.join(seq_dir, '*_R1*'))) > 0:
                for fname in glob.glob(os.path.join(seq_dir, '*_R1*.fastq.gz')):
                    fbasename = os.path.basename(fname)
                    seqname = os.path.split(fbasename)[1].split('_R')[0] #this is to get the seqid
                    forwrd = "{seq}_R1_001.fastq.gz".format(seq=seqname)
                    forwardname = os.path.join(seq_dir, forwrd)
                    reverse = "{seq}_R2_001.fastq.gz".format(seq=seqname)
                    reversename = os.path.join(seq_dir, reverse)
                    outdir = os.path.join(out_dir, seqname)
                    with open(runsnippy, 'a+') as file:
                        file.write("snippy --outdir '{ID}' --R1 '{R1}' --R2 '{R2}' --ref {ref}.fasta --cpus 20\n".format(ID=outdir,R1=forwardname,R2=reversename,ref=refloc))
                    with open(snippycore, 'a+') as file:
                        file.write('{ID} '.format(ID=seqname))


            #add in the se data if there are minion sequences
            if len(glob.glob(os.path.join(seq_dir, '*MIN*'))) > 0:
                for fname in glob.glob(os.path.join(seq_dir, '*MIN*.fastq.gz')):
                    fbasename = os.path.basename(fname)
                    seqname = os.path.split(fbasename)[1].split('.')[0]
                    outdir = os.path.join(out_dir, seqname)
                    with open(runsnippy, 'a+') as file:
                        file.write("snippy --outdir '{ID}' --se '{fnm}' --ref {ref}.fasta --cpus 20\n".format(ID=outdir,fnm=fname,ref=refloc))
                    with open(snippycore, 'a+') as file:
                        file.write('{ID} '.format(ID=seqname))

            #now add the snippycore command to the runfile
            with open(runsnippy, 'a+') as file:
                snipend = open(snippycore, 'r')
                file.write(snipend.read())

        #now have to use the script we just created above
        make_executable(runsnippy)

        #need to create a script that runs the script created above?
        activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/snippy'
        template2 = "#!/bin/bash\ncd {} && {} && {}".format(out_dir,activate, runsnippy)
        runsnippycluster = os.path.join(work_dir, 'runsnippy_on_cluster.sh')
        with open(runsnippycluster, 'w+') as file:
            file.write(template2)
        make_executable(runsnippycluster)
        #run the script
        os.system(runsnippycluster)

        #make a tree of the core alignment (not cleaned using gubbins)
        activategub = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/gubbins'
        iqtreecmd = 'iqtree -m GTR -s core.full.aln'
        templatetre = "#!/bin/bash\ncd {} && {} && {}".format(out_dir,activategub,iqtreecmd)
        coretreescript = os.path.join(work_dir, 'iqtree_full_alignment.sh')
        with open(coretreescript, 'w+') as file:
            file.write(templatetre)
        make_executable(coretreescript)
        # Run shell script
        os.system(coretreescript)

        #if the cleanup flag is included, then cleanup the alignment file and create a tree
        if argument_dict['cleanup_&_tree']:
            #cleanup the alignment file
            clncmd = 'snippy-clean_full_aln core.full.aln > clean.full.aln'
            templatecln = "#!/bin/bash\ncd {} && {} && {}".format(out_dir,activate, clncmd)
            snippyclean_script = os.path.join(work_dir, 'snippy_clean_alignment.sh')
            with open(snippyclean_script, 'w+') as file:
                file.write(templatecln)
            make_executable(snippyclean_script)
            # Run shell script
            os.system(snippyclean_script)

            #now run gubbins on the cleaned alignment
            activategub = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/gubbins'
            #gubbinscmd = 'run_gubbins.py -p gubbins.fasttree -t fasttree --model GTR clean.full.aln'
            gubbinscmd2 = 'run_gubbins.py -p gubbins.iqtree -t iqtree --model GTR clean.full.aln'
            templategub = "#!/bin/bash\ncd {} && {} && {}".format(out_dir,activategub,gubbinscmd2)
            gubbins_script = os.path.join(work_dir, 'run_gubbins.sh')
            with open(gubbins_script, 'w+') as file:
                file.write(templategub)
            make_executable(gubbins_script)
            # Run shell script
            os.system(gubbins_script)

            #now run snp-analysis
            snpcmd = 'snp-sites -c gubbins.iqtree.filtered_polymorphic_sites.fasta > clean.core.aln'
            templatesnps = "#!/bin/bash\ncd {} && {} && {}".format(out_dir,activate,snpcmd)
            snp_script = os.path.join(work_dir, 'run_snp-sites.sh')
            with open(snp_script, 'w+') as file:
                file.write(templatesnps)
            make_executable(snp_script)
            # Run shell script
            os.system(snp_script)

            #build a tree using iqtree
            treecmd = 'iqtree -m GTR -s clean.core.aln'
            templatetree = "#!/bin/bash\ncd {} && {} && {}".format(out_dir,activategub,treecmd)
            tree_script = os.path.join(work_dir, 'run_iqtree.sh')
            with open(tree_script, 'w+') as file:
                file.write(templatetree)
            make_executable(tree_script)
            # Run shell script
            os.system(tree_script)

            #make an output folder and move the new files
            cleantree = os.path.join(out_dir, 'cleanup_and_tree')
            os.makedirs(cleantree, exist_ok=True)
            output_files = os.listdir(out_dir)
            for file in output_files:
                if file.startswith("gubbins"):
                    shutil.move(os.path.join(out_dir,file), os.path.join(cleantree,file))
            for file in output_files:
                if file.startswith("clean."):
                    shutil.move(os.path.join(out_dir,file), os.path.join(cleantree,file))
            #move the treefile back to the base of the output directory
            treefiles = os.listdir(cleantree)
            for file in treefiles:
                if file == 'clean.core.aln.treefile':
                    shutil.move(os.path.join(cleantree,file),os.path.join(out_dir, '{}.snippy_final.clean.core.snp.aln.tree'.format(issue.id)))
            for file in treefiles:
                if file.startswith('gubbins.iqtree.node_labelled.final_tree'):
                    shutil.move(os.path.join(cleantree,file),os.path.join(out_dir, '{}.snippy.gubbins.iqtree.tree'.format(issue.id)))

        #rename the core.tab and core.txt files to include the redmine id
        output_files = os.listdir(out_dir)
        for file in output_files:
            if file == 'core.txt':
                shutil.move(os.path.join(out_dir,file),os.path.join(out_dir, '{}.snippy_core.txt'.format(issue.id)))
        for file in output_files:
            if file == 'core.tab':
                shutil.move(os.path.join(out_dir,file),os.path.join(out_dir, '{}.snippy_core.tab'.format(issue.id)))
        for file in output_files:
            if file == 'core.full.aln.treefile':
                shutil.move(os.path.join(out_dir,file),os.path.join(out_dir, '{}.snippy_core_alignment.iqtree.treefile'.format(issue.id)))

        # Zip output
        snippy_filename = 'snippy_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=out_dir,
                                  output_dir=work_dir,
                                  output_filename=snippy_filename)
        zip_filepath += '.zip'
        
                #upload some files to the redmine request
        output_list = list()

        if argument_dict['cleanup_&_tree']:
            #treefiles
            output_dict = dict()
            gubbinstree_name = '{}.snippy.gubbins.iqtree.tree'.format(issue.id)
            output_dict['path'] = os.path.join(work_dir, 'output', gubbinstree_name)
            output_dict['filename'] = gubbinstree_name
            output_list.append(output_dict)

            output_dict = dict()
            coresnptree_name = '{}.snippy_final.clean.core.snp.aln.tree'.format(issue.id)
            output_dict['path'] = os.path.join(work_dir, 'output', coresnptree_name)
            output_dict['filename'] = coresnptree_name
            output_list.append(output_dict)

        #core.txt summary file, and the .tab file
        output_dict = dict()
        coretxt = '{}.snippy_core.txt'.format(issue.id)
        output_dict['path'] = os.path.join(work_dir, 'output', coretxt)
        output_dict['filename'] = coretxt
        output_list.append(output_dict)

        output_dict = dict()
        coretab = '{}.snippy_core.tab'.format(issue.id)
        output_dict['path'] = os.path.join(work_dir, 'output', coretab)
        output_dict['filename'] = coretab
        output_list.append(output_dict)

        output_dict = dict()
        alntre = '{}.snippy_core_alignment.iqtree.treefile'.format(issue.id)
        output_dict['path'] = os.path.join(work_dir, 'output', alntre)
        output_dict['filename'] = alntre
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
                status_id=4,
                uploads=output_list,
                notes='Snippy analysis complete!\n\n'
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
        shutil.rmtree(ref_dir)
        # Remove the zip file
        os.remove(zip_filepath)

    except Exception as e:
        sentry_sdk.capture_exception(e)
        #trying to upload the core.txt summary file when/if there is an exception (if snippy fails)
        output_list = list()
        output_dict = dict()
        coretxt = '{}.snippy_core.txt'.format(issue.id)
        output_dict['path'] = os.path.join(work_dir, 'output', coretxt)
        output_dict['filename'] = coretxt
        output_list.append(output_dict)
        redmine_instance.issue.update(resource_id=issue.id,uploads=output_list,#added uploads to the exception
                                      notes='Something went wrong! Please check the core.txt file. '
                                            'If there are sequences with an Alignment of 0, remove from request '
                                            'and re-submit. Otherwise, please contact a bioinformatician '
                                            'to investigate.')




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
    snippy_redmine()
