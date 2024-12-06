import os
import glob
import click
import pickle
import shutil
import subprocess
from biotools import mash
from nastools.nastools import retrieve_nas_files
from externalretrieve import upload_to_ftp


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def pyseer_redmine(redmine_instance, issue, work_dir, description):
    try:
        # Unpickle Redmine objects
        redmine_instance = pickle.load(open(redmine_instance, 'rb'))
        issue = pickle.load(open(issue, 'rb'))
        description = pickle.load(open(description, 'rb'))

        #analyses supported by the automator
        analysistypes = [
            'roary', 'kmer', 'snv'
        ]

        distances = [
            'mash', 'bcgtree', 'custom'
        ]

        models = [
            'linearmixed', 'fixed'
        ]
        #dictionary of phenotype types argument flags to pass to the script

        # Variable to hold supplied arguments
        argument_dict = {
            'analysistype': str(),
            'dimensions': False,
            'model': 'linearmixed',
            'distance': 'mash',
            'bootstraps': 100,
            'min_proteomes': 2,
            'aa_substitution_model': 'AUTO',
            'continuous': False,
            'minsupp': 2,
            'maxsupp': 1000,
            'roaryissue': False,
            'maxdimensions': 2,
        }


        # Parse description to figure out what SEQIDs we need to run on.
        seqids = list()
        analysistype = "" # Just for the output at the end
        for item in description:
            item = item.rstrip()
            if 'analysistype' in item:
                argument_dict['analysistype'] = item.split('=')[1].lower()
                continue
            if 'dimensions' in item:
                argument_dict['dimensions'] = item.split('=')[1]
                continue
            if 'model' in item:
                argument_dict['model'] = item.split('=')[1].lower()
                continue
            if 'distance' in item:
                argument_dict['distance'] = item.split('=')[1].lower()
                continue
            if 'bootstraps' in item:
                argument_dict['bootstraps'] = item.split('=')[1]
                continue
            if 'min_proteomes' in item:
                argument_dict['min_proteomes'] = item.split('=')[1]
                continue
            if 'aa_substitution_model' in item:
                argument_dict['aa_substitution_model'] = item.split('=')[1].upper()
                continue
            if 'continuous' in item:
                argument_dict['continuous'] = True
                continue
            if 'minsupp' in item:
                argument_dict['minsupp'] = item.split('=')[1]
                continue
            if 'maxsupp' in item:
                argument_dict['maxsupp'] = item.split('=')[1]
                continue
            if 'roaryissue' in item:
                argument_dict['roaryissue'] = item.split('=')[1]
                continue
            if 'maxdimensions' in item:
                argument_dict['maxdimensions'] = item.split('=')[1]
                continue
            # Otherwise the item should be a SEQID
            seqids.append(item)

        # Ensure that the analysis type is provided
        if not argument_dict['analysistype']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not identify an analysis type. '
                                                'Please ensure that the first line of the issue contains one'
                                                ' of the following keywords: {ats}'.format(ats=', '.join(analysistypes)),
                                          status_id=4)
            return
        elif argument_dict['analysistype'] not in analysistypes:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied analysis type {at} current not in the supported '
                                                'list of analyses: {ats}'.format(at=argument_dict['analysistype'],
                                                                                 ats=', '.join(analysistypes)),
                                          status_id=4)
            return

        # Ensure that the distance type (phylogeny) is provided
        if not argument_dict['distance']:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not identify distance type. '
                                                'Please ensure that the first line of the issue contains one'
                                                ' of the following keywords: {ats}'.format(ats=', '.join(distances)),
                                          status_id=4)
            return
        elif argument_dict['distance'] not in distances:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied distance type {at} current not in the supported '
                                                'list of analyses: {ats}'.format(at=argument_dict['distance'],
                                                                                 ats=', '.join(distances)),
                                          status_id=4)
            return
        #ensure a proper model type is used. the default is linearmixed (as recommended for GWAS by the pyseer tutorial)
        if argument_dict['distance'] not in distances:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: supplied model type {at} current not in the supported '
                                                'list of models: {ats}'.format(at=argument_dict['model'],
                                                                                 ats=', '.join(models)),
                                          status_id=4)
            return  
        
        # No SEQIDs
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
        else:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='ERROR: Analysis requires an attached tsv file of '
                                                'traits. The automator could not find any attached files. '
                                                'Please create a new issue with the traits.tsv file attached and try '
                                                'again. (It must be named "traits.tsv")',
                                          status_id=4)
            return
        #make sure a traits.tsv file was uploaded with the correct name
        if os.path.isfile(os.path.join(work_dir,'traits.tsv')):
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Attached tsv file of traits found.',
                                          status_id=2)
        if argument_dict['distance'] == 'custom':
            if not os.path.isfile(os.path.join(work_dir,'tree.newick')):
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Custom distance/phylogeny analysis requires an attached tree '
                                                    'file in newick format. This file must be named "tree.newick". '
                                                    'The automator could not find an attached tree.newick file. '
                                                    'Please create a new issue with the tree.newick file attached and try '
                                                    'again.',
                                              status_id=4)

        # Create folder to drop FASTA files
        assemblies_folder = os.path.join(work_dir, 'assemblies')
        os.makedirs(assemblies_folder, exist_ok=True)

        # Create output folder for prokka
        output_folder = os.path.join(work_dir, 'output')
        os.makedirs(output_folder)

        # Create pyseer output folder
        pyseer_folder = os.path.join(work_dir, 'pyseer_output')
        os.makedirs(pyseer_folder)

        # Extract FASTA files.
        retrieve_nas_files(seqids=seqids, outdir=assemblies_folder, filetype='fasta', copyflag=False)
        missing_fastas = verify_fasta_files_present(seqids, assemblies_folder)
        if missing_fastas:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find the following requested SEQIDs on '
                                                'the OLC NAS: {}'.format(missing_fastas))

        #if roary is chosen, pull the gene_presence_absence.Rtab file from the previous redmine issue
        if argument_dict['analysistype'] == 'roary':
            if not argument_dict['roaryissue']:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='ERROR: Analysis type {at} requires a gene_presence_absence '
                                                    'file from a previous redmine analysis. Please create a new issue with '
                                                    'roaryissue=redmine_issue_number '
                                                    'included in the issue.'.format(at=argument_dict['analysistype']),
                                              status_id=4)
            else:
                #download the gene_presence_absence file from listed issue
                roary_dir = os.path.join('/mnt/nas2/redmine/bio_requests/', '{roary}/').format(roary=argument_dict['roaryissue'])
                Rtab_file = os.path.join(roary_dir, 'gene_presence_absence.Rtab')
                if os.path.isfile(Rtab_file):
                    shutil.copyfile(os.path.join(roary_dir,'gene_presence_absence.Rtab'), os.path.join(work_dir,'gene_presence_absence.Rtab'))
                if not os.path.isfile(Rtab_file):
                    redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                                  notes='ERROR: Could not find gene_presence_absence.Rtab for specified redmine issue.'
                                                        ' Redmine biorequest specified was {roary}'.format(roary=argument_dict['roaryissue']))

        #first need to make a distance matrix using a tree
        activatepyseer = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/pyseer'

        #check for attached phylogeny file, if distance=custom selected
        if argument_dict['distance'] == 'custom':
            phy_dist = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/phylogeny_distance.py'
            if argument_dict['model'] == 'linearmixed':
                phy_dist_cmd = 'python {phy_dist} --lmm tree.newick > distances_lmm.tsv'.format(phy_dist=phy_dist)
                scree_plot_cmd = 'scree_plot_pyseer distances_lmm.tsv'
            if argument_dict['model'] == 'fixed':
                phy_dist = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/phylogeny_distance.py'
                phy_dist_cmd = 'python {phy_dist} tree.newick > distances.tsv'.format(phy_dist=phy_dist)
                scree_plot_cmd = 'scree_plot_pyseer distances.tsv'
                # Create another shell script to execute within the pyseer conda environment
                templatephy = "#!/bin/bash\n{} && cd {} && {} && {}".format(activatepyseer, work_dir, phy_dist_cmd, scree_plot_cmd)
                phy_dist_script = os.path.join(work_dir, 'phylogeny_dist.sh')
                with open(phy_dist_script, 'w+') as file:
                    file.write(templatephy)
                make_executable(phy_dist_script)
                # Run shell script
                os.system(phy_dist_script)

        #will use mash sketch or bcgtree...
        if argument_dict['distance'] == 'mash':
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Running mash.',
                                          status_id=2)
            mash_sketch_cmd = 'mash sketch -s 10000 -o mash_sketch assemblies/*.fasta'
            mash_dist_cmd = 'mash dist mash_sketch.msh mash_sketch.msh| square_mash > distances.tsv'
            scree_plot_cmd = 'scree_plot_pyseer distances.tsv'

            # Create another shell script to execute within the pyseer conda environment
            templatemsh = "#!/bin/bash\n{activate} && cd {work_dir} && {mashsketch} && {mashdist}".format(activate=activatepyseer,
                                                                                                          work_dir=work_dir,
                                                                                                          mashsketch=mash_sketch_cmd,
                                                                                                          mashdist=mash_dist_cmd)
            mash_script = os.path.join(work_dir, 'run_mash_sketch.sh')
            with open(mash_script, 'w+') as file:
                file.write(templatemsh)
            make_executable(mash_script)

            # Run shell script
            os.system(mash_script)
            #create a mashtree phylogeny
            cmd = '/home/ubuntu/bin/mashtree --numcpus 24 --outtree {output_newick} {input_fastas}'\
                .format(output_newick=os.path.join(work_dir, 'mash.tree'),
                        input_fastas=os.path.join(work_dir, 'assemblies', '*.fasta'))
            returncode = subprocess.call(cmd, shell=True, env={'PERL5LIB': '$PERL5LIB:/home/ubuntu/lib/perl5'})
            if returncode != 0:
                raise Exception('Tree creation command mashtree for {} had return code {}'.format(issue.id, returncode))
            phy_dist = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/phylogeny_distance.py'
            if argument_dict['model'] == 'linearmixed':
                phy_dist_cmd = 'python {phy_dist} --lmm mash.tree > distances_lmm.tsv'.format(phy_dist=phy_dist)
                scree_plot_cmd = 'scree_plot_pyseer distances_lmm.tsv'
                # Create another shell script to execute within the pyseer conda environment
                templatephy = "#!/bin/bash\n{} && cd {} && {} && {}".format(activatepyseer, work_dir, phy_dist_cmd, scree_plot_cmd)
                phy_dist_script = os.path.join(work_dir, 'phylogeny_dist_lmm.sh')
                with open(phy_dist_script, 'w+') as file:
                    file.write(templatephy)
                make_executable(phy_dist_script)
                # Run shell script
                os.system(phy_dist_script)
#            if argument_dict['model'] == 'fixed':
#                phy_dist = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/phylogeny_distance.py'
#                phy_dist_cmd = 'python {phy_dist} mash.tree > distances.tsv'.format(phy_dist=phy_dist)
#                scree_plot_cmd = 'scree_plot_pyseer distances.tsv'
                # Create another shell script to execute within the pyseer conda environment
#                templatephy = "#!/bin/bash\n{} && cd {} && {} && {}".format(activatepyseer, work_dir, phy_dist_cmd, scree_plot_cmd)

#                phy_dist_script = os.path.join(work_dir, 'phylogeny_dist.sh')
#                with open(phy_dist_script, 'w+') as file:
#                    file.write(templatephy)
#                make_executable(phy_dist_script)

                # Run shell script
#                os.system(phy_dist_script)

        #if argument is for bcgtree (bacterial core genome tree), we need to create the tree first
        if argument_dict['distance'] == 'bcgtree':
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Running bcgTree.',
                                          status_id=2)
            # Create the folder in which .fa files are to be stored
            bcgtree_folder = os.path.join(work_dir, 'bcgtree')
            os.makedirs(bcgtree_folder)
            # These unfortunate hard coded paths appear to be necessary
            activatep = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/roary'
            prokka = '/mnt/nas2/virtual_environments/roary/bin/prokka'

            for assembly in glob.glob(os.path.join(assemblies_folder, '*.fasta')):
                seqid = os.path.split(assembly)[1].split('.')[0]
                # Prepare command
                prokkacmd = '{prokka} --outdir {output_folder} --prefix {seqid} {assembly}'.format(prokka=prokka,
                                                                                                   output_folder=os.path.join(output_folder, seqid),
                                                                                                   seqid=seqid,assembly=assembly)

                # Create another shell script to execute within the conda environment
                template = "#!/bin/bash\n{} && {}".format(activatep, prokkacmd)
                prokka_script = os.path.join(work_dir, 'run_prokka.sh')
                with open(prokka_script, 'w+') as file:
                    file.write(template)
                make_executable(prokka_script)
    
                # Run shell script
                os.system(prokka_script)
                faa_file = '{seqid}.faa'.format(seqid=seqid)
                output_faa = os.path.join(work_dir, 'bcgtree', faa_file)
                # Copy the .faa files to the output folder
                shutil.copyfile(os.path.join(output_folder, seqid, faa_file),
                                output_faa)

            #create config file to pass to bcgtree
            config_file = os.path.join(bcgtree_folder, 'config.txt')
            bootstraps = "--bootstraps={boots}".format(boots=argument_dict['bootstraps'])
            outdir = "--outdir={out}".format(out=bcgtree_folder)
            minproteomes = "--min-proteomes={minp}".format(minp=argument_dict['min_proteomes'])
            aasubmod = '--raxml-aa-substitution-model "{aamod}"'.format(aamod=argument_dict['aa_substitution_model'])
            configtemplate = "--threads=8\n{}\n{}\n{}\n{}".format(bootstraps, outdir, minproteomes, aasubmod)
            # with open(config_file, 'w+') as file:
            with open(config_file, 'w+') as file:
                file.write(configtemplate)
            for proteome in glob.glob(os.path.join(bcgtree_folder, '*.faa')):
                proteomeid = os.path.split(proteome)[1].split('.')[0]
                proteomefa = os.path.split(proteome)[1]
                # add a line to the config file for each proteome/seqid in the output folder            
                with open(config_file, 'a+') as file:
                    file.write("\n--proteome {proteomeid}={proteome} ".format(proteomeid=proteomeid, proteome=proteomefa))
            # Run bcgtree on the .faa files
            activateb = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/bcgTree'
            bcgtree = '/mnt/nas2/virtual_environments/bcgTree/bcgTree/bin/bcgTree.pl'

            #this command moves to the folder containing the .faa files and the config.txt file, then runs bcgTree
            bcgtree_cmd = 'cd {bcgtree_dir} && {bcgtree} @{config}'.format(bcgtree_dir=bcgtree_folder, 
                                                                           bcgtree=bcgtree,
                                                                           config=config_file)
            # Create another shell script to execute within the bcgtree conda environment
            templateb = "#!/bin/bash\n{} && {}".format(activateb, bcgtree_cmd)

            bcgtree_script = os.path.join(work_dir, 'run_bcgtree.sh')
            with open(bcgtree_script, 'w+') as file:
                file.write(templateb)
            make_executable(bcgtree_script)

            # Run shell script
            os.system(bcgtree_script)

            #move alignment files to new subfolder
            output_files = os.listdir(bcgtree_folder)
            alignments_folder = os.path.join(bcgtree_folder, 'alignment_files')
            os.makedirs(alignments_folder)        
            for file in output_files:
                if file.endswith(".aln"):
                    shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
            for file in output_files:
                if file.endswith(".aln-gb"):
                    shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
            for file in output_files:
                if file.endswith(".aln-gb.comp"):
                    shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
            for file in output_files:
                if file.endswith(".aln-gb.fa"):
                    shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
            for file in output_files:
                if file.endswith(".aln-gb.htm"):
                    shutil.move(os.path.join(bcgtree_folder,file), os.path.join(alignments_folder,file))
            #move id files to a new subfolder
            ids_folder = os.path.join(bcgtree_folder, 'ids_files')
            os.makedirs(ids_folder)        
            for file in output_files:
                if file.endswith(".ids"):
                    shutil.move(os.path.join(bcgtree_folder,file), os.path.join(ids_folder,file))
            #move protein fasta files to a new subfolder
            fa_files = os.listdir(bcgtree_folder)
            fastas_folder = os.path.join(bcgtree_folder, 'fasta_files')
            os.makedirs(fastas_folder)        
            for file in fa_files:
                if file.endswith(".fa"):
                    shutil.move(os.path.join(bcgtree_folder,file), os.path.join(fastas_folder,file))
            #move the tree file to the working directory
            shutil.move(os.path.join(bcgtree_folder,"RAxML_bestTree.final"), os.path.join(work_dir,"RAxML_bestTree.final")) 

            #now need to get distances from the bcgtree phylogeny
            if argument_dict['model'] == 'linearmixed':
                phy_dist = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/phylogeny_distance.py'
                phy_dist_cmd = 'python {phy_dist} --lmm RAxML_bestTree.final > distances_lmm.tsv'.format(phy_dist=phy_dist)
                scree_plot_cmd = 'scree_plot_pyseer distances_lmm.tsv'
                # Create another shell script to execute within the pyseer conda environment
                templatephy = "#!/bin/bash\n{} && cd {} && {} && {}".format(activatepyseer, work_dir, phy_dist_cmd, scree_plot_cmd)

                phy_dist_script = os.path.join(work_dir, 'phylogeny_dist_lmm.sh')
                with open(phy_dist_script, 'w+') as file:
                    file.write(templatephy)
                make_executable(phy_dist_script)

                # Run shell script
                os.system(phy_dist_script)
            if argument_dict['model'] == 'fixed':
                phy_dist = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/phylogeny_distance.py'
                phy_dist_cmd = 'python {phy_dist} RAxML_bestTree.final > distances.tsv'.format(phy_dist=phy_dist)
                scree_plot_cmd = 'scree_plot_pyseer distances.tsv'
                # Create another shell script to execute within the pyseer conda environment
                templatephy = "#!/bin/bash\n{} && cd {} && {} && {}".format(activatepyseer, work_dir, phy_dist_cmd, scree_plot_cmd)

                phy_dist_script = os.path.join(work_dir, 'phylogeny_dist.sh')
                with open(phy_dist_script, 'w+') as file:
                    file.write(templatephy)
                make_executable(phy_dist_script)

                # Run shell script
                os.system(phy_dist_script)


        #run analysis depending on analysistype
        #below is for the kmer analysistype
        if argument_dict['analysistype'] == 'kmer':
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='k-mer analysis chosen. '
                                                'Counting k-mers using fsm-lite.',
                                          status_id=2)
            nodes = str('# SBATCH -N 1')
            processors = str('#SBATCH --ntasks=25')
            memory = str('#SBATCH --mem=100000')
            #run fsm-lite to count kmers
            #create fsm_file_list file to pass to fsm-lite
            fsm_file = os.path.join(work_dir, 'fsm_file_list.txt')
            # with open(fsm_file, 'w+') as file:
            for fasta in glob.glob(os.path.join(assemblies_folder, '*.fasta')):
                fastaid = os.path.split(fasta)[1].split('.')[0]
                # add a line to the fsm file for each seqid in the assemblies folder            
                with open(fsm_file, 'a+') as file:
                    file.write("{fastaid}\t{fastaid}.fasta\n".format(fastaid=fastaid))
            fsm_cmd = 'fsm-lite -l {fsm_file} -s {minsp} -S 610 -v -t fsm_kmers | gzip -c - > {work}/fsm_kmers.txt.gz'.format(fsm_file=fsm_file,minsp=argument_dict['minsupp'],work=work_dir)
            templatefsm = "#!/bin/bash\n{}\n{}\n{}\n{} && cd {} && {}".format(nodes,processors,\
                                                                              memory,activatepyseer,\
                                                                              assemblies_folder,\
                                                                              fsm_cmd)
            fsm_script = os.path.join(work_dir, 'fsm-lite_script.sh')
            with open(fsm_script, 'w+') as file:
                file.write(templatefsm)
            make_executable(fsm_script)
            # Run shell script
            os.system(fsm_script)

            if argument_dict['model'] == 'linearmixed':
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='Running pyseer using linear mixed model. (This step may take a while)',
                                              status_id=2)
                #run pyseer with --lmm (I believe this is the longest step)
                pyseer_lmm_cmd = 'pyseer --lmm --phenotypes traits.tsv --kmers fsm_kmers.txt.gz --similarity distances_lmm.tsv --output-patterns kmer_patterns.txt --cpu 20'
                #if continuous data is givin, add the continuous flag to the command
                pyseer_lmm_cmd += ' --continuous' if argument_dict['continuous'] else ''
                #if going to allow additional flags, need to put the output as its own function
                pyseer_out = ' > pyseer_kmers.txt'
                # Create another shell script to execute within the pyseer conda environment
                templatepylmm = "#!/bin/bash\n{}\n{}\n{}\n{} && cd {} && {}{}".format(nodes,processors,\
                                                                                      memory,activatepyseer,work_dir,\
                                                                                      pyseer_lmm_cmd,pyseer_out)
                pylmm_script = os.path.join(work_dir, 'pyseer_lmm_script.sh')
                with open(pylmm_script, 'w+') as file:
                    file.write(templatepylmm)
                make_executable(pylmm_script)
                # Run shell script
                os.system(pylmm_script)
            if argument_dict['model'] == 'fixed':
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='Running pyseer using fixed effects model. (This step may take a while)',
                                              status_id=2)                
                #run pyseer with the other association model (fixed effects (SEER) generalized linear model) available with pyseer        
                pyseer_cmd = 'pyseer --phenotypes traits.tsv --kmers fsm_kmers.txt.gz --distances distances.tsv --output-patterns kmer_patterns.txt --cpu 20'
                #if continuous data is givin, add the continuous flag to the command
                pyseer_cmd += ' --continuous' if argument_dict['continuous'] else ''
                #if going to allow additional flags, need to put the output as its own function
                pyseer_out = ' > pyseer_kmers.txt'
                # Create another shell script to execute within the pyseer conda environment
                templateser = "#!/bin/bash\n{}\n{}\n{}\n{} && cd {} && {}{}".format(nodes,processors,memory,\
                                                                                    activatepyseer,\
                                                                                    work_dir,\
                                                                                    pyseer_cmd,pyseer_out)

                seer_script = os.path.join(work_dir, 'seer_script.sh')
                with open(seer_script, 'w+') as file:
                    file.write(templateser)
                make_executable(seer_script)

                # Run shell script
                os.system(seer_script)
            
            #now need find out what the significance threshold is using the number of unique kmer patterns
            sig_kmers = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/count_patterns.py'
            sig_kmers_cmd = 'python {count} kmer_patterns.txt > significant_kmers_threshold.txt'.format(count=sig_kmers)
            # Create another shell script to execute within the pyseer conda environment
            templatesig = "#!/bin/bash\n{} && cd {} && {}".format(activatepyseer,work_dir,sig_kmers_cmd)

            sig_kmers_script = os.path.join(work_dir, 'sig_kmers_script.sh')
            with open(sig_kmers_script, 'w+') as file:
                file.write(templatesig)
            make_executable(sig_kmers_script)
            # Run shell script
            os.system(sig_kmers_script)

            #now read the line from the significant kmers file to find the p-value (line 2)
            lines_to_read = [1]
            s_file = open(os.path.join(work_dir,'significant_kmers_threshold.txt'))
            for index, line in enumerate(s_file):
                if index in lines_to_read:
                    pvaluethreshold = line.split(':\t')[1]
            #output the qq-plot
            qqplot = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/qq_plot.py'
            qq_cmd = 'python {qqplot} pyseer_kmers.txt'.format(qqplot=qqplot)
            # Create another shell script to execute within the pyseer conda environment
            templateqq = "#!/bin/bash\n{} && cd {} && {}".format(activatepyseer,work_dir,qq_cmd)

            qqplot_script = os.path.join(work_dir, 'qqplot_script.sh')
            with open(qqplot_script, 'w+') as file:
                file.write(templateqq)
            make_executable(qqplot_script)
            # Run shell script
            os.system(qqplot_script)

            #now filter for the significant kmers using the threshold above
            print_text = str('{print $0}')
            filter_kmers = "cat <(head -1 pyseer_kmers.txt) <(awk '$4<{pvalue} {print_t}' pyseer_kmers.txt) "\
                           "> significant_kmers.txt".format(pvalue=pvaluethreshold, print_t=print_text)
            # Create another shell script to execute within the pyseer conda environment
            templatesigkmers = "#!/bin/bash\n{} && cd {} && {}".format(activatepyseer,work_dir,filter_kmers)

            sigkmers_script = os.path.join(work_dir, 'sigkmers_script.sh')
            with open(sigkmers_script, 'w+') as file:
                file.write(templatesigkmers)
            make_executable(sigkmers_script)
            # Run shell script
            os.system(sigkmers_script)
            
            #now need to annotate kmers
            #if a references.txt file was uploaded
            if os.path.isfile(os.path.join(work_dir,'references.txt')):
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='Annotating significant k-mers. '
                                                    'An attached references.txt file was found. '
                                                    'Attached file will be used for k-mer annotations',
                                              status_id=2)
            else:
                redmine_instance.issue.update(resource_id=issue.id,
                                              notes='Annotating significant k-mers. '
                                                    'No attached references.txt file, '
                                                    'a reference file will be created for k-mer annotations.',
                                              status_id=2)
                #create references file to pass to pyseer
                refs_file = os.path.join(work_dir, 'references.txt')
                for fasta in glob.glob(os.path.join(assemblies_folder, '*.fasta')):
                    fastaid = os.path.split(fasta)[1].split('.')[0]
                    # add a line to the fsm file for each seqid in the assemblies folder            
                    with open(refs_file, 'a+') as file:
                        file.write("{fastaid}.fasta\t{fastaid}.gff\tref\n".format(fastaid=fastaid))
            
            #now annotate the kmers
            #first have to move the gff files to the assemblies folder
            if argument_dict['distance'] == 'bcgtree':
                gff_file = '{seqid}.gff'.format(seqid=seqid)
                output_gff = os.path.join(assemblies_folder, gff_file)
            #I think I need to add something to make sure the gff file moves to the assembly folder
                shutil.copyfile(os.path.join(output_folder, seqid, gff_file),
                                output_gff) #changed to output_gff from assemblies_folder
            else:
                # run prokka if bcgtree was not used, as this means prokka was not run yet
                #These unfortunate hard coded paths appear to be necessary
                activatep = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/roary'
                prokka = '/mnt/nas2/virtual_environments/roary/bin/prokka'
                for assembly in glob.glob(os.path.join(assemblies_folder, '*.fasta')):
                    seqid = os.path.split(assembly)[1].split('.')[0]
                    # Prepare command
                    prokkacmd = '{prokka} --outdir {output_folder} --prefix {seqid} {assembly}'.format(prokka=prokka,
                                                                                                       output_folder=os.path.join(output_folder, seqid),
                                                                                                       seqid=seqid,assembly=assembly)

                    # Create another shell script to execute within the conda environment
                    templatep = "#!/bin/bash\n{} && {}".format(activatep, prokkacmd)
                    prokka_script = os.path.join(work_dir, 'run_prokka.sh')
                    with open(prokka_script, 'w+') as file:
                        file.write(templatep)
                    make_executable(prokka_script)
    
                    # Run shell script
                    os.system(prokka_script)
                    #move the files
                    gff_file = '{seqid}.gff'.format(seqid=seqid)
                    output_gff = os.path.join(assemblies_folder, gff_file)
                    shutil.copyfile(os.path.join(output_folder, seqid, gff_file),
                                    output_gff)

            #make a script to move to the assembly folder, then annotate the kmers
            annotate_cmd = 'annotate_hits_pyseer {work}/significant_kmers.txt {work}/references.txt {work}/annotated_kmers.txt'.format(work=work_dir)
            templateannotate = "#!/bin/bash\n{} && cd {} && {}".format(activatepyseer,assemblies_folder,annotate_cmd)

            annotate_script = os.path.join(work_dir, 'annotate_script.sh')
            with open(annotate_script, 'w+') as file:
                file.write(templateannotate)
            make_executable(annotate_script)
            # Run shell script
            os.system(annotate_script)

            #now going to summarise the annotations. the gene_hits output file can be used to create graphs in R
            summarise = '/mnt/nas2/virtual_environments/pyseer/pyseer/scripts/summarise_annotations.py'
            summarise_cmd = 'python {summarise} annotated_kmers.txt > gene_hits.txt'.format(summarise=summarise)
            templatesummarise = "#!/bin/bash\n{} && cd {} && {}".format(activatepyseer,work_dir,summarise_cmd)
            summarise_script = os.path.join(work_dir, 'summarise_script.sh')
            with open(summarise_script, 'w+') as file:
                file.write(templatesummarise)
            make_executable(summarise_script)
            # Run shell script
            os.system(summarise_script)

        #run if the analysis type is roary
        if argument_dict['analysistype'] == 'roary':
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='Roary analysis chosen. '
                                                'Running pyseer using gene_presence_absence file.',
                                          status_id=2)
            nodes = str('# SBATCH -N 1')
            processors = str('#SBATCH --ntasks=25')
            memory = str('#SBATCH --mem=100000')
            if argument_dict['model'] == 'fixed':
                pyseer_cmd = 'pyseer --phenotypes traits.tsv --pres gene_presence_absence.Rtab --distances distances.tsv --save-m {dist}_mds --max-dimensions {md}'.format(dist=argument_dict['distance'],
                                                                                                                                                                           md=argument_dict['maxdimensions'])
                #if continuous data is givin, add the continuous flag to the command
                pyseer_cmd += ' --continuous' if argument_dict['continuous'] else ''
                #if going to allow additional flags, need to put the output as its own function
                pyseer_out = ' > pyseer_COGs.txt'
                # Create another shell script to execute within the pyseer conda environment
                templateser = "#!/bin/bash\n{}\n{}\n{}\n{} && cd {} && {}{}".format(nodes,processors,memory,\
                                                                                    activatepyseer,\
                                                                                    work_dir,\
                                                                                    pyseer_cmd,pyseer_out)

                seer_script = os.path.join(work_dir, 'seer_script.sh')
                with open(seer_script, 'w+') as file:
                    file.write(templateser)
                make_executable(seer_script)

                # Run shell script
                os.system(seer_script)
            if argument_dict['model'] == 'linearmixed':
                pyseer_cmd = 'pyseer --lmm --phenotypes traits.tsv --pres gene_presence_absence.Rtab --similarity distances_lmm.tsv --output-patterns roary_patterns.txt '\
                             '--cpu 10'
                #if continuous data is givin, add the continuous flag to the command
                pyseer_cmd += ' --continuous' if argument_dict['continuous'] else ''
                #if going to allow additional flags, need to put the output as its own function
                pyseer_out = ' > pyseer_lmm_COGs.txt'
                # Create another shell script to execute within the pyseer conda environment
                templateser = "#!/bin/bash\n{}\n{}\n{}\n{} && cd {} && {}{}".format(nodes,processors,memory,\
                                                                                    activatepyseer,\
                                                                                    work_dir,\
                                                                                    pyseer_cmd,pyseer_out)

                seer_script = os.path.join(work_dir, 'seer_script.sh')
                with open(seer_script, 'w+') as file:
                    file.write(templateser)
                make_executable(seer_script)

                # Run shell script
                os.system(seer_script)

        # Zip prokka output
        output_filename = 'prokka_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=output_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'

        sas_url = upload_to_ftp(local_file=zip_filepath)
        # Prepare upload
        if sas_url:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='Prokka process complete!\n\n'
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
        # Remove the zip file
#        os.remove(zip_filepath)
        #remove the output folder
        shutil.rmtree(output_folder)

        # if bcgtree used, zip the bcgtree folder
        if argument_dict['distance'] == 'bcgtree':
            output_filename = 'bcgtree_output_{}'.format(issue.id)
            zip_filepath = zip_folder(results_path=bcgtree_folder,
                                      output_dir=work_dir,
                                      output_filename=output_filename)
            zip_filepath += '.zip'
            # Prepare upload
            # This file can get too big to upload to Redmine, so we should put it on the FTP.
            sas_url = upload_to_ftp(local_file=zip_filepath)
            # Prepare upload
            if sas_url:
            # Wrap up issue
                redmine_instance.issue.update(
                    resource_id=issue.id,
                    status_id=4,
                    notes='Analysis with bcgTree complete!\n\n'
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

        #move output files to pyseer folder
        #this probably should have been done during the pyseer process, but I don't want to go back and edit it for fear of ruining this working script
        output_files = os.listdir(work_dir)       
        for file in output_files:
            if file.endswith(".tsv"):
                shutil.move(os.path.join(work_dir,file), os.path.join(pyseer_folder,file))
        for file in output_files:
            if file.endswith(".msh"):
                shutil.move(os.path.join(work_dir,file), os.path.join(pyseer_folder,file))
        for file in output_files:
            if file.endswith(".txt"):
                shutil.move(os.path.join(work_dir,file), os.path.join(pyseer_folder,file))
        for file in output_files:
            if file.endswith(".png"):
                shutil.move(os.path.join(work_dir,file), os.path.join(pyseer_folder,file))
        for file in output_files:
            if file.endswith(".pkl"):
                shutil.move(os.path.join(work_dir,file), os.path.join(pyseer_folder,file))
            
        #move the tree file to the pyseer directory
        if os.path.isfile(os.path.join(work_dir,'tree.newick')):
            shutil.move(os.path.join(work_dir,'tree.newick'), os.path.join(pyseer_folder,'tree.newick'))
        if os.path.isfile(os.path.join(work_dir,'RAxML_bestTree.final')):
            shutil.move(os.path.join(work_dir,'RAxML_bestTree.final'), os.path.join(pyseer_folder,'RAxML_bestTree.final'))
        if os.path.isfile(os.path.join(work_dir,'mash.tree')):
            shutil.move(os.path.join(work_dir,'mash.tree'), os.path.join(pyseer_folder,'mash.tree'))
        #move other files to the pyseer directory
        if os.path.isfile(os.path.join(work_dir,'fsm_kmers.txt.gz')):
            shutil.move(os.path.join(work_dir,'fsm_kmers.txt.gz'), os.path.join(pyseer_folder,'fsm_kmers.txt.gz'))


        # Zip pyseer output
        pyseer_filename = 'pyseer_output_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=pyseer_folder,
                                  output_dir=work_dir,
                                  output_filename=pyseer_filename)
        zip_filepath += '.zip'

        sas_url = upload_to_ftp(local_file=zip_filepath)
        # Prepare upload
        if sas_url:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=4,
                notes='Pyseer process complete!\n\n'
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
        # Remove the zip file
#        os.remove(zip_filepath)
        #remove the pyseer folder
#        shutil.rmtree(pyseer_folder)

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


def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path


if __name__ == '__main__':
    pyseer_redmine()
