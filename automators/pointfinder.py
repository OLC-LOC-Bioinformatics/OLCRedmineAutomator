from accessoryFunctions.accessoryFunctions import make_path
from nastools.nastools import retrieve_nas_files
from biotools import mash
import pickle
import shutil
import click
import glob
import os
import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN


def write_report(summary_dict, seqid, genus, key):
    """
    Parse the PointFinder outputs, and write the summary report for the current analysis type
    :param summary_dict: nested dictionary containing data such as header strings, and paths to reports
    :param seqid: name of the strain,
    :param genus: MASH-calculated genus of current isolate
    :param key: current result type. Options are 'prediction', and 'results'
    """
    # Set the header string if the summary report doesn't already exist
    try:
        if not os.path.isfile(summary_dict[genus][key]['summary']):
            header_string = summary_dict[genus][key]['header']
        else:
            header_string = str()
    except KeyError:
        # If the genus is not in the summary_dict, default to Escherichia
        # 'Strain,Genus,Mutation,NucleotideChange,AminoAcidChange,Resistance,PMID,\n'
        if not os.path.isfile(summary_dict['Escherichia'][key]['summary']):
            header_string = summary_dict['Escherichia'][key]['header']
        else:
            header_string = str()
    summary_string = str()
    try:
        # Read in the predictions
        with open(summary_dict[genus][key]['output'], 'r') as outputs:
            # Skip the header
            next(outputs)
            for line in outputs:
                # Skip empty lines
                if line != '\n':
                    # When processing the results outputs, add the seqid to the summary string
                    if key == 'results':
                        summary_string += '{seq},{genus},'.format(seq=seqid,
                                                                  genus=genus)
                    # Clean up the string before adding it to the summary string - replace commas
                    # with semi-colons, and replace tabs with commas
                    summary_string += line.replace(',', ';').replace('\t', ',')
        # Ensure that there were results to report
        if summary_string:
            if not summary_string.endswith('\n'):
                summary_string += '\n'
        else:
            summary_string += '{seq},{genus}\n'.format(seq=seqid,
                                                       genus=genus)
        # Write the summaries to the summary file
        with open(summary_dict[genus][key]['summary'], 'a+') as summary:
            # Write the header if necessary
            if header_string:
                summary.write(header_string)
            summary.write(summary_string)
    # If the genus isn't one that is covered by the PointFinder database, still include the strain information
    except KeyError:
        summary_string += '{seq},{genus}\n'.format(seq=seqid,
                                                   genus=genus)
        # Write the summaries to the summary file
        with open(summary_dict['Escherichia'][key]['summary'], 'a+') as summary:
            # Write the header if necessary
            if header_string:
                summary.write(header_string)
            summary.write(summary_string)


def write_table_report(summary_dict, seqid, genus):
    """
    Parse the PointFinder table output, and write a summary report
    :param summary_dict: nested dictionary containing data such as header strings, and paths to reports
    :param seqid: name of the strain,
    :param genus: MASH-calculated genus of current isolate
    """
    # Set the header string if the summary report doesn't already exist
    if not os.path.isfile(summary_dict[genus]['table']['summary']):
            summary_string = summary_dict[genus]['table']['header']
    else:
        summary_string = str()
    summary_string += '{seq},'.format(seq=seqid)
    # Read in the predictions
    with open(summary_dict[genus]['table']['output'], 'r') as outputs:
        for line in outputs:
            #
            for header_value in summary_dict[genus]['table']['header'].split(',')[:-1]:
                if line.startswith('{hv}\n'.format(hv=header_value)):
                    # Iterate through the lines following the match
                    for subline in outputs:
                        if subline != '\n':
                            if subline.startswith('Mutation'):
                                for detailline in outputs:
                                    if detailline != '\n':
                                        summary_string += '{},'.format(detailline.split('\t')[0])
                                    else:
                                        break
                            else:
                                summary_string += '{},'.format(subline.replace(',', ';').replace('\t', ',').rstrip())
                                break
                        else:
                            break
                        break
    # Ensure that there were results to report
    if summary_string:
        if not summary_string.endswith('\n'):
            summary_string += '\n'
        # Write the summaries to the summary file
        with open(summary_dict[genus]['table']['summary'], 'a+') as summary:
            summary.write(summary_string)


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def pointfinder_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))
    # Parse description to get list of SeqIDs
    seqids = list()
    for i in range(0, len(description)):
        item = description[i]
        item = item.upper()
        # Minimal check to make sure IDs provided somewhat resemble a valid sample ID
        if item.isalpha():
            pass
        else:
            seqids.append(item)

    # Run Mash
    with open(os.path.join(work_dir, 'seqid.txt'), 'w') as f:
        for seqid in seqids:
            f.write(seqid + '\n')
    # Drop FASTA files into workdir
    retrieve_nas_files(seqids=seqids,
                       outdir=work_dir,
                       filetype='fasta',
                       copyflag=False)
    # Create output directory
    output_dir = os.path.join(work_dir, 'output')
    make_path(output_dir)
    # Get all of the FASTA files
    fasta_list = sorted(glob.glob(os.path.join(work_dir, '*.fasta')))
    # Set the folder to store all the PointFinder outputs
    pointfinder_output_dir = os.path.join(work_dir, 'pointfinder_outputs')
    # Initialise a dictionaries to store the mash-calculated, and pointfinder-formatted genus outputs for each strain
    genus_dict = dict()
    organism_dict = dict()
    # Create lists to store missing and unprocessed seqids
    unprocessed_seqs = list()
    missing_seqs = list()
    mash_fails = list()
    # Dictionary to convert the mash-calculated genus to the pointfinder format
    pointfinder_org_dict = {'Campylobacter': 'campylobacter',
                            'Escherichia': 'e.coli',
                            'Shigella': 'e.coli',
                            '‎Mycobacterium': 'tuberculosis',
                            'Neisseria': 'gonorrhoeae',
                            'Salmonella': 'salmonella'}
    # Reverse look-up dictionary
    rev_org_dict = {'campylobacter': 'Campylobacter',
                    'e.coli': 'Escherichia',
                    'tuberculosis': 'Mycobacterium',
                    'gonorrhoeae': 'Neisseria',
                    'salmonella': 'Salmonella'}
    summary_dict = {
        'Salmonella':
            {
                'prediction':
                    {
                        'header': 'Strain,Colitsin,Colistin,Spectinomycin,Quinolones,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'Salmonella_prediction_summary.csv')
                    },
                'table':
                    {
                        'header': 'Strain,parE,parC,gyrA,pmrB,pmrA,gyrB,16S_rrsD,23S,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'Salmonella_table_summary.csv')
                    },
                'results':
                    {
                        'header': 'Strain,Genus,Mutation,NucleotideChange,AminoAcidChange,Resistance,PMID,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'PointFinder_results_summary.csv')
                    }
            },
        'Escherichia':
            {
                'prediction':
                    {
                        'header': 'Strain,Colistin,GentamicinC,gentamicinC,Streptomycin,Macrolide,Sulfonamide,'
                                  'Tobramycin,Neomycin,Fluoroquinolones,Aminocoumarin,Tetracycline,KanamycinA,'
                                  'Spectinomycin,B-lactamResistance,Paromomycin,Kasugamicin,Quinolones,G418,'
                                  'QuinolonesAndfluoroquinolones,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'Escherichia_prediction_summary.csv')
                    },
                'table':
                    {
                        'header': 'Strain,parE,parC,folP,gyrA,pmrB,pmrA,16S_rrsB,16S_rrsH,gyrB,ampC,16S_rrsC,23S,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'Escherichia_table_summary.csv')
                    },
                'results':
                    {
                        'header': 'Strain,Genus,Mutation,NucleotideChange,AminoAcidChange,Resistance,PMID,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'PointFinder_results_summary.csv')
                    }
            },
        'Campylobacter':
            {

                'prediction':
                    {
                        'header': 'Strain,LowLevelIncreaseMIC,AssociatedWithT86Mutations,Macrolide,Quinolone,'
                                  'Streptinomycin,Erythromycin,IntermediateResistance,HighLevelResistance_'
                                  'nalidixic_and_ciprofloxacin,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'Campylobacter_prediction_summary.csv')
                    },
                'table':
                    {
                        'header': 'Strain,L22,rpsL,cmeR,gyrA,23S,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'Campylobacter_table_summary.csv')
                    },
                'results':
                    {
                        'header': 'Strain,Genus,Mutation,NucleotideChange,AminoAcidChange,Resistance,PMID,\n',
                        'output': str(),
                        'summary': os.path.join(pointfinder_output_dir, 'PointFinder_results_summary.csv')
                    }
            }
    }

    # Run mash screen on each of the assemblies
    for item in fasta_list:
        seqid = os.path.splitext(os.path.basename(item))[0]
        screen_file = os.path.join(output_dir, '{seqid}_screen.tab'.format(seqid=seqid))
        mash.screen('/mnt/nas2/databases/confindr/databases/refseq.msh',
                    item,
                    threads=8,
                    w='',
                    i='0.95',
                    output_file=screen_file,
                    returncmd=True)
        screen_output = mash.read_mash_screen(screen_file)
        # Determine the genus from the screen output file
        for screen in screen_output:
            # Extract the genus from the mash results
            mash_organism = screen.query_id.split('/')[-3]
            # Use the organism as a key in the pointfinder database name conversion dictionary
            try:
                mash_genus = pointfinder_org_dict[mash_organism]
            except KeyError:
                mash_genus = 'NA'
            # Populate the dictionaries with the seqid, and the calculated genus/pointfinder name
            genus_dict[seqid] = mash_genus
            organism_dict[seqid] = mash_organism
    # Delete all of the FASTA files
    for fasta in fasta_list:
        os.remove(fasta)
    # # Delete the output folder
    # shutil.rmtree(output_dir)

    # Pointfinder
    # These unfortunate hard coded paths appear to be necessary
    activate = 'source /home/ubuntu/miniconda3/bin/activate /mnt/nas2/virtual_environments/pointfinder'
    pointfinder_py = '/mnt/nas2/virtual_environments/pointfinder/pointfinder-3.0/pointfinder-3.0.py'
    # Database locations
    pointfinder_db = '/mnt/nas2/databases/assemblydatabases/0.3.4/pointfinder'
    # List of organisms in the pointfinder database
    pointfinder_list = ['campylobacter', 'e.coli', 'tuberculosis', 'gonorrhoeae', 'salmonella']
    os.makedirs(pointfinder_output_dir, exist_ok=True)
    # Pointfinder cannot handle an entire folder of sequences; each sample must be processed independently
    for seqid in sorted(seqids):
        # If the seqid isn't present in the dictionary, it is because the assembly could not be found - or because
        # MASH screen failed
        try:
            # Look up the PointFinder and the MASH-calculated genera
            pointfinder_genus = genus_dict[seqid]
            genus = rev_org_dict[pointfinder_genus]
            # If the genus isn't in the pointfinder database, do not attempt to process it
            if pointfinder_genus in pointfinder_list:
                # Create folder to drop FASTA files
                assembly_folder = os.path.join(work_dir, seqid)
                make_path(assembly_folder)
                # Extract FASTA files.
                retrieve_nas_files(seqids=[seqid],
                                   outdir=assembly_folder,
                                   filetype='fasta',
                                   copyflag=False)
                fasta = os.path.join(assembly_folder, '{seqid}.fasta'.format(seqid=seqid))
                # Prepare command
                cmd = 'python {py} -i {fasta} -s {orgn} -p {db} -o {output} -m blastn -m_p {blast_path}'\
                    .format(py=pointfinder_py,
                            fasta=fasta,
                            orgn=pointfinder_genus,
                            db=pointfinder_db,
                            output=pointfinder_output_dir,
                            blast_path='/mnt/nas2/virtual_environments/pointfinder/bin/blastn'
                            )
                # Create another shell script to execute within the PlasmidExtractor conda environment
                template = "#!/bin/bash\n{} && {}".format(activate, cmd)
                pointfinder_script = os.path.join(work_dir, 'run_pointfinder.sh')
                with open(pointfinder_script, 'w+') as file:
                    file.write(template)
                # Modify the permissions of the script to allow it to be run on the node
                make_executable(pointfinder_script)
                # Run shell script
                os.system(pointfinder_script)
                # Find the pointfinder outputs
                summary_dict[genus]['prediction']['output'] = \
                    glob.glob(os.path.join(pointfinder_output_dir, '{seq}*prediction.txt'.format(seq=seqid)))[0]
                summary_dict[genus]['table']['output'] = \
                    glob.glob(os.path.join(pointfinder_output_dir, '{seq}*table.txt'.format(seq=seqid)))[0]
                summary_dict[genus]['results']['output'] = \
                    glob.glob(os.path.join(pointfinder_output_dir, '{seq}*results.txt'.format(seq=seqid)))[0]
                # Process the predictions
                write_report(summary_dict=summary_dict,
                             seqid=seqid,
                             genus=genus,
                             key='prediction')
                # Process the results summary
                write_report(summary_dict=summary_dict,
                             seqid=seqid,
                             genus=genus,
                             key='results')
                # Process the table summary
                write_table_report(summary_dict=summary_dict,
                                   seqid=seqid,
                                   genus=genus)
            else:
                unprocessed_seqs.append(seqid)
        except KeyError:
            if not os.path.isfile(os.path.join(output_dir, '{seq}_screen.tab'.format(seq=seqid))):
                missing_seqs.append(seqid)
            else:
                mash_fails.append(seqid)
    # Attempt to clear out the tmp folder from the pointfinder_output_dir
    try:
        shutil.rmtree(os.path.join(pointfinder_output_dir, 'tmp'))
    except FileNotFoundError:
        pass
    # Zip output
    output_filename = 'pointfinder_output'
    zip_filepath = zip_folder(results_path=pointfinder_output_dir,
                              output_dir=work_dir,
                              output_filename=output_filename)
    zip_filepath += '.zip'
    # Prepare upload
    output_list = [
        {
            'filename': os.path.basename(zip_filepath),
            'path': zip_filepath
        }
    ]

    # Create a note to add to the updated Redmine issue
    notes = 'Pointfinder process complete!'
    # If there are missing, or unprocessed sequences, add details to the note
    if unprocessed_seqs:
        seq_list = list()
        for sequence in unprocessed_seqs:
            seq_list.append('{seqid} ({organism})'.format(seqid=sequence,
                                                          organism=organism_dict[sequence]))
        if len(unprocessed_seqs) > 1:
            notes += '\n The following sequences were not processed, as they were determined to be genera not ' \
                     'present in the pointfinder database: {seqs}'.format(seqs=', '.join(seq_list))
        else:
            notes += '\n The following sequence was not processed, as it was determined to be a genus not ' \
                     'present in the pointfinder database: {seqs}'.format(seqs=', '.join(seq_list))
    if missing_seqs:
        if len(missing_seqs) > 1:
            notes += '\n The following sequences were not processed, as they could not be located in the strain ' \
                     'database: {seqs}'.format(seqs=', '.join(missing_seqs))
        else:
            notes += '\n The following sequence was not processed, as it could not be located in the strain database:' \
                     ' {seqs}'.format(seqs=', '.join(missing_seqs))
    if mash_fails:
        if len(mash_fails) > 1:
            notes += '\n The following sequences could not be processed by MASH screen: {seqs}'\
                .format(seqs=', '.join(mash_fails))
        else:
            notes += '\n The following sequence could not be processed by MASH screen: {seqs}'\
                .format(seqs=', '.join(mash_fails))
    # Create a list of all the folders - will be used to clean up the working directory
    folders = glob.glob(os.path.join(work_dir, '*/'))
    # Remove all the folders
    for folder in folders:
        if os.path.isdir(folder):
            shutil.rmtree(folder)
    # Wrap up issue
    redmine_instance.issue.update(resource_id=issue.id,
                                  uploads=output_list,
                                  status_id=4,
                                  notes=notes)


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
    pointfinder_redmine()
