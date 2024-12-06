import os
import re
import glob
import click
import pickle
import shutil
import pandas as pd
from biotools import mash
from nastools.nastools import retrieve_nas_files


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def strainmash_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    # Reference path
    typestrain_db_path = '/mnt/nas2/databases/mashdatabases/RefSeqSketches.msh'

    # Parse description to get list of SeqIDs
    seqids = []
    for i in range(0, len(description)):
        item = description[i]
        item = item.upper()
        seqids.append(item)

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
    os.makedirs(output_dir, exist_ok=True)

    # Get all of the FASTA files
    fasta_list = glob.glob(os.path.join(work_dir, '*.fasta'))

    # Run mash_screen on everything
    for item in fasta_list:
        output_filepath = os.path.join(output_dir, (item.replace('.fasta', '') + '_strainmash.txt'))

        mash_screen(reference=typestrain_db_path,
                    queryfile=item,
                    outname=output_filepath)


    # Move all files to the actual output folder - this is a workaround to a weird bug
    output_files = glob.glob(os.path.join(work_dir, '*strainmash.txt'))
    for file in output_files:
        os.rename(file, os.path.join(output_dir, os.path.basename(file)))

    # Zip output folder
    shutil.make_archive(output_dir, 'zip', work_dir, 'output')

    # Glob zip
    zip_file = glob.glob(os.path.join(work_dir, '*.zip'))[0]

    # Output list containing dictionaries with file path as the key
    output_list = [
        {
            'path': os.path.join(work_dir, zip_file),
            'filename': os.path.basename(zip_file)
        }
    ]

    # Upload files, set status to Feedback
    redmine_instance.issue.update(resource_id=issue.id,
                                  uploads=output_list,
                                  status_id=4,
                                  notes='STRAINMASH complete. See attached file for results.')

    # Delete all of the FASTA files
    for fasta in fasta_list:
        os.remove(fasta)

    # Delete the output folder
    shutil.rmtree(output_dir)


def mash_screen(reference, queryfile, outname):
    mash.screen(reference, queryfile, output_file=outname, w='')
    mash2pandas(outname)


def mash2pandas(outname):
    # Pandas
    df = pd.read_csv(outname, sep='\t', header=None)
    try:
        df = df.drop(df.columns[[5]], axis=1)
    except:
        pass
    df = df.loc[df[0] >= 0.7]
    df = df.sort_values(df.columns[0], ascending=False)

    # Get rid of path from names in output (artifact of the straindb)
    df[4] = df[4].apply(lambda x: os.path.basename(x))

    # Grab only the accession with some regex
    df[4] = df[4].apply(lambda x: re.findall(r"^\w{3}_\d+\.\d", x)[0])

    # Reset the index
    df.reset_index(drop=True, inplace=True)

    # Take only the top 10 rows to significantly cut down on Entrez retrieval time
    df = df[:10]

    # Map the extract_entrez_data function on each accession
    df['Organism'] = df[4].map(extract_species)

    # Print out top hit
    tophit_name = df[4].iloc[0]
    tophit_score = df[1].iloc[0]
    print(df)
    print('\nTop Hit: {}\nScore: {}'.format(tophit_name, tophit_score))

    # Header names
    colnames = [
        'MashDistance',
        'NumMatchingHashes',
        'MedianMultiplicity',
        'Pvalue',
        'ReferenceStrain',
        'Organism'
    ]

    # Overwrite mash.screen output with a nice Pandas df
    try:
        df.to_csv(outname, sep='\t', header=colnames, index=False)
    except:
        df.to_csv(outname, sep='\t', header=None, index=False)


def extract_species(accession):
    try:
        df = pd.read_csv('/mnt/nas/Databases/RefSeq/assembly_summary_refseq.txt',
                         skiprows=1,
                         delimiter='\t')
        df_filtered = df[df.values == accession]
        organism = df_filtered['organism_name'].values[0]

        # Consistency
        if organism is 'na':
            organism = 'N/A'
    except:
        organism = 'N/A'

    return organism


if __name__ == '__main__':
    strainmash_redmine()
