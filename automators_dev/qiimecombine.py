import os
import glob
import click
import pickle
import zipfile
import datetime
import pandas as pd


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def qiimecombine_redmine(redmine_instance, issue, work_dir, description):
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    """
    Description is expected to be in the following format.
    level=taxlevel
    column_header=column_sample  # Allow searching on more than one column header by 
    optional! startdate-enddate. If not provided, assume we want to search all.
    
    So any example would be:
    level=5
    sample_type=meat
    180401-190408
    """
    # Parse description
    try:
        level, column_headers, column_contents, start_date, end_date = parse_description(description)
    except:  # Don't try to be too specific here, this could blow up in many ways.
        redmine_instance.issue.update(resource_id=issue.id,
                                      status_id=4,
                                      notes='Your description was not formatted correctly. See DOCS for a description '
                                            'of how to format the description.')
        return

    tax_barplots = glob.glob('/mnt/nas2/processed_sequence_data/miseq_assemblies/*/qiime2/taxonomy_barplot.qzv')
    # As it turns out, qiime2 qzv files are actually just zip files with a bunch of data/metadata.
    # https://github.com/joey711/phyloseq/issues/830
    # Getting at data seems to be easiest if we just unzip and read in the relevant csv files with pandas
    for tax_barplot in tax_barplots:
        # Check the date of the assembly to see if it's in our range.
        date_string = tax_barplot.split('/')[-3]
        run_date = string_to_year(date_string)
        if start_date < run_date < end_date:
            cmd = 'cp {} {}'.format(tax_barplot, work_dir)
            os.system(cmd)
            output_dir = os.path.join(work_dir, tax_barplot.split('/')[-3])
            with zipfile.ZipFile(os.path.join(work_dir, 'taxonomy_barplot.qzv'), 'r') as zipomatic:
                zipomatic.extractall(output_dir)

    # Now we should have folders for all of our QIIME2 runs (named with as run dates).
    # Grab the csv file for level of interest for all of them.
    dataframe_list = list()
    csv_files = glob.glob(os.path.join(work_dir, '*', '*', 'data', 'level-{}.csv'.format(level)))
    for csv_file in csv_files:
        df = pd.read_csv(csv_file)
        for column in df.columns:
            df.rename(columns={column: column.replace(' ', '_').upper()}, inplace=True)
        try:
            df = df.apply(lambda x: x.astype(str).str.upper())
            for j in range(len(column_headers)):
                column_header = column_headers[j]
                column_content = column_contents[j]
                df = df.loc[df[column_header] == column_content]
            dataframe_list.append(df)
        except KeyError:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find column {} in run {}'.format(column_header,
                                                                                                     csv_file.split('/')[-4]))
    output_list = list()
    output_dict = dict()
    output_dict['path'] = os.path.join(work_dir, 'results.csv')
    output_dict['filename'] = 'QIIMEcombine_results.csv'
    output_list.append(output_dict)
    result_df = pd.concat(dataframe_list, ignore_index=True, sort=False)
    result_df.fillna(0, inplace=True)
    result_df.to_csv(os.path.join(work_dir, 'results.csv'), index=False)

    # Clean up files.
    cmd = 'rm -r {}'.format(os.path.join(work_dir, '*/'))
    os.system(cmd)

    # Upload files, set status to Feedback
    redmine_instance.issue.update(resource_id=issue.id,
                                  status_id=4,
                                  uploads=output_list,
                                  notes='QIIMECombine complete. See attached file for results.')


def string_to_year(yymmdd_string):
    year = int('20' + yymmdd_string[:2])
    month = int(yymmdd_string[2:4])
    day = int(yymmdd_string[4:6])
    return datetime.datetime(year, month, day)


def parse_description(description):
    level = description[0].split('=')[1]
    column_headers = list()
    column_content = list()
    if '=' not in description[-1]:  # If last line of description doesn't have =, it contains a date.
        for i in range(1, len(description) - 1):
            column_headers.append(description[i].split('=')[0].upper().replace(' ', '_'))
            column_content.append(description[i].split('=')[1].upper())
        start_date = description[-1].split('-')[0]
        start_date = string_to_year(start_date)
        end_date = description[-1].split('-')[1]
        end_date = string_to_year(end_date)
    else:  # If no date range included in description, have date range include everything.
        for i in range(1, len(description)):
            column_headers.append(description[i].split('=')[0].upper().replace(' ', '_'))
            column_content.append(description[i].split('=')[1].upper())
        start_date = string_to_year('111111')
        end_date = string_to_year('991111')

    return level, column_headers, column_content, start_date, end_date


if __name__ == '__main__':
    qiimecombine_redmine()
