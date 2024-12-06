import os
import glob
import click
import ftplib
import pickle
import shutil
import sentry_sdk
import pandas
from amrsummary import before_send
from automator_settings import SENTRY_DSN
# Dropbox
from upload_to_dropbox import upload_to_dropbox
from tokens import (
    DROPBOX_ACCESS_TOKEN,
    DROPBOX_APP_KEY, 
    DROPBOX_APP_SECRET,
    DROPBOX_REFRESH_TOKEN
)


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def reportretrieve_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    print('External retrieving!')
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        #make a directory to hold the report folders
        out_dir = os.path.join(work_dir, str(issue.id))
        os.makedirs(out_dir, exist_ok=True)
        # Parse description to figure out what SEQIDs we need to run on.
        seqid_list = list()
        seqids = list() #Ashley added this so can create another list of seqids because subsetting wasn't working and I didnt want to mess with OG code
        for item in description:
            item = item.upper()
            if item != '':
                seqid_list.append(item)
            seqids.append(item)

        report_path_list = list()
        # Go through CombinedMetadata sheets to find which folders we need to copy to FTP.
        metadata_sheets = glob.glob('/mnt/nas2/processed_sequence_data/*/*/reports/combinedMetadata.csv')

        for metadata_sheet in metadata_sheets:
            with open(metadata_sheet) as csvfile:
                lines = csvfile.readlines()
                for i in range(1, len(lines)):
                    x = lines[i].split(',')
                    if x[0] in seqid_list:  # First entry in the row should be the SEQID.
                        report_path = os.path.abspath(metadata_sheet)
                        if report_path not in report_path_list:
                            report_path_list.append(report_path)
                        seqid_list.remove(x[0])

        # Warn the user if reports couldn't be found for some SEQIDs.
        if len(seqid_list) > 0:
            redmine_instance.issue.update(resource_id=issue.id,
                                          notes='WARNING: Could not find reports for the following SEQIDs: '
                                                '{}'.format(seqid_list))

        # Go through the report path list and copy reports folders, while renaming them.
        for report in report_path_list:
            complete_path = os.path.abspath(report)
            report_folder = os.path.split(complete_path)[:-1]
            new_folder_name = report_folder[0].split('/')[-2] + '_' + report_folder[0].split('/')[-1]
            if os.path.isdir(os.path.join(work_dir, str(issue.id), new_folder_name)):  # Very slim possiblity two folders
                # could have the same name. This takes care of that. No way there should ever be more than two.
                new_folder_name = new_folder_name + '_2'
            cmd = 'cp -r {report_folder} {new_folder}'.format(report_folder=report_folder[0],
                                                              new_folder=os.path.join(work_dir, str(issue.id), new_folder_name))
            os.system(cmd)

        #concatenate the reports from all folders into one new folder/report file
        #create a folder to hold the sequences
        combined_dir = os.path.join(work_dir,str(issue.id),'combined_reports')
        os.makedirs(combined_dir, exist_ok=True)

        # Now let's combine all of the reports from all of the subfolders in the output directory
        metadata_files = glob.glob("{}/{}/*/legacy_combinedMetadata.csv".format(work_dir, str(issue.id)))

        # Open a new outfile that will be the concatenation of the results files
        mdoutputfile = os.path.join(combined_dir, 'legacy_combinedMetadata.csv')
        with open(mdoutputfile, 'w') as outfile:
            for f in metadata_files:
                with open(f, 'r') as infile:
                    for line in infile:
                        outfile.write(line.rstrip('\n') + '\n')
        
        #read our combined dataframe back into a pandas dataframe we will call df
        print(seqids)
        df = pandas.read_csv(mdoutputfile)
        df2 = df[df['SeqID'].isin(seqids)]
        
        # Create a new file name for the subsetted dataframe
        mdfilename = 'legacy_combinedMetadata_redmine-{}.csv'.format(issue.id)
        mdfileloc = os.path.join(combined_dir,mdfilename)
        
        # Write the subsetted dataframe to a new file
        try:
            df2.to_csv(mdfileloc, index=False)
        # If there are any errors, just copy the original file to the new
        # location
        except CError:
            shutil.copy(mdoutputfile, mdfileloc)

        #upload the subsetted combined_metadata reports file to the redmine request
        output_list = list()
        output_dict = dict()
        output_dict['path'] = os.path.join(combined_dir,mdfilename)
        output_dict['filename'] = mdfilename
        output_list.append(output_dict)
        
        redmine_instance.issue.update(
            resource_id=issue.id,
            uploads=output_list,
            notes='Report Retrieve process complete!\n\n'
            'Results are attached.'
        )

        # Now make a zip folder that we'll upload to the FTP.
        shutil.make_archive(
            root_dir=os.path.join(work_dir, str(issue.id)),
            format='zip',
            base_name=os.path.join(work_dir, str(issue.id))
        )

        # Set the path of the zip file
        zip_filepath = os.path.join(work_dir, str(issue.id) + '.zip')
        
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
                status_id=2,
                notes='Report retrieve process complete!\n\n'
                      'Results are available at the following URL:\n'
                      '{url}'.format(url=download_link)
            )
        else:
            redmine_instance.issue.update(
                resource_id=issue.id,
                status_id=2,
                notes='Upload of reports was unsuccessful due to '
                'connectivity issues. Please try again later.'
            )

        # And finally, do some file cleanup.
        try:
            #shutil.rmtree(os.path.join(work_dir, str(issue.id)))
            os.remove(os.path.join(work_dir, str(issue.id) + '.zip'))
        except:
            pass

        #remove the folder now that it has uploaded files to correct locations
        shutil.rmtree(os.path.join(work_dir, str(issue.id)))

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon: {}'.format(e))


if __name__ == '__main__':
    reportretrieve_redmine()
