import os
import glob
import click
import ftplib
import pickle
import shutil
from automator_settings import FTP_USERNAME, FTP_PASSWORD


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def reportretrieve_redmine(redmine_instance, issue, work_dir, description):
    print('External retrieving!')
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        os.makedirs(os.path.join(work_dir, str(issue.id)))
        # Parse description to figure out what SEQIDs we need to run on.
        seqid_list = list()
        for item in description:
            item = item.upper()
            if item != '':
                seqid_list.append(item)

        report_path_list = list()
        # Go through CombinedMetadata sheets to find which folders we need to copy to FTP.
        metadata_sheets = glob.glob('/mnt/nas/WGSspades/*/reports/combinedMetadata.csv')
        metadata_sheets += glob.glob('/mnt/nas/External_WGSspades/*/*/reports/combinedMetadata.csv')

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
            cmd = 'cp {report_folder} {new_folder}'.format(report_folder=report_folder[0],
                                                           new_folder=os.path.join(work_dir, str(issue.id), new_folder_name))
            os.system(cmd)

        # Now make a zip folder that we'll upload to the FTP.
        shutil.make_archive(root_dir=os.path.join(work_dir, str(issue.id)),
                            format='zip',
                            base_name=os.path.join(work_dir, str(issue.id)))

        # Now need to login to the FTP to upload the zipped folder.
        s = ftplib.FTP('ftp.agr.gc.ca', user=FTP_USERNAME, passwd=FTP_PASSWORD)
        s.cwd('outgoing/cfia-ak')
        f = open(os.path.join(work_dir, str(issue.id) + '.zip'), 'rb')
        s.storbinary('STOR {}.zip'.format(str(issue.id)), f)
        f.close()
        s.quit()

        # And finally, do some file cleanup.
        try:
            shutil.rmtree(os.path.join(work_dir, str(issue.id)))
            os.remove(os.path.join(work_dir, str(issue.id) + '.zip'))
        except:
            pass

        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='Report Retrieve process complete!\n\n'
                                            'Results are available at the following FTP address:\n'
                                            'ftp://ftp.agr.gc.ca/outgoing/cfia-ak/{}'.format(str(issue.id) + '.zip'))
    except Exception as e:
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Send this error traceback to your friendly '
                                            'neighborhood bioinformatician: {}'.format(e))


if __name__ == '__main__':
    reportretrieve_redmine()
