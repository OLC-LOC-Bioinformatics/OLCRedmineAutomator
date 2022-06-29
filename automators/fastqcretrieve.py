import os
import glob
import click
import ftplib
import pickle
import shutil
import sentry_sdk
from amrsummary import before_send
from automator_settings import SENTRY_DSN
from externalretrieve import upload_to_ftp
from automator_settings import FTP_USERNAME, FTP_PASSWORD


@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def fastqcretrieve_redmine(redmine_instance, issue, work_dir, description):
    sentry_sdk.init(SENTRY_DSN, before_send=before_send)
    print('External retrieving!')
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
#        os.makedirs(os.path.join(work_dir, str(issue.id)))
        # Variable to hold supplied arguments
        #amino acid substitution models supported by automator
        platforms = [
            'miseq', 'nextseq'
        ]

        argument_dict = {
            'run': str(),
            'platform': 'miseq'
        }
        # Parse description to figure out what SEQIDs we need to run on.
#        seqrun = list()
        for item in description:
            item = item.upper()
            if 'RUN' in item:
                argument_dict['run'] = item.split('=')[1].upper()
                continue
            if 'PLATFORM' in item:
                argument_dict['platform'] = item.split('=')[1].lower()
                continue
#            seqrun.append(item)

        # Warn the user if no runid was provided.
        # make folder in the redmine biorequest working directory
        new_folder = os.path.join(work_dir, 'fastqc')
        os.makedirs(new_folder)
        raw_dir = os.path.join(new_folder, 'Raw')
        os.makedirs(raw_dir)
        trcr_dir = os.path.join(new_folder, 'TrimmedCorrected')
        os.makedirs(trcr_dir)
        # Go to the run, and move the fastqc files into a folder which we will then move to our working directory
        if argument_dict['platform'] == 'miseq':
            miseq = '/mnt/nas2/processed_sequence_data/miseq_assemblies/'
            run = '/mnt/nas2/processed_sequence_data/miseq_assemblies/{run}/'.format(run=argument_dict['run'])
            subdirectories = os.listdir(path='/mnt/nas2/processed_sequence_data/miseq_assemblies/{run}/'.format(run=argument_dict['run']))
        if argument_dict['platform'] == 'nextseq':
            miseq = '/mnt/nas2/processed_sequence_data/nextseq_assemblies/'
            run = '/mnt/nas2/processed_sequence_data/nextseq_assemblies/{run}/'.format(run=argument_dict['run'])
            subdirectories = os.listdir(path='/mnt/nas2/processed_sequence_data/nextseq_assemblies/{run}/'.format(run=argument_dict['run']))
        print(subdirectories)
#        fastqc_dir = os.path.join(run, 'fastqc')
#        os.makedirs(fastqc_dir)
#        raw_dir = os.path.join(fastqc_dir, 'Raw')
#        os.makedirs(raw_dir)
#        trcr_dir = os.path.join(fastqc_dir, 'TrimmedCorrected')
#        os.makedirs(trcr_dir)
        for subdirectorie in subdirectories:
            if subdirectorie != 'BestAssemblies' and subdirectorie != 'confindr' and subdirectorie != 'InterOp' and subdirectorie != 'raw_assemblies' and subdirectorie != 'reports' and not subdirectorie.endswith(".gz") and not subdirectorie.endswith(".txt") and not subdirectorie.endswith(".sh") and not subdirectorie.endswith(".xml") and not subdirectorie.endswith(".csv") and not subdirectorie.endswith(".out"):
                sub_raws = os.listdir(os.path.join(run, subdirectorie, 'fastqc/Raw'))
                print(sub_raws)
                sub_trcrs = os.listdir(os.path.join(run, subdirectorie, 'fastqc/trimmedcorrected'))
                print(sub_trcrs)
                for file in sub_raws:
                    if file.endswith(".html"):
                        shutil.copyfile(os.path.join(run, subdirectorie, 'fastqc/Raw', file), os.path.join(raw_dir, file))
                for file in sub_trcrs:
                    if file.endswith(".html"):
                        shutil.copyfile(os.path.join(run, subdirectorie, 'fastqc/trimmedcorrected', file), os.path.join(trcr_dir, file))
                        continue

        # Now make a zip folder that we'll upload to the FTP.
        output_filename = 'fastqc_{}'.format(issue.id)
        zip_filepath = zip_folder(results_path=new_folder,
                                  output_dir=work_dir,
                                  output_filename=output_filename)
        zip_filepath += '.zip'

        upload_successful = upload_to_ftp(local_file=zip_filepath)

        if upload_successful:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Fastqc Retrieve process complete!\n\n'
                                                'Results are available at the following FTP address:\n'
                                                'ftp://ftp.agr.gc.ca/outgoing/cfia-ac/{l}'.format(l=os.path.split(zip_filepath)[1]))
        else:
            redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                          notes='Upload of result files was unsuccessful due to FTP connectivity issues. '
                                                'Please try again later.')

        # And finally, do some file cleanup.
        shutil.rmtree(new_folder)
        os.remove(zip_filepath)

    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! We log this automatically and will look into the '
                                            'problem and get back to you with a fix soon.')

def zip_folder(results_path, output_dir, output_filename):
    output_path = os.path.join(output_dir, output_filename)
    shutil.make_archive(output_path, 'zip', results_path)
    return output_path

if __name__ == '__main__':
    fastqcretrieve_redmine()
