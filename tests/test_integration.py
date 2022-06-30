#!/usr/bin/env python

from validator_helper import validate
from redminelib import Redmine
import urllib.request
import sentry_sdk
import hashlib
import argparse
import tempfile
import logging
import zipfile
import time
import os

"""
Test script to make sure that adding an automator/changing some dependency somewhere doesn't break all the things.
"""

DEV_REDMINE_URL = 'http://192.168.1.2:8080'
REDMINE_PROJECT = 'test'
OUTGOING_FTP = 'ftp://ftp.agr.gc.ca/outgoing/cfia-ak'


def create_issue(redmine_instance, subject, description, attachments=None):
    if attachments is None:
        issue = redmine_instance.issue.create(project_id=REDMINE_PROJECT,
                                              subject=subject,
                                              description=description)
    else:
        issue = redmine_instance.issue.create(project_id=REDMINE_PROJECT,
                                              subject=subject,
                                              description=description,
                                              uploads=attachments)

    return issue.id


def create_test_issues(redmine):

    """
    Creates issues to be tested.
    :param redmine: instantiated redmine instance
    :return: dictionary with issue subjects as keys and redmine issue IDs as values.
    """

    logging.info('Creating Issues!')
    # CREATE ALL THE ISSUES - store the issue ID for each of the issues created in the dictionary.
    issues_created = dict()

    issues_created['amrsummary'] = create_issue(redmine_instance=redmine,
                                                subject='AMR Summary',
                                                description='2014-SEQ-0276')
    issues_created['autoclark'] = create_issue(redmine_instance=redmine,
                                               subject='autoclark',
                                               description='fasta\n2014-SEQ-0276')
    issues_created['ecgf'] = create_issue(redmine_instance=redmine,
                                          subject='ecgf',
                                          description='2018-LET-0106')
    issues_created['ec_typer'] = create_issue(redmine_instance=redmine,
                                              subject='ec_typer',
                                              description='2014-SEQ-0276')
    issues_created['externalretrieve'] = create_issue(redmine_instance=redmine,
                                                      subject='External Retrieve',
                                                      description='fasta\n2014-SEQ-0276')
    issues_created['diversitree'] = create_issue(redmine_instance=redmine,
                                                 subject='diversitree',
                                                 description='2\n2014-SEQ-0275\n2014-SEQ-0276\n2014-SEQ-0277\n2014-SEQ-0278')
    issues_created['geneseekr'] = create_issue(redmine_instance=redmine,
                                               subject='geneseekr',
                                               description='analysis=sixteens\n2014-SEQ-0276')
    issues_created['pointfinder'] = create_issue(redmine_instance=redmine,
                                                 subject='pointfinder',
                                                 description='2014-SEQ-0276')
    issues_created['prokka'] = create_issue(redmine_instance=redmine,
                                            subject='prokka',
                                            description='2014-SEQ-0276')
    issues_created['resfinder'] = create_issue(redmine_instance=redmine,
                                               subject='resfinder',
                                               description='2014-SEQ-0276')
    issues_created['sipprverse'] = create_issue(redmine_instance=redmine,
                                                subject='sipprverse',
                                                description='analysis=mlst\n2014-SEQ-0276')
    issues_created['staramr'] = create_issue(redmine_instance=redmine,
                                             subject='staramr',
                                             description='2014-SEQ-0282')
    logging.info('Issue creation complete.')
    return issues_created


def validate_csv_content(query_csv, ref_csv):
    if query_csv.endswith('.tsv'):
        column_list = validate.find_all_columns(csv_file=ref_csv, columns_to_exclude=[], range_fraction=0.05, separator='\t')
    else:
        column_list = validate.find_all_columns(csv_file=ref_csv, columns_to_exclude=[], range_fraction=0.05)
    if 'ec_typer_report.tsv' in ref_csv:  # EC typer report doesn't play nicely with validator helper, so just take an md5sum
        if md5(ref_csv) == md5(query_csv):
            return True
        else:
            return False
    validator = validate.Validator(reference_csv=ref_csv,
                                   test_csv=query_csv,
                                   column_list=column_list,
                                   identifying_column='Strain')
    if validator.same_columns_in_ref_and_test() and validator.all_test_columns_in_ref_and_test() and validator.check_samples_present() and validator.check_columns_match():
        return True
    else:
        return False


def validate_attachments(issue, file_to_find, validate_content=False):
    validation_status = 'Unknown'
    if len(issue.attachments) > 0:
        attachment_found = False
        for attachment in issue.attachments:
            if attachment.filename == file_to_find:
                attachment_found = True
                with tempfile.TemporaryDirectory() as tmpdir:
                    attachment.download(savepath=tmpdir, filename=file_to_find)
                    if os.path.getsize(os.path.join(tmpdir, file_to_find)) > 0:
                        if validate_content is True:
                            ref_csv = 'tests/ref_csvs/{}'.format(file_to_find)
                            content_ok = validate_csv_content(query_csv=os.path.join(tmpdir, file_to_find),
                                                              ref_csv='tests/ref_csvs/{}'.format(file_to_find))
                            if content_ok is True:
                                validation_status = 'Validated'
                            else:
                                validation_status = 'Reference content does not match query.'
                        else:
                                validation_status = 'Validated'
                    else:
                        validation_status = 'Uploaded file has zero size.'
        if attachment_found is False:
            validation_status = 'Could not find {}'.format(file_to_find)
    else:
        validation_status = 'No attachments, validation failed.'
    return validation_status


def md5(fname):
    """
    Gets MD5sum for a file - use for quick and dirty comparisons of files.
    Hooray for stackoverflow: https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    """
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def validate_csv_in_zip(issue, zip_file, report_file, ref_csv):
    """
    Crude validation that checks that CSV contents are the exact same between some folder in a zip
    file and a reference that we keep on hand.
    :param issue:
    :param zip_file:
    :param report_file:
    :param ref_csv:
    :return:
    """
    validation_status = 'Unknown'
    if len(issue.attachments) > 0:
        attachment_found = False
        for attachment in issue.attachments:
            if attachment.filename == zip_file:
                attachment_found = True
                with tempfile.TemporaryDirectory() as tmpdir:
                    downloaded_zip = os.path.join(tmpdir, zip_file)
                    attachment.download(savepath=tmpdir, filename=zip_file)
                    zipped_archive = zipfile.ZipFile(downloaded_zip)
                    zipped_archive.extractall(path=tmpdir)
                    if md5(ref_csv) == md5(os.path.join(tmpdir, report_file)):
                        validation_status = 'Validated'
                    else:
                        validation_status = 'Contents not identical. Check on what changed.'
        if attachment_found is False:
            validation_status = 'Could not find {}'.format(zip_file)
    else:
        validation_status = 'No attachments, validation failed.'
    return validation_status


def validate_ftp_upload(ftp_file, min_file_size=0):
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            urllib.request.urlretrieve(os.path.join(OUTGOING_FTP, ftp_file),
                                       os.path.join(tmpdir, ftp_file))
            if os.path.getsize(os.path.join(tmpdir, ftp_file)) > min_file_size:
                validation_status = 'Validated'
            else:
                validation_status = 'Uploaded file has zero size.'
    except:  # If file doesn't exist, something went wrong. Set validation status accordingly.
        validation_status = 'No uploaded file found on FTP'
    return validation_status


def monitor_issues(redmine, issue_dict, timeout):
    """
    Monitors issues for completion.
    :param redmine: Instantiated redmine instance.
    :param issue_dict: Dictionary created by create_test_issues (issue subjects as keys, issue IDs as values)
    :param timeout: Number of seconds to wait for all issues to finish running
    :return: Dictionary with issue subjects as keys, and some sort of status message for each as values
    """
    all_complete = False
    total_time_taken = 0
    increment = 30  # Check on issues every increment seconds
    # Create dictionary to track progress of each issue created.
    issues_validated = dict()
    for issue_subject in issue_dict:
        issues_validated[issue_subject] = 'Unknown'

    # Now loop through all the things!
    while all_complete is False and total_time_taken < timeout:
        logging.info('Checking tasks for completion')
        all_complete = True  # Set flag to true - if any of our issues are not finished, this gets set back to False
        # Check all issues for completion.
        for issue_subject in issue_dict:
            issue_id = issue_dict[issue_subject]
            issue = redmine.issue.get(issue_id, include=['attachments'])

            # AMR Summary #####
            if issue_subject == 'amrsummary' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:  # 4 corresponds to issue being complete.
                    issues_validated[issue_subject] = validate_attachments(issue, 'amr_summary.csv',
                                                                           validate_content=True)
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:  # If hasn't finished yet, just wait!
                    all_complete = False

            # CLARK #####
            elif issue_subject == 'autoclark' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    # TODO: Validate xlsx files?
                    issues_validated[issue_subject] = validate_attachments(issue, 'abundance.xlsx')
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False

            # ECGF #####
            elif issue_subject == 'ecgf' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_csv_in_zip(issue=issue,
                                                                          zip_file='eCGF_output_{}.zip'.format(issue.id),
                                                                          report_file='{}_summary_report.csv'.format(issue.id),
                                                                          ref_csv='/mnt/nas2/redmine/applications/OLCRedmineAutomator/tests/ref_csvs/124_summary_report.csv')  # Hardcoded paths stink, fix me!
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False

            # EC TYPER #####
            elif issue_subject == 'ec_typer' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_attachments(issue, 'ec_typer_report.tsv',
                                                                           validate_content=True)
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False

            # EXTERNAL RETRIEVE #####
            elif issue_subject == 'externalretrieve' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_ftp_upload('{}.zip'.format(issue.id),
                                                                          min_file_size=1000)
                else:
                    all_complete = False

            # DIVERSITREE #####
            elif issue_subject == 'diversitree' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    # TODO: HTML validation - should just be able to md5sum or something on file contents.
                    issues_validated[issue_subject] = validate_attachments(issue, 'diversitree_report.html')
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False

            # GENESEEKR #####
            elif issue_subject == 'geneseekr' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_csv_in_zip(issue=issue,
                                                                          zip_file='geneseekr_output.zip',
                                                                          report_file='2014-SEQ-0276_blastn_sixteens_full.tsv',
                                                                          ref_csv='/mnt/nas2/redmine/applications/OLCRedmineAutomator/tests/ref_csvs/2014-SEQ-0276_blastn_sixteens_full.tsv')  # Hardcoded paths stink, fix me!
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False

            # POINTFINDER #####
            elif issue_subject == 'pointfinder' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_csv_in_zip(issue=issue,
                                                                          zip_file='pointfinder_output.zip',
                                                                          report_file='PointFinder_results_summary.csv',
                                                                          ref_csv='/mnt/nas2/redmine/applications/OLCRedmineAutomator/tests/ref_csvs/PointFinder_results_summary.csv')  # Hardcoded paths stink, fix me!
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False

            # PROKKA #####
            elif issue_subject == 'prokka' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_ftp_upload('prokka_output_{}.zip'.format(issue.id),
                                                                          min_file_size=1000)
                else:
                    all_complete = False

            # RESFINDER #####
            elif issue_subject == 'resfinder' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_attachments(issue, 'resfinder_blastn.xlsx')
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False

            # SIPPRVERSE #####
            elif issue_subject == 'sipprverse' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_csv_in_zip(issue=issue,
                                                                          zip_file='sipprverse_output.zip',
                                                                          report_file='mlst.csv',
                                                                          ref_csv='/mnt/nas2/redmine/applications/OLCRedmineAutomator/tests/ref_csvs/mlst.csv')  # Hardcoded paths stink, fix me!
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False

            # STARAMR #####
            elif issue_subject == 'staramr' and issues_validated[issue_subject] == 'Unknown':
                if issue.status.id == 4:
                    issues_validated[issue_subject] = validate_csv_in_zip(issue=issue,
                                                                          zip_file='staramr_output.zip',
                                                                          report_file='salmonella/summary.tsv',
                                                                          ref_csv='/mnt/nas2/redmine/applications/OLCRedmineAutomator/tests/ref_csvs/summary.tsv')  # Hardcoded paths stink, fix me!
                    logging.info('{} is complete, status is {}'.format(issue_subject, issues_validated[issue_subject]))
                else:
                    all_complete = False
        total_time_taken += increment
        time.sleep(increment)

    return issues_validated


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sentry_dsn',
                        type=str,
                        required=True,
                        help='DSN for sentry. Can be made if you sign up at sentry.io and create a python project.')
    parser.add_argument('-k', '--api_key',
                        type=str,
                        required=True,
                        help='API Key for Redmine dev environment.')
    parser.add_argument('-l', '--long_test',
                        default=False,
                        action='store_true',
                        help='Activate this flag to do a longer test that includes automators that take a few '
                             'hours to run (NOT ACTUALLY IMPLEMENTED YET).')
    args = parser.parse_args()
    sentry_dsn = args.sentry_dsn
    api_key = args.api_key
    sentry_sdk.init(sentry_dsn)
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    redmine = Redmine(DEV_REDMINE_URL,
                      key=api_key,
                      requests={
                          'verify': False,
                          'timeout': 10,
                      })

    """
    Need to: create an issue for each of the automators we actually care about testing, and then check on outputs
    for each of them after a reasonable timeframe.
    Automators to test - more to be added as we go: 
    AMRsummary
    AutoCLARK
    eCGF
    ECTyper
    External Retrieve
    DiversiTree
    GeneSeekr
    PointFinder
    Prokka
    ResFinder
    sipprverse
    staramr
    
    Automators that take a long time (roughly 1 hour or more) to run and will only go when the --long
    flag is passed as an argument:
    cowsnphr_src
    snvphyl
    
    """

    issues_created = create_test_issues(redmine)

    issues_validated = monitor_issues(redmine=redmine,
                                      issue_dict=issues_created,
                                      timeout=3600)
    for issue_subject in issues_validated:
        if issues_validated[issue_subject] != 'Validated':
            sentry_sdk.capture_message('Redmine validation failing - automator {} with id '
                                       '{} had status {}'.format(issue_subject,
                                                                 issues_created[issue_subject],
                                                                 issues_validated[issue_subject]))
        print('{}: {}'.format(issue_subject, issues_validated[issue_subject]))

