# begin-doc-include

import os
import re
import click
import pickle
import pylatex as pl
import autoroga_extract_report_data as extract_report_data

from datetime import datetime
from pylatex.utils import bold, italic
from autoroga_database import update_db

"""
This script is as per request from Ray Allain, who always types autotoga in the subject instead of autoroga...

"""

@click.command()
@click.option('--redmine_instance', help='Path to pickled Redmine API instance')
@click.option('--issue', help='Path to pickled Redmine issue')
@click.option('--work_dir', help='Path to Redmine issue work directory')
@click.option('--description', help='Path to pickled Redmine description')
def redmine_toga(redmine_instance, issue, work_dir, description):
    """
    Main method for generating autoTOGA
    """
    # Unpickle Redmine objects
    redmine_instance = pickle.load(open(redmine_instance, 'rb'))
    issue = pickle.load(open(issue, 'rb'))
    description = pickle.load(open(description, 'rb'))

    try:
        # Parse description to figure out what SEQIDs we need to run on.
#        seqids = list()
#        for item in description:
#            item = item.upper().rstrip()
            # Otherwise the item should be a SEQID
#            seqids.append(item)

        redmine_instance.issue.update(resource_id=issue.id, status_id=4,
                                      notes='      ////\\\\\n'
                                            '      |      |\n'
                                            '     @  O  O  @\n'
                                            '      |  ~   |         \__\n'
                                            '       \ -- /          |\ |\n'
                                            '     ___|  |___        | \|\n'
                                            '    /      /  /\      /|__|\n'
                                            '   /   ___/  /  \    / /\n'
                                            '  /  /|      |\  \  / /\n'
                                            ' /  / |      | \  \/ /\n'
                                            '<  <  |      |  \   /\n'
                                            ' \  \ |      |   \_/\n'
                                            '  \  \|------|\n'
                                            '    \_|------|\n'
                                            '      |      |\n'
                                            '      |      |\n'
                                            '      |      |\n'
                                            '      |      |\n'
                                            '      |      |\n'
                                            '      |      |\n'
                                            '      |______|\n'
                                            '     _|  |  |\n'
                                            ' cccC_Cccc___)\n'
                                            '\n'
                                            'Did you actually intend to request an auto ROGA? Please re-submit with the correct subject.')


    except Exception as e:
        sentry_sdk.capture_exception(e)
        redmine_instance.issue.update(resource_id=issue.id,
                                      notes='Something went wrong! Please contact a bioinformatician '
                                            'to investigate. (Did you mean to submit an auto ROGA instead?)')

if __name__ == '__main__':
    redmine_toga()
