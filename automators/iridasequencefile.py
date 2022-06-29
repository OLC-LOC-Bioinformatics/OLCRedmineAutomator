import os


class SequencePair:

    csv_sheet = 'SampleSheet.csv'
    bak_sheet = 'SampleSheet.bak2'

    def __init__(self, sequence_info):

        self.seqid_paths = list()  # List has a maximum size of 2 SEQ-IDs (R1 & R2 pair)
        self.seqid_info = sequence_info  # Info given through the text document on Redmine
        self.nas_sample_sheet_path = None  # Path to the sample sheet on the nas
        self.both_exist = False  # Flag when a pair of SEQ-ID are found
        self.run_date = "01012000"  # Generic run date
        self.csv_file = True

        # File paths for the files to be put on the drive
        self.seqid_mounted_folder = None
        self.sample_sheet_mounted_path = None

    def add_nas_seqid_path(self, seqid_path):
        """
        This will add the path of the current SEQ-ID in the nas to the seqid_path list.
        This method ensures that no more than 2 sequence paths are ever added to the list and moved to the drive.
        If more than 2 files are found, throw an error since this should never happen and there is a major problem.
        :param seqid_path: Path in the nas where the matching SEQ-ID name was found 
        """

        if len(self.seqid_paths) is 1:
            self.both_exist = True  # both found, no need to continue searching
            self.seqid_paths.append(seqid_path)

        elif len(self.seqid_paths) is 0:
            self.seqid_paths.append(seqid_path)
        else:
             raise ImportWarning("More than 2 files found that have the SEQ-ID: %s", self.seqid_info.seq_id)

    def add_sample_sheet(self, nas_sheet_dir):
        """
        Add the path of the SampleSheet in the nas for the SEQ-ID pair
        :param nas_sheet_dir: the directory of the SampleSheet in the nas
        """
        # Create the paths for both a .csv and .bak SampleSheet
        csv_sheet_path = os.path.join(nas_sheet_dir, SequencePair.csv_sheet)
        bak_sheet_path = os.path.join(nas_sheet_dir, SequencePair.bak_sheet)

        # Check if the SampleSheet is a .csv file as is the default
        if os.path.exists(csv_sheet_path):
            self.nas_sample_sheet_path = csv_sheet_path

        # Check if the SampleSheet is a .bak file, this is uncommon but does occur
        # Set the csv_file flag to false so that it can be dealt with appropriately
        elif os.path.exists(bak_sheet_path):
            self.nas_sample_sheet_path = bak_sheet_path
            self.csv_file = False


class SequenceInfo(object):
    def __init__(self, text_line):
        """
        This separates a inputted line from the text file into an object, the line uses \t decimeters
        :param text_line: this should have the format 'SampleName   SampleId    SampleProject' 
        """
        input_list = text_line.split()
        self.sample_name = str(input_list[0]).rstrip()
        self.sample_id = str(input_list[1]).rstrip()
        self.sample_project = str(input_list[2]).rstrip()
