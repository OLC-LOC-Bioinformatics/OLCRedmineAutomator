import os
import csv
import glob
import gzip
import shutil
import zipfile
from iridasequencefile import SequencePair


class MassExtractor(object):

    def __init__(self, nas_mnt):
        self.missing = list()
        self.low_quality = list()
        self.nas_mnt = nas_mnt
        self.seqid_rows = list()
        self.generic_sample_sheet_path = ""
        self.seqid_mounted_path = ""

    def move_files(self, sequences_info, outputfolder):
        """
        The main method that processes the inputs from Redmine and completes the task of moving the files to the mounted
        drive by running a series of tasks
        :param sequences_info: List of inputted Sequence_Info objects that will be used inplace of a text line
        :param outputfolder: The folder path on the drive where the folder from the nas must be moved to
        :return: Any files that could not be found on the nas 
        """

        if sequences_info is None:
            raise ValueError('No input files were found.')

        # The SEQ-IDs need to be places in a specific directory for the irida uploader to work
        self.seqid_mounted_path = os.path.join(outputfolder, "Data/Intensities/BaseCalls")
        # Create the directory if it does not exist
        if not os.path.exists(self.seqid_mounted_path):
            os.makedirs(self.seqid_mounted_path)
        completed_counter = 0
        for sequence in sequences_info:
            completed_counter += 1
            print('Moving {} of {} sets of files to external drive.'.format(completed_counter, len(sequences_info)))
            # Create a SequencePair object that will store all relevant paths
            sequence_pair = SequencePair(sequence_info=sequence)
            # Set the path to check depending on the SEQ-ID abbreviation used to speed up the process
            if 'MER' in sequence.sample_name:
                path_to_check = os.path.join('/mnt/nas2/raw_sequence_data/merged_sequences/*.fastq.gz')
            else:
                path_to_check = os.path.join('/mnt/nas2/raw_sequence_data/*/*/*.fastq.gz')
            # else:
            #     path_to_check = os.path.join(self.nas_mnt, 'External_MiSeq_Backup', '*', '*', '*.fastq.gz')

            for path in glob.iglob(path_to_check):
                if sequence.sample_name in path:

                    # Add a the path of the sequence to the object
                    sequence_pair.add_nas_seqid_path(seqid_path=path)

                    # If no SampleSheet associated with the fastq pair then add one
                    if sequence_pair.nas_sample_sheet_path is None:
                        sequence_pair.add_sample_sheet(path.split(sequence.sample_name)[0])
                    # If both of the pair are found then exit the search
                    if sequence_pair.both_exist:
                        break
            # Move the files for the sequence pair to the drive and add their SampleSheet data to list
            self.mount_seqid_files(sequence_pair)
            self.add_seqid_csv_data(sequence_pair)
        # Mount the generic SampleSheet and then append it will the new information
        print(self.seqid_rows)
        self.mount_generic_samplesheet(outputfolder)
        self.append_generic_csv(self.generic_sample_sheet_path)

        return self.missing, self.low_quality

    def add_seqid_csv_data(self, sequence_pair):
        """
        For each pair of sequence pair, add their data into the generic SampleSheet that will be moved to the drive 
        and used for the upload. For regular pairs, open their SampleSheet row that is logged on the nas and change it
        according to the inputted parameters. Next append it to the list that will be used to populate the 
        generic sheet.
        For MER sequence a custom row must be created from scratch since there are no SampleSheets for them on the nas
        :param sequence_pair: All information for a pair of SEQ-IDs
        """

        if "MER" in sequence_pair.seqid_info.sample_name:
            # Use custom parameters to be put in the SampleSheet for MER sequences
            self.seqid_rows.append(self.get_default_merge_sequence_row(sequence_pair))
        else:
            nas_csv_samplesheet = sequence_pair.nas_sample_sheet_path
            delimiter = ','
            # Open the SampleSheet for the sequence on the nas
            with open(nas_csv_samplesheet, 'r') as input_file:
                reader = csv.reader(input_file, delimiter=delimiter)
                for row in reader:
                    if len(row) > 8:  # incase of improper formatted row
                        # Search through the rows until the sample name is found in the first row of the SampleSheet
                        if sequence_pair.seqid_info.sample_name in row[0]:
                            row[0] = sequence_pair.seqid_info.sample_id   # Change the Sample_Name in the csv to the Sample ID
                            row[1] = sequence_pair.seqid_info.sample_id  # Change the Sample_ID in the csv to the input Sample ID
                            row[8] = sequence_pair.seqid_info.sample_project  # Change Sample_Project in the csv to the input
                            row[9] = sequence_pair.seqid_info.sample_name  # Change the description in the csv to the Sample Name

                            # If the length of the row is longer than the template, delete the extra columns
                            if len(row) > 10:
                                i = 10 - len(row)
                                del row[i:]
                            # Add the row to the list - if the sequences had high enough average quality.
                            if sequence_pair.seqid_info.sample_id not in self.low_quality:
                                self.seqid_rows.append(row)
                            break

    @ staticmethod
    def get_default_merge_sequence_row(sequence_pair):
        """
        Return the default row of data to be inputted into the data sheet for all merge type sequences 
        :param sequence_pair: All information for a pair of SEQ-IDs
        """
        return [sequence_pair.seqid_info.sample_id,  # Sample ID
                      sequence_pair.seqid_info.sample_id,  # Sample Name
                      "",  # Sample Plate
                      "",  # Sample Well
                      "na",  # I7 Index ID
                      "na",  # index
                      "na",  # I5 Index ID
                      "na",  # index2
                      sequence_pair.seqid_info.sample_project,  # Sample Project
                      sequence_pair.seqid_info.sample_name]  # Description

    def append_generic_csv(self, sample_sheet_path):
        """
        Add the rows from the different SampleSheets on the nas to the generic SampleSheet
        :param sample_sheet_path: Path of the default SampleSheet in the local folder
        """
        delimiter = ','
        with open(sample_sheet_path, 'a') as output_file:
            append = csv.writer(output_file, delimiter=delimiter)
            for row in self.seqid_rows:
                append.writerow(row)

    def mount_seqid_files(self, sequence_pair):
        """
        Put the sequence pair referenced in the sequence_files onto the drive - with specific path 
        "drivepath/Data/Intensities/BaseCalls"
        :param sequence_pair: All information for a pair of SEQ-IDs
        """
        # Figure out if which file is forward/reverse reads.
        if len(sequence_pair.seqid_paths) != 2:
            self.missing.append(sequence_pair.seqid_info.sample_name)
            return
        if '_R1' in sequence_pair.seqid_paths[0]:
            forward_reads = sequence_pair.seqid_paths[0]
            reverse_reads = sequence_pair.seqid_paths[1]
        else:
            forward_reads = sequence_pair.seqid_paths[1]
            reverse_reads = sequence_pair.seqid_paths[0]

        forward_out = os.path.join(self.seqid_mounted_path.replace(' ', '\\ '), sequence_pair.seqid_info.sample_id + '_S1_L001_R1_001.fastq.gz')
        reverse_out = os.path.join(self.seqid_mounted_path.replace(' ', '\\ '), sequence_pair.seqid_info.sample_id + '_S1_L001_R2_001.fastq.gz')
        if not os.path.isfile(forward_out):
            knowndepth = check_depth(forward_reads)
            if knowndepth > 200:
                # Check genome size - used to downsample extremely high coverage stuff to 200x coverage.
                genome_size = check_genome_size(forward_reads, reverse_reads)
                # Now call reformat.sh, set samplebasestarget to 200X coverage. If reads have less than that, this just
                # acts as a copy, otherwise, will downsample to 200X.
                # Very occasionally genome_size will get set to None. If that's the case, just assume a 5MB genome.
                if genome_size is None:
                    genome_size = 5000000
                samplereadstarget = int(float(genome_size) * (1.0 / 3.0))
                sqtk = os.path.relpath("/mnt/nas/users/julie/seqtk/seqtk")
                cmd = '{sqtk} sample -s100 {forward_reads} {samplereadstarget} | gzip > {forward_out} && {sqtk} sample -s100 {reverse_reads} {samplereadstarget} | gzip > {reverse_out}'.format(sqtk=sqtk, forward_reads=forward_reads, reverse_reads=reverse_reads, forward_out=forward_out, reverse_out=reverse_out, samplereadstarget=samplereadstarget)
            else:
                cmd = 'cp {forward_reads} {forward_out} && cp {reverse_reads} {reverse_out}'.format(forward_reads=forward_reads, forward_out=forward_out, reverse_reads=reverse_reads, reverse_out=reverse_out)
            print(cmd)
            os.system(cmd)
        average_forward_score, average_reverse_score = find_average_qscore(forward_reads, reverse_reads)
        if average_forward_score < 30.0:
            self.low_quality.append(sequence_pair.seqid_info.sample_id)
            # os.remove(forward_out)
            # os.remove(reverse_out)

    def mount_generic_samplesheet(self, outputfolder):
        """
        Put the generic SampleSheet on the drive in the root folder
        :param outputfolder: Root folder on the drive where the Redmine request is to be placed
        """
        # Create the directory and the make path to place the SampleSample sheet on the drive
        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)
        self.generic_sample_sheet_path = os.path.join(outputfolder, 'SampleSheet.csv')

        # Re-create the path for the local generic SampleSheet
        local_dir_path = os.path.dirname(os.path.realpath(__file__))
        local_generic_samplesheet_path = os.path.join(local_dir_path, 'iridadummysheet.csv')

        # Copy the local SampleSheet onto the drive in the new location
        shutil.copy(local_generic_samplesheet_path, self.generic_sample_sheet_path)


def check_depth(forward_reads):
    if forward_reads.startswith("/mnt/nas2/raw_sequence_data/miseq/"):
        pathcomp = forward_reads.split("/")
        metafile = "/mnt/nas2/processed_sequence_data/miseq_assemblies/" + pathcomp[5] + "/reports/combinedMetadata.csv"
        if os.path.isfile(metafile):
            samplename = pathcomp[6].split("_", 1)[0]
            with open(metafile) as f:
                lines = f.readlines()
            try:
                cols = lines[0].split(",")
                colofinterest = [c for c in range(len(cols) - 1) if cols[c] == "AverageCoverageDepth"][0]
                for line in lines:
                    if line.startswith(samplename + ","):
                        return float(line.split(",")[colofinterest])
                        break
            except:
                return 9001
        else:
            return 9002
    elif forward_reads.startswith("/mnt/nas2/raw_sequence_data/merged_sequences"):
        samplename = forward_reads.split("/")[5].split("_", 1)[0]
        files_to_check = os.path.join('/mnt/nas2/processed_sequence_data/merged_assemblies/*/reports/combinedMetadata.csv*')
        for metafile in glob.iglob(files_to_check):
            if os.path.isfile(metafile):
                with open(metafile) as f:
                    lines = f.readlines()
                try:
                    cols = lines[0].split(",")
                    colofinterest = [c for c in range(len(cols) - 1) if cols[c] == "AverageCoverageDepth"][0]
                    for line in lines:
                        if line.startswith(samplename + ","):
                            return float(line.split(",")[colofinterest])
                            break
                except:
                    return 9001
    else:
        return 9003 



def check_genome_size(forward_reads, reverse_reads):
    genome_size = None
    if forward_reads.startswith("/mnt/nas2/raw_sequence_data/miseq/"):
        pathcomp = forward_reads.split("/")
        metafile = "/mnt/nas2/processed_sequence_data/miseq_assemblies/" + pathcomp[5] + "/reports/combinedMetadata.csv"
        if os.path.isfile(metafile):
            try:
                samplename = pathcomp[6].split("_", 1)[0]
                with open(metafile) as f:
                   lines = f.readlines()
                cols = lines[0].split(",")
                colofinterest = [c for c in range(len(cols) - 1) if cols[c] == "TotalLength"][0]
                for line in lines:
                    if line.startswith(samplename + ","):
                        genome_size = int(line.split(",")[colofinterest])
                        if genome_size > 1000000:
                            return genome_size
                            break
            except:
                print("Something went wrong with finding a genome size.")
                return None
    elif forward_reads.startswith("/mnt/nas2/raw_sequence_data/merged_sequences/"):
            files_to_check = os.path.join('/mnt/nas2/processed_sequence_data/merged_assemblies/*/reports/combinedMetadata.csv*')
            samplename = forward_reads.split("/")[5].split("_", 1)[0]
            for metafile in glob.iglob(files_to_check):
                if os.path.isfile(metafile):
                    try:
                        with open(metafile) as f:
                            lines = f.readlines()
                        cols = lines[0].split(",")
                        colofinterest = [c for c in range(len(cols) - 1) if cols[c] == "TotalLength"][0]
                        for line in lines:
                            if line.startswith(samplename + ","):
                                genome_size = int(line.split(",")[colofinterest])
                                if genome_size > 1000000:
                                    return genome_size
                                    break
                    except:
                        break
 

    # if genome_size doesn't pass sanity check, just return None.
    print("Could not find genome size.")
    return None


# Heng Li readfq - super fast!
def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs);  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break

def qscore_from_fastqc(fqcpath, fqc):
    reading=False
    qsum=0
    count=0
    with zipfile.ZipFile(fqcpath) as z:
        with z.open(fqc) as f:
            for l in f:
                line = l.decode("UTF-8")
                if reading:
                    if line.startswith(">"):
                        return qsum / count
                    else:
                        splitline = line.split("\t")
                        qsum += float(splitline[0]) * float(splitline[1])
                        count += float(splitline[1])
                elif line == "#Quality\tCount\n":
                    reading=True 


def find_average_qscore(forward_reads, reverse_reads):
    # these methods seem to get slightly different results, but
    # I'm not super confident that the old way's better, and the new way's so much faster...
    try:
        pathcomp = forward_reads.split("/")
        fqdir = "/mnt/nas2/processed_sequence_data/miseq_assemblies/" + pathcomp[5] + "/" + pathcomp[6].split("_")[0] + "/fastqc/Raw"
        n = pathcomp[6].split(".")[0]
        average_forward_score = qscore_from_fastqc(fqdir + "/" + n + "_fastqc.zip", n + "_fastqc/fastqc_data.txt")
        n = reverse_reads.split("/")[6].split(".")[0]
        average_reverse_score = qscore_from_fastqc(fqdir + "/" + n + "_fastqc.zip", n + "_fastqc/fastqc_data.txt")
    except:
        forward_read_qualities = list()
        reverse_read_qualities = list()

        if forward_reads.endswith('.gz'):
            for name, seq, qual in readfq(gzip.open(forward_reads, 'rt')):
                for q in qual:
                    forward_read_qualities.append((ord(q) - 33))
        else:
            for name, seq, qual in readfq(open(forward_reads)):
                for q in qual:
                    forward_read_qualities.append((ord(q) - 33))

        if reverse_reads.endswith('.gz'):
            for name, seq, qual in readfq(gzip.open(reverse_reads, 'rt')):
                for q in qual:
                    reverse_read_qualities.append((ord(q) - 33))
        else:
            for name, seq, qual in readfq(open(reverse_reads)):
                for q in qual:
                    reverse_read_qualities.append((ord(q) - 33))

        average_forward_score = sum(forward_read_qualities)/len(forward_read_qualities)
        average_reverse_score = sum(reverse_read_qualities)/len(reverse_read_qualities)

    return average_forward_score, average_reverse_score
