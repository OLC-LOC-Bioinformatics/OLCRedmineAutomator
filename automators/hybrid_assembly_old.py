#!/usr/bin/env python

import os
import sys
import glob
import shutil
import logging
import argparse
import datetime
import subprocess
import multiprocessing


def run_cmd(cmd, logfile):
    if logfile is None:
        subprocess.call(cmd, shell=True)
    else:
        with open(logfile, 'a+') as f:
            f.write(cmd)
            subprocess.call(cmd, shell=True, stdout=f, stderr=f)


def trim_reads(forward_in, reverse_in, forward_trimmed, reverse_trimmed, threads=8, logfile=None):
    cmd = 'bbduk.sh in={forward_in} in2={reverse_in} out={forward_trimmed} out2={reverse_trimmed} ' \
          'ref=adapters trimq=10 qtrim=w minlength=50 threads={threads}'.format(forward_in=forward_in,
                                                                                reverse_in=reverse_in,
                                                                                forward_trimmed=forward_trimmed,
                                                                                reverse_trimmed=reverse_trimmed,
                                                                                threads=threads)
    run_cmd(cmd, logfile)


def correct_reads(forward_in, reverse_in, forward_corrected, reverse_corrected, threads=8, logfile=None):
    cmd = 'tadpole.sh in={forward_in} in2={reverse_in} out={forward_corrected} out2={reverse_corrected} ' \
          'mode=correct threads={threads}'.format(forward_in=forward_in,
                                                  reverse_in=reverse_in,
                                                  forward_corrected=forward_corrected,
                                                  reverse_corrected=reverse_corrected,
                                                  threads=threads)
    run_cmd(cmd, logfile)


def assemble(forward_illumina, reverse_illumina, pacbio, output_dir, threads=8, logfile=None):
    cmd = 'unicycler -1 {forward_illumina} -2 {reverse_illumina} -l {pacbio} -o {output_dir} ' \
          '--no_correct -t {threads}'.format(forward_illumina=forward_illumina,
                                             reverse_illumina=reverse_illumina,
                                             pacbio=pacbio,
                                             threads=threads,
                                             output_dir=output_dir)
    run_cmd(cmd, logfile)


def file_check(file_path):
    file_is_present = True
    if not os.path.isfile(file_path):
        logging.error('Your specified input file(s) could not be found! '
                      'Missing file(s) were: {}'.format(file_path))
        file_is_present = False
    return file_is_present


def main(forward_illumina, reverse_illumina, pacbio, output_dir, assembly_name, threads=8, logfile=None, keep_files=False):
    # Get logger setup.
    logging.basicConfig(format='\033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%H:%M:%S')
    # Start of by doing some checks that files specified actually exist. Boot user if
    file_presence = list()
    file_presence.append(file_check(forward_illumina))
    file_presence.append(file_check(reverse_illumina))
    file_presence.append(file_check(pacbio))
    if False in file_presence:
        sys.exit(1)
    # With input files verified, make the output directory that we want to use.
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    # Now trim and correct illumina reads before we pass them to unicycler so we don't have to use SPAdes error
    # correction, which is real slow.
    trimmed_forward = os.path.join(output_dir, assembly_name + '_trimmed_R1.fastq.gz')
    trimmed_reverse = os.path.join(output_dir, assembly_name + '_trimmed_R2.fastq.gz')
    corrected_forward = os.path.join(output_dir, assembly_name + '_corrected_R1.fastq.gz')
    corrected_reverse = os.path.join(output_dir, assembly_name + '_corrected_R2.fastq.gz')
    logging.info('Trimming Illumina reads...')
    trim_reads(forward_in=forward_illumina,
               reverse_in=reverse_illumina,
               forward_trimmed=trimmed_forward,
               reverse_trimmed=trimmed_reverse,
               threads=threads,
               logfile=logfile)
    logging.info('Correcting Illumina reads...')
    correct_reads(forward_in=trimmed_forward,
                  reverse_in=trimmed_reverse,
                  forward_corrected=corrected_forward,
                  reverse_corrected=corrected_reverse,
                  threads=threads,
                  logfile=logfile)
    # Once Illumina reads are ready to go, it's time to assemble!
    logging.info('Assembling using Unicycler. This will take quite a while...')
    assemble(forward_illumina=corrected_forward,
             reverse_illumina=corrected_reverse,
             pacbio=pacbio,
             output_dir=output_dir,
             threads=threads,
             logfile=logfile)

    # With assembly done, clean up all the files we don't really care about from Unicycler (so all except the
    # assembly FASTA, unless user said to keep them), and rename the assembly to something better.
    shutil.move(os.path.join(output_dir, 'assembly.fasta'), os.path.join(output_dir, assembly_name + '.fasta'))
    if not keep_files:
        logging.info('Cleaning up Unicycler output files...')
        files_to_remove = glob.glob(os.path.join(output_dir, '*.gfa'))
        files_to_remove += glob.glob(os.path.join(output_dir, '*.fastq.gz'))
        for item in files_to_remove:
            os.remove(item)
        shutil.rmtree(os.path.join(output_dir, 'pilon_polish'))
    logging.info('Hybrid Assembly complete. Assembly file is {}'.format(os.path.join(output_dir, assembly_name + '.fasta')))


if __name__ == '__main__':
    cpu_count = multiprocessing.cpu_count()
    parser = argparse.ArgumentParser(description='A wrapper for Unicycler that performs read trimming and correction'
                                                 ' on Illumina reads before using Unicycler. Only works for hybrid'
                                                 ' assemblies, and designed for PacBio reads. Will update for Nanopore'
                                                 ' if/when we actually get a MinION up and running.')
    parser.add_argument('-p', '--pacbio',
                        type=str,
                        required=True,
                        help='Path to PacBio reads - recommended to use the filtered subreads file.')
    parser.add_argument('-1', '--forward',
                        type=str,
                        required=True,
                        help='Path to Illumina forward reads.')
    parser.add_argument('-2', '--reverse',
                        type=str,
                        required=True,
                        help='Path to Illumina reverse reads.')
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        required=True,
                        help='Directory where you want to create your output. Will be created - DO NOT USE an existing'
                             ' directory, as files may be deleted from it.')
    parser.add_argument('-n', '--name',
                        type=str,
                        default='hybrid_assembly',
                        help='Name for your output assembly. Defaults to assembly.fasta - to change, use only the'
                             ' base filename, no extension (the .fasta is added automatically)')
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=cpu_count,
                        help='Number of threads to use. Defaults to all on your machine.')
    parser.add_argument('-l', '--logfile',
                        type=str,
                        default=datetime.datetime.now().strftime('%Y%m%d-%H%m%S') + '_log.txt',
                        help='Name of logfile to use. Defaults to YYYYMMDD-HHMMSS_log.txt in your current working'
                             ' directory.')
    parser.add_argument('-k', '--keep_files',
                        default=False,
                        action='store_true',
                        help='By default, all files created by Unicycle except the final assembly will be deleted.'
                             ' If for some reason you want to keep them, set this flag.')
    args = parser.parse_args()

    pacbio = args.pacbio
    forward_illumina = args.forward
    reverse_illumina = args.reverse
    output_dir = args.output_dir
    assembly_name = args.name
    threads = args.threads
    logfile = args.logfile
    keep_files = args.keep_files

    main(forward_illumina=forward_illumina,
         reverse_illumina=reverse_illumina,
         pacbio=pacbio,
         output_dir=output_dir,
         assembly_name=assembly_name,
         threads=threads,
         logfile=logfile,
         keep_files=keep_files)
