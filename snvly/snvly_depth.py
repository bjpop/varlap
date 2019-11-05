'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 16 Oct 2019 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import pysam
import csv
import os.path
import pathlib


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "snvly_depth"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Compute various bits of information about variants in one or more BAM files'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        'bams', nargs='+', metavar='BAM', type=str, help='Filepaths of BAM files')
    parser.add_argument(
        '--labels', required=False, nargs='*', metavar='LABEL', type=str, help='Labels for BAM files')
    parser.add_argument(
        '--chroms', required=True, nargs='+', metavar='CHROM', type=str, help='Chromosomes to consider')
    parser.add_argument(
        '--step', required=True, metavar='STEP', type=int, help='Step size of chromosome position')
    parser.add_argument(
        '--sample', default='', required=False, metavar='SAMPLE', type=str, help='Sample identifier')
    parser.add_argument(
        '--noheader', action='store_true', help='Suppress output header row')
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    return parser.parse_args()



def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))



def write_header(options, bam_labels):
    header_general = ["chrom", "pos", "pos normalised", "sample"]
    bam_headers = [label for label in bam_labels]
    if not options.noheader:
        header_regions = []
        header = header_general + bam_headers 
        print(",".join(header))


# give each input BAM file a label, either from the labels given
# on the command line, or, if there are too few of those, from
# the name of the file
def get_bam_labels(labels, bam_files):
    result = []
    for index, bam_filename in enumerate(bam_files):
        if labels is not None and index <= len(labels):
            result.append(labels[index])
        else:
            name = pathlib.PurePath(bam_filename).name
            result.append(name)
    return result


def get_chrom_size(readers, chrom):
    if len(readers) > 0:
        first_reader = readers[0]
        return first_reader.get_reference_length(chrom)
    else:
        return None


def process_bams(options):
    bam_labels = get_bam_labels(options.labels, options.bams)
    bam_readers = [BamReader(filepath) for filepath in options.bams]
    write_header(options, bam_labels)
    for chrom in options.chroms:
        chrom_size = get_chrom_size(bam_readers, chrom)
        if chrom_size is not None:
            for pos in range(1, chrom_size, options.step):
                pos_normalised = pos / chrom_size
                bams_depths = [reader.depth(chrom, pos) for reader in bam_readers]
                write_output_row(chrom, pos, pos_normalised, options.sample, bams_depths)
    for reader in bam_readers:
        reader.close()


def write_output_row(chrom, pos, pos_normalised, sample, bams_depths):
    row_bams = [str(depth) for depth in bams_depths] 
    row = [chrom, str(pos), str(pos_normalised), sample] + row_bams 
    print(",".join(row))
 

MAX_PILEUP_DEPTH = 1000000000

class BamReader(object):
    def __init__(self, filepath):
        logging.info(f"Reading BAM file: {filepath}")
        resolved_path = os.path.realpath(filepath)
        logging.info(f"Resolved filepath to: {resolved_path}")
        self.samfile = pysam.AlignmentFile(resolved_path, "rb")

    # pos is expected to be 1-based 
    def depth(self, chrom, pos):
        zero_based_pos = pos - 1
        num_reads = 0
        for pileupcolumn in self.samfile.pileup(chrom, zero_based_pos, zero_based_pos+1,
                                               truncate=True, stepper='samtools',
                                               ignore_overlaps=False, ignore_orphans=True,
                                               max_depth=MAX_PILEUP_DEPTH): 
            for read in pileupcolumn.pileups:
                if not read.is_del and not read.is_refskip:
                    num_reads += 1
        return num_reads 

    def get_reference_length(self, chrom):
        return self.samfile.get_reference_length(chrom)

    def close(self):
        self.samfile.close()


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    process_bams(options)
    logging.info("Completed")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
