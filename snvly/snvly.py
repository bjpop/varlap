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


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "snvly"


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
    description = 'Read one or more FASTA files, compute simple stats for each file'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        '--tumour', metavar='BAM', type=str, help='Filepath of tumour BAM file')
    parser.add_argument(
        '--normal', metavar='BAM', type=str, help='Filepath of normal BAM file')
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


class Counts(object):
    def __init__(self):
        self.A = 0
        self.T = 0
        self.G = 0
        self.C = 0

    def increment_base_count(self, base):
        if base == 'A':
            self.A += 1
        elif base == 'T':
            self.T += 1
        elif base == 'G':
            self.G += 1
        elif base == 'C':
            self.C += 1
        else:
            exit(f"Unrecognised base: {base}")

    def get_count(self, base):
        if base == 'A':
            return self.A
        elif base == 'T':
            return self.T
        elif base == 'G':
            return self.G
        elif base == 'C':
            return self.C

    def __str__(self):
        return f"A:{self.A}, T:{self.T}, G:{self.G}, C:{self.C}"


def get_variants():
    result = set()
    for line in sys.stdin:
        if line.startswith('#'):
            continue
        fields = line.strip().split()
        if len(fields) >= 5:
            chrom, pos, _id, ref, alt = fields[:5]
            pos = int(pos)
            result.add((chrom, pos, ref, alt))
    return result

header = ["chrom", "pos", "ref", "alt",
           "tumour depth", "tumour A", "tumour T", "tumour G", "tumour C", "tumour ref count", "tumour alt count", "tumour avg nm", "tumour avg base qual", "tumour avg map qual", "tumour avg align len",
           "normal depth", "normal A", "normal T", "normal G", "normal C", "normal ref count", "normal alt count", "normal avg nm", "normal avg base qual", "normal avg map qual", "normal avg align len"]

          

def process_variants_bams(variants, tumour_filepath, normal_filepath):
    tumour_data = process_bam(variants, tumour_filepath)
    normal_data = process_bam(variants, normal_filepath)
    writer = csv.DictWriter(sys.stdout, fieldnames=header)
    writer.writeheader()
    for variant in sorted(variants):
        chrom, pos, ref, alt = variant
        this_tumour = tumour_data[variant]
        this_normal = normal_data[variant]
        row = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
               "tumour depth": this_tumour["depth"], "tumour A": this_tumour["A"],
               "tumour T": this_tumour["T"], "tumour G": this_tumour["G"],
               "tumour C": this_tumour["C"], "tumour ref count": this_tumour["ref count"],
               "tumour alt count": this_tumour["alt count"],
               "tumour avg nm": this_tumour["avg nm"], "tumour avg base qual": this_tumour["avg base qual"],
               "tumour avg map qual": this_tumour["avg map qual"], "tumour avg align len": this_tumour["avg align len"],
               "normal depth": this_normal["depth"], "normal A": this_normal["A"],
               "normal T": this_normal["T"], "normal G": this_normal["G"],
               "normal C": this_normal["C"], "normal ref count": this_normal["ref count"],
               "normal alt count": this_normal["alt count"],
               "normal avg nm": this_normal["avg nm"], "normal avg base qual": this_normal["avg base qual"],
               "normal avg map qual": this_normal["avg map qual"], "normal avg align len": this_normal["avg align len"]}
        writer.writerow(row)
 

VALID_DNA_BASES = "ATGC"

def process_bam(variants, filepath):
    result = {}
    samfile = pysam.AlignmentFile(filepath, "rb")
    for (chrom, pos, ref, alt) in variants:
        zero_based_pos = pos - 1
        counts = Counts()
        coverage = 0
        num_considered_reads = 0
        this_nm = 0
        this_base_qual = 0
        this_map_qual = 0
        this_align_len = 0
        for pileupcolumn in samfile.pileup(chrom, zero_based_pos, zero_based_pos+1, truncate=True, stepper='samtools',
                                           ignore_overlaps=False, ignore_orphans=True,
                                           max_depth=1000000000): 
            coverage = pileupcolumn.nsegments
            for pileupread in pileupcolumn.pileups:
                this_alignment = pileupread.alignment
                mapping_quality = this_alignment.mapping_quality 
                alignment_length = this_alignment.query_alignment_length
                nm_count = 0
                cigar_stats = this_alignment.get_cigar_stats()[0]
                if len(cigar_stats) == 11:
                    nm_count = cigar_stats[10]
                if not pileupread.is_del and not pileupread.is_refskip:
                    this_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                    base_qual = pileupread.alignment.query_qualities[pileupread.query_position]
                    if this_base in VALID_DNA_BASES:
                        num_considered_reads += 1
                        this_nm += nm_count
                        this_base_qual += base_qual
                        this_map_qual += mapping_quality 
                        this_align_len += alignment_length
                        counts.increment_base_count(this_base)
        alt_count = counts.get_count(alt)
        ref_count = counts.get_count(ref)
        average_nm = this_nm / num_considered_reads
        average_base_qual = this_base_qual / num_considered_reads
        average_map_qual = this_map_qual / num_considered_reads
        average_align_len = this_align_len / num_considered_reads
        result[(chrom, pos, ref, alt)] = {"depth": coverage, "A": counts.A, "T": counts.T,
                                          "G": counts.G, "C": counts.C, "ref count": ref_count,
                                          "alt count": alt_count,
                                          "avg nm": average_nm,
                                          "avg base qual": average_base_qual,
                                          "avg map qual": average_map_qual,
                                          "avg align len": average_align_len}
    return result



def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    variants = get_variants()
    process_variants_bams(variants, options.tumour, options.normal)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
