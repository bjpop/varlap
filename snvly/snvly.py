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
from intervaltree import Interval, IntervalTree
import os.path


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
    description = 'Compute various bits of information about somatic variants in tumour and normal BAM files'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        '--tumour', metavar='BAM', type=str, help='Filepath of tumour BAM file')
    parser.add_argument(
        '--normal', metavar='BAM', type=str, help='Filepath of normal BAM file')
    parser.add_argument(
        '--sample', required=True, metavar='SAMPLE', type=str, help='Sample identifier')
    parser.add_argument(
        '--regions', required=False, metavar='REGIONS', type=str, help='Filepath of genomic regions-of-interest in BED format')
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
    return sorted(result)



def make_output_headers(regions):
    header_general = ["chrom", "pos", "ref", "alt", "sample"]
    header_bam = ["depth", "A", "T", "G", "C", "N", "ref count", "alt count", "alt VAF", "avg NM", "avg INDEL", "avg clipping", "avg base qual", "avg map qual", "avg align len"]
    header_tumour = ["tumour " + h for h in header_bam]
    header_normal = ["normal " + h for h in header_bam]
    header_regions = []
    if regions is not None:
        header_regions = sorted(regions.region_labels)
    return header_general, header_regions, header_bam, header_tumour,  header_normal


def get_variant_region_intersections(regions, chrom, pos):
    result = {}
    if regions is not None:
        for label in regions.region_labels:
            result[label] = False 
        for _start, _end, overlapped_label in regions.get_overlapping_regions(chrom, pos, pos+1):
            result[overlapped_label] = True 
    return result
 

def process_variants_bams(options, regions, variants):
    tumour_data = process_bam(options.tumour, variants)
    normal_data = process_bam(options.normal, variants)
    header_general, header_regions, header_bam, header_tumour, header_normal = make_output_headers(regions)
    output_header = header_general + header_regions + header_tumour + header_normal
    writer = csv.DictWriter(sys.stdout, fieldnames=output_header)
    if not options.noheader:
        writer.writeheader()
    for variant in variants:
        chrom, pos, ref, alt = variant
        variant_region_counts = get_variant_region_intersections(regions, chrom, pos)
        this_tumour = tumour_data[variant]
        this_normal = normal_data[variant]
        row_general = [("chrom", chrom), ("pos", pos), ("ref", ref), ("alt", alt), ("sample", options.sample)]
        row_regions = list(variant_region_counts.items())
        row_tumour = [("tumour " + h, this_tumour[h]) for h in header_bam]
        row_normal = [("normal " + h, this_normal[h]) for h in header_bam] 
        row = dict(row_general + row_regions + row_tumour + row_normal)
        writer.writerow(row)
 

VALID_DNA_BASES = "ATGCN"

class Features(object):
    def __init__(self):
        self.nm_sum = 0
        self.base_qual_sum = 0
        self.map_qual_sum = 0
        self.align_len_sum = 0
        self.clipping_sum = 0
        self.indel_sum = 0
        self.base_counts = { base: 0 for base in VALID_DNA_BASES }
        

    # See: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
    def alignment_features(self, alignment):
        self.map_qual_sum += alignment.mapping_quality 
        self.align_len_sum += alignment.query_alignment_length
        cigar_stats = alignment.get_cigar_stats()[0]
        if len(cigar_stats) == 11:
            self.nm_sum += cigar_stats[10]
            # count soft and hard clips together
            self.clipping_sum += cigar_stats[4] + cigar_stats[5]
            self.indel_sum += cigar_stats[1] + cigar_stats[2]


    def base_features(self, read):
        base = read.alignment.query_sequence[read.query_position].upper()
        if base in VALID_DNA_BASES:
            self.base_counts[base] += 1
        else:
            exit(f"Unrecognised base: {base}")
        self.base_qual_sum += read.alignment.query_qualities[read.query_position]




def process_bam(filepath, variants):
    result = {}
    logging.info(f"Reading BAM file: {filepath}")
    resolved_path = os.path.realpath(filepath)
    logging.info(f"Resolved filepath to: {resolved_path}")
    samfile = pysam.AlignmentFile(resolved_path, "rb")
    for (chrom, pos, ref, alt) in variants:
        zero_based_pos = pos - 1
        features = Features()
        num_considered_reads = 0
        for pileupcolumn in samfile.pileup(chrom, zero_based_pos, zero_based_pos+1,
                truncate=True, stepper='samtools',
                ignore_overlaps=False, ignore_orphans=True,
                max_depth=1000000000): 
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    this_alignment = pileupread.alignment
                    num_considered_reads += 1
                    features.alignment_features(this_alignment)
                    features.base_features(pileupread)
        alt_count = features.base_counts[alt]
        ref_count = features.base_counts[ref]
        if num_considered_reads > 0:
            average_nm = features.nm_sum / num_considered_reads
            average_base_qual = features.base_qual_sum / num_considered_reads
            average_map_qual = features.map_qual_sum / num_considered_reads
            average_align_len = features.align_len_sum / num_considered_reads
            alt_vaf = alt_count / num_considered_reads 
            average_clipping = features.clipping_sum / num_considered_reads
            average_indels = features.indel_sum / num_considered_reads
        else:
            average_nm = average_base_qual = average_map_qual = average_align_len = alt_vaf = average_clipping = average_indels = ''
        result[(chrom, pos, ref, alt)] = {
            "depth": num_considered_reads,
            "A": features.base_counts['A'],
            "T": features.base_counts['T'], 
            "G": features.base_counts['G'], 
            "C": features.base_counts['C'], 
            "N": features.base_counts['N'], 
            "ref count": ref_count,
            "alt count": alt_count,
            "alt VAF": alt_vaf,
            "avg NM": average_nm,
            "avg base qual": average_base_qual,
            "avg map qual": average_map_qual,
            "avg align len": average_align_len,
            "avg clipping": average_clipping,
            "avg INDEL": average_indels }
    return result


class Regions(object):
    def __init__(self):
        self.chrom_regions = {}
        self.region_labels = set()

    def add_region(self, chrom, one_based_inclusive_start, one_based_exclusive_end, label):
        if chrom not in self.chrom_regions:
            self.chrom_regions[chrom] = IntervalTree()
        self.chrom_regions[chrom][one_based_inclusive_start:one_based_exclusive_end] = label
        self.region_labels.add(label)

    def get_overlapping_regions(self, chrom, one_based_inclusive_start, one_based_exclusive_end):
        if chrom in self.chrom_regions:
            this_regions = self.chrom_regions[chrom]
            for interval in this_regions[one_based_inclusive_start: one_based_exclusive_end]:
                start, end, label = interval
                yield start, end, label

'''
Regions are represented in BED format as:
chrom	start	end 	name

where:
start and end are zero based, and the region is defined by the half-closed interval:
[start, end)

The result is a dictionary with chromosomes as keys, and interval trees as values, and
a set of all the region names from the file.
'''
NUM_FEATURE_FIELDS = 4

def get_regions(filepath):
    result = None
    if filepath is not None:
        result = Regions()
        logging.info(f"Reading regions from {filepath}")
        num_regions = 0
        with open(filepath) as file:
            for row in file:
                if row.startswith("#"):
                    logging.info(f"Skipping likely header row from regions file: {row.strip()}")
                    continue
                else:
                    fields = row.strip().split()
                    if len(fields) >= NUM_FEATURE_FIELDS:
                        try:
                            chrom, start, end, label = fields[:NUM_FEATURE_FIELDS]
                            one_based_start = int(start) + 1
                            one_based_end = int(end) + 1
                        except:
                            logging.info(f"Skipping likely malformed (could not parse) region row: {row}")
                        else:
                            result.add_region(chrom, one_based_start, one_based_end, label)
                            num_regions += 1
                    else:
                        logging.info(f"Skipping likely malformed (too short) region row: {row}")
        num_labels = len(result.region_labels)
        logging.info(f"Read {num_regions} regions with {num_labels} labels")
    return result


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    # variants are sorted
    variants = get_variants()
    regions = get_regions(options.regions)
    process_variants_bams(options, regions, variants)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
