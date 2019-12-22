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
import pathlib


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "varlap"


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
        '--sample', default='', required=False, metavar='SAMPLE', type=str, help='Sample identifier')
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
            if ref in VALID_DNA_BASES and alt in VALID_DNA_BASES:
                pos = int(pos)
                result.add((chrom, pos, ref, alt))
            else:
                logging.warning(f"Skipping variant: {chrom}:{pos},{ref}>{alt}")
    logging.info(f"Read {len(result)} variants")
    return sorted(result)


def write_header(options, bam_labels, regions):
    header_general = ["chrom", "pos", "ref", "alt", "pos normalised", "sample"]
    bam_headers = [label + " " + field for label in bam_labels for field in LocusFeatures.fields]
    if not options.noheader:
        header_regions = []
        if regions is not None:
            header_regions = ["region " + label for label in regions.get_labels()]
        header = header_general + header_regions + bam_headers 
        print(",".join(header))


def get_variant_region_intersections(regions, chrom, pos):
    result = []
    if regions is not None:
        overlaps = {}
        for _start, _end, label in regions.get_overlapping_regions(chrom, pos, pos+1):
            overlaps[label] = True
        for output_label in regions.get_labels():
            if output_label in overlaps:
                result.append(str(overlaps[output_label]))
            else:
                result.append('False')
    return result
 

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


# Return the variant position as a fraction of the entire chromosome
# length. 
def get_chrom_pos_fraction(readers, chrom, pos):
    # We assume (require) that all BAM files are aligned
    # to the same reference genome, so it doesn't matter which
    # one we use to decide where in the chromosome a position appears.
    # Thus it is safe to use the first one.
    if len(readers) > 0:
        first_reader = readers[0]
        this_chrom_length = first_reader.get_reference_length(chrom)
        return pos / this_chrom_length
    else:
        return ''


def process_variants_bams(options, regions, variants):
    bam_labels = get_bam_labels(options.labels, options.bams)
    bam_readers = [BamReader(filepath) for filepath in options.bams]
    write_header(options, bam_labels, regions)
    for chrom, pos, ref, alt in variants:
        pos_normalised = get_chrom_pos_fraction(bam_readers, chrom, pos)
        region_counts = get_variant_region_intersections(regions, chrom, pos)
        bams_features = [reader.variant_features(chrom, pos, ref, alt) for reader in bam_readers]
        write_output_row(chrom, pos, ref, alt, pos_normalised, options.sample, region_counts, bams_features)
    for reader in bam_readers:
        reader.close()


def write_output_row(chrom, pos, ref, alt,  pos_normalised, sample, regions, bam_features):
    row_bams = [str(x) for bam in bam_features for x in bam.as_list()]
    row = [chrom, str(pos), ref, alt, str(pos_normalised), sample] + regions + row_bams 
    print(",".join(row))
 

VALID_DNA_BASES = "ATGCN"

class ReadFeatures(object):

    fields = ["avg NM", "avg base qual", "avg map qual",
              "avg align len", "avg clipped bases", "avg indel bases",
              "fwd strand", "rev strand", "normalised read position"]

    def __init__(self):
        self.nm = 0
        self.base_qual = 0
        self.map_qual = 0
        self.align_len = 0
        self.clipping = 0
        self.indel = 0
        self.forward_strand = 0
        self.reverse_strand = 0
        self.normalised_read_position = 0 
        self.num_reads = 0

    # See: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_cigar_stats
    def count(self, read):
        self.num_reads += 1
        query_position = read.query_position
        query_length = read.alignment.query_length
        if query_position is not None and query_length > 0:
            self.normalised_read_position += query_position / read.alignment.query_length
        if not read.is_del and not read.is_refskip:
            self.base_qual += read.alignment.query_qualities[query_position]
        alignment = read.alignment
        self.map_qual += alignment.mapping_quality 
        self.align_len += alignment.query_alignment_length
        cigar_stats = alignment.get_cigar_stats()[0]
        if len(cigar_stats) == 11:
            self.nm += cigar_stats[10]
            # count soft and hard clips together
            self.clipping += cigar_stats[4] + cigar_stats[5]
            self.indel += cigar_stats[1] + cigar_stats[2]
        if alignment.is_reverse:
            self.reverse_strand += 1
        else:
            self.forward_strand += 1

    def as_list(self):
        # normalise the results to the number of observed reads in total
        if self.num_reads > 0:
            result = [x / self.num_reads for x in [
                self.nm, self.base_qual, self.map_qual,
                self.align_len, self.clipping, self.indel,
                self.forward_strand, self.reverse_strand,
                self.normalised_read_position]]
        else:
            result = ['' for _ in ReadFeatures.fields]
        return result 

class BaseCounts(object):
    fields = ["depth", "A", "T", "G", "C", "N", "ref", "alt", "alt vaf"]

    def __init__(self, ref, alt):
        self.ref = ref
        self.alt = alt
        self.A = 0
        self.T = 0
        self.G = 0
        self.C = 0
        self.N = 0

    def count(self, base):
        if base == 'A':
            self.A += 1
        elif base == 'T':
            self.T += 1
        elif base == 'G':
            self.G += 1
        elif base == 'C':
            self.C += 1
        elif base == 'N':
            self.N += 1

    def as_list(self):
        total_depth = self.A + self.T + self.G + self.C + self.N
        ref_count = getattr(self, self.ref)
        alt_count = getattr(self, self.alt)
        if total_depth > 0:
            alt_vaf = alt_count / total_depth
        else:
            alt_vaf = ''
        return [total_depth, self.A, self.T, self.G, self.C, self.N, ref_count, alt_count, alt_vaf]


class LocusFeatures(object):
    fields = BaseCounts.fields + \
             ["ref " + x for x in ReadFeatures.fields] + \
             ["alt " + x for x in ReadFeatures.fields] + \
             ["all " + x for x in ReadFeatures.fields]

    def __init__(self, ref, alt):
        self.ref = ref
        self.alt = alt
        # counts of DNA bases at this pileup position
        self.base_counts = BaseCounts(ref, alt) 
        # features where the read contains the reference base at this position
        self.ref_read_features = ReadFeatures()
        # features where the read contains the alternative base at this position
        self.alt_read_features = ReadFeatures()
        # features for all reads that overlap this position, regardless of the base
        self.all_read_features = ReadFeatures()

    def count(self, read):
        # Get the DNA base from the current read if we can
        if not read.is_del and not read.is_refskip:
            base = read.alignment.query_sequence[read.query_position].upper()
            self.base_counts.count(base)
            if base == self.ref:
                # only count features of reads containing the reference base
                self.ref_read_features.count(read)
            elif base == self.alt:
                # only count features of reads containing the alternative base
                self.alt_read_features.count(read)
        # count features of reads regardless of the base
        self.all_read_features.count(read)

    def as_list(self):
        return self.base_counts.as_list() + \
               self.ref_read_features.as_list() + \
               self.alt_read_features.as_list() + \
               self.all_read_features.as_list()


MAX_PILEUP_DEPTH = 1000000000

class BamReader(object):
    def __init__(self, filepath):
        logging.info(f"Reading BAM file: {filepath}")
        resolved_path = os.path.realpath(filepath)
        logging.info(f"Resolved filepath to: {resolved_path}")
        self.samfile = pysam.AlignmentFile(resolved_path, "rb")

    # pos is expected to be 1-based 
    def variant_features(self, chrom, pos, ref, alt):
        zero_based_pos = pos - 1
        features = LocusFeatures(ref, alt)
        for pileupcolumn in self.samfile.pileup(chrom, zero_based_pos, zero_based_pos+1,
                                               truncate=True, stepper='samtools',
                                               ignore_overlaps=False, ignore_orphans=True,
                                               max_depth=MAX_PILEUP_DEPTH): 
            for read in pileupcolumn.pileups:
                features.count(read)
        return features 

    def get_reference_length(self, chrom):
        return self.samfile.get_reference_length(chrom)

    def close(self):
        self.samfile.close()


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

    def get_labels(self):
        return sorted(list(self.region_labels))


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
    logging.info("Completed")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
