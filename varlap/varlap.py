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
from collections import namedtuple
from intervaltree import Interval, IntervalTree
import os.path
import pathlib
from copy import copy

'''
Note: on the command line we allow the user to specify a desired "varclass"
to indicate the class of variant they want to consider. This can be "SNV" or
"INDEL" at the moment, but may be extended to other classes in the future.

However, within the program we break the classes of variants down further into
"vartype"s, which currently can be "SNV", "INS", "DEL". In this case "INS" and
"DEL" vartypes belong to the varclass "INDEL".
'''


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_BAD_FILE_FORMAT = 3
PROGRAM_NAME = "varlap"
VALID_VARIANT_TYPES = ["SNV", "INDEL"]
DEFAULT_INPUT_FILE_FORMAT = "CSV"
VALID_INPUT_FILE_FORMATS = ["CSV", "TSV", "VCF"]


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
        '--format', required=False, metavar='FORMAT', type=str, default=DEFAULT_INPUT_FILE_FORMAT, choices=VALID_INPUT_FILE_FORMATS,
        help='File format of input variants file. Options: %(choices)s. Default: %(default)s')
    parser.add_argument(
        '--maxindelsize', required=False, metavar='SIZE', type=int,
        help='Maximum size of indels to consider. Indels above this size will be excluded from analysis and output.')
    parser.add_argument(
        '--noheader', action='store_true', help='Suppress output header row')
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--varclass', required=True, metavar='TYPE', type=str,
        choices=VALID_VARIANT_TYPES, help="Type of variants to consider. Options: %(choices)s")
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


'''
VCF docs on CHROM, POS, REF, ALT:

CHROM - chromosome: An identifier from the reference genome or an angle-bracketed ID String (“<ID>”)
pointing to a contig in the assembly file (cf. the ##assembly line in the header). All entries for a specific
CHROM should form a contiguous block within the VCF file. The colon symbol (:) must be absent from all
chromosome names to avoid parsing errors when dealing with breakends. (String, no white-space permitted,
Required).

POS - position: The reference position, with the 1st base having position 1. Positions are sorted numerically,
in increasing order, within each reference sequence CHROM. It is permitted to have multiple records with the
same POS. Telomeres are indicated by using positions 0 or N+1, where N is the length of the corresponding
chromosome or contig. (Integer, Required)

REF - reference base(s): Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted.
The value in the POS field refers to the position of the first base in the String. For simple insertions and
deletions in which either the REF or one of the ALT alleles would otherwise be null/empty, the REF and ALT
Strings must include the base before the event (which must be reflected in the POS field), unless the event
occurs at position 1 on the contig in which case it must include the base after the event; this padding base is
not required (although it is permitted) for e.g. complex substitutions or other events where all alleles have at
least one base represented in their Strings. If any of the ALT alleles is a symbolic allele (an angle-bracketed
ID String “<ID>”) then the padding base is required and POS denotes the coordinate of the base preceding
the polymorphism. Tools processing VCF files are not required to preserve case in the allele Strings. (String,
Required).

ALT - alternate base(s): Comma separated list of alternate non-reference alleles. These alleles do not have to
be called in any of the samples. Options are base Strings made up of the bases A,C,G,T,N,*, (case insensitive)
or an angle-bracketed ID String (“<ID>”) or a breakend replacement string as described in the section on
breakends. The ‘*’ allele is reserved to indicate that the allele is missing due to a upstream deletion. If there
are no alternative alleles, then the missing value should be used. Tools processing VCF files are not required
to preserve case in the allele String, except for IDs, which are case sensitive. (String; no whitespace, commas,
or angle-brackets are permitted in the ID String itself)

'''

# XXX support for DELINS

class VariantReader(object):
    def __init__(self, file, file_format, varclass, max_indel_size=None):
        self.file = file
        self.file_format = file_format
        self.varclass = varclass
        self.max_indel_size = max_indel_size
        if file_format == "VCF":
            self.reader = vcf_reader(file)
            self.fieldnames = ["chrom", "pos", "ref", "alt", "vartype"]
        elif file_format == "CSV":
            self.reader = csv_tsv_reader(file, "CSV")
            self.fieldnames = self.reader.fieldnames + ["vartype"]
        elif file_format == "TSV":
            self.reader = csv_tsv_reader(file, "TSV")
            self.fieldnames = self.reader.fieldnames + ["vartype"]
        else:
            exit_with_error(f"Unknown file format: {format}", EXIT_BAD_FILE_FORMAT)
        self.total_variants_in_input = 0
        self.num_variants_analysed = 0
        self.num_variants_skipped = 0


    def get_variants(self):
        '''Read variants from input VCF file, yield one at a time'''
        for input_row in self.reader: 
            self.total_variants_in_input += 1
            if is_valid_input_row(input_row):
                this_ref = input_row["ref"]
                # allow possibly multiple alts in the same variant, split them into separate alleles
                alts = input_row["alt"].split(",")
                for this_alt in alts:
                    this_var_type = get_var_type(this_ref, this_alt)
                    if is_acceptable_variant(dict(input_row), self.varclass, this_var_type, this_ref, this_alt, self.max_indel_size):
                        output_row = copy(input_row)
                        output_row["alt"] = this_alt
                        output_row["pos"] = int(input_row["pos"])
                        output_row["vartype"] = this_var_type
                        self.num_variants_analysed += 1
                        yield output_row
            else:
                logging.warning(f"Skipping invalid input row: {dict(input_row)}")

    def log_totals(self):
        self.num_variants_skipped = self.total_variants_in_input - self.num_variants_analysed
        logging.info(f"Total variants in input: {self.total_variants_in_input}")
        logging.info(f"Num variants kept for analysis: {self.num_variants_analysed}")
        logging.info(f"Num variants skipped: {self.num_variants_skipped}")


# Check if an input variant is acceptable for this analysis and valid
# The purpose of this method is to reject and thus skup any input variants
# that we cannot intepret within the current analysis
def is_acceptable_variant(row, varclass, vartype, ref, alt, max_indel_size):
    if not is_only_DNA_bases(ref) or not is_only_DNA_bases(alt):
        # ref and alt must only consist of DNA bases
        logging.info(f"Skipping variant with non-DNA bases: {row}")
        return False
    if not is_desired_type(varclass, vartype):
        # this particular variant must be of a type that is compatible
        # with the varclass specified on the command line, (SNV, INDEL(
        logging.info(f"Skipping variant of unwanted type: {vartype} {row}")
        return False
    if not is_within_max_size(varclass, max_indel_size, ref, alt):
        # The variant must not be longer than the maximum size, if specified
        logging.info(f"Skipping variant that is too large {row}")
        return False
    if varclass == "INDEL" and not is_valid_indel(ref, alt):
        # Check INDELs for appropriate formatting
        logging.info(f"Skipping invalid INDEL: {row}")
        return False
    return True
    
VALID_DNA_BASES = set("ATGC")
    
def is_only_DNA_bases(sequence):
    return set(sequence.upper()).issubset(VALID_DNA_BASES) 


# Check that an input variant is within some size bound. This is normally
# only relevant for INDELs where the max_indel_size can be set on
# the command line. However, for completeness we also check that SNVs
# are indeed the expected size of 1 DNA base.
def is_within_max_size(varclass, max_indel_size, ref, alt):
    if varclass == "SNV":
        return len(ref) == 1 and len(alt) == 1
    elif varclass == "INDEL" and max_indel_size is not None:
        this_indel_size = abs(len(ref) - len(alt))
        return this_indel_size <= max_indel_size
    return True

# Check that the INDEL is specified in a way that we can interpet:
# there must be at least 1 context base, the shortest of ref and alt
# must be a prefix of the other. They must not have the same length.
def is_valid_indel(ref, alt):
    ref_len = len(ref)
    alt_len = len(alt)
    if ref_len < alt_len:
        shortest = ref
        longest = alt
    elif alt_len < ref_len:
        shortest = alt
        longest = ref
    else:
        # ref and alt are the same length, not currently supported 
        return False
    return ref_len >= 1 and alt_len >= 1 and longest.startswith(shortest)

REQUIRED_INPUT_VARIANT_FIELDS = set(["chrom", "pos", "ref", "alt"])

def is_valid_input_row(row):
    return set(row.keys()).issuperset(REQUIRED_INPUT_VARIANT_FIELDS)

def vcf_reader(file):
    for line in file:
        if line.startswith('#'):
            continue
        fields = line.strip().split()
        # Technically VCF requires the first 8 fields to be defined, but we want to be as liberal
        # as possible in accepting inputs.
        if len(fields) >= 5:
            chrom, pos, _id, ref, alt = fields[:5]
            yield {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
        else:
            logging.warning(f"Skipping input row: {line}")


def csv_tsv_reader(file, format="CSV"):
    if format == "CSV":
        reader = csv.DictReader(file)
    elif format == "TSV":
        reader = csv.DictReader(file, delimiter="\t")
    else:
        exit_with_error(f"Unknown file format: {format}", EXIT_BAD_FILE_FORMAT)
    return reader

def is_desired_type(varclass, vartype):
    return (varclass == "SNV" and vartype == "SNV" or
            varclass == "INDEL" and vartype in ["INS", "DEL"])


def get_var_type(ref, alt):
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    elif len(ref) > len(alt):
        return "DEL"
    elif len(alt) > len(ref):
        return "INS"
    else:
        logging.warning(f"Cannot determine the type of variant with ref: {ref} and alt: {alt}")
        return "UNKNOWN"



def write_header(options, variant_fieldnames, bam_labels, regions):
    header_general = variant_fieldnames + ["pos normalised", "sample"]
    if options.varclass == "SNV":
        bam_headers = [label + " " + field for label in bam_labels for field in LocusFeaturesSNV.fields]
    elif options.varclass == "INDEL":
        bam_headers = [label + " " + field for label in bam_labels for field in LocusFeaturesINDEL.fields]
    if not options.noheader:
        header_regions = []
        if regions is not None:
            header_regions = ["region " + label for label in regions.get_labels()]
        header = header_general + header_regions + bam_headers 
        print(",".join(header))


def get_variant_region_intersections(regions, variant):
    chrom = variant["chrom"]
    pos = variant["pos"]
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
def get_chrom_pos_fraction(readers, variant):
    chrom = variant["chrom"]
    pos = variant["pos"]
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

def variant_as_list(variant, fieldnames):
    return [str(variant[f]) for f in fieldnames]

def process_variants_bams(options, regions):
    bam_labels = get_bam_labels(options.labels, options.bams)
    bam_readers = [BamReader(filepath, options.varclass) for filepath in options.bams]
    variant_reader = VariantReader(sys.stdin, options.format, options.varclass, options.maxindelsize)
    write_header(options, variant_reader.fieldnames, bam_labels, regions)
    for variant in variant_reader.get_variants():
        pos_normalised = get_chrom_pos_fraction(bam_readers, variant)
        region_counts = get_variant_region_intersections(regions, variant)
        bams_features = [bam_reader.variant_features(variant) for bam_reader in bam_readers]
        output_variant = variant_as_list(variant, variant_reader.fieldnames)
        write_output_row(output_variant, pos_normalised, options.sample, region_counts, bams_features)
    variant_reader.log_totals()
    for reader in bam_readers:
        reader.close()


def write_output_row(variant, pos_normalised, sample, regions, bam_features):
    row_bams = [str(x) for bam in bam_features for x in bam.as_list()]
    row = variant + [str(pos_normalised), sample] + regions + row_bams 
    print(",".join(row))
 


class ReadFeatures(object):

    fields = ["avg NM", "avg base qual", "avg map qual",
              "avg align len", "avg clipped bases", "avg indel bases",
              "fwd strand", "rev strand", "supplementary", "normalised read position"]

    def __init__(self):
        self.nm = 0
        self.base_qual = 0
        self.map_qual = 0
        self.align_len = 0
        self.clipping = 0
        self.indel = 0
        self.forward_strand = 0
        self.reverse_strand = 0
        self.supplementary = 0
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
        if alignment.is_supplementary:
            self.supplementary += 1

    def as_list(self):
        # normalise the results to the number of observed reads in total
        if self.num_reads > 0:
            result = [x / self.num_reads for x in [
                self.nm, self.base_qual, self.map_qual,
                self.align_len, self.clipping, self.indel,
                self.forward_strand, self.reverse_strand, self.supplementary,
                self.normalised_read_position]]
        else:
            result = ['' for _ in ReadFeatures.fields]
        return result 

# XXX this should be SNVAlleleCounts, so we can have a similar INDELAlleleCounts

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



# True if the intervals of two indels overlap
def interval_overlaps(start1, end1, start2, end2):
    return not((end1 < start2) or (start1 > end2))

'''
Cigar operations in BAM file

M	BAM_CMATCH	0
I	BAM_CINS	1
D	BAM_CDEL	2
N	BAM_CREF_SKIP	3
S	BAM_CSOFT_CLIP	4
H	BAM_CHARD_CLIP	5
P	BAM_CPAD	6
=	BAM_CEQUAL	7
X	BAM_CDIFF	8
B	BAM_CBACK	9

cigartuples is returned as a list of tuples of (operation, length).
'''


IndelEvent = namedtuple("IndelEvent", ["indel_type", "start", "end", "bases"])

# Determine the allele in the read at the locus of an INDEL variant
def indels_overlapping_variant(read, var_start, var_end):
    alignment = read.alignment
    cigar = alignment.cigartuples
    # this is the position in the reference where the first non-clipped
    # base in the read aligns to
    ref_pos = alignment.reference_start
    read_pos = 0 
    result = []

    for operation, length in cigar:
        if operation == 0:
            # match
            ref_pos += length
            read_pos += length
        elif operation == 1:
            # INS
            this_start = ref_pos
            this_end = ref_pos + length - 1
            if interval_overlaps(var_start, var_end, this_start, this_end):
                inserted_bases = alignment.query_sequence[read_pos: read_pos + length].upper()
                result.append(IndelEvent(indel_type="INS", start=this_start, end=this_end, bases=inserted_bases))
            read_pos += length
        elif operation == 2:
            # DEL
            this_start = ref_pos
            this_end = ref_pos + length - 1
            if interval_overlaps(var_start, var_end, this_start, this_end):
                result.append(IndelEvent(indel_type="DEL", start=this_start, end=this_end, bases=""))
            ref_pos += length
        elif operation == 3:
            # REF_SKIP
            read_pos += length
        elif operation == 4:
            # SOFT_CLIP
            # XXX handle soft clips as potential insertions
            read_pos += length
        elif operation == 5:
            # HARD_CLIP
            pass
        elif operation == 6:
            # PAD
            pass
        elif operation == 7:
            ref_pos += length
            read_pos += length
            # EQUAL
        elif operation == 8:
            # DIFF
            ref_pos += length
            read_pos += length 
        elif operation == 9:
            # BACK, XXX not sure what this does
            pass

    return result 
        

# get the genome coordinates of where an INDEL variant will actually
# start, as opposed to the location of where the variant is reported
# the starting position must take into account the context bases that
# are given when the variant is reported in the VCF file.
def get_indel_start_coord(pos, ref, alt):
    len_ref = len(ref)
    len_alt = len(alt)
    shortest_len = min(len_ref, len_alt)
    return pos + shortest_len
    

class LocusFeaturesINDEL(object):

    fields = ["depth", "ref_allele", "alt_allele", "other_allele", "alt_vaf", "overlapping_indels"] + \
             ["ref " + x for x in ReadFeatures.fields] + \
             ["alt " + x for x in ReadFeatures.fields] + \
             ["all " + x for x in ReadFeatures.fields]

    def __init__(self, pos, ref, alt):
        self.ref = ref
        self.alt = alt
        # NOTE: INDELS are reported as one base to the left of the actual variant
        # to allow for a context base, we retain pos for the purpose and use 
        # start as the zero-based position of the first base of the actual variant
        self.pos = pos
        self.start = get_indel_start_coord(pos, ref, alt) 
        # XXX this size needs for be fixed for DELINS
        self.size = abs(len(ref) - len(alt))
        self.end = self.start + self.size - 1
        # XXX this size needs for be fixed for DELINS
        if len(ref) > len(alt):
            self.indel_type = "DEL"
            self.bases = ""
        elif len(alt) > len(ref):
            self.indel_type = "INS"
            self.bases = alt[1:]
        else:
            self.indel_type = "UNKNOWN"
        # features where the read contains the reference allele at this position
        self.ref_read_features = ReadFeatures()
        # features where the read contains the alternative allele at this position
        self.alt_read_features = ReadFeatures()
        # features for all reads that overlap this position, regardless of the base
        self.all_read_features = ReadFeatures()
        self.ref_count = 0
        self.alt_count = 0
        self.other_count = 0 
        self.depth = 0
        self.overlapping_indels_count = 0

    def count(self, read):
        self.depth += 1
        # get all the INDELs in the read that overlap this particular variant
        overlapping_indels = indels_overlapping_variant(read, self.start, self.end)
        self.overlapping_indels_count += len(overlapping_indels)
        # Check if any of the (possibly empty) overlapping INDELs support the ALT allele
        read_supports_alt = False
        for event in overlapping_indels:
            if event.indel_type == self.indel_type:
                if event.start == self.start and event.end == self.end:
                    if self.indel_type == "DEL" or (self.indel_type == "INS" and event.bases == self.bases):
                        read_supports_alt = True
                        break
        read_supports_ref = False
        if not overlapping_indels and read.query_position is not None:
            if self.indel_type == "INS":
                read_bases = read.alignment.query_sequence[read.query_position].upper()
            elif self.indel_type == "DEL":
                read_bases = read.alignment.query_sequence[read.query_position: read.query_position + self.size + 1].upper()
            if read_bases == self.ref:
                read_supports_ref = True
        if read_supports_ref:
            self.ref_count += 1
            self.ref_read_features.count(read)
        elif read_supports_alt:
            self.alt_count += 1
            self.alt_read_features.count(read)
        else:
            self.other_count += 1 
        self.all_read_features.count(read)

    def as_list(self):
        alt_vaf = ''
        if self.depth > 0:
            alt_vaf = self.alt_count / self.depth
        return [self.depth, self.ref_count, self.alt_count, self.other_count, alt_vaf, self.overlapping_indels_count] + \
               self.ref_read_features.as_list() + \
               self.alt_read_features.as_list() + \
               self.all_read_features.as_list()


class LocusFeaturesSNV(object):
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
    def __init__(self, filepath, varclass):
        logging.info(f"Reading BAM file: {filepath}")
        resolved_path = os.path.realpath(filepath)
        logging.info(f"Resolved filepath to: {resolved_path}")
        self.samfile = pysam.AlignmentFile(resolved_path, "rb")
        self.varclass = varclass 

    # pos is expected to be 1-based 
    def variant_features(self, variant):
        chrom = variant["chrom"]
        pos = variant["pos"]
        ref = variant["ref"]
        alt = variant["alt"]
        vartype = variant["vartype"] 
        zero_based_pos = pos - 1
        if vartype == "SNV" and self.varclass == "SNV":
            features = LocusFeaturesSNV(ref, alt)
        elif vartype in ["INS", "DEL"] and self.varclass == "INDEL":
            # NOTE: indel variants are reported one base to the left of the actual variant
            # to allow for the reference context base 
            features = LocusFeaturesINDEL(zero_based_pos, ref, alt)
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
    regions = get_regions(options.regions)
    process_variants_bams(options, regions)
    logging.info("Completed")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
