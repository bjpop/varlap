# Overview 

Snvly is designed to compute various bits of information about SNVs from corresponding BAM files.

Common use cases are to consider somatic variants in the context of tumour and normal alignments, or germline variants against normal alignments. 
However, snvly is quite flexible and allows the use of any number of BAM files as input. 

It is particularly intended for quality control checking of SNVs. For example, you can use the output of snvly
to check for outliers in the data that might be indicative of errors. Therefore, many of the bits of information collected
by snvly are often signs of things going wrong in the sequencing and alignment.

These are some examples of the kinds of information collected per variant:
* Base quality.
* Mapping quality of overlapping reads.
* Number of mismatches in overlapping reads (NM, "edit distance").
* Alignment length of overlapping reads.
* Counts of DNA bases at locus of variant.
* Number of INDELs in overlapping reads.
* Number of clipped bases in overlapping reads.
* Number of overlapping reads on forward and reverse strands.
* Optionally, whether the variant overlaps a genomic region of interest (such as a repeat region) - the user must supply a bed file of regions.

These features are computed for all reads overlapping the variant locus, and for those reads that contain only the reference and alt alleles.

Snvly is also accompanied by some helpful auxilliary programs that can:
* Produce summary plots of the computed data.
* Detect outliers in selected features in the data.

In the examples below, `$` indicates the command line prompt.

# Assumptions

Snvly does not do any filtering or interpretation of the variants in the input VCF file. For example, it does not assume the variants
originate from any particular sample, it does not check any genotype information in the VCF file, and it does not filter the variants
in any way. It simply checks each variant locus in each of the supplied BAM files. The only information it uses from the input VCF file
is the chromosome, position, reference allele and alternative allele. Currently it only supports variants with one alternative allele
in a given row in the input file.

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bjpop/snvly/master/LICENSE).

# Installing

Clone this repository: 
```
$ git clone https://github.com/bjpop/snvly
```

Move into the repository directory:
```
$ cd snvly
```

Python 3 is required for this software.

Snvly can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
$ python3 -m venv snvly_dev
$ source snvly_dev/bin/activate
$ pip install -U /path/to/snvly
```
2. Into the global package database for all users:
```
$ pip install -U /path/to/snvly
```
3. Into the user package database (for the current user only):
```
$ pip install -U --user /path/to/snvly
```

# Output

Snvly generates a CSV file on the standard output device.

The following column headings are always present for every output row:

* chrom: the chromosome on which the variant appears (same as input VCF)
* pos: the genomic position of the variant on the chromosomes (1-based coordinates, same as input VCF)
* ref: the reference allele of the variant (same as input VCF)
* alt: the alternative allele of the variant (same as input VCF)
* sample: the sample identifier if provided on the command line with `--sample`, otherwise this is empty. This feature is useful if you intend to analyse many samples separately and then combine their output CSV files together, for a group analysis.

If a `--regions` command line argument is given, snvly will produce a column for every unique feature in the input regions BED file. For example, if the input regions file contained information about repeats, and there were 6 different kinds of repeats in the file, then there would be 6 columns in the output, one for each feature label. For a given region column, each row in the output would contain either True or False indicating whether the variant intersected with the given feature.

For each BAM file in the input, snvly generates the following output columns:

* depth: total depth of sequencing observed at this locus. Specifically this is the sum of the number of As, Ts, Gs, Cs and Ns at this locus in reads from the BAM file that overlap the locus 
* A: number of As at this locus in reads from the BAM file that overlap the locus
* T: number of Ts at this locus in reads from the BAM file that overlap the locus
* G: number of Gs at this locus in reads from the BAM file that overlap the locus
* C: number of Cs at this locus in reads from the BAM file that overlap the locus
* N: number of Ns at this locus in reads from the BAM file that overlap the locus
* ref: number of reference allele bases at this locus in reads from the BAM file that overlap the locus
* alt: number of alternative allele bases at this locus in reads from the BAM file that overlap the locus
* alt vaf: variant allele fraction (between 0 and 1) at this locus. Specifically this is the number of alternative allele bases divided by the total number of bases observed at this locus.

The following information is collected for reads overlapping the locus in 3 ways: 1) only reads that contain the reference allele at the locus, 2) only reads that contain the alternative allele at the locus, and 3) all reads overlapping the locus regardless of the base at the locus:

* avg NM: average edit distance in reads overlapping the locus 
* avg base qual: average base quality of refererence bases at this locus
* avg map qual: average mapping quality of reads overlapping the locus
* avg align len: average alignment length of reads overlapping the locus
* avg clipped bases: average number of soft and hard clipped bases in reads overlapping the locus
* avg indel bases: average number of INDELs in reads overlapping the locus
* fwd strand: average number of reads overlapping the locus that map to the forward strand
* rev strand: average number of reads overlapping the locus that map to the reverse strand

# General behaviour

Get help:
```
snvly -h
usage: snvly [-h] [--labels [LABEL [LABEL ...]]] [--sample SAMPLE]
             [--regions REGIONS] [--noheader] [--version] [--log LOG_FILE]
             BAM [BAM ...]

Compute various bits of information about variants in one or more BAM files

positional arguments:
  BAM                   Filepaths of BAM files

optional arguments:
  -h, --help            show this help message and exit
  --labels [LABEL [LABEL ...]]
                        Labels for BAM files
  --sample SAMPLE       Sample identifier
  --regions REGIONS     Filepath of genomic regions-of-interest in BED format
  --noheader            Suppress output header row
  --version             show program's version number and exit
  --log LOG_FILE        record program progress in LOG_FILE
```

# Optional regions BED file

Snvly permits the use of an optional input "regions" BED file that specifies genomic regions that may be of interest in your analysis. A typical example might be to supply information about repeat regions, such as those computed by [RepeatMasker](http://www.repeatmasker.org).

The BED file must contain 4 columns: chromosome, start, end, label. The start and end coordinates follow the BED convention of being zero-based and semi-closed. The label can be any string that you like. Each unique label in this file will result in a column in the output, therefore it is adviseable to avoid having too many unique labels.

Below is an example of the first few lines of a possible region file for humans representing repeats from RepeatMasker:

```
#genoName genoStart genoEnd repClass
1 16777160 16777470 SINE
1 25165800 25166089 SINE
1 33553606 33554646 LINE
1 58720067 58720973 LINE
1 83886030 83886750 LTR
```

# Example usage

Consider somatic variants in the context of tumour and normal BAM files, using a regions bed file
```
snvly --sample sample_id --labels tumour normal --regions regions.bed -- tumour.bam normal.bam < variants.vcf > variants.snvly.csv
```

Consider germline variants in the context of a normal BAM:
```
snvly --sample sample_id --labels normal -- normal.bam < variants.vcf > variants.snvly.csv
```

# Computing outliers

The `snvly_outliers` program reads the output of `snvly` and annotates features on variants that are statistically outliers in the context of the entire data set.

You must specify what features you would like to consider (column headings) and what chromosomes to consider. Note that the selection of chromosomes may
influence what data points are considered outliers, due to underlying biological differences. For example, in human data, it is adviseable to consider 
the autosomes separately from the X and Y chromosomes. 

Outliers are computed based on the interquartile range of the data, using the so-called (Tukey's fences)[https://en.wikipedia.org/wiki/Outlier#Tukey's_fences]. Data points
outside the below range are considered outliers:

```
    [Q1 - k(Q3 - Q1), Q3 + k(Q3 - Q1)]

    where Q1 = the first quartile
          Q3 = the second quartile
          k = a scaling factor, determined by the command line parameter --stringency
```
Suggested values for k are 1.5 for outliers and 3 for "far" outliers. Some experimentation with this parameter may be useful for a given set of data.

```
snvly_outliers -h
usage: snvly_outliers [-h] [--stringency FLOAT] --chroms CHROM [CHROM ...]
                      --features FEATURES [FEATURES ...] [--noheader]
                      [--version] [--log LOG_FILE]
                      DATA

Compute outliers in snvly outputs

positional arguments:
  DATA                  Filepaths of snvly CSV results file

optional arguments:
  -h, --help            show this help message and exit
  --stringency FLOAT    Stringency factor for detecting outliers, default=1.5
  --chroms CHROM [CHROM ...]
                        Consider variants on these chromosomes
  --features FEATURES [FEATURES ...]
                        Features to consider for outlier values
  --noheader            Suppress output header row
  --version             show program's version number and exit
  --log LOG_FILE        record program progress in LOG_FILE
```

Example usage:

```
snvly_outliers --stringency 3 --chroms '1' '2' '3' --features 'tumour all avg NM' 'tumour alt vaf' 'normal alt vaf' variants.snvly.csv > variants.snvly.outliers.csv
```


# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[snvly issue tracker](https://github.com/bjpop/snvly/issues)
