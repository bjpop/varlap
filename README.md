# Overview 

Varlap is primarily a quality control tool for genetic variants arising from high throughput DNA sequencing, where the
variants have been called by aligning DNA sequencing reads to a reference genome. It takes as input a set of DNA variants
and one or more BAM files. Varlap considers the genomic locus of each variant in each of the supplied
BAM files and records information about the corresponding alignment context at that locus. For example, one of the metrics
it calculates is the average edit distance of reads overlapping the variant locus. This can be a useful metric because
regions with significantly higher average edit distance are more likely to contain erroneous variant calls. 
Varlap outputs a CSV file containing one row per input variant, with columns recording the various computed metrics about that variant.
Subsequent analysis of this output (such as outlier detection) can be used to identify potentially problematic variants and
samples.

Common use cases are to consider somatic variants in the context of tumour and normal alignments, or germline variants against normal alignments. 
However, varlap is quite flexible and allows the use of any number of BAM files as input. 

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

Varlap is also accompanied by some helpful auxilliary programs that can:
* Detect outliers in selected features in the data.

In the examples below, `$` indicates the command line prompt.

# Assumptions

Varlap does not do any filtering or interpretation of the variants in the input VCF file. For example, it does not assume the variants
originate from any particular sample, it does not check any genotype information in the VCF file, and it does not filter the variants
in any way. It simply checks each variant locus in each of the supplied BAM files. The only information it uses from the input VCF file
is the chromosome, position, reference allele and alternative allele. Currently it only supports variants with one alternative allele
in a given row in the input file.

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/bjpop/varlap/master/LICENSE).

# Installing

## Installing from source

Clone this repository: 
```
$ git clone https://github.com/bjpop/varlap
```

Move into the repository directory:
```
$ cd varlap
```

Python 3 (>= 3.6) is required for this software.

Varlap can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
$ python3 -m venv varlap_dev
$ source varlap_dev/bin/activate
$ pip install -U /path/to/varlap
```
2. Into the global package database for all users:
```
$ pip install -U /path/to/varlap
```
3. Into the user package database (for the current user only):
```
$ pip install -U --user /path/to/varlap
```

## Docker container 

### Building the Docker container from source

The file `Dockerfile` contains instructions for building a Docker container for varlap, assuming you have a copy of the varlap source code on your computer.

If you have Docker installed on your computer you can build the container by running the following command within the top directory of the source code: 
```
$ docker build -t varlap .
```
See below for information about running varlap within the Docker container.

### Pulling the Docker container from Docker Hub 

Alternatively, you can pull latest version of the Docker container from Docker Hub like so:
```
docker pull bjpop/varlap:latest
```
Or, if you are using Singularity, then you can pull the docker container like so:
```
singularity pull docker://bjpop/varlap:latest
```

# Output

Varlap generates a CSV file on the standard output device.

The following column headings are always present for every output row:

* chrom: the chromosome on which the variant appears (same as input VCF)
* pos: the genomic position of the variant on the chromosomes (1-based coordinates, same as input VCF)
* ref: the reference allele of the variant (same as input VCF)
* alt: the alternative allele of the variant (same as input VCF)
* sample: the sample identifier if provided on the command line with `--sample`, otherwise this is empty. This feature is useful if you intend to analyse many samples separately and then combine their output CSV files together, for a group analysis.

If a `--regions` command line argument is given, varlap will produce a column for every unique feature in the input regions BED file. For example, if the input regions file contained information about repeats, and there were 6 different kinds of repeats in the file, then there would be 6 columns in the output, one for each feature label. For a given region column, each row in the output would contain either True or False indicating whether the variant intersected with the given feature.

For each BAM file in the input, varlap generates the following output columns:

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
$ varlap -h
usage: varlap [-h] [--labels [LABEL [LABEL ...]]] [--sample SAMPLE] [--regions REGIONS] [--format FORMAT] [--maxindelsize SIZE] [--noheader] [--version] --varclass TYPE
              [--log LOG_FILE]
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
  --format FORMAT       File format of input variants file. Options: CSV, TSV, VCF. Default: CSV
  --maxindelsize SIZE   Maximum size of indels to consider. Indels above this size will be excluded from analysis and output.
  --noheader            Suppress output header row
  --version             show program's version number and exit
  --varclass TYPE       Type of variants to consider. Options: SNV, INDEL
  --log LOG_FILE        record program progress in LOG_FILE
```

# Optional regions BED file

Varlap permits the use of an optional input "regions" BED file that specifies genomic regions that may be of interest in your analysis. A typical example might be to supply information about repeat regions, such as those computed by [RepeatMasker](http://www.repeatmasker.org).

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
varlap --sample sample_id --labels tumour normal --regions regions.bed -- tumour.bam normal.bam < variants.vcf > variants.varlap.csv
```

Consider germline variants in the context of a normal BAM:
```
varlap --sample sample_id --labels normal -- normal.bam < variants.vcf > variants.varlap.csv
```

# Computing outliers

The `varlap_outliers` program reads the output of `varlap` and annotates features on variants that are statistically outliers in the context of the entire data set.

You must specify what features you would like to consider (column headings) and what chromosomes to consider. Note that the selection of chromosomes may
influence what data points are considered outliers, due to underlying biological differences. For example, in human data, it is adviseable to consider 
the autosomes separately from the X and Y chromosomes. Note that each feature is considered independently for outliers, and thus no testing is done for outliers in mutliple dimensions simultaneously. 

The output will contain all the columns from the input data plus an additional N boolean columns, one for each feature specified on the command line. If a variant is an outlier with respect
to some feature, it will contain True in the correspond column, otherwise False.

Outliers are computed based on the interquartile range of the data, using the so-called [Tukey's fences](https://en.wikipedia.org/wiki/Outlier#Tukey's_fences). Data points
outside the below range are considered outliers:

```
    [Q1 - k(Q3 - Q1), Q3 + k(Q3 - Q1)]

    where Q1 = the first quartile
          Q3 = the second quartile
          k = a scaling factor, determined by the command line parameter --stringency
```
Suggested values for k are 1.5 for outliers and 3 for "far" outliers. Some experimentation with this parameter may be useful for a given set of data.

```
$ varlap_outliers -h
usage: varlap_outliers [-h] [--stringency FLOAT] --chroms CHROM [CHROM ...]
                      --features FEATURES [FEATURES ...] [--noheader]
                      [--version] [--log LOG_FILE]
                      DATA

Compute outliers in varlap outputs

positional arguments:
  DATA                  Filepaths of varlap CSV results file

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
$ varlap_outliers --stringency 3 --chroms '1' '2' '3' --features 'tumour all avg NM' 'tumour alt vaf' 'normal alt vaf' -- variants.varlap.csv > variants.varlap.outliers.csv
```

# Running within the Docker container

The following section describes how to run varlap within the Docker container. It assumes you have Docker (or Singularity) installed on your computer and have built (or pulled) the container as described above. 
The container behaves in the same way as the normal version of varlap, however there are some Docker-specific details that you must be aware of.

The general syntax for running varlap within Docker is as follows:
```
$ docker run -i varlap CMD
```
where CMD should be replaced by the specific command line invocation of varlap. Specific examples are below.

Display the help message:
```
$ docker run -i varlap varlap -h
```
Note: it may seem strange that `varlap` is mentioned twice in the command. The first instance is the name of the Docker container and the second instance is the name of the varlap executable that you want to run inside the container.

Display the version number:
```
$ docker run -i varlap varlap --version
```

Read from multuple input BAM files named on the command line, where all the files are in the same directory. You must replace `DATA` with the absolute file path of the directory containing the BAM files:  
```
$ docker run -i -v DATA:/in varlap varlap --labels tumour normal -- /in/tumour.bam /in/normal.bam < sample.vcf > sample.varlap.csv
```
The argument `DATA:/in` maps the directory called DATA on your local machine into the `/in` directory within the Docker container.

Logging progress to a file in the directory OUT: 
```
$ docker run -i -v DATA:/in -v OUT:/out varlap varlap --log /out/logfile.txt --labels tumour normal -- /in/tumour.bam /in/normal.bam < sample.vcf > sample.varlap.csv
```
Replace `OUT` with the absolute path of the directory to write the log file. For example, if you want the log file written to the current working directory, replace `OUT` with `$PWD`.
As above, you will also need to replace `DATA` with the absolite path to the directory containing your input BAM files.

## Usage with Singularity

Singularity can be used to run Docker containers. This can be useful in some environments where Docker is not available (e.g. High Performance Computing systems).

The principles are similar to using Docker, though some of the command line syntax is different. The example below shows how to run varlap on a tumour normal pair, assuming that the Docker container has been imported as `varlap_latest.sif`. Instructions for pulling the container using Singularity are provided above.

```
singularity exec --containall -B DATA:/bam varlap_latest.sif varlap --labels tumour normal -- /bam/tumour.bam /bam/normal.bam < sample.vcf > sample.varlap.csv
```

As with the Docker example, `DATA` is the path to a directory on your local machine that contains the input BAM files.

# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[varlap issue tracker](https://github.com/bjpop/varlap/issues)
