# Overview 

Snvly is designed to compute various bits of information about somatic SNVs from corresponding tumour and normal BAM files.

It is particularly intended for quality control checking of somatic SNVs. For example, you can use the output of snvly
to check for outliers in the data that might be indicative of errors. Therefore, many of the bits of information collected
by snvly are often signs of things going wrong in the sequencing and alignment.

These are some examples of the kinds of information collected per variant:
* Base quality.
* Mapping quality of overlapping reads.
* Number of mismatches in overlapping reads (NM, "edit distance").
* Alignment length of overlapping reads.
* Counts of DNA bases at locus of variant.

In the examples below, `$` indicates the command line prompt.

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

# General behaviour

Get help:
```
$ snvly -h
usage: snvly [-h] [--tumour BAM] [--normal BAM] --sample SAMPLE [--version]
             [--log LOG_FILE]

Compute various bits of information about somatic variants in tumour and
normal BAM files

optional arguments:
  -h, --help       show this help message and exit
  --tumour BAM     Filepath of tumour BAM file
  --normal BAM     Filepath of normal BAM file
  --sample SAMPLE  Sample identifier
  --version        show program's version number and exit
  --log LOG_FILE   record program progress in LOG_FILE

```

Run on a somatic VCF file supplying tumour and normal bams for the same donor:
```
snvly --sample sample_id --tumour tumour.bam --normal normal.bam < variants.vcf
```

Or from a compressed VCF file:
```
zcat variants.vcf.gz | snvly --sample sample_id --tumour tumour.bam --normal normal.bam
```

# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[snvly issue tracker](https://github.com/bjpop/snvly/issues)
