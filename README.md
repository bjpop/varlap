# Overview 

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
snvly -h
usage: snvly [-h] [--tumour BAM] [--normal BAM] [--version] [--log LOG_FILE]

Read one or more FASTA files, compute simple stats for each file

optional arguments:
  -h, --help      show this help message and exit
  --tumour BAM    Filepath of tumour BAM file
  --normal BAM    Filepath of normal BAM file
  --version       show program's version number and exit
  --log LOG_FILE  record program progress in LOG_FILE

```

Run on a somatic VCF file supplying tumour and normal bams for the same donor:
```
snvly --tumour sample.tumour.bam --normal sample.normal.bam < sample.vcf
```

Or from a compressed VCF file:
```
zcat sample.vcf.gz | snvly --tumour sample.tumour.bam --normal sample.normal.bam
```

# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[snvly issue tracker](https://github.com/bjpop/snvly/issues)
