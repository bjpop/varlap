#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''
Read a VCF file of somatic variants and check
their pileup in the corresponding BAM file for the
same donor. The purpose is to see how well supported
a somatic SNV is or isn't supported by the germline
evidence. This might be useful for ruling out
false positives.
'''


setup(
    name='snvly',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['snvly'],
    package_dir={'snvly': 'snvly'},
    entry_points={
        'console_scripts': ['snvly = snvly.snvly:main',
                           'snvly_outliers = snvly.snvly_outliers:main',
                           'snvly_dist_plots = snvly.snvly_dist_plots:main', 
                           'snvly_scatter_plots = snvly.snvly_scatter_plots:main']
    },
    url='https://github.com/bjpop/snvly',
    license='LICENSE',
    description=('Check somatic SNVs against their germline support'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["pysam", "intervaltree", "pandas", "seaborn"],
)
