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
    name='varlap',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['varlap'],
    package_dir={'varlap': 'varlap'},
    entry_points={
        'console_scripts': ['varlap = varlap.varlap:main',
                           'varlap_outliers = varlap.varlap_outliers:main',
                           'varlap_dist_plots = varlap.varlap_dist_plots:main', 
                           'varlap_scatter_plots = varlap.varlap_scatter_plots:main',
                           'varlap_depth = varlap.varlap_depth:main',
                           'varlap_depth_plots = varlap.varlap_depth_plots:main']
    },
    url='https://github.com/bjpop/varlap',
    license='LICENSE',
    description=('Check somatic SNVs against their germline support'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["pysam==0.16.0.1", "intervaltree==3.0.2", "pandas==0.25.2"],
)
