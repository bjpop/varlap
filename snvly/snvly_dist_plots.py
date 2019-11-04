'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 16 Oct 2019 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX


Generate distribution plots of snvly features 
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "snvly_plots"
DEFAULT_OUT_DIR = "plots"

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
    description = 'Generate plots of snvly features'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        'data',  metavar='DATA', type=str, help='Filepaths of snvly CSV results file')
    parser.add_argument(
        '--outdir',  metavar='DIR', default=DEFAULT_OUT_DIR, type=str,
        help=f'Name of output directory. Default={DEFAULT_OUT_DIR}')
    parser.add_argument(
        '--features',  metavar='FEATURE', nargs="+", required=True, type=str,
        help=f'Features to plot')
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


def plots(options):
    df = pd.read_csv(options.data, sep=",", dtype={'chrom': str})
    plot_distributions(df, options)
    if 'chrom' in df.columns:
        plot_distributions_by(df, options, 'chrom')
    else:
        logging.warn(f"Feature: chrom does not exist in data, skipping")
    if 'sample' in df.columns:
        plot_distributions_by(df, options, 'sample')
    else:
        logging.warn(f"Feature: sample does not exist in data, skipping")

def plot_distributions(df, options):
    for feature in options.features:
        if feature in df.columns:
            plt.clf()
            plt.suptitle('')
            fig, ax = plt.subplots(figsize=(10,8))
            sns.distplot(df[feature], kde=False, bins=100) 
            filename = os.path.join(options.outdir, 'dist', feature.replace(' ', '_') + '.dist.png')
            ax.set(xlabel=feature, ylabel='num variants')
            plt.tight_layout()
            plt.savefig(filename)
            plt.close()
        else:
            logging.warn(f"Feature: {feature} does not exist in data, skipping")


def plot_distributions_by(df, options, key):
    for feature in options.features:
        if feature in df.columns:
            plt.clf()
            plt.suptitle('')
            fig, ax = plt.subplots(figsize=(10,8))
            sns.boxplot(data=df, x=key, y=feature) 
            filename = os.path.join(options.outdir, "dist" + "_" + key,
                                    feature.replace(' ', '_') + '.dist_' + key + '.png')
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.set(xlabel=key, ylabel=feature)
            plt.tight_layout()
            plt.savefig(filename)
            plt.close()
        else:
            logging.warn(f"Feature: {feature} does not exist in data, skipping")


PLOT_DIRS = ["dist", "dist_chrom", "dist_sample"]

def make_output_directories(outdir):
    for subdir in PLOT_DIRS:
        os.makedirs(os.path.join(outdir, subdir), exist_ok=True)


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    make_output_directories(options.outdir)
    plots(options)
    logging.info("Completed")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
