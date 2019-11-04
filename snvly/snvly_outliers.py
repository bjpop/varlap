'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 16 Oct 2019 
License     : MIT 
Maintainer  : bjpope@unimelb.edu.au 
Portability : POSIX

Outlier detection for snvly
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import pandas as pd


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
PROGRAM_NAME = "snvly_outliers"
DEFAULT_STRINGENCY = 1.5


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
    description = 'Compute outliers in snvly outputs'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        'data',  metavar='DATA', type=str, help='Filepaths of snvly CSV results file')
    parser.add_argument(
        '--stringency',  default=DEFAULT_STRINGENCY, metavar='FLOAT', type=float,
        help=f'Stringency factor for detecting outliers, default={DEFAULT_STRINGENCY}')
    parser.add_argument(
        '--chroms', required=True, nargs='+', metavar='CHROM', type=str,
        help=f'Consider variants on these chromosomes')
    parser.add_argument(
        '--features', nargs='+', required=True, metavar='FEATURES', type=str,
        help=f'Features to consider for outlier values')
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


def outliers(options):
    df = pd.read_csv(options.data, sep=",", dtype={'chrom': str})
    # only consider the selected chromosomes
    df = df[df['chrom'].isin(options.chroms)]
    outlier_feature_labels = []
    for feature in options.features:
        new_outlier_feature = feature + ' outlier'
        outlier_feature_labels.append(new_outlier_feature)
        df[new_outlier_feature] = feature_outliers(df[feature], options.stringency)

    df['num outliers'] = df[outlier_feature_labels].sum(axis=1)
    df.drop(df[df['num outliers'] == 0].index, inplace=True)
    df.to_csv(sys.stdout, sep=",")


def feature_outliers(feature_values, stringency):
    Q1 = feature_values.quantile(0.25)
    Q3 = feature_values.quantile(0.75)
    IQR = Q3 - Q1
    return (feature_values < (Q1 - stringency * IQR)) | (feature_values > (Q3 + stringency * IQR))

def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    outliers(options)
    logging.info("Completed")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
