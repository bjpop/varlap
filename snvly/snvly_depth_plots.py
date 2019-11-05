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
PROGRAM_NAME = "snvly_depth_plots"
DEFAULT_OUT_DIR = "depth_plots"
DEFAULT_ALPHA = 0.3
DEFAULT_LINEWIDTH = 0

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
    description = 'Generate plots of scaled depth of coverage'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        'data',  metavar='DATA', type=str, help='Filepaths of snvly CSV results file')
    parser.add_argument(
        '--labels', required=True, nargs='+', metavar='LABEL', type=str, help='columns to plot')
    parser.add_argument(
        '--outdir',  metavar='DIR', default=DEFAULT_OUT_DIR, type=str,
        help=f'Name of output directory. Default={DEFAULT_OUT_DIR}')
    parser.add_argument(
        '--nolegend', action='store_true', 
        help=f'Turn off the legend in the plot')
    parser.add_argument(
        '--alpha',  metavar='ALPHA', type=float, default=DEFAULT_ALPHA,
        help=f'Alpha value for plotting points (default: {DEFAULT_ALPHA})')
    parser.add_argument(
        '--linewidth',  metavar='WIDTH', type=int, default=DEFAULT_LINEWIDTH,
        help=f'Line width value for plotting points (default: {DEFAULT_LINEWIDTH})')
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
    df_main = pd.read_csv(options.data, sep=",", dtype={'chrom': str})

    for label in options.labels:
        if label in df_main.columns:
            # Ignore samples with depth == 0, otherwise they tend to skew the mean depth
            df = df_main[df_main[label] > 0]
            # Get the mean depth of coverage per sample
            df_sample_mean = df.groupby("sample").agg("mean")
            df = df_main.join(df_sample_mean[label], how="left", on="sample", rsuffix=' mean')
            df['scaled depth'] = df[label] / df[label + " mean"]

            plt.clf()
            plt.suptitle('')
            fig, ax = plt.subplots(figsize=(10,8))
            g=sns.scatterplot(data=df, x='pos normalised', y='scaled depth', hue='sample', alpha=options.alpha, linewidth=options.linewidth) 
            if options.nolegend:
                g.legend_.remove()
            filename = os.path.join(options.outdir, label + '.depth.png')
            ax.set(xlabel='pos normalised', ylabel='depth of coverage fold change')
            plt.tight_layout()
            plt.savefig(filename)
            plt.close()
        else:
            logging.warn(f"Label {label} not in data, skipping")

def make_output_directories(outdir):
    os.makedirs(outdir, exist_ok=True)


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
