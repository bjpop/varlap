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
PROGRAM_NAME = "snvly_scatter_plots"
DEFAULT_OUT_DIR = "scatter_plots"
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
    description = 'Generate plots of snvly features'
    parser = ArgumentParser(description=description)
    parser.add_argument(
        'data',  metavar='DATA', type=str, help='Filepaths of snvly CSV results file')
    parser.add_argument(
        '--outdir',  metavar='DIR', default=DEFAULT_OUT_DIR, type=str,
        help=f'Name of output directory. Default={DEFAULT_OUT_DIR}')
    parser.add_argument(
        '--nolegend', action='store_true', 
        help=f'Turn off the legend in the plot')
    parser.add_argument(
        '--hues',  metavar='FEATURE', type=str, nargs='+',
        help=f'Name of features to use for colouring dots (e.g. "sample" "chrom")')
    parser.add_argument(
        '--alpha',  metavar='ALPHA', type=float, default=DEFAULT_ALPHA,
        help=f'Alpha value for plotting points (default: {DEFAULT_ALPHA})')
    parser.add_argument(
        '--linewidth',  metavar='WIDTH', type=int, default=DEFAULT_LINEWIDTH,
        help=f'Line width value for plotting points (default: {DEFAULT_LINEWIDTH})')
    parser.add_argument(
        '--features',  metavar='FEATURE', nargs="+", required=True, type=str,
        help=f'Features to plot, format: feature1,feature2 (no spaces between feature names, e.g. "pos normalised","tumour depth")')
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
    for feature_pair in options.features:
        fields = feature_pair.split(",")
        if len(fields) == 2:
            feature1, feature2 = fields
            if {feature1, feature2}.issubset(df.columns):
                plt.clf()
                plt.suptitle('')
                fig, ax = plt.subplots(figsize=(10,8))
                if options.hues:
                    for hue in options.hues:
                        scatter_plot(options, df, feature1, feature2, hue, options.alpha, options.linewidth)
                else:
                    scatter_plot(options, df, feature1, feature2, None, options.alpha, options.linewidth)
            else:
                logging.warn(f"One or both of these features was not in the data: '{feature1}' '{feature2}', skipping")
        else:
            logging.warn(f"Badly formatted feature pair {feature_pair}, skipping")


def scatter_plot(options, df, feature1, feature2, hue, alpha, linewidth):
    plt.clf()
    plt.suptitle('')
    fig, ax = plt.subplots(figsize=(10,8))
    g=sns.scatterplot(data=df, x=feature1, y=feature2, hue=hue, alpha=options.alpha, linewidth=options.linewidth) 
    g.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    if options.nolegend:
        g.legend_.remove()
    feature1_str = feature1.replace(' ', '_')
    feature2_str = feature2.replace(' ', '_')
    hue_str = '' if hue is None else hue
    filename = os.path.join(options.outdir, 'scatter',
        feature1_str + "-" + feature2_str + '-' + hue_str + '.scatter.png')
    ax.set(xlabel=feature1, ylabel=feature2)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

PLOT_DIRS = ["scatter"]

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
