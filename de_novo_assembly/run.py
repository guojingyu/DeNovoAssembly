"""
This piece of script contains the entry method to run the package as a
DNA sequence assembler

Author Jingyu Guo
"""

import argparse
import logging
import datetime
import sys


# Configure logging
timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
logging.basicConfig(filename="../" + timestamp + "_denovo_assembly.log",
                    stream=sys.stdout,
                    level=logging.DEBUG)
stderrLogger=logging.StreamHandler()
stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
logger = logging.getLogger(__name__)
logger.addHandler(stderrLogger)

def main():
    """
    main method to run the DNA sequence assembly
    :return:
    """
    # get input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--fastafile',
                        help="fasta formatted text file containing the sequences for assembly",
                        default='./data/coding_challenge_data_set.fasta')
    parser.add_argument('--include_reverse_complement', action='store_false',
                        help="include to the reverse complement sequences of the fasta file sequence")
    parser.add_argument('-k', '--kmerlength',
                        help="parameter k to control the length of kmer as the edge of De Bruijn Graph",
                        default=3)
    parser.add_argument('--graph', action='store_false',
                        help="plot the De Bruijn Graph constructed from kmer of the input sequence")
    parser.add_argument('--outputfile', action='store_false',
                        help="output assembled DNA sequences into file with input name and time stamp appended")
    parser.add_argument('--output_longest_assembly', action='store_true',
                        help="output the longest assembled DNA sequence. If "
                             "more than one same longest DNA sequence found, randomly return one")

    args = parser.parse_args()

    logger.info(datetime.datetime.now().strftime("%H:%M:%S") +
                " Begin assemble DNA sequences ... ")
    # read in fasta file

    # build the De Bruijn Graph

    # Euler Walk to find the path

    # output to file if needed or

    # print to screen


if __name__ == "__main__":
    main()