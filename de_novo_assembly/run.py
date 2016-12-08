#!/usr/bin/env python

"""
This piece of script contains the entry method to run the package as a
DNA sequence assembler

Author Jingyu Guo
"""

import argparse
import logging
import datetime
import sys
from utils import get_sequences_from_fasta_file, include_reverse_complement,\
    output_assembly, clock_now
from de_bruijn_graph import DeBruijnGraph, visualize_graph
from eulerian import eulerian_random_walk


# Configure logging
timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
logging.basicConfig(filename="../" + timestamp + "_denovo_assembly.log",
                    stream=sys.stdout,
                    level=logging.DEBUG)
stderrLogger=logging.StreamHandler()
stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
logger = logging.getLogger(__name__)
logger.addHandler(stderrLogger)

def run(inputfastafile, outputfile, kmerlength,
        reverse_complement = False, graph = False,
        longest_assembly_only = False, print_to_console =
        True):
    """
    main method to run the DNA sequence assembly
    :return:
    """
    logger.info(clock_now() +
                " Beginning: Begin assemble DNA sequences ... ")
    # read in fasta file
    try:
        seq_dict = get_sequences_from_fasta_file(inputfastafile)
    except IOError as err:
        logger.error(clock_now() + " Exit: " + err.message + ". Quitted.")
        sys.exit(1)
    except RuntimeError as err:
        logger.error(clock_now() + " Exit: " + err.message + ". Quitted.")
        sys.exit(1)

    # if needed to include the reverse complement sequences
    if reverse_complement:
        seq_dict = include_reverse_complement(seq_dict)

    # build the De Bruijn Graph
    try:
        dbg = DeBruijnGraph(seq_dict,int(kmerlength))
        if graph:
            visualize_graph(dbg.G)
    except ValueError as err:
        print err.message
        logger.error(clock_now() + " Exit: Value Error for input k. Quitted.")
        sys.exit(1)
    except RuntimeError as err:
        print err.message
        logger.error(clock_now() + " Exit: Other errors. Quitted.")
        sys.exit(1)

    # Euler Walk to find the paths
    # the return should be a list of sequences
    try:
        assembly = eulerian_random_walk(dbg)
    except RuntimeError as err:
        print "Exit: fatal errors during random walk. Quitted."
        logger.error(clock_now() + " Exit: fatal errors during random walk. "
                                   "Quitted.")
        sys.exit(1)
    except ValueError as err:
        print "Eulerian: errors during random walk. no output can be generated."
        logger.error(clock_now() + " Eulerian: errors during random walk. no output can be generated.")
        logger.error(err.message)
        sys.exit(1)


    # output to file if needed or
    if len(assembly) != 0:
        output_assembly(assembly,outputfile,longest_assembly_only)
    else:
        print "Exit: no assembled DNA returned. Quitted."
        logger.warn(clock_now() + " Exit: no assembled DNA returned. Quitted.")
        sys.exit(1)

    # print to screen
    if print_to_console:
        if longest_assembly_only:
            longest_seq = max(assembly, key=len)
            print(longest_seq)
            print(len(longest_seq))
        else:
            for sequence in sorted(assembly,key=len,reverse=True):
                print("Assembly Results with length :")
                print(sequence)
                print(len(sequence))
    # log the termination
    logger.info(clock_now() +
                " Done.")

    return assembly

if __name__ == "__main__":
    # get input parameters
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfastafile',
                        help="fasta formatted text file containing the sequences for assembly",
                        default='../data/dummy_data.fasta')
    parser.add_argument('--reverse_complement', action='store_true',
                        help="include or not. to include to the reverse "
                             "complement sequences of the fasta file sequence",
                        default=False)
    parser.add_argument('-k', '--kmerlength',
                        help="parameter k to control the length of kmer as the edge of De Bruijn Graph",
                        default=6)
    parser.add_argument('--graph', action='store_true',
                        help="include or not. plot the De Bruijn Graph constructed from kmer of the input sequence",
                        default=False)
    parser.add_argument('-o', '--outputfile',
                        help="output assembled DNA sequences into file with "
                             "input name and time stamp appended",
                        default='./output_DNA_assembly.txt')
    parser.add_argument('--output_longest_assembly', action='store_true',
                        help="include or not. output the longest assembled "
                             "DNA sequence. If more than one same longest DNA sequence found, "
                             "return the first one in rank",
                        default=False)
    parser.add_argument('--print_to_console', action='store_true',
                        help="include or not. to print assembled DNA sequence to "
                             "console",
                        default=False)

    args = parser.parse_args()

    assembly = run(args.inputfastafile, args.outputfile, args.kmerlength,
                   args.reverse_complement, args.graph,
                   args.output_longest_assembly, args.print_to_console)

    # print seq and length
    # for sequence in sorted(assembly, key=len, reverse=True):
    #     print(sequence + '\n')
    #     print(len(sequence))