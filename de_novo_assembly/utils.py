"""
This script is to provide utility functions, suchas functions for text
based file IO


Author : Jingyu Guo
"""

from Bio import SeqIO
import networkx as nx
import matplotlib.pyplot as plt
import os.path
import sys
import logging
import datetime

# Configure logging
timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
logging.basicConfig(filename="../logs/" + timestamp + "_denovo_assembly.log",
                    stream=sys.stdout,
                    level=logging.DEBUG)
stderrLogger=logging.StreamHandler()
stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
logger = logging.getLogger(__name__)
logger.addHandler(stderrLogger)

def get_sequences_from_fasta_file(file, include_reverse_complement= True):
    """
    read fasta file and return all relevant sequences
    :param file: fasta file
    :param include_reverse_complement: if to include the reverse complement
    sequences
    :return: a dictionary of Seq objects (biopython)
    """
    if os.path.isfile(file):
        try:
            record_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
            logger.info(datetime.datetime.now().strftime("%H:%M:%S") +
                " Read the input file: " + file)
            logger.info(datetime.datetime.now().strftime("%H:%M:%S") +
                        " A number of " + str(len(record_dict.keys())) +
                        " records read in.")
            if include_reverse_complement:
                reversed_record_dict = {}
                for seq_id,seq in record_dict.iteritems():
                    reversed_record_dict[seq_id+"_rev_com"] = seq.reverse_complement()
                logger.info(datetime.datetime.now().strftime("%H:%M:%S") +
                            " A number of " + str(len(reversed_record_dict.keys())) +
                        " reverse complement records added.")
                record_dict.update(reversed_record_dict)
            return record_dict
        except IOError as err:
            print(err.args)
            print("Reading fasta file failed. Please check format :" + file)
            logger.error(datetime.datetime.now().strftime("%H:%M:%S") + " FileRead: " + err)
            sys.exit(1)

    else:
        print("The provided file not exist :" + file)
        logger.error(datetime.datetime.now().strftime("%H:%M:%S") + " The provided file not exist :" + file)
        sys.exit(1)


def visualize_graph(graph):
    """

    :param G:
    :return:
    """
    nx.draw(graph, cmap=plt.get_cmap('jet'))
    plt.show()
    pass

if __name__ == "__main__":
    # run the file reading
    fa_test_file = "../data/coding_challenge_data_set.fasta"
    get_sequences_from_fasta_file(fa_test_file)

    # run the graph plot