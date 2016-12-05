"""
This script contains a class for De Bruijn Graph
DeBruijn(k-mers)
    from a node for each (k-1)-mer from k-mers
    for each k-mer in k-mers
        connect its prefix node with its suffix node by an edge

EulerianCycle(BalancedGraph)
    form a cycle by randomly walking in BalancedGraph (avoiding already
    visited edges)
    while cycle is not Eularian
        select a node newStart in Cycle with still unexplored outgoing edges
        form a Cycle' by traversing Cycle from newStart and randomly walking
        afterwards Cycle <- Cycle'
    return Cycle

PairedDeBruijnGraphs(paired k-mers):
    represent every paired k-mer as an edge between its paired prefix and
    paired suffix
    glue all nodes with identical labels


Author Jingyu Guo
"""
import networkx as nx
import logging
import datetime
import sys
from eulerian import is_euler, has_euler_path, has_euler_circuit

# Configure logging
timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
logging.basicConfig(filename=timestamp+"_denovo_assembly.log",
                    stream=sys.stdout,
                    level=logging.DEBUG)
stderrLogger=logging.StreamHandler()
stderrLogger.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
logger = logging.getLogger(__name__)
logger.addHandler(stderrLogger)


class Kmer():
    def __init__(self, kmer_as_a_str):
        """
        constructor for class Kmer to adapt the string
        representation into an edge presentation to efficiently convert
        the
        :param kmer: a str of the kmer
        """
        self.kmer = kmer_as_a_str
        self.k = len(self.kmer)
        # create the (k-1)-mers as De Bruijn node
        self.l_node = self.kmer[:-1]
        self.r_node = self.kmer[1:]


class DeBruijnGraph():
    """A de Bruijn Graph is a multigraph built from a collection of strings
    feed from input and k-mer length k. Nodes of DBG are
     (k-1)-mers and edges join a left k-1-mer to a right k-1-mer is the
     k-mer itself. """
    def __init__(self, sequence_dict, k = 3, *args, **kwargs):
        """
        constructor
        :param sequence_dict: sequence dict
        :param k: k for kmer length
        :param args:
        :param kwargs:
        """
        self.sequences = sequence_dict.values()
        min_len = len(min(self.sequences, key=len))
        if k < 2: # if k is smaller than 2 the assembly does not work
            logger.info("DBG : k is too small : " + str(k))
            raise Exception('DBG : k is too small : ' + str(k))
        # or if k is larger than the min length of given sequences
        elif k >= min_len:
            logger.info("DBG : k is larger than the min length of provide "
                        "sequence : " + str(k))
            raise Exception("DBG : k is larger than the min length of provide "
                        "sequence : " + str(k))
        else:
            self.k = k

        self.G = nx.MultiDiGraph()

        # construct the graph
        self.build_DBG()

        # check for subgraphs
        self.subgraphs = list(nx.connected_component_subgraphs(self.G))

        # TODO check subgraph features


    def convert_read_to_kmer(self,seq):
        """
        method to convert the reads into kmers
        Using generator here to save memory
        :return: an iterator
        """
        # convert sequence into kmers
        seq = DeBruijnGraph.legit_DNA_base_filter(seq)
        yield (Kmer(seq[i:i + self.k]) for i in xrange(
                len(seq) - self.k + 1))


    @staticmethod
    def legit_DNA_base_filter(seq):
        """
        For the current implementation, this assembly only solve the DNA sequence with ATCG included.
        The ambiguity base such as N will not be used
        :param seq: seq of the DNA
        :return: filter seq of the DNA with only ATCG and to uppercase
        """
        return "".join([base for base in seq.upper() if base in "ATCG"])

    def build_DBG(self):
        """
        This is the essential method to build the DBG.
        The build is assuming an exact matching rule -- referring to "perfect
        sequence". For real NGS project, such assumption does not hold due
        to reasons like sequencing error or polyploidy and the matching
        criteria has to be loosened/updated accordingly.
        One kmer: one edge and two nodes
        The logic here can be described as:
        1. if a l node (k-1 prefix of kmer) is in G, don't add it but get
        the node instead for later
        2. if not, create the node, but also get the node for later
        3. repeat above for the r node (k-1 suffix of kmer)
        4. Add an edge with both l and r node
        So each kmer has one edge and as a multi graph, it is allowed to
        have mroe than one edge between any two nodes
        Please be noted that two attributes start and end, were appended to
        each node added to the DBG -- these were to track the beginning and
        end of any Euler path if exists.

        This methods also handles the merge of nodes of the graph G into one
        new_node, and in the same time layout all the edges that pointed
        to or from one of these nodes.

        :param kmer: a kmer typed as Kmer class
        :return:
        """
        for seq in self.sequences:
            for kmer in self.convert_read_to_kmer(seq):
                L, R = None, None
                if kmer.l_node not in self.G:
                    self.G.add_node(kmer.l_node, start=False, end=False)
                L = kmer.l_node
                if kmer.r_node not in self.G:
                    self.G.add_node(kmer.r_node, start=False, end=False)
                R = kmer.r_node

                # then add the edge
                self.G.add_edge((L, R), str=kmer.kmer_str)

        # check if there is Euler Path exists in the graph
        if self.mark_euler_path():
            self.euler_path_exists = True
        else:
            self.euler_path_exists = False
        assert self.euler_path_exists is not None # should never fail this one


