"""
Test De Bruijn Graph build

Author : Jingyu Guo
"""

import unittest
from Bio.SeqRecord import SeqRecord
from de_novo_assembly.de_bruijn_graph import DeBruijnGraph, Kmer

class KmerTests(unittest.TestCase):
    def setUp(self):
        self.ref_seq_1_2mers = ["AT","TT","TA","AG","GA","AC","CC","CT",
                                    "TG"]
        self.ref_seq_1_3mers = ["ATT", "TTA", "TAG", "AGA", "GAC", "ACC",
                                    "CCT","CTG"]
        self.ref_seq_1_9mers = ["ATTAGACCT", "TTAGACCTG"]
        self.test_seq_1_2mers = [Kmer(kmer) for kmer in self.ref_seq_1_2mers]
        self.test_seq_1_3mers = [Kmer(kmer) for kmer in self.ref_seq_1_3mers]
        self.test_seq_1_9mers = [Kmer(kmer) for kmer in self.ref_seq_1_9mers]

    def test_kmer(self):
        assert self.test_seq_1_2mers[3].kmer_str == 'AG'
        assert self.test_seq_1_2mers[3].l_node == 'A'
        assert self.test_seq_1_2mers[3].r_node == 'G'
        assert self.test_seq_1_2mers[3].k == 2

        assert self.test_seq_1_3mers[5].kmer_str == 'ACC'
        assert self.test_seq_1_3mers[2].l_node == 'TA'
        assert self.test_seq_1_3mers[-1].r_node == 'TG'

        assert self.test_seq_1_9mers[0].kmer_str == 'ATTAGACCT'
        assert self.test_seq_1_9mers[1].l_node == 'TTAGACCT'
        assert self.test_seq_1_9mers[-1].r_node == 'TAGACCTG'


class DeBruijnGraphTests(unittest.TestCase):

    def setUp(self):
        self.sequence_1 = "ATTAGACCTG"

        self.sequence_2 = "ATTTAGACCCTG"
        self.sequence_3 = "AGACCCTGAGTCG"

        self.test_seq_1 = {'seq_1': SeqRecord(self.sequence_1)}
        self.test_seq_2 = {'seq_2': SeqRecord(self.sequence_2),
                           'seq_3': SeqRecord(self.sequence_3)}

        self.dbg_1 = DeBruijnGraph(self.test_seq_1,k=4)

        self.dbg_2 = DeBruijnGraph(self.test_seq_2,k=5)


    def test_graph_node_edge_number_1(self):
        """
        Test graph 1 nodes
        :return:
        """
        assert len(self.dbg_1.G.nodes()) == 8
        assert len(self.dbg_1.G.edges()) == 7 # as a path


    def test_has_euler_circuit(self):
        """
        outdegree of DBG node == indegree of DBG node if only Eulerian
        circuit exists.
        This is to ensure that the established DBG satisfy Eulerian graph
        requirement.
        :return:
        """
        # now link the dbg_1 to make a circle
        self.dbg_1.G.add_edge("CTG","ATT")
        for node in self.dbg_1.G.nodes_iter():
            if self.dbg_1.G.in_degree(node) != self.dbg_1.G.out_degree(node):
                assert False
        assert True


    def test_has_euler_path(self):
        """
        outdegree of DBG node == indegree of DBG node
        But allow two nodes that has 1 as in or out in case only Eulerian
        path exists.
        This is to ensure that the established DBG satisfy Eulerian graph
        requirement.
        :return:
        """
        counter_unbalanced_node = 0
        for node in self.dbg_1.G.nodes_iter():
            if self.dbg_1.G.in_degree(node) != self.dbg_1.G.out_degree(node):
                counter_unbalanced_node += 1
        if counter_unbalanced_node != 2:
            assert False
        assert True


if __name__ == '__main__':
    unittest.main()