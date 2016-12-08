"""
Test De Bruijn Graph build

Author : Jingyu Guo
"""

import unittest
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
        self.ref_seq_1_2mers = ["AT", "TT", "TA", "AG", "GA", "AC", "CC", "CT",
                                "TG"]
        self.ref_seq_1_3mers = ["ATT", "TTA", "TAG", "AGA", "GAC", "ACC",
                                "CCT", "CTG"]
        self.ref_seq_1_9mers = ["ATTAGACCT", "TTAGACCTG"]

        self.sequence_2 = "ATTTAGACCCTG"
        self.ref_seq_2_2mers = ["AT","TT","TA","AG","GA","AC","CC","CT",
                                    "TG"]
        self.ref_seq_2_3mers = ["ATT", "TTA", "TAG", "AGA", "GAC", "ACC",
                                    "CCT","CTG"]
        self.ref_seq_2_9mers = ["ATTAGACCT", "TTAGACCTG"]

        #self.dbg = DeBruijnGraph()


    def test_balanced_graph(self):
        """
        outdegree of DBG node == indegree of DBG node
        But allow two nodes that has 1 or each
        This is to ensure that the established DBG satisfy Eulerian graph
        requirement.
        :return:
        """
        pass


if __name__ == '__main__':
    unittest.main()