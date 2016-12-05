"""
Test file reading

Author : Jingyu Guo
"""

import unittest
from de_novo_assembly.de_bruijn_graph import DeBruijnGraph, KmerNode


class DeBruijnGraphTests(unittest.TestCase):

    def setUp(self):
        self.sequence_1 = "ATTAGACCTG"
        self.ref_seq_1_2mers = ["AT","TT","TA","AG","GA","AC","CC","CT",
                                    "TG"]
        self.ref_seq_1_3mers = ["ATT", "TTA", "TAG", "AGA", "GAC", "ACC",
                                    "CCT","CTG"]
        self.ref_seq_1_9mers = ["ATTAGACCT", "TTAGACCTG"]
        self.test_seq_1_2mers = KmerNode.create_kmer(self.sequence_1,k=2)
        self.test_seq_1_3mers = KmerNode.create_kmer(self.sequence_1, k=2)
        self.test_seq_1_9mers = KmerNode.create_kmer(self.sequence_1, k=2)


    def test_balanced_graph(self):
        """
        outdegree of DBG == indegree of DBG
        This is to ensure that the established DBG satisfy Eulerian graph
        requirement.
        :return:
        """
        assert(len(self.fasta_record_dict.keys()) == 50,
               "Reading fasta: number of reads not matched 50: ", len(self.fasta_record_dict.keys()))

    def test_2mer(self):
        # assert
        pass

    def test_3mer(self):
        # assert
        pass

    def test_9mer(self):
        # assert
        pass


if __name__ == '__main__':
    unittest.main()