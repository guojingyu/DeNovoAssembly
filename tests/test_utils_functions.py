"""
Test major utils functions : file reading, reverse_complement

Author : Jingyu Guo
"""

import unittest
from de_novo_assembly.utils import get_sequences_from_fasta_file, include_reverse_complement


class FastaFileReadTests(unittest.TestCase):

    def setUp(self):
        self.fasta_file_1 = "../data/coding_challenge_data_set.fasta"
        self.fasta_record_dict = read_fasta_file(self.fasta_file_1)

    def test_number_of_reads(self):
        assert(len(self.fasta_record_dict.keys()) == 50,
               "Reading fasta: number of reads not matched 50: ", len(self.fasta_record_dict.keys()))

    def test_fasta_read_seq_ex1(self):
        # assert
        pass

    def test_fasta_read_seq_ex2(self):
        # assert
        pass

    def test_fasta_read_seq_ex2(self):
        # assert
        pass


class ReverseComplementReadTests(unittest.TestCase):

    def setUp(self):
        self.fasta_file_1 = "../data/coding_challenge_data_set.fasta"
        self.fasta_record_dict = read_fasta_file(self.fasta_file_1)

    def test_number_of_reads(self):
        assert(len(self.fasta_record_dict.keys()) == 50,
               "Reading fasta: number of reads not matched 50: ", len(self.fasta_record_dict.keys()))

    def test_fasta_read_seq_ex1(self):
        # assert
        pass

    def test_fasta_read_seq_ex2(self):
        # assert
        pass

    def test_fasta_read_seq_ex2(self):
        # assert
        pass

if __name__ == '__main__':
    unittest.main()