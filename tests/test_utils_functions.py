#!/usr/bin/env python

"""
Test major utils functions : file reading, reverse_complement

Author : Jingyu Guo
"""

import unittest
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from de_novo_assembly.utils import get_sequences_from_fasta_file, include_reverse_complement


class FastaFileReadTests(unittest.TestCase):

    def setUp(self):
        self.fasta_file_1 = "../data/coding_challenge_data_set.fasta"
        self.fasta_record_dict = get_sequences_from_fasta_file(self.fasta_file_1)

    def test_number_of_reads(self):
        """Reading fasta: test number of reads not matched 50"""
        assert len(self.fasta_record_dict.keys()) == 50

    def test_fasta_read_seq_1(self):
        """Reading fasta: test length read Rosalind_1836 not matched to ref
        length 994"""
        assert len(self.fasta_record_dict['Rosalind_1836'].seq) == 994

    def test_fasta_read_seq_2(self):
        """Reading fasta: test error in reads"""
        assert self.fasta_record_dict['Rosalind_1836'].seq[0:5] == 'GCGCC'

    def test_fasta_read_seq_3(self):
        """Reading fasta: test error in reads"""
        assert self.fasta_record_dict['Rosalind_1836'].seq[-5:] == 'AGAGG'


class ReverseComplementReadTests(unittest.TestCase):

    def setUp(self):
        self.fasta_file_1 = "../data/coding_challenge_data_set.fasta"
        self.sequence_1 = "ATCGGTCTGA"
        self.reverse_complement_1 = "TCAGACCGAT"

        self.sequence_2 = "A"
        self.reverse_complement_2 = "T"

        self.tests = {'test_1': Seq(self.sequence_1, generic_dna),
                      'test_2': Seq(self.sequence_2, generic_dna)}
        self.original_seq_record = len(self.tests.keys())

        # feed into include_reverse_complement
        self.included_reverse_complement = include_reverse_complement(self.tests)
        print self.tests.keys()

    def test_number_of_reads(self):
        """Reverse Complement: test output record number mismatched"""
        assert self.original_seq_record * 2 == len(self.included_reverse_complement.keys())

    def test_reverse_complement_read_1(self):
        """Reverse Complement: test sequence test_1 failed reverse complement."""
        assert str(self.included_reverse_complement['test_1_rev_com']) == self.reverse_complement_1

    def test_reverse_complement_read_2(self):
        """Reverse Complement: test sequence test_2 failed reverse complement."""
        assert str(self.included_reverse_complement['test_2_rev_com']) == self.reverse_complement_2


if __name__ == '__main__':
    unittest.main()