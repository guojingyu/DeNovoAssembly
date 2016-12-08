#!/usr/bin/env python

"""
Test eulerian functions including the random walk

Author : Jingyu Guo
"""

import unittest
from de_novo_assembly.de_bruijn_graph import DeBruijnGraph, Kmer
import networkx
from de_novo_assembly.eulerian import has_euler_path, has_euler_circuit, make_contig_from_path


class EulerianTests(unittest.TestCase):
    def setUp(self):
        ## Build two DeBruijnGraph
        self.




    def test_has_euler_circuit(self):
        pass

    def test_has_euler_path(self):
        pass

    def test_eulerian_random_walk(self):
        pass





class PathToSequenceTests(unittest.TestCase):
    def setUp(self):
        self.path_1 = [('1','2'),('2','3'),('3','4')]
        self.seq_1 = '1234'
        self.path_2 = [('1','2'),('2','3'),('3','4'),('4','1')]
        self.seq_2 = '12341'
        self.path_3 = [('1', '2'), ('2', '1')]
        self.seq_3 = '121'
        self.path_4 = [('1','2'),('2','3'),('3','1')]
        self.seq_4 = '1231'
        self.path_4 = [('1', '2'), ('2', '3'), ('4', '1')] # path 4 is not
        # eulerian and would expect to generate an error

    def test_path_to_sequence(self):
        assert self.seq_1 == make_contig_from_path(self.path_1)
        assert self.seq_2 == make_contig_from_path(self.path_2)
        assert self.seq_3 == make_contig_from_path(self.path_3)

    def test_path_to_sequence_throw_error(self):
        try:

            _ = make_contig_from_path(self.path_4)
            print "If saw this line print out, something is wrong"
            assert False
        except ValueError as err:
            # "expected to be here"
            assert True

