#!/usr/bin/env python

"""
Test eulerian functions including the random walk

Author : Jingyu Guo
"""

import unittest
from de_novo_assembly.de_bruijn_graph import DeBruijnGraph
from Bio.SeqRecord import SeqRecord
from de_novo_assembly.eulerian import has_euler_path, has_euler_circuit, \
    make_contig_from_path, eulerian_random_walk


class EulerianTests(unittest.TestCase):
    def setUp(self):
        self.sequence_1 = "ATTAGACCTG"

        self.sequence_2 = "ATTTAGACCCTG"
        self.sequence_3 = "AGACCCTGAGTCG"

        self.test_seq_1 = {'seq_1': SeqRecord(self.sequence_1)}
        self.test_seq_2 = {'seq_2': SeqRecord(self.sequence_2),
                           'seq_3': SeqRecord(self.sequence_3)}

        self.dbg_1 = DeBruijnGraph(self.test_seq_1,k=4)
        # now link the dbg_1 to make a Eulerian circle
        self.dbg_1.G.add_edge("CTG", "ATT")

        self.dbg_2 = DeBruijnGraph(self.test_seq_2,k=6)
        self.reference_eulerian_path_2 = [('ATTT', 'TTTA'), ('TTTA', 'TTAG'),
                                          ('TTAG', 'TAGA'), ('TAGA', 'AGAC'),
                                          ('AGAC', 'GACC'), ('GACC', 'ACCC'),
                                          ('ACCC', 'CCCT'), ('CCCT', 'CCTG'),
                                          ('CCTG', 'CTGA'), ('CTGA', 'TGAG'),
                                          ('TGAG', 'GAGT'), ('GAGT', 'AGTC'),
                                          ('AGTC', 'GTCG')]

    def test_has_euler_path_function(self):
        assert has_euler_circuit(self.dbg_1.G) == True
        flag,_,_ = has_euler_path(self.dbg_1.G)
        assert flag == False

    def test_has_euler_circuit_function(self):
        assert has_euler_circuit(self.dbg_2.G) == False
        flag, _, _ = has_euler_path(self.dbg_2.G)
        assert flag == True

    def test_eulerian_random_walk(self):
        print eulerian_random_walk(self.dbg_2)


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


    def test_path_to_sequence(self):
        assert self.seq_1 == make_contig_from_path(self.path_1)
        assert self.seq_2 == make_contig_from_path(self.path_2)
        assert self.seq_3 == make_contig_from_path(self.path_3)


