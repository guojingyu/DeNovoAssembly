"""
Test eulerian functions including the random walk

Author : Jingyu Guo
"""

import unittest
from de_novo_assembly.de_bruijn_graph import DeBruijnGraph, Kmer
import networkx
from de_novo_assembly.eulerian import has_euler_path, has_euler_circuit, \
    eulerian_random_walk, make_contig_from_path


class EulerianTests(unittest.TestCase):
    pass