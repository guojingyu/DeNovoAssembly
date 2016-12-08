#!/usr/bin/env python

"""
Test sample dummy data

Author : Jingyu Guo
"""

import unittest
from de_novo_assembly.run import run

class DummySampleTests(unittest.TestCase):
    """

    """
    def setUp(self):
        self.reference_assembly = "ATTAGACCTGCCGGAATAC"

        self.inputfastafile = "../data/dummy_data.fasta"
        self.outputfile = "../dummy_sample_test_output.txt"
        self.kmerlength = 7
        self.print_to_console = False

        self.assembly = run(self.inputfastafile, self.outputfile,
                            self.kmerlength, self.print_to_console)

    def test_dummy_assembly_result(self):
        assert self.assembly[0] == self.reference_assembly

if __name__ == '__main__':
    unittest.main()