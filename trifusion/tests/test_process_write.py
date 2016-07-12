#!/usr/bin/python2

import os
import shutil
import unittest
from data_files import *

from trifusion.process.sequence import AlignmentList
from trifusion.process.base import Base

x = Base()


class ProcessWriteTest(unittest.TestCase):

    def setUp(self):

        aln_obj = AlignmentList(dna_data_fas)
        self.aln_obj = aln_obj.concatenate(remove_temp=True,
                                           alignment_name="test",
                                           dest=".")
        aln_obj.clear_alignments()
        os.makedirs("output")
        self.output_file = os.path.join("output", "test")

    def tearDown(self):

        self.aln_obj._clear_alignment_temp()
        shutil.rmtree("test")
        shutil.rmtree("output")

    def test_write_fasta(self):

        self.aln_obj.write_to_file(["fasta"],
                                   self.output_file)

        if self.assertTrue(x.check_sizes(self.aln_obj.alignment,
                                         self.output_file + ".fas")):
            self.assertEqual(x.autofinder(self.output_file + ".fas")[0],
                             "fasta")

    def test_write_nexus(self):

        self.aln_obj.write_to_file(["nexus"],
                                   self.output_file)

        if self.assertTrue(x.check_sizes(self.aln_obj.alignment,
                                         self.output_file + ".nex")):
            self.assertEqual(x.autofinder(self.output_file + ".nex")[0],
                             "nexus")

    def test_write_mcmctree(self):

        self.aln_obj.write_to_file(["mcmctree"],
                                   self.output_file)

        if self.assertTrue(x.check_sizes(self.aln_obj.alignment,
                                         self.output_file + "_mcmctree.phy")):
            self.assertEqual(x.autofinder(
                self.output_file + "_mcmctree.phy")[0], "phylip")

    def test_write_phy(self):

        self.aln_obj.write_to_file(["phylip"],
                                   self.output_file)

        if self.assertTrue(x.check_sizes(self.aln_obj.alignment,
                                         self.output_file + ".phy")):
            self.assertEqual(x.autofinder(self.output_file + ".phy")[0],
                             "phylip")

    def test_write_stockholm(self):

        self.aln_obj.write_to_file(["stockholm"],
                                   self.output_file)

        if self.assertTrue(x.check_sizes(self.aln_obj.alignment,
                                         self.output_file + ".stockholm")):
            self.assertEqual(x.autofinder(self.output_file + ".stockholm")[0],
                             "stockholm")

    def test_write_gphocs(self):

        self.aln_obj.write_to_file(["gphocs"],
                                   self.output_file)

    def test_write_interleave(self):

        self.aln_obj.write_to_file(["phylip"],
                                   self.output_file,
                                   interleave=True)

    def test_write_gap(self):

        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   self.output_file,
                                   gap="?")

    def test_write_model_phy(self):

        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   self.output_file,
                                   model_phylip="LG")

    def test_write_outgoup_list(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   self.output_file,
                                   outgroup_list=["spa", "spb"])

    def test_write_use_charset(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   self.output_file,
                                   use_charset=False)

    def test_write_partition_file(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   self.output_file,
                                   partition_file=False)

    def test_write_output_dir(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   "teste",
                                   output_dir="test2")
        shutil.rmtree("test2")

    def test_write_ldhat(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   self.output_file,
                                   ld_hat=True)

if __name__ == "__main__":
    unittest.main()
