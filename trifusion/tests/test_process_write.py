#!/usr/bin/python2

import os
import shutil
import unittest
from data_files import *

from trifusion.process.sequence import AlignmentList, Alignment
from trifusion.process.base import Base

x = Base()

sql_db = "sequencedb"


class ProcessWriteSinglesTest(unittest.TestCase):

    def setUp(self):

        aln_obj = AlignmentList([], sql_db=sql_db)
        self.con = aln_obj.con
        self.aln_obj = Alignment(dna_data_fas[0], sql_con=aln_obj.cur)

        if not os.path.exists("output"):
            os.makedirs("output")

        self.output_file = os.path.join("output", "test")

    def tearDown(self):

        self.aln_obj = None
        self.con.close()
        os.remove(sql_db)
        shutil.rmtree("output")

    def test_write_gphocs(self):

        self.aln_obj.write_to_file(["gphocs"],
                                   self.output_file)

    def test_write_mcmctree(self):

        self.aln_obj.write_to_file(["mcmctree"],
                                    self.output_file)

    def test_write_nexus(self):

        self.aln_obj.write_to_file(["nexus"],
                                   self.output_file)


class ProcessWriteTest(unittest.TestCase):

    def setUp(self):

        aln_obj = AlignmentList(dna_data_fas, sql_db=sql_db)
        self.con = aln_obj.con
        self.aln_obj = aln_obj.concatenate(alignment_name="test")
        os.makedirs("output")
        self.output_file = os.path.join("output", "test")

    def tearDown(self):

        shutil.rmtree("output")
        self.con.close()
        os.remove(sql_db)

    def test_write_fasta(self):

        self.aln_obj.write_to_file(["fasta"],
                                   self.output_file)

        self.assertEqual(x.autofinder(self.output_file + ".fas")[0],
                         "fasta")

    def test_write_fasta_interleave(self):

        self.aln_obj.write_to_file(["fasta"], self.output_file,
                                   interleave=True)

    def test_write_nexus(self):

        self.aln_obj.write_to_file(["nexus"],
                                   self.output_file)

        self.assertEqual(x.autofinder(self.output_file + ".nex")[0],
                         "nexus")

    def test_write_nexus_interleave(self):

        self.aln_obj.write_to_file(["nexus"], self.output_file,
                                   interleave=True)

    def test_write_mcmctree(self):

        self.aln_obj.write_to_file(["mcmctree"],
                                   self.output_file)

        self.assertEqual(x.autofinder(
            self.output_file + "_mcmctree.phy")[0], "phylip")

    def test_write_phy(self):

        self.aln_obj.write_to_file(["phylip"],
                                   self.output_file)

        self.assertEqual(x.autofinder(self.output_file + ".phy")[0],
                         "phylip")

    def test_write_stockholm(self):

        self.aln_obj.write_to_file(["stockholm"],
                                   self.output_file)

        self.assertEqual(x.autofinder(self.output_file + ".stockholm")[0],
                         "stockholm")

    def test_write_gphocs(self):

        self.aln_obj.write_to_file(["gphocs"],
                                   self.output_file)

    def test_write_ima2(self):

        ima2_params = [ima2_pop_file,
                       "(1,2):3)4:5",
                       "IS",
                       "1"]

        self.aln_obj.write_to_file(["ima2"], self.output_file,
                                   ima2_params=ima2_params)

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
