#!/usr/bin/python2

import os
import shutil
import unittest
from data_files import *

from trifusion.process.sequence import AlignmentList, Alignment
from trifusion.process.base import Base

x = Base()

temp_dir = ".temp"
sql_db = ".temp/sequencedb"


class ProcessWriteSinglesTest(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList([dna_data_fas[0]], sql_db=sql_db)

        if not os.path.exists("output"):
            os.makedirs("output")

        self.output_file = os.path.join("output", "test")

    def tearDown(self):

        self.aln_obj.con.close()
        os.remove(sql_db)
        shutil.rmtree("output")
        shutil.rmtree(temp_dir)

    def test_write_gphocs(self):

        self.aln_obj.write_to_file(["gphocs"],
                                   output_file=self.output_file)

    def test_write_mcmctree(self):

        self.aln_obj.write_to_file(["mcmctree"],
                                   output_file=self.output_file)

    def test_write_nexus(self):

        self.aln_obj.write_to_file(["nexus"],
                                   output_file=self.output_file)


class ProcessWriteMultisTest(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList([dna_data_fas[0]], sql_db=sql_db)
        os.makedirs("output")
        self.output_dir = os.path.join("output")

    def tearDown(self):

        shutil.rmtree("output")
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

    def test_custom_taxaset_nexus(self):
        """
        Test explicitly for the head of the nexus, which should update the
        ntax parameter when changing the active taxa set
        """

        self.aln_obj.update_taxa_names(["spa", "spb", "spc"])

        self.aln_obj.write_to_file(["nexus"], output_dir=self.output_dir)

        # Get the specific line with the ntax parameter
        header_line = ""
        with open(os.path.join(self.output_dir, "BaseConc1.nex")) as fh:
            while not header_line.startswith("dimensions"):
                header_line = next(fh).strip()

        ref_header = "dimensions ntax=3 nchar=85 ;"
        self.assertEqual(header_line, ref_header)


class ProcessWriteTest(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList(dna_data_fas, sql_db=sql_db)
        self.aln_obj.concatenate()
        os.makedirs("output")
        self.output_file = os.path.join("output", "test")

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        shutil.rmtree("output")
        shutil.rmtree(temp_dir)

    def test_write_fasta(self):

        self.aln_obj.write_to_file(["fasta"],
                                   output_file=self.output_file)

        self.assertEqual(x.autofinder(self.output_file + ".fas")[0],
                         "fasta")

    def test_write_fasta_interleave(self):

        self.aln_obj.write_to_file(["fasta"],
                                   output_file=self.output_file,
                                   interleave=True)

    def test_write_nexus(self):

        self.aln_obj.write_to_file(["nexus"],
                                   output_file=self.output_file)

        self.assertEqual(x.autofinder(self.output_file + ".nex")[0],
                         "nexus")

    def test_write_nexus_interleave(self):

        self.aln_obj.write_to_file(["nexus"], output_file=self.output_file,
                                   interleave=True)

    def test_write_mcmctree(self):

        self.aln_obj.write_to_file(["mcmctree"],
                                   output_file=self.output_file)

        self.assertEqual(x.autofinder(
            self.output_file + "_mcmctree.phy")[0], "phylip")

    def test_write_phy(self):

        self.aln_obj.write_to_file(["phylip"],
                                   output_file=self.output_file)

        self.assertEqual(x.autofinder(self.output_file + ".phy")[0],
                         "phylip")

    def test_write_stockholm(self):

        self.aln_obj.write_to_file(["stockholm"],
                                   output_file=self.output_file)

        self.assertEqual(x.autofinder(self.output_file + ".stockholm")[0],
                         "stockholm")

    def test_write_gphocs(self):

        self.aln_obj.write_to_file(["gphocs"],
                                   output_file=self.output_file)

    def test_write_ima2(self):

        ima2_params = [ima2_pop_file,
                       "(1,2):3)4:5",
                       "IS",
                       "1"]

        self.aln_obj.write_to_file(["ima2"], output_file=self.output_file,
                                   ima2_params=ima2_params)

    def test_write_interleave(self):

        self.aln_obj.write_to_file(["phylip"],
                                   output_file=self.output_file,
                                   interleave=True)

    def test_write_upper_case_phy(self):

        self.aln_obj.write_to_file(["phylip"], output_file=self.output_file,
                                   upper_case=True)

        flag = True
        with open(self.output_file + ".phy") as fh:
            next(fh)
            for line in fh:
                seq = line.strip().split()[1]
                if not seq.isupper():
                    flag = False

        self.assertTrue(flag, True)

    def test_write_upper_case_fasta(self):

        self.aln_obj.write_to_file(["fasta"], output_file=self.output_file,
                                   upper_case=True)

        flag = True
        with open(self.output_file + ".fas") as fh:
            for line in fh:
                if line.startswith(">") or line.strip() == "":
                    continue
                else:
                    if not line.strip().isupper():
                        flag = False

        self.assertTrue(flag, True)

    def test_write_gap(self):

        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   output_file=self.output_file,
                                   gap="?")

    def test_write_model_phy(self):

        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   output_file=self.output_file,
                                   model_phylip="LG")

    def test_write_outgoup_list(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   output_file=self.output_file,
                                   outgroup_list=["spa", "spb"])

    def test_write_use_charset(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   output_file=self.output_file,
                                   use_charset=False)

    def test_write_partition_file(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   output_file=self.output_file,
                                   partition_file=False)

    def test_write_output_dir(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   output_dir="test2")
        shutil.rmtree("test2")

    def test_write_ldhat(self):
        self.aln_obj.write_to_file(["fasta", "phylip", "nexus", "mcmctree",
                                    "stockholm", "gphocs"],
                                   output_file=self.output_file,
                                   ld_hat=True)

    def test_write_snapp(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        os.remove(sql_db)
        self.aln_obj = AlignmentList(variable_data,  sql_db=sql_db)
        self.aln_obj.concatenate()
        self.aln_obj.write_to_file(["snapp"], output_file=self.output_file)

        with open(self.output_file + "_snapp.nex") as fh:
            res = sorted(fh.readlines())

        with open(snapp_output[0]) as fh:
            ref = sorted(fh.readlines())

        self.assertEqual(res, ref)

    def test_get_non_contiguous_partitions(self):

        self.aln_obj.partitions.merge_partitions(["BaseConc1.fas", "BaseConc3.fas",
                                       "BaseConc7.fas"], "non_contiguous")

        self.aln_obj.write_to_file(["mcmctree", "stockholm", "gphocs",
                                    "snapp"], output_file=self.output_file)

    def test_write_non_contiguous_partitions(self):

        self.aln_obj.partitions.merge_partitions(
            ["BaseConc1.fas", "BaseConc3.fas",
             "BaseConc7.fas"], "non_contiguous")

        self.aln_obj.write_to_file(["phylip", "nexus"],
                                   output_file=self.output_file)

if __name__ == "__main__":
    unittest.main()
