#!/usr/bin/python2

import os
import shutil
import unittest
from collections import OrderedDict
from data_files import *

try:
    from process.sequence import AlignmentList, Alignment
except ImportError:
    from trifusion.process.sequence import AlignmentList, Alignment


def compare_inst(inst1, inst2, blacklist=None):
    """
    Compares attributes between two instancies and returns equality test
    :param inst1: instance 1
    :param inst2: instance 2
    :param blacklist: list of attributes that should be ignored in comparison
    """

    d1 = dict((x, y) for x, y in inst1.__dict__.items() if x not in blacklist)
    d2 = dict((x, y) for x, y in inst2.__dict__.items() if x not in blacklist)

    return d1 == d2


class LoadAlignmentsTest(unittest.TestCase):

    def setUp(self):

        self.aln_obj = AlignmentList([])

    def tearDown(self):

        self.aln_obj.clear_alignments()

    def test_class_instance(self):

        self.aln_obj = AlignmentList([])
        self.assertIsInstance(self.aln_obj.alignments, OrderedDict)

    def test_load_fas(self):

        self.aln_obj = AlignmentList(dna_data_fas)

    def test_load_phy(self):

        self.aln_obj = AlignmentList(dna_data_phy)

    def test_load_nex(self):

        self.aln_obj = AlignmentList(dna_data_nex)

    def test_load_stc(self):

        self.aln_obj = AlignmentList(dna_data_stc)

    def test_load_loci(self):

        self.aln_obj = AlignmentList(dna_data_loci)

    def test_load_nexus_par(self):

        self.aln_obj = AlignmentList(concatenated_medium_nexus)
        self.assertTrue(self.aln_obj.partitions.partitions)

    def test_load_wrong_type(self):

        self.aln_obj = AlignmentList(bad_file)
        self.assertTrue(self.aln_obj.bad_alignments)
        shutil.rmtree("bad_file")

    def test_duplicate_files(self):

        self.aln_obj = AlignmentList(dna_data_loci + dna_data_loci)
        self.assertTrue(self.aln_obj.duplicate_alignments)

    def test_unequal_length(self):

        self.aln_obj = AlignmentList(unequal_file)
        self.assertTrue(self.aln_obj.non_alignments)
        shutil.rmtree("unequal_length")

    def test_load_no_data(self):

        self.aln_obj = AlignmentList(no_data)
        shutil.rmtree("no_data")


class AlignmentManipulationTest(unittest.TestCase):

    def setUp(self):

        self.aln_obj = AlignmentList(dna_data_fas)

    def tearDown(self):

        self.aln_obj.clear_alignments()

    def test_clear_alns(self):

        self.aln_obj.clear_alignments()
        aln = AlignmentList([])

        self.assertTrue(compare_inst(self.aln_obj, aln, ["log_progression",
                                                         "locus_length",
                                                         "partitions"]))

    def test_update_act_anls(self):

        self.aln_obj.update_active_alignments(["BaseConc1.fas",
                                               "BaseConc2.fas"])

        self.assertEqual(list(self.aln_obj.alignments.keys()),
                         ["BaseConc1.fas", "BaseConc2.fas"])

    def test_update_act_alns_err(self):

        self.aln_obj.update_active_alignments(["BaseConc1.fas",
                                               "BaseConc2.fas",
                                               "Wrong_name"])

        self.assertEqual(list(self.aln_obj.alignments.keys()),
                         ["BaseConc1.fas", "BaseConc2.fas"])

    def test_update_aln_shelve(self):

        self.aln_obj.update_active_alignment("BaseConc1.fas", "shelve")

        self.assertEqual(list(self.aln_obj.alignments.keys()),
                         ["BaseConc2.fas", "BaseConc3.fas", "BaseConc4.fas",
                          "BaseConc5.fas", "BaseConc6.fas", "BaseConc7.fas"])

    def test_update_aln_act(self):

        self.aln_obj.update_active_alignments([])
        self.aln_obj.update_active_alignment("BaseConc1.fas", "active")

        self.assertEqual(list(self.aln_obj.alignments.keys()),
                         ["BaseConc1.fas"])

    def test_add_aln_obj(self):

        fl = self.aln_obj.alignments.keys()

        aln = Alignment(dna_data_loci[0])

        self.aln_obj.add_alignments([aln])

        self.assertEqual(self.aln_obj.alignments.keys(),
                         fl + ["c97d5m4p2.loci"])

    def test_retrieve_alignment(self):

        aln = self.aln_obj.retrieve_alignment("BaseConc1.fas")

        aln2 = Alignment(dna_data_fas[0])

        self.assertTrue(compare_inst(aln, aln2,
                                     ["log_progression", "locus_length",
                                      "partitions"]))

    def test_concatenation(self):

        aln = self.aln_obj.concatenate(remove_temp=True, alignment_name="test",
                                       dest=".")
        aln.write_to_file(["fasta"], "test")

        with open("trifusion/tests/data/BaseConcatenation.fas") as fh1, \
                open("test.fas") as fh2:
            self.assertEqual(fh1.read(), fh2.read())

        os.remove("test.fas")
        shutil.rmtree("test")


if __name__ == "__main__":
    unittest.main()
