#!/usr/bin/python2

import os
import shutil
import unittest
from os.path import join
from collections import OrderedDict
from data_files import *

try:
    from process.sequence import AlignmentList, Alignment
except ImportError:
    from trifusion.process.sequence import AlignmentList, Alignment

sql_db = "sequencedb"

data_path = join("trifusion/tests/data/")


def compare_inst(inst1, inst2, blacklist=None):
    """
    Compares attributes between two instancies and returns equality test
    :param inst1: instance 1
    :param inst2: instance 2
    :param blacklist: list of attributes that should be ignored in comparison
    """

    d1 = dict((x, y) for x, y in inst1.__dict__.items() if x not in blacklist)
    d2 = dict((x, y) for x, y in inst2.__dict__.items() if x not in blacklist)

    print(d1)
    print(d2)

    try:
        return d1 == d2
    except ValueError:
        return any(d1) == any(d2)


class LoadAlignmentsTest(unittest.TestCase):

    def setUp(self):

        self.aln_obj = AlignmentList([], sql_db=sql_db)

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        os.remove(sql_db)

    def test_class_instance(self):

        self.aln_obj = AlignmentList([], sql_db=sql_db)
        self.assertIsInstance(self.aln_obj.alignments, OrderedDict)

    def test_load_fas(self):

        self.aln_obj = AlignmentList(dna_data_fas, sql_db=sql_db)

    def test_load_single_fas(self):

        single_aln = Alignment(dna_data_fas[0], sql_cursor=self.aln_obj.cur)

    def test_load_phy(self):

        self.aln_obj = AlignmentList(dna_data_phy, sql_db=sql_db)

    def test_load_single_phy(self):

        single_aln = Alignment(dna_data_phy[0], sql_cursor=self.aln_obj.cur)

    def test_load_single_interleave_phy(self):

        single_aln = Alignment(phylip_interleave[0],
                               sql_cursor=self.aln_obj.cur)

    def test_load_nex(self):

        self.aln_obj = AlignmentList(dna_data_nex, sql_db=sql_db)

    def test_load_single_nex(self):

        single_aln = Alignment(dna_data_nex[0], sql_cursor=self.aln_obj.cur)

    def test_load_stc(self):

        self.aln_obj = AlignmentList(dna_data_stc, sql_db=sql_db)

    def test_load_single_stc(self):

        single_aln = Alignment(dna_data_stc[0], sql_cursor=self.aln_obj.cur)

    def test_load_loci(self):

        self.aln_obj = AlignmentList(dna_data_loci, sql_db=sql_db)

    def test_load_single_loci(self):

        single_aln = Alignment(dna_data_loci[0], sql_cursor=self.aln_obj.cur)

    def test_load_nexus_par(self):

        self.aln_obj = AlignmentList(concatenated_medium_nexus, sql_db=sql_db)
        self.assertTrue(self.aln_obj.partitions.partitions)

    def test_load_wrong_type(self):

        self.aln_obj = AlignmentList(bad_file, sql_db=sql_db)
        self.assertTrue(self.aln_obj.bad_alignments)

    def test_duplicate_files(self):

        self.aln_obj = AlignmentList(dna_data_loci + dna_data_loci,
                                     sql_db=sql_db)
        self.assertTrue(self.aln_obj.duplicate_alignments)

    def test_unequal_length(self):

        self.aln_obj = AlignmentList(unequal_file, sql_db=sql_db)
        self.assertTrue(self.aln_obj.non_alignments)

    def test_load_no_data(self):

        self.aln_obj = AlignmentList(no_data, sql_db=sql_db)


class AlignmentManipulationTest(unittest.TestCase):

    def setUp(self):

        self.aln_obj = AlignmentList(dna_data_fas, sql_db=sql_db)

    def tearDown(self):

        try:
            self.aln_obj.clear_alignments()
        except:
            pass
        self.aln_obj.con.close()
        os.remove(sql_db)

    def test_clear_alns(self):

        self.aln_obj.clear_alignments()
        aln = AlignmentList([], sql_db=sql_db)

        self.assertTrue(compare_inst(self.aln_obj, aln, ["log_progression",
                                                         "locus_length",
                                                         "partitions",
                                                         "cur",
                                                         "con"]))

    def test_update_act_anls(self):

        self.aln_obj.update_active_alignments([join(data_path,
                                                    "BaseConc1.fas"),
                                               join(data_path,
                                                    "BaseConc2.fas")])

        self.assertEqual(list(self.aln_obj.alignments.keys()),
                         [join(data_path, "BaseConc1.fas"),
                          join(data_path, "BaseConc2.fas")])

    def test_update_act_alns_err(self):

        self.aln_obj.update_active_alignments([join(data_path,
                                                    "BaseConc1.fas"),
                                               join(data_path,
                                                    "BaseConc2.fas"),
                                               join(data_path,
                                                    "Wrong_name")])

        self.assertEqual(list(self.aln_obj.alignments.keys()),
                         [join(data_path, "BaseConc1.fas"),
                          join(data_path, "BaseConc2.fas")])

    def test_update_aln_shelve(self):

        self.aln_obj.update_active_alignment(join(data_path, "BaseConc1.fas"),
                                             "shelve")

        self.assertEqual(list(self.aln_obj.alignments.keys()),
                         [join(data_path, "BaseConc2.fas"),
                          join(data_path, "BaseConc3.fas"),
                          join(data_path, "BaseConc4.fas"),
                          join(data_path, "BaseConc5.fas"),
                          join(data_path, "BaseConc6.fas"),
                          join(data_path, "BaseConc7.fas")])

    def test_update_aln_act(self):

        self.aln_obj.update_active_alignments([])
        self.aln_obj.update_active_alignment(join(data_path, "BaseConc1.fas"),
                                             "active")

        self.assertEqual(list(self.aln_obj.alignments.keys()),
                         [join(data_path, "BaseConc1.fas")])

    def test_add_aln_obj(self):

        fl = self.aln_obj.alignments.keys()

        aln = Alignment(dna_data_loci[0], sql_cursor=self.aln_obj.cur)

        self.aln_obj.add_alignments([aln])

        self.assertEqual(self.aln_obj.alignments.keys(),
                         fl + [join(data_path, "c97d5m4p2.loci")])

    def test_remove_taxa_from_list(self):

        taxa_list = [
            "1285_RAD_original",
            "130a_RAD_original",
            "137a_RAD_original",
            "1427_RAD_original",
            "167a_RAD_original"
        ]

        expected_taxa = [tx for tx in self.aln_obj.taxa_names if
                         tx not in taxa_list]

        self.aln_obj.remove_taxa(taxa_list)

        self.assertEqual(self.aln_obj.taxa_names, expected_taxa)

    def test_remove_taxa_from_file(self):

        taxa_list = [
            "1285_RAD_original",
            "130a_RAD_original",
            "137a_RAD_original",
            "1427_RAD_original",
            "167a_RAD_original"
        ]

        expected_taxa = [tx for tx in self.aln_obj.taxa_names if
                         tx not in taxa_list]

        self.aln_obj.remove_taxa(taxa_to_remove)

        self.assertEqual(self.aln_obj.taxa_names, expected_taxa)

    def test_remove_taxa_from_list_inverse(self):

        taxa_list = [
            "1285_RAD_original",
            "130a_RAD_original",
            "137a_RAD_original",
            "1427_RAD_original",
            "167a_RAD_original"
        ]

        expected_taxa = [tx for tx in self.aln_obj.taxa_names if
                         tx not in taxa_list]

        self.aln_obj.remove_taxa(taxa_list, mode="inverse")

        self.assertEqual(self.aln_obj.taxa_names, taxa_list)

    #
    # def test_retrieve_alignment(self):
    #
    #     aln = self.aln_obj.retrieve_alignment("BaseConc1.fas")
    #
    #     aln2 = Alignment(dna_data_fas[0], dest="new_one")
    #
    #     self.assertTrue(compare_inst(aln, aln2,
    #                                  ["log_progression", "locus_length",
    #                                   "partitions"]))

    def test_concatenation(self):

        aln = self.aln_obj.concatenate(alignment_name="test")
        aln.write_to_file(["fasta"], "test")

        with open("trifusion/tests/data/BaseConcatenation.fas") as fh1, \
                open("test.fas") as fh2:
            self.assertEqual(fh1.read(), fh2.read())

        os.remove("test.fas")

if __name__ == "__main__":
    unittest.main()
