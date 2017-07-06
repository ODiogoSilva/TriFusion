#!/usr/bin/python2
# -*- coding: utf-8 -*-

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

temp_dir = ".temp"
sql_db = ".temp/sequencedb"

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

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList([], sql_db=sql_db)

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

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

    def test_load_interleave_nex(self):

        single_aln = Alignment(concatenated_interleave_nexus[0],
                               sql_cursor=self.aln_obj.cur)

    def test_load_stc(self):

        self.aln_obj = AlignmentList(dna_data_stc, sql_db=sql_db)

    def test_load_single_stc(self):

        single_aln = Alignment(dna_data_stc[0], sql_cursor=self.aln_obj.cur,
                               db_idx=self.aln_obj._idx + 1)

    def test_load_loci(self):

        self.aln_obj = AlignmentList(dna_data_loci, sql_db=sql_db)

    def test_load_single_loci(self):

        single_aln = Alignment(dna_data_loci[0], sql_cursor=self.aln_obj.cur,
                               db_idx=self.aln_obj._idx + 1)

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

    def test_alternative_missing(self):

        self.aln_obj = AlignmentList(alternative_missing, sql_db=sql_db)

        self.assertEqual(self.aln_obj.sequence_code[1], "?")

    def test_dna_missing_default(self):

        self.aln_obj = AlignmentList(single_dna, sql_db=sql_db)

        self.assertEqual(self.aln_obj.sequence_code[1], "n")

    def test_protein_missing_default(self):

        self.aln_obj = AlignmentList(protein_no_missing, sql_db=sql_db)

        self.assertEqual(self.aln_obj.sequence_code[1], "x")

    def test_dna_missing_eval(self):

        self.aln_obj = AlignmentList(concatenated_medium_nexus, sql_db=sql_db)

        self.assertEqual(self.aln_obj.sequence_code[1], "n")

    def test_protein_missing_eval(self):

        self.aln_obj = AlignmentList(protein_normal_missing, sql_db=sql_db)

        self.assertEqual(self.aln_obj.sequence_code[1], "x")

    def test_non_ascii_taxon_names(self):

        self.aln_obj = AlignmentList(non_ascii, sql_db=sql_db)

        non_ascii_tx = [x for x in self.aln_obj.taxa_names if
                        x == '\xc3\xa9!"#$%&/=?\'~\xc2\xba\xc2\xaa^"><']
        self.assertEqual(len(non_ascii_tx), 1)

    def test_non_ascii_iteration(self):

        self.aln_obj = AlignmentList(non_ascii, sql_db=sql_db)

        non_ascii_tx = []

        for tx, _, _ in self.aln_obj.iter_alignments():

            if tx == '\xc3\xa9!"#$%&/=?\'~\xc2\xba\xc2\xaa^"><':
                non_ascii_tx.append(tx)

        self.assertEqual(len(non_ascii_tx), 1)

    def test_non_ascii_get_taxaidx(self):

        self.aln_obj = AlignmentList(non_ascii, sql_db=sql_db)

        aln = self.aln_obj.alignments.values()[0]

        non_ascii_tx = [x for x in aln.taxa_idx
                        if x == '\xc3\xa9!"#$%&/=?\'~\xc2\xba\xc2\xaa^"><']

        self.assertEqual(len(non_ascii_tx), 1)

    def test_non_ascii_iter_columns(self):

        self.aln_obj = AlignmentList(non_ascii, sql_db=sql_db)

        tx_list, _, _ = next(self.aln_obj.iter_columns(include_taxa=True))

        non_ascii_tx = [x for x in tx_list
                        if x == '\xc3\xa9!"#$%&/=?\'~\xc2\xba\xc2\xaa^"><']

        self.assertEqual(len(non_ascii_tx), 1)


class AlignmentManipulationTest(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList(dna_data_fas, sql_db=sql_db)

    def tearDown(self):

        try:
            self.aln_obj.clear_alignments()
        except:
            pass
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

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

        aln = Alignment(dna_data_loci[0], sql_cursor=self.aln_obj.cur,
                        sql_con=self.aln_obj.con,
                        db_idx=self.aln_obj._idx + 1, temp_dir=temp_dir)

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
    #                                   "_partitions"]))

    def test_concatenation(self):

        self.aln_obj.concatenate()
        self.aln_obj.write_to_file(["fasta"], output_file="test")

        with open("trifusion/tests/data/BaseConcatenation.fas") as fh1, \
                open("test.fas") as fh2:
            self.assertEqual(sorted(fh1.readlines()), sorted(fh2.readlines()))

        os.remove("test.fas")


class LoadBadAlignmentsTest(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = None

    def test_bad_extra_space_nexus_interleave(self):

        self.aln_obj = AlignmentList(bad_extraspace_interleave, sql_db=sql_db)

        self.assertEqual(self.aln_obj.non_alignments,
                         bad_extraspace_interleave)

    def test_no_final_colon_interleave(self):

        self.aln_obj = AlignmentList(bad_no_colon_interleave, sql_db=sql_db)

        aln_obj = self.aln_obj.alignments.values()[0]
        data = [aln_obj.name,
                aln_obj.locus_length,
                len(aln_obj.taxa_idx)]

        self.assertEqual(data, ["bad_no_colon_interleave.nex",
                                898,
                                12])

    def test_no_end_colon_interleave(self):

        self.aln_obj = AlignmentList(bad_no_end_interleave, sql_db=sql_db)

        aln_obj = self.aln_obj.alignments.values()[0]
        data = [aln_obj.name,
                aln_obj.locus_length,
                len(aln_obj.taxa_idx)]

        self.assertEqual(data, ["bad_no_end_interleave.nex",
                                898,
                                12])

    def test_bad_no_colon_nexus(self):

        self.aln_obj = AlignmentList(bad_no_colon, sql_db=sql_db)

        aln_obj = self.aln_obj.alignments.values()[0]
        data = [aln_obj.name,
                aln_obj.locus_length,
                len(aln_obj.taxa_idx)]

        self.assertEqual(data, ["bad_no_colon.nex",
                                898,
                                12])

    def test_bad_no_end_nexus(self):

        self.aln_obj = AlignmentList(bad_no_end, sql_db=sql_db)

        aln_obj = self.aln_obj.alignments.values()[0]
        data = [aln_obj.name,
                aln_obj.locus_length,
                len(aln_obj.taxa_idx)]

        self.assertEqual(data, ["bad_no_end.nex",
                                898,
                                12])

    def test_bad_no_header_nexus(self):

        self.aln_obj = AlignmentList(bad_no_header, sql_db=sql_db)

        self.assertEqual(self.aln_obj.bad_alignments, bad_no_header)

    def test_bad_no_matrix_nexus(self):

        self.aln_obj = AlignmentList(bad_no_matrix, sql_db=sql_db)

        self.assertEqual(self.aln_obj.bad_alignments, bad_no_matrix)

    def test_bad_no_format_line_nexus(self):

        self.aln_obj = AlignmentList(bad_no_format_line, sql_db=sql_db)

        self.assertEqual(self.aln_obj.bad_alignments, bad_no_format_line)

    def test_bad_space_in_middle_nexus(self):

        self.aln_obj = AlignmentList(bad_space_in_middle, sql_db=sql_db)

        aln_obj = self.aln_obj.alignments.values()[0]
        data = [aln_obj.name,
                aln_obj.locus_length,
                len(aln_obj.taxa_idx)]

        self.assertEqual(data, ["bad_space_in_middle.nex",
                                898,
                                12])

    def test_bad_wrong_dimensions_nexus(self):

        self.aln_obj = AlignmentList(bad_wrong_dimensions, sql_db=sql_db)

        self.assertEqual(self.aln_obj.bad_alignments, bad_wrong_dimensions)

    def test_bad_wrong_size_nexus(self):

        self.aln_obj = AlignmentList(bad_wrong_size, sql_db=sql_db)

        aln_obj = self.aln_obj.alignments.values()[0]
        data = [aln_obj.name,
                aln_obj.locus_length,
                len(aln_obj.taxa_idx)]

        self.assertEqual(data, ["bad_wrong_size.nex",
                                898,
                                12])

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

if __name__ == "__main__":
    unittest.main()
