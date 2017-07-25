#!/usr/bin/python2

import shutil
import unittest
from collections import OrderedDict
from os.path import join
import os
from data_files import *

try:
    from process.sequence import AlignmentList
    from process.error_handling import *
except ImportError:
    from trifusion.process.sequence import AlignmentList
    from trifusion.process.error_handling import *

temp_dir = ".temp"
sql_db = ".temp/sequencedb"

data_path = join("trifusion/tests/data/")


class AlignmentMissingFiltersTest(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList([], sql_db=sql_db)

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

    def test_filter_default(self):

        self.aln_obj.add_alignment_files(
            ["trifusion/tests/data/missing_data.phy",
             "trifusion/tests/data/missing_data2.phy"]
        )
        self.aln_obj.filter_missing_data(25, 50)

        s = []
        for aln in self.aln_obj:
            s.append(aln.locus_length)

        self.assertEqual(s, [42, 43])

    def test_filter_and_concat(self):

        self.aln_obj.add_alignment_files(
            ["trifusion/tests/data/missing_data.phy",
             "trifusion/tests/data/missing_data2.phy"]
        )

        self.aln_obj.filter_missing_data(25, 50, table_out="master_out")

        self.aln_obj.concatenate(table_in="master_out")

        self.assertEqual(self.aln_obj.size, 85)

    def test_no_filters(self):

        self.aln_obj.add_alignment_files(
            ["trifusion/tests/data/missing_data.phy",
             "trifusion/tests/data/missing_data2.phy"]
        )

        self.aln_obj.filter_missing_data(100, 100)

        s = []
        for aln in self.aln_obj:
            s.append(aln.locus_length)

        self.assertEqual(s, [50, 50])

    def test_no_missing(self):

        self.aln_obj.add_alignment_files(
            ["trifusion/tests/data/missing_data.phy",
             "trifusion/tests/data/missing_data2.phy"]
        )
        self.aln_obj.filter_missing_data(0, 0)

        s = []
        for aln in self.aln_obj:
            s.append(aln.locus_length)

        self.assertEqual(s, [0, 19])

    def test_no_data_aln_default_filters(self):

        self.aln_obj.add_alignment_files(
            ["trifusion/tests/data/missing_data3.phy"]
        )

        self.aln_obj.filter_missing_data(25, 50)

        s = None
        for aln in self.aln_obj:
            s = aln.locus_length

        self.assertEqual(s, 0)

    def test_no_data_aln_no_filters(self):

        self.aln_obj.add_alignment_files(
            ["trifusion/tests/data/missing_data3.phy"]
        )

        self.aln_obj.filter_missing_data(100, 100)

        s = None
        for aln in self.aln_obj:
            s = aln.locus_length

        self.assertEqual(s, 50)


class AlignmentTaxaFilters(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList([], sql_db=sql_db)

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

    def test_filter_min_taxa(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.filter_min_taxa(50)

        self.assertEqual(len(self.aln_obj.alignments), 5)

    def test_filter_min_taxa_max(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.filter_min_taxa(100)

        self.assertEqual(len(self.aln_obj.alignments), 1)

    def test_filter_min_taxa_min(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.filter_min_taxa(0)

        self.assertEqual(len(self.aln_obj.alignments), 7)

    def test_filter_by_taxa_include(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.filter_by_taxa(["spa", "spb", "spc", "spd"], "Contain")

        self.assertEqual(len(self.aln_obj.alignments), 2)

    def test_filter_by_taxa_exclude(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.filter_by_taxa(["spa", "spb", "spc", "spd"], "Exclude")

        self.assertEqual(len(self.aln_obj.alignments), 5)

    def test_filter_by_taxa_all(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.filter_by_taxa(["no_taxa"], "Contain")

        self.assertEqual(len(self.aln_obj.alignments), 0)

    def test_filter_by_taxa_from_file(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.filter_by_taxa("trifusion/tests/data/filter_taxa.txt",
                                    "Contain")

        self.assertEqual(len(self.aln_obj.alignments), 2)


class AlignmentCodonFilters(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList([],sql_db=sql_db)

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

    def test_codon_filter_pos1(self):

        self.aln_obj.add_alignment_files(codon_filter)

        self.aln_obj.filter_codon_positions([True, False, False],
                                            table_out="master_out")

        s = []
        for _, seq, _ in self.aln_obj.iter_alignments("master_out"):
            s.append(seq)

        self.assertEqual(s, ["a" * 16] * 10)

    def test_codon_filter_pos2(self):

        self.aln_obj.add_alignment_files(codon_filter)

        self.aln_obj.filter_codon_positions([False, True, False],
                                            table_out="master_out")

        s = []
        for _, seq, _ in self.aln_obj.iter_alignments("master_out"):
            s.append(seq)

        self.assertEqual(s, ["t" * 16] * 10)

    def test_codon_filter_pos3(self):

        self.aln_obj.add_alignment_files(codon_filter)

        self.aln_obj.filter_codon_positions([False, False, True],
                                            table_out="master_out")

        s = []
        for _, seq, _ in self.aln_obj.iter_alignments("master_out"):
            s.append(seq)

        self.assertEqual(s, ["g" * 16] * 10)

    def test_codon_filter_pos12(self):

        self.aln_obj.add_alignment_files(codon_filter)

        self.aln_obj.filter_codon_positions([True, True, False],
                                            table_out="master_out")

        s = []
        for _, seq, _ in self.aln_obj.iter_alignments("master_out"):
            s.append(seq)

        self.assertEqual(s, ["at" * 16] * 10)

    def test_codon_filter_pos13(self):

        self.aln_obj.add_alignment_files(codon_filter)

        self.aln_obj.filter_codon_positions([True, False, True],
                                            table_out="master_out")

        s = []
        for _, seq, _ in self.aln_obj.iter_alignments("master_out"):
            s.append(seq)

        self.assertEqual(s, ["ag" * 16] * 10)

    def test_codon_filter_all(self):

        self.aln_obj.add_alignment_files(codon_filter)

        self.aln_obj.filter_codon_positions([True, True, True])

        s = []
        for _, seq, _ in self.aln_obj.iter_alignments():
            s.append(seq)

        self.assertEqual(s, ["atg" * 16] * 10)


class AlignmentVariationFilters(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList([], sql_db=sql_db)

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

    def test_variation_filter_min(self):

        self.aln_obj.add_alignment_files(variable_data)

        self.aln_obj.filter_segregating_sites(None, None)

        self.assertEqual(len(self.aln_obj.alignments), 3)

    def test_variation_var_sites(self):

        self.aln_obj.add_alignment_files(variable_data)

        self.aln_obj.filter_segregating_sites(1, 2)

        self.assertEqual(len(self.aln_obj.alignments), 0)

    def test_variation_var_sites2(self):
        self.aln_obj.add_alignment_files(variable_data)

        self.aln_obj.filter_segregating_sites(1, 3)

        self.assertEqual(len(self.aln_obj.alignments), 1)

    def test_variation_inf_min(self):
        self.aln_obj.add_alignment_files(variable_data)

        self.aln_obj.filter_informative_sites(None, None)

        self.assertEqual(len(self.aln_obj.alignments), 3)

    def test_variation_inf_sites(self):

        self.aln_obj.add_alignment_files(variable_data)

        self.aln_obj.filter_informative_sites(1, 4)

        self.assertEqual(len(self.aln_obj.alignments), 1)

    def test_variation_inf_sites2(self):

        self.aln_obj.add_alignment_files(variable_data)

        self.aln_obj.filter_informative_sites(1, 1)

        print(self.aln_obj.alignments)

        self.assertEqual(len(self.aln_obj.alignments), 1)


if __name__ == "__main__":
    unittest.main()
