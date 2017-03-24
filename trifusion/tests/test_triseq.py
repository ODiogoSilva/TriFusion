#!/usr/bin/python

import os
import sys
import logging
import unittest
from os.path import join
from data_files import *
import shutil

from argparse import ArgumentTypeError, ArgumentError

try:
    from process.sequence import AlignmentList, Alignment
    from TriSeq import get_args, main, triseq_arg_check
    from base.sanity import triseq_arg_check
except ImportError:
    from trifusion.process.sequence import AlignmentList, Alignment
    from trifusion.TriSeq import get_args, main
    from trifusion.base.sanity import triseq_arg_check

output_dir = "triseq_test"
data_path = join("trifusion/tests/data/")


class TriSeqTest(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    def tearDown(self):
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

    def test_simple_concatenation(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["-of", "fasta",
                         "-o", join(output_dir, "teste"),
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        with open(join(data_path, "BaseConcatenation.fas")) as fh1, \
                open(join(output_dir, "teste.fas")) as fh2:
            self.assertEqual(fh1.read(), fh2.read())

    def test_simple_conversion(self):

        args = get_args(["-in", dna_data_fas[0],
                         "-of", "fasta",
                         "-c",
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        with open(join(data_path, "BaseConc1.fas")) as fh1, \
                open("BaseConc1.fas") as fh2:
            self.assertEqual(fh1.read().rstrip(), fh2.read().rstrip())

        os.remove("BaseConc1.fas")

    def test_reverse_concatenate(self):

        args = get_args(["-in", concatenated_small_phy[0],
                         "-r", concatenated_small_par[0],
                         "-of", "fasta",
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        exp = []
        data = []

        for fl in sorted([x for x in os.listdir(".")
                          if x.startswith("BaseConc")]):
            with open(fl) as fh:
                data.append(fh.read().rstrip())
            os.remove(fl)

        for fl in dna_data_fas:
            with open(fl) as fh:
                exp.append(fh.read().rstrip())

        self.assertEqual(exp, data)

    def test_convert_par2raxml(self):

        args = get_args(["-p", concatenated_small_par[0],
                         "-o", join(output_dir, "test"),
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        with open(concatenated_small_parNex[0]) as fh1, \
                open(join(output_dir, "test.charset")) as fh2:
            self.assertEqual(fh1.read(), fh2.read())

    def test_convert_raxml2nex(self):

        args = get_args(["-p", concatenated_small_parNex[0],
                         "-o", join(output_dir, "test"),
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        with open(concatenated_small_par[0]) as fh1, \
                open(join(output_dir, "test.part.File")) as fh2:
            self.assertEqual(fh1.read(), fh2.read().replace("LG", "DNA"))

    def test_convert_raxml2nex_with_model(self):

        args = get_args(["-p", concatenated_small_parNex[0],
                         "-o", join(output_dir, "test"),
                         "--model", "WAG",
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        with open(concatenated_small_par[0]) as fh1, \
                open(join(output_dir, "test.part.File")) as fh2:
            self.assertEqual(fh1.read(), fh2.read().replace("WAG", "DNA"))

    def test_select_alignments(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["-s", "spa",
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        fls = sorted(os.listdir("Taxa_selection"))

        shutil.rmtree("Taxa_selection")

        self.assertEqual(fls, ["BaseConc1.fas", "BaseConc7.fas"])

    def test_collapse(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--collapse",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        fls = sorted(os.listdir(output_dir))

        self.assertEqual(fls, ['teste.fas', 'teste.haplotypes'])

    def test_code_gaps(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--code-gaps",
                         "-o", join(output_dir, "teste"),
                         "-of", "nexus",
                         "-quiet"])
        triseq_arg_check(args)
        main(args)


    def test_consensus(self):

        args = get_args(["-in", dna_data_fas[0],
                        "--consensus", "Soft mask",
                         "-of", "fasta",
                         "-c",
                         "-quiet"])
        triseq_arg_check(args)
        main(args)

        with open("BaseConc1.fas") as fh:
            self.assertTrue(fh.read().startswith(">consensus"))

        os.remove("BaseConc1.fas")

    def test_consensus_singlefile(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--consensus", "Soft mask",
                         "--consensus-single-file",
                         "-of", "fasta",
                         "-c",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)

        self.assertTrue(os.path.exists("consensus.fas"))
        os.remove("consensus.fas")

    def test_missing_filter_perc_val(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--missing-filter", "10", "50",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_missing_filter_prop_val(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--missing-filter", "0.1", "0.5",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_min_taxa_perc_val(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--min-taxa", "50",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_min_taxa_prop_val(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--min-taxa", "0.5",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_contain_taxa(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--contain-taxa", "spa",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_exclude_taxa(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--exclude-taxa", "spa",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_codon_filter(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--codon-filter", "1", "2",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_variable_filter(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--variable-filter", "1", "2",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_informative_filter(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--informative-filter", "1", "2",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_interleave(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--interleave",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.fas")))

    def test_ima2_output(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--ima2-params", ima2_pop_file, "(1,2):3)4:5", "IS",
                         "1",
                         "-o", join(output_dir, "teste"),
                         "-of", "ima2",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        self.assertTrue(os.path.exists(join(output_dir, "teste.txt")))

    def test_remove_taxa(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["-rm", "spa",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        with open(join(output_dir, "teste.fas")) as fh:
            self.assertEqual(fh.read().count("spa"), 0)

    def test_grep_taxa(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["-grep", "spa",
                         "-o", join(output_dir, "teste"),
                         "-of", "fasta",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)
        with open(join(output_dir, "teste.fas")) as fh:
            self.assertEqual(fh.read().count("spa"), 1)

    def test_get_taxa(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["--get-taxa",
                         "-quiet"])

        triseq_arg_check(args)
        main(args)

        self.assertTrue(os.path.exists("Taxa_list.csv"))
        os.remove("Taxa_list.csv")


if __name__ == "__main__":
    unittest.main()
