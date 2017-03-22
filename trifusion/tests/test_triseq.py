#!/usr/bin/python

import os
import sys
import logging
import unittest
from os.path import join
from data_files import *
import shutil

try:
    from process.sequence import AlignmentList, Alignment
    from TriSeq import get_args, main
except ImportError:
    from trifusion.process.sequence import AlignmentList, Alignment
    from trifusion.TriSeq import get_args, main

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

        main(args)

        with open(join(data_path, "BaseConcatenation.fas")) as fh1, \
                open(join(output_dir, "teste.fas")) as fh2:
            self.assertEqual(fh1.read(), fh2.read())

    def test_simple_conversion(self):

        args = get_args(["-in", dna_data_fas[0],
                        "-of", "fasta",
                         "-c",
                         "-quiet"])

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

        main(args)

        with open(concatenated_small_parNex[0]) as fh1, \
                open(join(output_dir, "test.charset")) as fh2:
            self.assertEqual(fh1.read(), fh2.read())

    def test_convert_raxml2nex(self):

        args = get_args(["-p", concatenated_small_parNex[0],
                         "-o", join(output_dir, "test"),
                         "-quiet"])

        main(args)

        with open(concatenated_small_par[0]) as fh1, \
                open(join(output_dir, "test.part.File")) as fh2:
            self.assertEqual(fh1.read(), fh2.read().replace("LG", "DNA"))

    def test_select_alignments(self):

        args = get_args(["-in"] + dna_data_fas +
                        ["-s", "spa",
                         "-quiet"])

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

        main(args)

        fls = sorted(os.listdir(output_dir))

        self.assertEqual(fls, ['teste.fas', 'teste.haplotypes'])


if __name__ == "__main__":
    unittest.main()
