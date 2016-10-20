#!/usr/bin/python2

import os
import shutil
import unittest
from data_files import *

from trifusion.process.sequence import AlignmentList
from trifusion.process.data import Partitions, Zorro


class SeconaryOpsTest(unittest.TestCase):

    def setUp(self):

        self.aln_obj = AlignmentList([])

    def tearDown(self):

        self.aln_obj.clear_alignments()

    def test_collapse_single(self):

        self.aln_obj.add_alignment_files([variable_data[0]])

        self.aln_obj.collapse(haplotype_name="Testing",
                              haplotypes_file="teste")

        aln = self.aln_obj.alignments.values()[0]

        self.assertEqual(len(aln.alignment), 1)
        aln._clear_alignment_temp()

    def test_collapse_with_variation(self):

        self.aln_obj.add_alignment_files([variable_data[1]])

        self.aln_obj.collapse(haplotype_name="Testing",
                              haplotypes_file="teste")

        aln = self.aln_obj.alignments.values()[0]

        self.assertEqual(len(aln.alignment), 4)
        aln._clear_alignment_temp()

    def test_collapse_after_concatenation(self):

        self.aln_obj.add_alignment_files(variable_data)

        aln = self.aln_obj.concatenate(remove_temp=True, alignment_name="test",
                                       dest=".")
        aln.collapse(haplotype_name="Testing", haplotypes_file="teste",
                     dest="test")

        self.assertEqual(len(aln.alignment), 7)
        aln._clear_alignment_temp()
        shutil.rmtree("test")

    def test_gcoder(self):

        self.aln_obj.add_alignment_files(gcode_data)

        self.aln_obj.code_gaps()

        s = []
        for aln in self.aln_obj:
            for seq in aln.sequences():
                s.append(seq)

        res = [
            "aaaaaaaa-aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa10000",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa---aaaaaaaaaaa01000",
            "aaaaaaaaaaaa--aaaaaaaaaaaaaaaaaaaaa---aaaaaaaaaaa01100",
            "aaaaaaaaaaaa--aaaaaaaaaaaaaaaaaaaaa---aaaaaaaaaaa01100",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa---aaaaaaaaaaa01000",
            "aaaaaaaaaaaaaaaaaaaaaa----aaaaaaaaa---aaaaaaaaaaa01010",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa-aaa---aaaaaaaaaaa01001",
            "aaaaaaaaaaaaaaaaaaaaaa----aaaaaaaaa---aaaaaaaaaaa01010",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa---aaaaaaaaaaa01000",
            "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa---aaaaaaaaaaa01000"
        ]

        self.assertEqual(s, res)

    def test_consensus_multi_file(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.consensus("IUPAC")

        s = []
        for aln in self.aln_obj:
            s.append(len(aln.alignment))

        self.assertEqual(s, [1] * 7)

    def test_consensus_single_file(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        aln = self.aln_obj.consensus("IUPAC", True)

        self.assertEqual(len(aln.alignment), 7)

    def test_consensus_soft_mask(self):

        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.consensus("Soft mask")

        s = []
        for aln in self.aln_obj:
            s.append(len(aln.alignment))

        self.assertEqual(s, [1] * 7)

    def test_consensus_remove(self):
        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.consensus("Remove")

        s = []
        for aln in self.aln_obj:
            s.append(len(aln.alignment))

        self.assertEqual(s, [1] * 7)

    def test_consensus_first_seq(self):
        self.aln_obj.add_alignment_files(dna_data_fas)

        self.aln_obj.consensus("First sequence")

        s = []
        for aln in self.aln_obj:
            s.append(len(aln.alignment))

        self.assertEqual(s, [1] * 7)

    def test_reverse_concatenate(self):

        self.aln_obj.add_alignment_files(concatenated_small_phy)

        partition_obj = Partitions()
        # In case the partitions file is badly formatted or invalid, the
        # exception will be returned by the read_from_file method.
        partition_obj.read_from_file(concatenated_small_par[0])
        aln = self.aln_obj.concatenate(remove_temp=True, alignment_name="test",
                                       dest=".")

        aln.set_partitions(partition_obj)

        alns = aln.reverse_concatenate(dest="test_conc")

        self.assertEqual(len(alns.alignments), 7)

        shutil.rmtree("test_conc")

    def test_zorro(self):

        self.aln_obj.add_alignment_files(zorro_data_fas)

        # Generate zorro output
        zorro_data = Zorro(self.aln_obj, "_zorro")
        zorro_data.write_to_file("test")

        # Read zorro and reference files
        zorro_content = open("test_zorro.out").read()
        reference = open(zorro_out).read()

        self.assertEqual(zorro_content, reference)

        os.remove("test_zorro.out")

    def test_zorro_with_dir(self):

        self.aln_obj.add_alignment_files(zorro_data_fas)

        # Generate zorro output
        zorro_data = Zorro(self.aln_obj, "_zorro", "trifusion/tests/data/")
        zorro_data.write_to_file("test")

        # Read zorro and reference files
        zorro_content = open("test_zorro.out").read()
        reference = open(zorro_out).read()

        self.assertEqual(zorro_content, reference)

        os.remove("test_zorro.out")


# class MultipleSeconaryOpsTest(unittest.TestCase):
#
#     def test_sequential_secondary_operations_concat(self):
#
#         aln_obj = AlignmentList(dna_data_fas)
#
#         aln_obj.filter_min_taxa(0)
#
#         aln_obj.filter_by_taxa("Exclude", "NaN")
#
#         aln_obj.filter_codon_positions([True, True, True])
#
#         aln_obj.filter_missing_data(100, 100)
#
#         aln_obj.filter_segregating_sites(None, None)
#
#         aln_obj.filter_informative_sites(None, None)
#
#         aln_obj = aln_obj.concatenate(remove_temp=True, alignment_name="test",
#                                        dest=".")
#
#         aln_obj.collapse()
#
#         aln_obj.consensus("IUPAC")
#
#         shutil.rmtree("test_conc")

if __name__ == "__main__":
    unittest.main()
