#!/usr/bin/python2

import os
import unittest
from data_files import *
from os.path import join
import shutil

try:
    from process.sequence import AlignmentList
    from process.error_handling import *
    from process.data import Partitions
except ImportError:
    from trifusion.process.sequence import AlignmentList
    from trifusion.process.error_handling import *
    from trifusion.process.data import Partitions

temp_dir = ".temp"
sql_db = ".temp/sequencedb"

data_path = join("trifusion/tests/data/")


class SeconaryOpsTest(unittest.TestCase):

    def setUp(self):

        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)

        self.aln_obj = AlignmentList(dna_data_fas, sql_db=sql_db)

    def tearDown(self):

        self.aln_obj.clear_alignments()
        self.aln_obj.con.close()
        shutil.rmtree(temp_dir)

    def test_summary_stats_all(self):

        sum_table, table_data = self.aln_obj.get_summary_stats()

        self.assertEqual([sum_table, table_data],
                         [{'missing': '5 (0.04%)', 'taxa': 24, 'genes': 7,
                           'informative': '0 (0.0%)', 'gaps': '0 (0.0%)',
                           'avg_gaps': 0.0, 'avg_missing': 1.0, 'variable': '7 (1.18%)',
                           'seq_len': 595, 'avg_var': 1.0, 'avg_inf': 0.0},
                          [['Genes', 'Taxa', 'Alignment length', 'Gaps',
                            'Gaps per gene', 'Missing data',
                            'Missing data per gene', 'Variable sites',
                            'Variable sites per gene', 'Informative sites',
                            'Informative sites per gene'],
                           [7, 24, 595, '0 (0.0%)', 0.0, '5 (0.04%)', 1.0,
                            '7 (1.18%)', 1.0, '0 (0.0%)', 0.0]]])

    def test_summary_stats_one_active(self):

        sum_table, table_data = self.aln_obj.get_summary_stats([
            join(data_path, "BaseConc1.fas")])

        self.assertEqual([sum_table, table_data],
                         [{'missing': '1 (0.05%)', 'taxa': 24, 'genes': 1,
                           'informative': '0 (0.0%)', 'gaps': '0 (0.0%)',
                           'avg_gaps': 0.0, 'avg_missing': 1.0, 'variable': '1 (1.18%)',
                           'seq_len': 85, 'avg_var': 1.0, 'avg_inf': 0.0},
                          [['Genes', 'Taxa', 'Alignment length', 'Gaps',
                            'Gaps per gene', 'Missing data',
                            'Missing data per gene', 'Variable sites',
                            'Variable sites per gene', 'Informative sites',
                            'Informative sites per gene'],
                           [1, 24, 85, '0 (0.0%)', 0.0, '1 (0.05%)', 1.0,
                            '1 (1.18%)', 1.0, '0 (0.0%)', 0.0]]])

    def test_single_aln_outlier_mdata(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.outlier_missing_data(),
                         {"exception": "single_alignment"})

    def test_single_aln_outlier_mdata_sp(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        print(self.aln_obj.alignments)

        self.assertEqual(self.aln_obj.outlier_missing_data_sp(),
                         {"exception": "single_alignment"})

    def test_single_aln_outlier_seg(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.outlier_segregating(),
                         {"exception": "single_alignment"})

    def test_single_aln_outlier_seg_sp(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        print(self.aln_obj.alignments)

        self.assertEqual(self.aln_obj.outlier_segregating_sp(),
                         {"exception": "single_alignment"})

    def test_single_aln_outlier_seqsize(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.outlier_sequence_size(),
                         {"exception": "single_alignment"})

    def test_single_aln_outlier_seqsize_sp(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.outlier_sequence_size_sp(),
                         {"exception": "single_alignment"})

    def test_single_aln_average_seqsize_per_species(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.average_seqsize_per_species(),
                         {"exception": "single_alignment"})

    def test_single_aln_average_seqsize(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.average_seqsize(),
                         {"exception": "single_alignment"})

    def test_single_aln_sequence_similarity(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.sequence_similarity(),
                         {"exception": "single_alignment"})

    def test_single_aln_sequence_segregation(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.sequence_segregation(),
                         {"exception": "single_alignment"})

    def test_single_aln_length_polymorphism_correlation(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.length_polymorphism_correlation(),
                         {"exception": "single_alignment"})

    def test_single_aln_taxa_distribution(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.taxa_distribution(),
                         {"exception": "single_alignment"})

    def test_single_aln_cumulative_missing_genes(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.cumulative_missing_genes(),
                         {"exception": "single_alignment"})

    def test_single_aln_gene_occupancy(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.gene_occupancy(),
                         {"exception": "single_alignment"})

    def test_single_aln_missing_data_distribution(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.missing_data_distribution(),
                         {"exception": "single_alignment"})

    def test_single_aln_missing_genes_average(self):

        self.aln_obj.update_active_alignments([dna_data_fas[0]])

        self.assertEqual(self.aln_obj.missing_genes_average(),
                         {"exception": "single_alignment"})

    def test_no_data(self):

        self.aln_obj = AlignmentList([], sql_db=sql_db)

        self.assertEqual(self.aln_obj.gene_occupancy(),
                         {'exception': "empty_data"})

    def test_gene_occupancy(self):

        self.assertTrue(self.aln_obj.gene_occupancy())

    def test_missing_data_distribution(self):

        self.assertTrue(self.aln_obj.missing_data_distribution())

    def test_missing_data_per_species(self):

        self.assertTrue(self.aln_obj.missing_data_per_species())

    def test_missing_genes_per_species(self):

        self.assertTrue(self.aln_obj.missing_genes_per_species())

    def test_missing_genes_average(self):

        self.assertTrue(self.aln_obj.missing_genes_average())

    def test_average_seqsize_per_species(self):

        self.assertTrue(self.aln_obj.average_seqsize_per_species())

    def test_average_seqsize(self):

        self.assertTrue(self.aln_obj.average_seqsize())

    def test_characters_proportion(self):

        self.assertTrue(self.aln_obj.characters_proportion())

    def test_characters_proportion_per_species(self):

        self.assertTrue(self.aln_obj.characters_proportion_per_species())

    def test_characters_proportion_gene(self):

        self.assertTrue(self.aln_obj.characters_proportion_gene(
            join(data_path, "BaseConc1.fas"), 10
        ))

    def test_sequence_similarity(self):

        self.assertTrue(self.aln_obj.sequence_similarity())

    def test_sequence_similarity_per_species(self):

        self.assertTrue(self.aln_obj.sequence_similarity_per_species())

    def test_sequence_similarity_gene(self):

        self.assertTrue(self.aln_obj.sequence_similarity_gene(
            join(data_path, "BaseConc1.fas"), 10))

    def test_sequence_conservation(self):

        self.assertTrue(self.aln_obj.sequence_conservation_gnp(
            join(data_path, "BaseConc1.fas"), 10
        ))

    def test_sequence_segregation(self):

        self.assertTrue(self.aln_obj.sequence_segregation())

    def test_sequence_segregation_per_species(self):

        self.assertTrue(self.aln_obj.sequence_segregation_per_species())

    def test_sequence_segregation_gene(self):

        self.assertTrue(self.aln_obj.sequence_segregation_gene(
            join(data_path, "BaseConc1.fas"), 10))

    def test_length_polymorphism_correlation(self):

        self.assertTrue(self.aln_obj.length_polymorphism_correlation())

    def test_allele_frequency_spectrum(self):

        self.assertTrue(self.aln_obj.allele_frequency_spectrum())

    def test_allele_frequency_spectrum_gene(self):

        self.assertTrue(self.aln_obj.allele_frequency_spectrum_gene(
            join(data_path, "BaseConc1.fas"), None))

    def test_taxa_distribution(self):

        self.assertTrue(self.aln_obj.taxa_distribution())

    def test_cumulative_missing_genes(self):

        self.assertTrue(self.aln_obj.cumulative_missing_genes())

    def test_outlier_missing_data(self):

        self.assertTrue(self.aln_obj.outlier_missing_data())

    def test_outlier_missing_data_sp(self):

        self.assertTrue(self.aln_obj.outlier_missing_data_sp())

    def test_outlier_segregating(self):

        self.assertTrue(self.aln_obj.outlier_segregating())

    def test_outlier_segregating_sp(self):

        self.assertTrue(self.aln_obj.outlier_segregating_sp())

    def test_outlier_sequence_size(self):

        self.assertTrue(self.aln_obj.outlier_sequence_size())

    def test_outlier_sequence_size_sp(self):

        self.assertTrue(self.aln_obj.outlier_sequence_size_sp())


if __name__ == "__main__":
    unittest.main()
