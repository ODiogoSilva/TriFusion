#!/usr/bin/env python2
#
#  Copyright 2012 Unknown <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.


import time
import argparse
import configparser
from process.sequence import *
from base.plotter import *
from TriSeq import CleanUp
from process.base import print_col, GREEN, RED, YELLOW

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Command line interface for "
                                     "TriFusion Statistics module")

    # Main execution
    main_exec = parser.add_argument_group("Main execution")
    main_exec.add_argument("-in", dest="infile", nargs="+",
                           help="Provide the input files.")
    main_exec.add_argument("-o", dest="project_name",
                           help="Name of the output directory")
    main_exec.add_argument("-cfg", dest="config_file",
                           help="Name of the configuration file with the "
                           "statistical analyses to be executed")
    main_exec.add_argument("--generate-cfg", dest="generate_cfg",
                           action="store_const", const=True,
                           help="Generates a configuration template file")

    arg = parser.parse_args()


def generate_cfg_template():

    template_fh = open("stats_template.ini", "w")

    template_fh.write("""
# Configuration template file for statistical analyses that can be passed to
# TriStats.py using the -cfg option. Section and option names mimic the names
# Present in the Statistics side panel.

# For each option, three values are possible but not necessarily all used.
# The allowed values are provided above each option.
# The values are:
#   - species: Creates a plot focusing on the taxa
#   - average: Creates a plot that averages over all alignments
#   - gene: Creates a plot for a single gene. Can only be used with single
#           input files

[General Information]
# Options available: species average
distribution_sequence_size: species average
# Options available: species average
proportion_nucleotides_residues: species average Ai
# Options available: average
distribution_taxa_frequency: average

[Polymorphism and Variation]
# Options available: gene species average
sequence_similarity: species average
# Options available: gene species average
segregating_sites: species average
# Options available: average
alignment_pol_correlation: average

[Missing Data]
# Options available: average
gene_occupancy: average
# Options available: species average
distribution_missing_genes: species average
# Options available: species average
distribution_missing_data: species average
# Options available: average
cumulative_distribution_missing_genes: average

[Outlier Detection]
# Options available: species average
missing_data_outliers: species average
# Options available: species average
segregating_sites_outliers: species average
    """)

    template_fh.close()


@CleanUp
def main():

    print_col("Executing TriStats module at %s %s" % (
        time.strftime("%d/%m/%Y"), time.strftime("%I:%M:%S")), GREEN, 2)

    if arg.generate_cfg:
        print_col("Generating configuration template file", GREEN, 2)
        return generate_cfg_template()

    # Arguments
    input_files = arg.infile
    output_dir = arg.project_name
    config_file = arg.config_file

    # Read configuration file
    print_col("Reading configuration file", GREEN, 2)
    settings = configparser.ConfigParser()
    settings.read(config_file)

    # Parse alignments
    print_col("Parsing %s alignments" % len(input_files), GREEN, 2)
    alignments = AlignmentList(input_files, dest=".tmp/")

    # Create output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Variable mapping each available option with the appropriate statistics
    # and plotting methods
    func_map = {
        ("General Information", "distribution_sequence_size", "species"):
            [alignments.average_seqsize_per_species,
             (box_plot, "avg_seqsize_species.png")],

        ("General Information", "distribution_sequence_size", "average"):
            [alignments.average_seqsize,
             (histogram_plot, "avg_seqsize.png")],

        ("General Information", "proportion_nucleotides_residues", "species"):
            [alignments.characters_proportion_per_species,
             (stacked_bar_plot, "char_proportions_sp.png")],

        ("General Information", "proportion_nucleotides_residues", "average"):
            [alignments.characters_proportion,
             (bar_plot, "char_proportions.png")],

        ("General Information", "distribution_taxa_frequency", "average"):
            [alignments.taxa_distribution,
             (histogram_plot, "distribution_taxa_frequency.png")],

        ("Polymorphism and Variation", "sequence_similarity", "species"):
            [alignments.sequence_similarity_per_species,
             (triangular_heat, "similarity_distribution_sp.png")],

        ("Polymorphism and Variation", "sequence_similarity", "average"):
            [alignments.sequence_similarity,
             (histogram_plot, "similarity_distribution.png")],

        ("Polymorphism and Variation", "sequence_similarity", "gene"):
            [alignments.sequence_similarity_gene,
             (sliding_window, "similarity_distribution_gn.png")],

        ("Polymorphism and Variation", "segregating_sites", "species"):
            [alignments.sequence_segregation_per_species,
             (triangular_heat, "segregating_sites_sp.png")],

        ("Polymorphism and Variation", "segregating_sites", "average"):
            [alignments.sequence_segregation,
             (histogram_plot, "segregating_sites.png")],

        ("Polymorphism and Variation", "segregating_sites", "gene"):
            [alignments.sequence_segregation_gene,
             (sliding_window, "segregating_sites_gn.png")],

        ("Polymorphism and Variation", "alignment_pol_correlation", "average"):
            [alignments.length_polymorphism_correlation,
             (scatter_plot, "length_polymorphism_correlation.png")],

        ("Missing Data", "gene_occupancy", "average"):
            [alignments.gene_occupancy,
             (interpolation_plot, "gene_occupancy.png")],

        ("Missing Data", "distribution_missing_genes", "species"):
            [alignments.missing_genes_per_species,
             (bar_plot, "missing_gene_distribution.png")],

        ("Missing Data", "distribution_missing_genes", "average"):
            [alignments.missing_genes_average,
             (histogram_plot, "missing_gene_distribution_avg.png")],

        ("Missing Data", "distribution_missing_data", "species"):
            [alignments.missing_data_per_species,
             (stacked_bar_plot, "missing_data_distribution_sp.png")],

        ("Missing Data", "distribution_missing_data", "average"):
            [alignments.missing_data_distribution,
             (histogram_smooth, "missing_data_distribution.png")],

        ("Missing Data", "cumulative_distribution_missing_genes", "average"):
            [alignments.cumulative_missing_genes,
             (bar_plot, "cumulative_distribution_missing_genes.png")],

        ("Outlier Detection", "missing_data_outliers", "species"):
            [alignments.outlier_missing_data_sp,
             (outlier_densisty_dist, "Missing_data_outliers_sp.png")],

        ("Outlier Detection", "missing_data_outliers", "average"):
            [alignments.outlier_missing_data,
             (outlier_densisty_dist, "Missing_data_outliers.png")],

        ("Outlier Detection", "segregating_sites_outliers", "species"):
            [alignments.outlier_segregating_sp,
             (outlier_densisty_dist, "Segregating_sites_outliers_sp.png")],

        ("Outlier Detection", "segregating_sites_outliers", "average"):
            [alignments.outlier_segregating,
             (outlier_densisty_dist, "Segregating_sites_outliers.png")]
    }

    print_col("Parsing configuation file options", GREEN, 2)

    # Iterate over each individual option
    for section in settings.sections():
        for option, val in settings.items(section):
            for i in val.split():

                # Check if current option is available or suported
                if (section, option, i) in func_map:
                    print_col("Generating plot for option: %s - %s - %s" %
                              (section, option, i), GREEN, 2)
                    # Get appropriate method list
                    funcs = func_map[(section, option, i)]
                    # Retrieve plot data using statistics method
                    plot_data = funcs[0]()
                    # Generate plot object
                    plot_obj, _, lgd = funcs[1][0](**plot_data)
                    plot_obj.tight_layout()

                    # Save plot to file, including the legend object, if
                    # available
                    if lgd:
                        plot_obj.savefig(join(output_dir, funcs[1][1]),
                                         bbox_extra_artists=(lgd,), dpi=200)
                    else:
                        plot_obj.savefig(join(output_dir, funcs[1][1]), dpi=200)
                else:
                    print_col("Invald option: %s - %s - %s. Skipping." %
                              (section, option, i), YELLOW, 2)

if __name__ == '__main__':
    main()


__author__ = "Diogo N. Silva"