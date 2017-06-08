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
#

import time
import argparse
from os.path import join
import os

try:
    from process.base import print_col, RED, GREEN, YELLOW
    from TriSeq import CleanUp
    from ortho import OrthomclToolbox as OT
    from ortho import protein2dna
except ImportError:
    from trifusion.process.base import print_col, RED, GREEN, YELLOW
    from trifusion.TriSeq import CleanUp
    from trifusion.ortho import OrthomclToolbox as OT
    from trifusion.ortho import protein2dna


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Command line interface for "
                                     "the TriFusion Orthology explore module.")

    # Main execution
    main_exec = parser.add_argument_group("Main execution")
    main_exec.add_argument("-in", dest="infile", nargs="*",
                           help="Provide the OrthoMCL group file(s). Multiple "
                           "group files are allowed but with restricted "
                           "access to some options.")
    main_exec.add_argument("-o", dest="output_dir", required=True,
                           help="Output directory.")

    # Filter options
    filter_opts = parser.add_argument_group("Filter options")
    filter_opts.add_argument("-t1", "--gene-threshold", dest="gn_threshold",
                             type=int, help="Provide the "
                             "maximum threshold of gene copy numbers per "
                             "species that should be allowed. If the option "
                             "is not specified all gene copy values will be "
                             "allowed.")
    filter_opts.add_argument("-t2", "--species-threshold", dest="sp_threshold",
                             help="Provide the minimum "
                             "number of species that should be allowed. If "
                             "the option is not specified, there will be no "
                             "lower limit.")

    # Conversion options
    convert_opts = parser.add_argument_group("Conversion options")
    convert_opts.add_argument("-g2f", dest="groups2fasta",
                              action="store_const", const=True,
                              help="Retrieve the"
                              " protein sequences of each ortholog cluster to a"
                              "single file per cluster. This retrieval is "
                              "performed after filtering ortholog using "
                              "options '-t1' and '-t2'. The protein database "
                              "used in the ortholog search must be provided "
                              "using the --protein-db.")
    convert_opts.add_argument("--protein-db", dest="protein_db",
                              help="Path to the protein database in fasta "
                              "format. Can be either a single file containing"
                              " sequences for all taxa or individual files "
                              "for each taxon.")
    convert_opts.add_argument("-p2d", dest="protein2dna", help="Convert "
                              "unaligned protein sequence files into their "
                              "corresponding nucleotide sequences. A CDS "
                              "database must be provided using the --cds-db "
                              "option and the protein files must be provided "
                              "using the --protein-db option. A nucleotide "
                              "sequence file will be generated for each "
                              "protein sequence file provided in --protein-db.")
    convert_opts.add_argument("--cds-db", dest="dna_db", help="Path to the "
                              "CDS database in fasta format. Can be either a "
                              "single file containing sequences for all taxa "
                              "or individual files for each taxon.")

    # Plotting options
    plot_opts = parser.add_argument_group("Plotting options")
    plot_opts.add_argument("-s", dest="plots", nargs="*",
                           choices=["1", "2", "3", "4", "5"], help="Use the "
                           "available choices to perform statistical analyses "
                           "on group files. All analyses are performed on the"
                           " provided group files, before applying any filters."
                           "Note that most analyses require a"
                           " single group file as input. Choice description:"
                           "\n\t1: Compares the number of total and filtered "
                           "ortholog clusters between multiple group files.\n"
                           "\t2: [Single group required] Distribution of the "
                           "number of taxa across orthologs.\n\t3: [Single "
                           "group required] Plots the available and missing "
                           "data proportions per taxon.\n\t4: [Single group "
                           "required] Plots the cumulative number of multiple "
                           "gene copies per taxon.\n\t5: [Single group "
                           "required] Distribution of gene copy numbers across"
                           "orthologs.")

    output_opts = parser.add_argument_group("Output options")
    output_opts.add_argument("-e", dest="export", action="store_const",
                             const=True, help="Exports the filtered groups "
                             "into a new file.")

    arg = parser.parse_args()


def main_check():
    """
    Performs sanity checks to argument combinations
    """

    if arg.protein2dna and arg.infile:
        print_col("Group file operations are ignored when specifying "
                  "conversion options.", YELLOW, 3)

    # Check if protein and cds data bases are provided when required
    if arg.groups2fasta and not arg.protein_db:
        print_col("A protein database must be provided to convert group "
                  "files into sequence files using the --protein-db option. "
                  "Exiting.", RED, 3)

    if arg.protein2dna and (not arg.protein_db and not arg.dna_db):
        print_col("A CDS data base and protein sequence files must be "
                  "provided to convert protein sequences into nucleotide "
                  "sequences using the --cds-db and --protein-db options, "
                  "respectively. Exiting.", RED, 3)

    # Print warnings when trying to execute options that are not available to
    # multiple input group files
    if len(arg.infile) > 1:
        if arg.groups2fasta or arg.protein2dna:
            print_col("Conversion options are only available for single "
                      "group files input.", YELLOW, 3)

    if arg.groups2fasta and not arg.gn_threshold and not arg.sp_threshold:
        print_col("No filters have been specified for the conversion of "
                  "group files into protein sequences. This may result in a "
                  "very large number of output files.", YELLOW, 3)


@CleanUp
def main():

    print_col("Executing TriOrtho module at %s %s" % (
        time.strftime("%d/%m/%Y"), time.strftime("%I:%M:%S")), GREEN, 3)

    # Create tmp dir
    os.makedirs(".tmp")

    # Arguments
    groups_file = arg.infile
    output_dir = arg.output_dir

    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if arg.protein2dna:
        print_col("Converting protein sequences into nucleotide sequences",
                  GREEN, 3)
        # Create database
        print_col("Creating database", GREEN, 3)
        id_db = protein2dna.create_db(arg.dna_db, ".tmp")
        # Create query for USEARCH
        print_col("Creating query", GREEN, 3)
        query_db = protein2dna.create_query(arg.protein_db, ".tmp")
        # Execute search
        print_col("Executing search", GREEN, 3)
        protein2dna.pair_search(".tmp")
        pair_db = protein2dna.get_pairs(".tmp")
        # Convert files
        print_col("Converting files", GREEN, 3)
        protein2dna.convert_protein_file(pair_db, query_db, id_db, output_dir)
        return print_col("Protein to nucleotide conversion complete", GREEN, 3)

    gene_threshold = arg.gn_threshold
    species_threshold = arg.sp_threshold
    protein_db = arg.protein_db

    if len(groups_file) == 1:

        print_col("Parsing group file", GREEN, 3)
        group_file = groups_file[0]
        group_object = OT.GroupLight(group_file, gene_threshold,
                                     species_threshold)

        # Check for plotting options
        if arg.plots:

            plt_methods = {"2": [group_object.bar_species_distribution,
                                 "Species distribution"],
                           "3": [group_object.bar_species_coverage,
                                 "Species data coverage"],
                           "4": [group_object.bar_genecopy_per_species,
                                 "Gene copies per species"],
                           "5": [group_object.bar_genecopy_distribution,
                                 "Gene copy distribution"]}

            for i in arg.plots:
                if i == "1":
                    print_col("Plotting option 1 requires multiple group "
                              "files. Skipping.", YELLOW, 3)
                    continue

                # Generate plot data and file
                print_col("Generating plot for %s" % plt_methods[i][1],
                          GREEN, 3)
                plot_obj, _, table = plt_methods[i][0](dest=output_dir)

        # Export filtered group file
        if arg.export:
            print_col("Exporting filtered group file using %s maximum gene "
                      "copies and %s minimum taxa representation" %
                      (gene_threshold, species_threshold), GREEN, 3)
            group_object.export_filtered_group(dest=output_dir)
            print_col("Filtering complete.\nTotal orthologs: %s;\nAfter gene "
                      "filter: %s;\nAfter species filter: %s;\nAfter both "
                      "filters: %s" % (len(group_object.species_frequency),
                                       group_object.num_gene_compliant,
                                       group_object.num_species_compliant,
                                       group_object.all_compliant), GREEN, 3)

        if arg.groups2fasta:
            print_col("Exporting group file as protein sequence files",
                      GREEN, 3)
            # Set sqlite file
            sqldb = join(".tmp", "group2protein.db")
            group_object.retrieve_sequences(sqldb, protein_db, output_dir)

    else:
        print_col("Parsing %s group files" % len(groups_file), GREEN, 3)
        multiple_groups_object = OT.MultiGroupsLight(".tmp",
                                                     groups_file,
                                                     gene_threshold,
                                                     species_threshold)

        if arg.plots:

            for i in arg.plots:

                if i != "1":
                    print_col("Plotting option %s requires a single group "
                              "file as input. Skipping." % str(i), YELLOW, 3)
                    continue

                print_col("Generating plot for Multiple group comparison",
                          GREEN, 3)
                multiple_groups_object.update_filters(gene_threshold,
                                                      species_threshold)
                multiple_groups_object.bar_orthologs(dest=output_dir)

        if arg.export:

            for gname, gobj in multiple_groups_object:

                gname = os.path.basename(gname)
                output_file = os.path.splitext(gname)[0] + "_filtered.txt"

                print_col("Exporting group file %s using %s maximum gene "
                          "copies and %s minimum taxa representation" %
                          (gname, gene_threshold, species_threshold), GREEN, 3)
                gobj.export_filtered_group(output_file_name=output_file,
                                           dest=output_dir)
                print_col("Filtering complete for group file %s.\nTotal "
                          "orthologs: %s;\nAfter gene filter: %s;\nAfter "
                          "species filter: %s;\nAfter both filters: %s" %
                          (gname,
                           len(gobj.species_frequency),
                           gobj.num_gene_compliant,
                           gobj.num_species_compliant,
                           gobj.all_compliant), GREEN, 3)

main_check()
main()


__author__ = "Diogo N. Silva"
