#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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

import warnings

# Suppress import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    import os
    import sys
    import shutil
    import time
    import argparse
    from glob import glob

    try:
        from process.base import print_col, RED, GREEN, YELLOW, CleanUp
        from process import sequence as seqset
        from process import data
        from process.error_handling import *
        from base.sanity import triseq_arg_check, mfilters, post_aln_checks, \
            check_infile_list
        from progressbar import ProgressBar, Timer, Bar, Percentage, \
            SimpleProgress
    except ImportError:
        from trifusion.process.base import print_col, RED, GREEN, YELLOW,\
            CleanUp
        from trifusion.process import sequence as seqset
        from trifusion.process import data
        from trifusion.process.error_handling import *
        from trifusion.base.sanity import triseq_arg_check, mfilters, \
            post_aln_checks, check_infile_list
        from trifusion.progressbar import ProgressBar, Timer, Bar,\
            Percentage, SimpleProgress


def gen_wgt(msg):

    bar_wdg = [
        "( ", SimpleProgress(), " ) ",
        Bar(),
        Percentage(),
        " [", Timer(), "] ",
    ]

    return bar_wdg

@CleanUp
def main_parser(arg, alignment_list):
    """ Function with the main operations of TriSeq """

    print_col("Executing TriSeq module at %s %s" % (
        time.strftime("%d/%m/%Y"), time.strftime("%I:%M:%S")), GREEN,
              quiet=arg.quiet)

    # Create temp directory
    tmp_dir = ".trifusion-temp"
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # Set path to temporary sqlite database
    sql_db = os.path.join(tmp_dir, "trifusion.db")

    # If database already exists, erase it. Make sure we start fresh.
    if os.path.exists(sql_db):
        os.remove(sql_db)

    # Defining main variables
    conversion = arg.conversion
    output_format = arg.output_format
    outfile = arg.outfile
    interleave = arg.interleave
    model_phy = arg.model_phy
    # outgroup_taxa = arg.outgroup_taxa

    # Defining output file name
    if conversion is None and arg.outfile is not None:
        outfile = "".join(arg.outfile)
    elif conversion is None and arg.outfile is not None:
        outfile = "".join(arg.outfile)
    elif arg.consensus and arg.consensus_single and not arg.outfile:
        outfile = "consensus"

    # The input file at this stage is not necessary
    # If just converting the partition file format do this and exit
    if arg.partition_file is not None:
        # Initializing Partitions instance and reading partitions file
        partition = data.Partitions()
        partition.read_from_file(arg.partition_file)
        if partition.partition_format == "nexus":
            partition.write_to_file("raxml", outfile, model_phy)
        else:
            partition.write_to_file("nexus", outfile)
        return 0

    # Support wildcars as arguments for windows
    fl = []
    if sys.platform in ["win32", "cygwin"]:
        for p in alignment_list:
            fl += glob(p)
        alignment_list = fl

    # Check input files for directories
    alignment_list, dirs, lost = check_infile_list(alignment_list)

    if dirs:
        print_col("Ignoring input files pointing to a directory: {}".format(
            " ".join(dirs)), YELLOW)
    if lost:
        print_col("Ignoring input files that do not exist: {}".format(
            " ".join(lost)), YELLOW)
    if not alignment_list:
        print_col("No valid input files have been provided. Terminating...",
                  RED)

    # Input alignments are mandatory from now on
    if not arg.quiet:
        pbar = ProgressBar(max_value=len(alignment_list), widgets=gen_wgt(""))
    else:
        pbar = None

    print_col("Parsing %s alignments" % len(alignment_list), GREEN,
              quiet=arg.quiet)
    alignments = seqset.AlignmentList(alignment_list, sql_db=sql_db,
                                      pbar=pbar)

    post_aln_checks(arg, alignments)

    # ################################ Utilities ##############################
    # Return a file with taxa list and exit
    if arg.get_taxa is True:
        print_col("Writing taxa to new file", GREEN, quiet=arg.quiet)
        alignments.write_taxa_to_file()
        return 0

    # Remove taxa
    if arg.remove:
        print_col("Removing taxa", GREEN, quiet=arg.quiet)
        alignments.remove_taxa(arg.remove)

    # Grep taxa
    if arg.grep:
        print_col("Grepping taxa", GREEN, quiet=arg.quiet)
        alignments.remove_taxa(arg.grep, mode="inverse")

    # Select alignments
    if arg.select:
        print_col("Selecting alignments", GREEN, quiet=arg.quiet)
        if not os.path.exists("Taxa_selection"):
            os.makedirs("Taxa_selection")

        # Check if any of the provided taxa is absent from the alignments
        absent_taxa = [x for x in arg.select if x not in alignments.taxa_names]
        if absent_taxa:
            print_col("The following taxa were not found in any alignment and"
                      " will be ignored: {}".format(" ".join(absent_taxa)),
                      YELLOW, quiet=arg.quiet)

        selected_alignments = alignments.select_by_taxa(arg.select,
                                                        mode="relaxed")
        for aln in selected_alignments:
            alignment_file = aln.path
            shutil.copy(alignment_file, "Taxa_selection")

        return

    # ############################# Main operations ###########################
    # Reverse concatenation
    if arg.reverse is not None:
        print_col("Reverse concatenating", GREEN, quiet=arg.quiet)
        if len(alignment_list) > 1:
            raise ArgumentError("Only one input file allowed for reverse "
                                "concatenation")
        if arg.reverse:
            aln = alignments.alignments.values()[0]
            # Initializing and reading partition file
            partition = data.Partitions()
            partition.read_from_file(arg.reverse[0])
            # Updating alignment _partitions
            aln.set_partitions(partition)
            alignments.set_partition_from_alignment(aln, reset=True)

        alignments.reverse_concatenate(pbar=pbar)

    # Filtering
    # Filter by minimum taxa
    if arg.min_taxa:
        print_col("Filtering by minimum taxa", GREEN, quiet=arg.quiet)
        alignments.filter_min_taxa(arg.min_taxa, pbar=pbar)

    # Filter by alignments that contain taxa
    if arg.contain_filter:
        print_col("Filtering alignment(s) including a taxa group", GREEN,
                  quiet=arg.quiet)
        alignments.filter_by_taxa(arg.contain_filter, "Contain", pbar=pbar)

    # Filter by alignments that exclude taxa
    if arg.exclude_filter:
        print_col("Filtering alignments excluding a taxa group", GREEN,
                  quiet=arg.quiet)
        alignments.filter_by_taxa(arg.exclude_filter, "Exclude", pbar=pbar)

    # Filter by codon position
    if arg.codon_filter:
        print_col("Filtering by codon positions", GREEN, quiet=arg.quiet)
        if alignments.sequence_code[0] == "DNA":
            codon_settings = [True if str(x) in arg.codon_filter else False
                              for x in range(1, 4)]
            alignments.filter_codon_positions(codon_settings, pbar=pbar)

    # Filter by missing data
    if arg.m_filter:
        print_col("Filtering by missing data", GREEN, quiet=arg.quiet)
        alignments.filter_missing_data(arg.m_filter[0], arg.m_filter[1],
                                       pbar=pbar, use_main_table=True)

    # Filtering by variable sites
    if arg.var_filter:
        print_col("Filtering by variable sites", GREEN, quiet=arg.quiet)
        alignments.filter_segregating_sites(arg.var_filter[0],
                                            arg.var_filter[1],
                                            pbar=pbar)

    # Filtering by informative sites
    if arg.inf_filter:
        print_col("Filtering by variable sites", GREEN, quiet=arg.quiet)
        alignments.filter_informative_sites(arg.inf_filter[0],
                                            arg.inf_filter[1],
                                            pbar=pbar)

    # Concatenation
    if not arg.conversion and arg.reverse is None and not arg.consensus:
        print_col("Concatenating", GREEN, quiet=arg.quiet)
        alignments.concatenate(pbar=pbar)

        # Concatenate zorro files
        if arg.zorro:
            zorro = data.Zorro(alignment_list, arg.zorro)
            zorro.write_to_file(outfile)

    # Collapsing
    if arg.collapse:
        print_col("Collapsing", GREEN, quiet=arg.quiet)
        alignments.collapse(use_main_table=True, pbar=pbar,
                            haplotypes_file=outfile)

    # Gcoder
    if arg.gcoder:
        print_col("Coding gaps", GREEN, quiet=arg.quiet)
        if output_format == ["nexus"]:
            alignments.code_gaps(use_main_table=True, pbar=pbar)

    # Consensus
    if arg.consensus:
        consensus_type = arg.consensus[0]
        print_col("Creating consensus sequences", GREEN, quiet=arg.quiet)
        alignments.consensus(consensus_type, single_file=arg.consensus_single,
                             pbar=pbar)

    # Write output
    print_col("Writing output", GREEN, quiet=arg.quiet)
    alignments.write_to_file(output_format, output_file=outfile,
                             interleave=interleave,
                             ima2_params=arg.ima2_params,
                             partition_file=True,
                             use_charset=True,
                             pbar=pbar)


def get_args(arg_list=None, unittest=False):

    # The inclusion of the argument definition in main, makes it possible to
    # import this file as a module and not triggering argparse. The
    # alternative of using a if __name__ == "__main__" statement does not
    # work well with the entry_points parameter of setup.py, since they call
    # the main function but do nothing inside said statement.
    parser = argparse.ArgumentParser(description="Command line interface for "
                                                 "TriFusion Process module")

    # Main execution
    main_exec = parser.add_argument_group("Main execution")
    main_exec.add_argument("-in", dest="infile", nargs="+", help="Provide the "
                           "input file name. If multiple files are provided"
                           ", please separated the names with spaces")
    main_exec.add_argument("-of", dest="output_format", nargs="+",
                           default=["nexus"],
                           choices=["nexus", "phylip", "fasta", "mcmctree",
                                    "ima2", "stockholm", "gphocs"],
                           help="Format of the output file(s). You may select "
                           "multiple output formats simultaneously (default is"
                           " '%(default)s')")
    main_exec.add_argument("-o", dest="outfile", help="Name of the output file")

    # Alternative modes
    alternative = parser.add_argument_group("Alternative execution modes")
    alternative.add_argument("-c", dest="conversion", action="store_const",
                             const=True, help="Used for conversion of the "
                             "input files passed as arguments with the -in "
                             "option. This flag precludes the usage of the -o "
                             "option, as the output file name is automatically "
                             "generated based on the input file name")
    alternative.add_argument("-r", dest="reverse", help="Reverse a concatenated"
                             " file into its original single locus alignments."
                             " A partition file similar to the one read by "
                             "RAxML must be provided", nargs="*")
    alternative.add_argument("-z", "--zorro-suffix", dest="zorro", type=str,
                             help="Use this option if you wish to concatenate "
                             "auxiliary Zorro files associated with each "
                             "alignment. Provide the sufix for the concatenated"
                             " zorro file")
    alternative.add_argument("-p", "--partition-file", dest="partition_file",
                             type=str, help="Using this option and providing "
                             "the partition file will convert it between a "
                             "RAxML or Nexus format")
    alternative.add_argument("-s", dest="select", nargs="*", help="Selects "
                             "alignments containing the provided taxa "
                             "(separate multiple taxa with whitespace)")
    alternative.add_argument("--collapse", dest="collapse",
                             action="store_const",
                             const=True, default=False, help="Use this flag if "
                             "you would like to collapse the input alignment(s)"
                             " into unique haplotypes")
    alternative.add_argument("--code-gaps", dest="gcoder", action="store_const",
                             const=True, default=False, help="Use this flag to "
                             "code the gaps of the alignment into a binary "
                             "state  matrix that is appended to the end of "
                             "the alignment")
    alternative.add_argument("--consensus", dest="consensus", nargs=1,
                             choices=["First sequence", "IUPAC", "Soft mask",
                                      "Remove"],
                             help="Creates a consensus of the final alignments"
                             " specifying how variation is handled")
    alternative.add_argument("--consensus-single-file",
                             dest="consensus_single", action="store_const",
                             const=True, default=False, help="Merges "
                             "consensus sequences in a single file")

    filter_g = parser.add_argument_group("Filter options")
    filter_g.add_argument("--missing-filter", dest="m_filter", nargs=2,
                          help="Use this  option if you wish to filter the"
                          " alignment's missing data. Along with this "
                          "option provide the threshold percentages for "
                          "gap and missing data, respectively (e.g. -filter "
                          "50 75 - filters alignments columns with more "
                          "than 50%% of gap+missing data and columns with "
                          "more than 75%% of true missing data)",
                          type=mfilters)
    filter_g.add_argument("--min-taxa", dest="min_taxa",
                          type=mfilters, help="Set the minimum percentage "
                          "of taxa that needs to be present in an alignment")
    filter_g.add_argument("--contain-taxa", dest="contain_filter",
                          nargs="*", help="Only "
                          "processes alignments that contain the specified "
                          "taxa")
    filter_g.add_argument("--exclude-taxa", dest="exclude_filter", help="Only "
                          "process alignments that do NOT contain the "
                          "specified taxa")
    filter_g.add_argument("--codon-filter", dest="codon_filter", nargs="*",
                          choices=["1", "2", "3"], help="Include only the "
                          "codon positions specified by this option (DNA only)")
    filter_g.add_argument("--variable-filter", dest="var_filter", nargs=2,
                          help="Provide minimum and maximum values of "
                          "variable sites for each alignment. Filters "
                          "alignments with a number of variable sites outside "
                          "the specified range", type=int)
    filter_g.add_argument("--informative-filter", dest="inf_filter", nargs=2,
                          help="Provide minimum and maximum values of "
                          "informative sites for each alignment. Filters "
                          "alignments with a number of informative sites "
                          "outside the specified range", type=int)

    # Formatting options
    formatting = parser.add_argument_group("Formatting options")
    formatting.add_argument("--model", dest="model_phy", default="LG",
                            choices=["DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG",
                                     "RTREV", "CPREV", "VT", "BLOSUM62",
                                     "MTMAM", "LG"],
                            help="This option only applies for the "
                            "concatenation  of protein data into phylip "
                            "format. Specify the model for all partitions "
                            "defined in the partition file  (default is '%("
                            "default)s')")
    formatting.add_argument("--interleave", dest="interleave",
                            action="store_const", const="interleave",
                            default=False, help="Specify this  option to "
                            "write output files in interleave format (currently"
                            " only supported for nexus files")
    formatting.add_argument("--ima2-params", dest="ima2_params", nargs=4,
                            help="Provide 4 additional arguments needed to "
                            "write the output in a format compliant with "
                            "IMa2. The order of the required arguments "
                            "(separated by whitespace is as follows: "
                            "[(1) File name of population mapping]"
                            "[(2) Population tree]"
                            "[(3) Mutational model]"
                            "[(4) Inheritance Scalar]. "
                            "Additional notes: (1) The  population mapping "
                            "file is a simple .csv file  containing two "
                            "columns separated by a semi-colon, "
                            "in which the first column contains the taxon "
                            "name and the second column contains the "
                            "corresponding population name; (2) The order of "
                            "the population names in the population tree "
                            "must be the same as the order in the file with "
                            "the population mapping")

    # Data manipulation
    manipulation = parser.add_argument_group("Data manipulation")
    manipulation.add_argument("-rm", dest="remove", nargs="*", help="Removes "
                              "the  specified taxa from the final alignment. "
                              "Unwanted taxa my be provided in a csv file "
                              "containing 1 column with a species name in "
                              "each line or they may be specified in the "
                              "command line and separated by whitespace")
    manipulation.add_argument("-grep", dest="grep", nargs="*", help="The "
                              "inverse  of the -rm command. It removes all "
                              "taxa from the alignment except for the ones "
                              "specified with this option. Taxa names may be "
                              "specified in a csv file containing 1 column "
                              "with a species name in each line or in the "
                              "command line separated by whitespace")
    # manipulation.add_argument("-outgroup", dest="outgroup_taxa", nargs="*",
    #                           help="Provide taxon names/number for the "
    #                           "outgroup  (This option is only supported for "
    #                           "NEXUS output format files)")

    utilities = parser.add_argument_group("Utilities")
    utilities.add_argument("--get-taxa", dest="get_taxa",
                           action="store_const", const=True,
                           default=False, help="Writes all taxa names into a"
                           " .csv file")

    miscellaneous = parser.add_argument_group("Miscellaneous")
    miscellaneous.add_argument("-quiet", dest="quiet", action="store_const",
                               const=True, default=False, help="Removes all "
                               "terminal output")

    args = parser.parse_args(arg_list)

    # Print help when no arguments are provided
    if len(sys.argv) == 1 and not unittest:
        parser.print_help()
        sys.exit(1)

    return args


def main():
    arguments = get_args()
    triseq_arg_check(arguments)
    main_parser(arguments, arguments.infile)


if __name__ == "__main__":

    main()


__author__ = "Diogo N. Silva"
