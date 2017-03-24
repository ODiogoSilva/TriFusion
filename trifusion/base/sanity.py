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

from argparse import ArgumentTypeError
import os

try:
    from process.base import print_col, RED, YELLOW
except ImportError:
    from trifusion.process.base import print_col, RED, YELLOW


def triseq_arg_check(arg):

    if arg.gcoder and "nexus" not in arg.output_format:
        print_col("Gap coding can only be performed for Nexus output format.",
                  RED)

    if arg.gcoder and "nexus" in arg.output_format and \
            arg.output_format != ["nexus"]:
        print_col("Gap coding can only be performed for Nexus output format."
                  " This operation will be ignored for other output formats.",
                  YELLOW, quiet=arg.quiet)

    if arg.conversion and arg.reverse:
        print_col("Ignoring conversion flag (-c) when specifying reverse"
                  " concatenation (-r)", YELLOW, quiet=arg.quiet)

    if arg.outfile and arg.reverse:
        print_col("Ignoring output file option (-o) when specifying reverse"
                  " concatenation (-r)", YELLOW, quiet=arg.quiet)

    if arg.partition_file is not None and arg.outfile is None:
        print_col("An output file must be provided with option '-o'", RED)

    if "ima2" in arg.output_format and arg.ima2_params is None:
        print_col("Additional arguments must be provided with the"
                  " option --ima2-params when selecting ima2 output format",
                  RED)

    if "ima2" in arg.output_format and len(arg.ima2_params) != 4:
        print_col("Four additional arguments must be provided with"
                  " option --ima2-params when selecting the "
                  "ima2 output format. %s were given" %
                  (len(arg.ima2_params)), RED)

    if arg.partition_file is not None:
        return 0

    if arg.conversion is None and arg.outfile is None and arg.reverse is None\
            and arg.select is None and arg.get_taxa is False:

        print_col(
            "If you wish to concatenate provide the output file name using "
            "the '-o' option. If you wish to convert a "
            "file, specify it using the '-c' option", RED)

    if len(arg.infile) == 1 and arg.conversion is None and arg.reverse is None\
            and arg.collapse is None:

        print_col(
            "Cannot perform concatenation of a single file. Please provide"
            " additional files to concatenate, or specify the conversion "
            "'-c' option", RED)

    if arg.zorro is not None and len(arg.infile) == 1:
        print_col(
            "The '-z' option cannot be invoked when only a single input "
            "file is provided. This option is reserved for"
            " concatenation of multiple alignment files", RED)

    if arg.consensus and arg.output_format != ["fasta"]:
        print_col("Output format must be only Fasta when using the "
                  "consensus option", RED)

    if not arg.consensus and arg.consensus_single:
        print_col("Ignoring consensus single file option (--consensus-single-"
                  "file) when the consensus operation is not specified",
                  YELLOW, quiet=arg.quiet)

    else:
        return 0


def post_aln_checks(arg, aln_obj):

    if arg.consensus == ["IUPAC"] and aln_obj.sequence_code[0] != "DNA":
        print_col("'IUPAC' option of the consensus operation can "
                  "only be performed on nucleotide alignments.", RED)
    if arg.codon_filter and aln_obj.sequence_code[0] != "DNA":
        print_col("The codon filter option (--codon-filter) can only be"
                  " performed on nucleotide alignments.", RED)
    if aln_obj.bad_alignments:
        print_col("The following input files could not be read or are empty"
                  ": {}".format(" ".join(aln_obj.bad_alignments)), YELLOW)
    if aln_obj.non_alignments:
        print_col("The following input files have alignments of unequal "
                  "length: {}".format(" ".join(aln_obj.non_alignments)),
                  YELLOW)

    else:
        return 0


def mfilters(filt):
    """
    Checks the type of the some filter options. Assures that the values
    are within the accepted boundaries of [0-100]. If the provided filter
    is a decimal, it converts to percentage automatically
    :param filt: (string) The value provided with the --missing-filter
    :return:
    """

    # Check if filt is convertable to float
    try:
        filt = float(filt)
    except ValueError:
        raise ArgumentTypeError("The value '{}' is not a "
                                "number between 0 and 100.".format(filt))

    # Check if filt is within acceptable range
    if not (0 <= filt <= 100):
        raise ArgumentTypeError("The value '{}' must be a number between "
                                "0 and 100.".format(int(filt)))

    # If filt is a float between 0 and 1, convert to percentage
    if 0 <= filt < 1:
        filt = int(filt * 100)
    else:
        filt = int(filt)

    return filt


def check_infile_list(infiles):

    dirs = []
    lost = []
    good_files = []

    for fpath in infiles:

        if not os.path.exists(fpath):
            lost.append(fpath)

        elif os.path.isdir(fpath):
            dirs.append(fpath)

        else:
            good_files.append(fpath)

    return good_files, dirs, lost
