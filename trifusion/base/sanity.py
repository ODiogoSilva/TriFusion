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

try:
    from process.base import print_col, RED
except ImportError:
    from trifusion.process.base import print_col, RED


def triseq_arg_check(arg):
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