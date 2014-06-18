#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
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
#  Author: Diogo N. Silva
#  Version: 0.1
#  Last update: 11/02/14

## TODO [URGENT]: This was originally in python2. I'll need to port to python3

import argparse
from stats import stats

parser = argparse.ArgumentParser(description="Filters alignment files and creates statistics and graphics for "
											 "alignments")
parser.add_argument("-in", dest="infile", nargs="+", required=True, help="Provide the input files")
parser.add_argument("-o", dest="project_name", required=True, help="Provide a name for the project")

modes = parser.add_argument_group("Report options")
modes.add_argument("-f", dest="full_report", action="store_const", const=True, help="Generate full report")
modes.add_argument("-m", dest="mode", nargs="+", choices=["1"], help="Specify which report(s): \n\t\t1: Basic "
																	"phylogenetic information")

arg = parser.parse_args()

def main():

	# Arguments
	input_list = arg.infile
	project = arg.project_name

	# Determine the mode
	if arg.full_report:
		mode = "all"
	else:
		mode = arg.mode

	if "1" in mode or mode == "all":

		report = stats.MultiReport(input_list)

		# Defining output file name based on project
		output_file = "%s_basic_table.csv" % project

		report.report_table(output_file)


if __name__ == '__main__':
	main()

