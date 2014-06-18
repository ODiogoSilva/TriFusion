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
from base import html_creator

parser = argparse.ArgumentParser(description="Filters alignment files and creates statistics and graphics for "
											 "alignments")
parser.add_argument("-in", dest="infile", nargs="+", required=True, help="Provide the input files")
parser.add_argument("-o", dest="project_name", required=True, help="Provide a name for the project")

modes = parser.add_argument_group("Report options")
modes.add_argument("-f", dest="full_report", action="store_const", const=True, help="Generate full report")
modes.add_argument("-m", dest="mode", choices=["1", "2", "2a", "2b"], help="Specify which report(s): "
																			   "\n\t\t1: Basic phylogenetic statistics;"
																			   "\n\t\t2: Species specific missing data "
																			   "plot and table;\n\t\t2a: Species "
																			   "specific missing data plot;\n\n\t2b: "
																			   "Species specific missing data table"
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

	# Initialize report and html creator instances
	report = stats.MultiReport(input_list)
	html_instance = html_creator.HTML_template()
	html_instance.addTitle(project)

	if "1" in mode or mode == "all":

		# Defining output file name based on project
		file_name = "%s_basic_table.csv" % project

		report.report_table(file_name)

	if "2a" in mode or "2" in mode or mode == "all":

		# Defining output file name based on project
		file_name = "%s_species_missing_data.svg" % project

		plot_object = report.species_missing_data(plot=True)
		plot_object.render_to_file(file_name)

		# Adding plot to html template object
		html_instance.addSinglePlot("Species specific missing data", file_name, heading_level=3)

	if "2b" in mode or "2" in mode or mode == "all":

		# Defining output file name based on project
		file_name = "%s_species_missing_data.csv" % project

		report.species_missing_data(table=True, output_file=file_name)

	# Finally, write the html file
	html_instance.write_file(project + "_report")


if __name__ == '__main__':
	main()

