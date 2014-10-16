#!/usr/bin/env python3
# -*- coding: utf-8 -*-#
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
#  Version:
#  Last update:
#

__author__ = 'DiogoNSilva'

import argparse
from ortho import OrthomclToolbox as OT

parser = argparse.ArgumentParser(description="Toolbox to analyse and filter OrthoMCL group files")

parser.add_argument("-in", dest="infile", nargs="+", help="Provide the OrthoMCL group file(s)")
parser.add_argument("-t1", "--gene-threshold", dest="gene_threshold", nargs=1, default=1, help="Provide the maximum "
					"threshold of gene copy numbers per species that should be allowed")
parser.add_argument("-t2", "--species-threshold", dest="species_threshold", nargs=1, default=200, help="Provide the "
					"minimum number of species that should be allowed.")
parser.add_argument("-e", dest="export", action="store_const", const=True, help="Exports the filtered groups into a "
					"new file")
parser.add_argument("-s", dest="stats", choices=["1"], help="Use the available choices to perform statistical "
					"analyses or summaries of the groups file. 1: Generates basic summary statistics on the number of "
					"clusters and sequences of the groups file.")
parser.add_argument("-g2f", dest="groups2fasta", help="Retrieves the sequences of each cluster to a single file per "
					"cluster. The BLAST database must be provided with this option")
parser.add_argument("-p", dest="pipeline", nargs="*", choices=["1"], help="TriOrtho can be used in pipeline format to "
					"perform sequential steps that users most often use. 1: Filter the original groups file and "
					"retrieve the Fasta sequences. ")

arg = parser.parse_args()

### Execution


def main():

	# Arguments
	groups_file = arg.infile
	gene_threshold = int(arg.gene_threshold[0])
	species_threshold = int(arg.species_threshold[0])
	pipeline_mode = arg.pipeline

	if len(groups_file) == 1:

		group_file = groups_file[0]
		group_object = OT.Group(group_file, gene_threshold, species_threshold)

		if pipeline_mode:

			if "1" in pipeline_mode:
				group_object.export_filtered_group()

		else:

			if arg.export:
				group_object.export_filtered_group()

			#if arg.stats:
				#if "1" in arg.stats:

			if arg.groups2fasta:
				database = arg.groups2fasta
				group_object.retrieve_fasta(database)

	else:
		multiple_groups_object = OT.MultiGroups(groups_file, gene_threshold, species_threshold)

		if arg.stats:
			if "1" in arg.stats:
				multiple_groups_object.basic_multigroup_statistics()


main()