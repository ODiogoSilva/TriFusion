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


from process.alignment import Alignment
from collections import OrderedDict

class UniReport(Alignment):
	""" This will create a report object for single alignments and will inherit from the Alignment class. As with its
	base class, it will be possible to create a UniReport object directly from an external file or from an ordered
	dictionary object
	"""

	def base_statistics(self):
		"""
		returns a list with the basic phylogenetic statistics, such as number of indels, overall percentage of
		indels variable sites,parsimony informative sites.
		"""

		from collections import Counter

		parsimonious = 0
		variable = 0
		indel = 0

		for column_position in range(self.locus_length):

			column = [char[column_position].lower() for char in self.alignment.values()]

			# Check for variable sites
			column_set = len(set(column) - set(["n",self.sequence_code[1]]))
			if column_set > 1:
				variable += 1

				# Check for parsimonious informative sites
				if len([val for char, val in Counter(column).items() if x > 1 and y not in "-" and y not in "n"]) > \
						1:  # Filters the Counter dictionary with key:val pairs with vals that are not missing
						# data or gap and higher than 1
					parsimonious += 1

			# Check for indels
			if self.sequence_code[1] in column:
				indel += 1

		cumulative_gap = 0
		for sequence in self.alignment.values():
			cumulative_gap += sequence.count(self.sequence_code[1]) + sequence.count("-")

		indel_perc = float(cumulative_gap) / (self.locus_length*len(self.alignment))

		return indel, indel_perc, variable, parsimonious

class MultiReport():
	""" This will create a report object for multiple alignments. Its most basic instance will be a list of UniReport
	objects
	"""
	def __init__(self, alignment_list):

		self.report_list = []

		for report in alignment_list:
			report_object = UniReport(report)

			self.report_list.append(report_object)


	def report_table (self, output_file):

		output_handle = open(output_file,"w")
		output_handle.write("Gene; Number of sites; Number of indels; Percentage of missing data; Variable sites; "
							"Parsimonious sites\n")

		table_contents = OrderedDict()

		for report in self.report_list:
			name = report.input_alignment
			sites = report.locus_length
			indel, indel_perc, variable, parsimonious = report.base_statistics()

			# In case I'll need this in anything other than a table
			table_contents[name] = [sites, indel, indel_perc, variable, parsimonious]

		for gene, vals in table_contents.items():

			output_handle.write("%s; %s; %s; %s; %s; %s\n" % (gene, vals[0], vals[1], vals[2], vals[3], vals[4]))

		output_handle.close()

__author__ = "Diogo N. Silva"
__copyright__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Diogo N. Silva"
__email__ = "o.diogosilva@gmail.com"
__status__ = "Prototype"
