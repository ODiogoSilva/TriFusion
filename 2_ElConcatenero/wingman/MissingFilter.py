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
#  

class MissingFilter ():
	""" Contains several methods used to trim and filter missing data from alignments. It's mainly used for inheritance """

	def __init__ (self, alignment_dict, gap_threshold=50, missing_threshold=75, gap_symbol="-", missing_symbol="n"):
		""" the gap_threshold variable is a cut-off to total_missing_proportion and missing_threshold in a cut-off to missing_proportion """

		self.alignment = alignment_dict
		self.gap = gap_symbol
		self.missing = missing_symbol

		# Definig thresholds
		self.gap_threshold = gap_threshold
		self.missing_threshold = missing_threshold

		# Basic filter
		self.filter_terminals()
		self.filter_columns()

	def filter_terminals (self):
		""" Given an alignment, this will replace the gaps in the extremities of the alignment with missing data """

		for taxa,seq in self.alignment.items():

			trim_seq = list(seq)
			counter, reverse_counter = 0, -1

			while trim_seq[counter] == self.gap:
				trim_seq[counter] = self.missing
				counter += 1

			while trim_seq[reverse_counter] == self.gap:
				trim_seq[reverse_counter] = self.missing
				reverse_counter -= 1

			seq = "".join(trim_seq)

			self.alignment[taxa] = seq

	def filter_columns (self,verbose=False):
		""" Here several missing data metrics are calculated, and based on some user defined thresholds, columns with inappropriate missing data are removed """

		taxa_number = len(self.alignment)
		self.old_locus_length = len(list(self.alignment.values())[0])

		filtered_alignment = dict((taxa, list(seq)) for taxa, seq in self.alignment.items())

		# Creating the column list variable
		for column_position in range(self.old_locus_length-1, -1, -1): # The reverse iteration over the sequences is necessary to maintain the column numbers when removing them

			if verbose == True:
				print ("\rFiltering alignment column %s out of %s" % (column_position+1, self.old_locus_length+1), end="")

			column = [char[column_position] for char in filtered_alignment.values()]

			# Calculating metrics
			gap_proportion = (float(column.count(self.gap))/float(taxa_number))*float(100)
			missing_proportion = (float(column.count(self.missing))/float(taxa_number))*float(100)
			total_missing_proportion = gap_proportion+missing_proportion

			if total_missing_proportion > float(self.gap_threshold):

				list(map ((lambda seq: seq.pop(column_position)), filtered_alignment.values()))

			elif missing_proportion > float(self.missing_threshold):

				list(map ((lambda seq: seq.pop(column_position)), filtered_alignment.values()))

		self.alignment = dict((taxa, "".join(seq)) for taxa,seq in filtered_alignment.items())
		self.locus_length = len(list(self.alignment.values())[0])

