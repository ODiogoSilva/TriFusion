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

from collections import OrderedDict


class Group ():
	""" This represents the main object of the orthomcl toolbox module. It is initialized with a file name of a
	orthomcl groups file and provides several methods that act on that group file. To process multiple Group objects,
	see MultiGroups object """

	def __init__(self, groups_file):

		# Initialize groups attribute
		self.groups = OrderedDict()
		# Parse groups file and populate groups attribute
		self.parse_groups(groups_file)

	def parse_groups(self, groups_file):
		"""
		Parses the ortholog clusters in the groups file and creates an ordered dictionary attributed containing the
		group number as key and the sequence references as values in list mode (e.x., {group1:[seq_spA, seq_spB (...)]})
		:param groups_file: File name for the orthomcl groups file
		:return: populates the groups attribute
		"""

		groups_file_handle = open(groups_file)

		for line in groups_file_handle:
			group_key = line.split(":")[0].strip()
			group_vals = line.split(":")[1].strip().split()

			self.groups[group_key] = group_vals
