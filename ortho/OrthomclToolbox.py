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

class Cluster():
	""" Object for clusters of the OrhtoMCL groups file. It is usefull to set a number of attributes that will make
	subsequent filtration and processing much easier """


class Group ():
	""" This represents the main object of the orthomcl toolbox module. It is initialized with a file name of a
	orthomcl groups file and provides several methods that act on that group file. To process multiple Group objects,
	see MultiGroups object """

	def __init__(self, groups_file, project_prefix="MyGroups"):

		# Initialize the project prefix for possible ouput files
		self.prefix = project_prefix
		# Initialize attributes for the parser_groups method
		self.groups = OrderedDict()
		self.name = None
		# Parse groups file and populate groups attribute
		self.__parse_groups(groups_file)

	def __parse_groups(self, groups_file):
		"""
		Parses the ortholog clusters in the groups file and populates the self.groups ordered dictionary containing the
		group number as key and the sequence references as values in list mode.
		For each group, it also creates a dictionary containing the gene frequency of each species. This dictionary
		is added as the second elements of the group's dictionary value.
		A final self.groups dictionary should be like: {groups1: [[seq_spA, seq_spB, seq_spC], {spA:1, spB:1, spC:1}]}
		:param groups_file: File name for the orthomcl groups file
		:return: populates the groups attribute
		"""

		def species_frequency(group_cluster):
			"""
			:param group_cluster: List containing the sequence references for each ortholog cluster
			:return: a dictionary containing the frequency of each species in this cluster
			"""

			species_list = set([field.split("|")[1] for field in group_cluster])
			species_frequency_dictionary = dict((species, frequency) for species, frequency in zip(species_list,
											map(lambda species: str(group_cluster).count(species), species_list)))

			return species_frequency_dictionary

		self.name = groups_file
		groups_file_handle = open(groups_file)

		for line in groups_file_handle:
			group_key = line.split(":")[0].strip()
			group_vals = line.split(":")[1].strip().split()

			group_species_frequency = species_frequency(group_vals)

			self.groups[group_key] = [group_vals, group_species_frequency]

	def basic_group_statistics(self, gene_threshold, species_threshold):
		"""
		This method creates a basic table in list format containing basic information of the groups file (total
		number of clusters, total number of sequences, number of clusters below the gene threshold, number of
		clusters below the species threshold and number of clusters below the gene AND species threshold
		:param gene_threshold: Integer with the maximum number of gene copies per species
		:param species_threshold: Integer with the minimum number of species per cluster
		:return: List containing number of [total clusters, total sequences, clusters above gene threshold,
		clusters above species threshold, clusters above gene and species threshold]
		"""
		# Total number of clusters
		total_cluster_num = len(self.groups)

		# Remaining counters
		total_sequence_num = 0
		clusters_gene_threshold = 0
		clusters_species_threshold = 0
		clusters_all_threshold = 0
		for group_vals in self.groups.values():
			# For total number of sequences
			sequence_num = len(group_vals[0])
			total_sequence_num += sequence_num

			species_freq = group_vals[1]
			# For clusters above species threshold
			if len(species_freq) >= species_threshold:
				clusters_species_threshold += 1

			# For clusters below gene threshold
			if max(species_freq.values()) <= gene_threshold:
				clusters_gene_threshold += 1

			if len(species_freq) >= species_threshold and max(species_freq.values()) <= gene_threshold:
				clusters_all_threshold += 1

		statistics = [total_cluster_num, total_sequence_num, clusters_species_threshold, clusters_gene_threshold,
					 clusters_all_threshold]

		return statistics


class MultiGroups ():
	""" Creates an object composed of multiple Group objects """

	def __init__(self, groups_files, project_prefix="MyGroups"):
		"""
		:param groups_files: A list containing the file names of the multiple group files
		:return: Populates the self.multiple_groups attribute
		"""

		self.prefix = project_prefix

		self.multiple_groups = []

		for group_file in groups_files:

			group_object = Group(group_file)
			self.multiple_groups.append(group_object)

	def basic_multigroup_statistics(self, gene_threshold, species_threshold,
									output_file_name="multigroup_base_statistics.csv"):
		"""
		:param gene_threshold: Integer with the maximum number of gene copies per species
		:param species_threshold:
		:param output_file_name:
		:return:
		"""

		# Creates the storage for the statistics of the several files
		statistics_storage = OrderedDict()

		for group in self.multiple_groups:
			group_statistics = group.basic_group_statistics(gene_threshold, species_threshold)
			statistics_storage[group.name] = group_statistics

		output_handle = open(self.prefix + "." + output_file_name, "w")
		output_handle.write("Group file; Total clusters; Total sequences; Clusters below gene threshold; Clusters "
							"above species threshold; Clusters below gene and above species thresholds\n")

		for group, vals in statistics_storage.items():
			output_handle.write("%s; %s\n" % (group, ";".join([str(x) for x in vals])))

		output_handle.close()