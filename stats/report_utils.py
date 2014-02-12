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

import pygal
import numpy as np

class Report ():
	def __init__ (self, alignment_dictionaries):
		self.alignment_list = alignment_dictionaries

		self.taxaNames = self.getTaxaNames()

	def getTaxaNames (self):
		""" This function is only used to retrieve all taxa names into a list, to be used by other methods """

		taxaNames = []

		for alignment in self.alignment_list:
			for taxon in alignment:
				if taxon in taxaNames:
					pass
				elif taxon not in taxaNames:
					taxaNames.append(taxon)

		return taxaNames

	def getSpeciesMissingData (self,table=False,plot=False,missing_data_symbol="x",gap_symbol="-"):
		""" Using the alignment_list structure, this returns a dictionary with taxa as keys and a list as value, containing ["Missing genes", "Missing data", "Gaps", "Total missing data", "Total genes", "Total characters"]. If an output_file is provided (other than None) this will write the resulting dictionary into a table. """

		species_intermediate_info = dict((taxon,[]) for taxon in self.taxaNames) # The value table will be [missing genes, missing characters, gaps]
		total_characters = 0

		# Iteration over all taxa for each alignment. Populares the species_intermediate_info dictionary with missing data information for each species for all alignments
		for alignment in self.alignment_list:
			total_characters += len(alignment.values()[0])
			for taxon in self.taxaNames:
				if species_intermediate_info[taxon] == []:
					species_intermediate_info[taxon].extend([0, 0,0])
				if taxon not in alignment:
					missing_gene = 1
					missing_characters = len(alignment.values()[0]) # Assuming all taxa have the same number of characters, as in a well formatted alignment
					species_intermediate_info[taxon][0] += missing_gene
					species_intermediate_info[taxon][1] += missing_characters
				elif taxon in alignment:
					missing_characters = alignment[taxon].count(missing_data_symbol)
					gaps = alignment[taxon].count("-")
					species_intermediate_info[taxon][1] += missing_characters
					species_intermediate_info[taxon][2] += gaps

		# Creation of the final list with missing and total data values for each species
		final_list = [(taxon, [values[0],len(self.alignment_list),values[1],values[2],(values[1]+values[2]),total_characters]) for taxon, values in species_intermediate_info.items()]

		# Calculating mean missing and gap data
		mean_missing = round(np.mean([float(missing_value[1][2])/float(missing_value[1][5]) for missing_value in final_list])*100,2)
		mean_gap = round(np.mean([float(missing_value[1][3])/float(missing_value[1][5]) for missing_value in final_list])*100,2)
		mean_total_missing = round(np.mean([float(missing_value[1][4])/float(missing_value[1][5]) for missing_value in final_list])*100,2)

		# Order species by missing data 
		final_list_sorted = sorted(final_list, key=lambda x: x[1][4],reverse=True)

		# If a table with the values is to be produced
		if table != False:
			output_handle = open("Report_results/"+table+".csv","w")

			output_handle.write("Taxon; Missing genes; Total genes; Missing data; Gaps; Total missing data; Total characters\n")

			for taxon, values in final_list_sorted:
				values = [str(x) for x in values]
				output_handle.write("%s; %s\n" % (taxon, ";".join(values)))

			output_handle.close()


		if plot == True:

			#### Plotting by characters ####

			line_chart = pygal.StackedBar(x_label_rotation=45,width=1200, legend_at_bottom=True, height=800,label_font_size=8,legend_font_size=20,margin=50,major_label_font_size=10)
			line_chart.title = "Character missing data per species (Mean missing data: %s%%; Mean gap data: %s%%; Total average missing data %s%%)" % (mean_missing, mean_gap, mean_total_missing)
			line_chart.x_labels = [taxon[0] for taxon in final_list_sorted]


			# Get list of gaps, missing data, and total characters proportions
			gaps_list, missing_list, characters_list, genes_list = [], [], [], []
			for taxon,value in final_list_sorted:
				gap_dic = {"value": float(value[3])/float(value[-1]), "label": str(value[3])}
				missing_dic = {"value": float(value[2])/float(value[-1]), "label":str(value[2])}
				characters_dic =  {"value":((float(value[5])-float(value[4]))/float(value[-1])), "label":str(value[5]-value[4])}
				genes_dic = {"value":float(value[0])/float(value[1]), "label": str(value[0])+" out of %s" % (value[1])}
				gaps_list.append(gap_dic)
				missing_list.append(missing_dic)
				characters_list.append(characters_dic)
				genes_list.append(genes_dic)

			line_chart.add("Gaps", gaps_list)
			line_chart.add("Missing Data", missing_list)
			line_chart.add("Effective characters", characters_list)

			#### Plotting by genes ####

			speciesGene_chart = pygal.StackedBar(x_label_rotation=45,legend_at_bottom=True,width=1200, height=800,label_font_size=8,legend_font_size=20,margin=50)
			speciesGene_chart.title = "Gene missing data per species"
			speciesGene_chart.x_labels = [taxon[0] for taxon in final_list_sorted]

			speciesGene_chart.add("Missing genes", genes_list)

			return line_chart, speciesGene_chart

		else:
			return 0

	def getSizeDistribution (self, table=False, plot=True):
		""" Function returns the distribution of average sequence sizes """

		sequence_sizes = []

		# Get average sequence length for each alignment
		for alignment in self.alignment_list:
			seq_sizes = [len(seq.replace("-","").replace("X","")) for seq in alignment.values()]
			avg_seq_size = np.mean(seq_sizes)
			sequence_sizes.append(avg_seq_size)

		sizes_frequency, size_labels = self.getDistribution(sequence_sizes, 30)

		if table != False:

			output_handle = open("Report_results/"+table+"_sizeDistribution.csv","w")
			output_handle.write("Range; Frequency\n")

			for i in range (len(size_labels)):
				output_handle.write("%s; %s\n" % (size_labels[i], sizes_frequency[i]))

			output_handle.close()

		if plot == True:
			barDistribution = pygal.Bar(x_label_rotation=45,width=1200, legend_at_bottom=True, height=800,label_font_size=8,legend_font_size=20,margin=50,major_label_font_size=10,show_legend=False)
			barDistribution.title = "Average sequence size distribution (Mean: %s; Median %s)" % (round(np.mean(sequence_sizes),2), round(np.median (sequence_sizes),2))
			barDistribution.x_labels = size_labels

			barDistribution.add("sequence size (aa)", sizes_frequency)

			return barDistribution

		else:

			return 0

	def getSimilarityDistribution (self, table=False, plot=True):
		""" Function returns the distribution of similarity among sequences """

		similarity_values = []

		# Get average similarity for each alignment
		for alignment in self.alignment_list:
			similarity_values.append(self.getAlignmnentSimilarity(alignment))

		similarity_frequency, similarity_labels = self.getDistribution(similarity_values, 30)

		if table != False:

			output_handle = open("Report_results/"+table+"_similarityDistribution.csv","w")
			output_handle.write("Range; Frequency\n")

			for i in range (len(similarity_labels)):
				output_handle.write("%s; %s\n" % (similarity_labels[i], similarity_frequency[i]))

			output_handle.close()

		if plot == True:
			barDistribution = pygal.Bar(x_label_rotation=45,width=1200, legend_at_bottom=True, height=800,label_font_size=8,legend_font_size=20,margin=50,major_label_font_size=10,show_legend=False)
			barDistribution.title = "Average sequence similarity distribution (Mean: %s; Median %s)" % (round(np.mean(similarity_values),2), round(np.median (similarity_values),2))
			barDistribution.x_labels = similarity_labels

			barDistribution.add("sequence similarity", similarity_frequency)

			return barDistribution

		else:

			return 0

	def getAlignmnentSimilarity (self, alignment):
		""" Calculates the average similarity of an alignment """

		def getSimilarity (seq1, seq2):
			"""" Calculates the pairwise sequence similarity """

			similarity = 0

			for char1, char2 in zip(*[seq1, seq2]):
				if char1 == "-" or char2 == "-":
					pass
				elif char1 == char2:
					similarity += 1

			similarity_percentage = float(similarity)/float(len(seq1))*100

			return similarity_percentage

		sequence_list = list(alignment.values())

		pairwise_similarities = []

		for sequence1 in alignment.values():
			sequence_list.remove(sequence1)
			pairwise_similarities.extend([getSimilarity(sequence1,sequence2) for sequence2 in sequence_list])

		avg_pairwise_similarities = np.mean(pairwise_similarities)

		return avg_pairwise_similarities


		# 	pairwise_comparisons.extend([(sequence1, sequence2) for sequence2 in sequence_list])

		# for pair in pairwise_comparisons:
		# 	pairwise_similarities.append(getSimilarity(pair[0], pair[1]))




	def getDistribution (self, values, bins):
		""" Function returns the distribution of average sequence sizes. It requires the number of bins on the plot """

		# Sort sequence sizes
		values = sorted(values)

		# Determine bin range
		bin_range = (max(values)-min(values))/bins
		bin_rest = (max(values)-min(values))%bins

		# Distribute sequence sizes over a number of categories determined by the bin range variable
		size_distribution = []
		for i in range (int(min(values)), int(max(values)), int(bin_range)):
			size_distribution.append([size for size in values if (i+bin_range)-1 > size > i ])
			last = i
		else:
			size_distribution.append([size for size in values if last+bin_rest > size > last])
		
		frequency_distribution = [len(sizes) for sizes in size_distribution]

		# Define labels based on the minimum and maximum size of each categorie
		labels = [str(i)+"-"+str(i+int(bin_range-1)) for i in range(int(min(values)), int(max(values)), int(bin_range))]

		return frequency_distribution, labels

	# def variation_profile (self, mode="all"):
	# 	""" Creates an amino acid variation profiles of a or a set of alignments """

	# 	aminoacid_table ={"A":["Alanine","nonpolar","neutral"],
	# "R":["Arginine","polar","positive"],
	# "N":["Asparagine","polar","neutral"],
	# "D":["Aspartate","polar","negative"],
	# "C":["Cysteine","nonpolar","neutral"],
	# "E":["Glutamate","polar","negative"],
	# "Q":["Glutamine","polar","neutral"],
	# "G":["Glycine","nonpolar","neutral"],
	# "H":["Histidine","polar","neutral"],
	# "I":["Isoleucine","nonpolar","neutral"],
	# "L":["Leucine","nonpolar","neutral"],
	# "K":["Lysine","polar","positive"],
	# "M":["Methionine","nonpolar","neutral"],
	# "F":["Phenylalanine","nonpolar","neutral"],
	# "P":["Proline","nonpolar","neutral"],
	# "S":["Serine","polar","neutral"],
	# "T":["Threonine","polar","neutral"],
	# "W":["Tryptophan","nonpolar","neutral"],
	# "Y":["Tyrosine","polar","neutral"],
	# "V":["Valine","nonpolar","neutral"],
	# "U":["Selenocysteine","",""],
	# "O":["Pyrrolysine","",""],
	# "X":["Missing","",""]	
	# }

	# 	profile_storage, variable_sites, parsimony_informative = [], 0, 0

	# 	aminoacid_frequency = dict((code,0) for code in aminoacid_table.keys())

	# 	for alignment in self.alignment_list:
	# 		for column in zip(*alignment.values()): # Iterate over each alignment column
	# 			column = column.upper().remove("X").remove("-") # Remove missing and gap data
	# 			profile_storage.append(len(column))

	# 			if len(column) > 1:
	# 				variable_sites += 1

	# 			if len(column) > 2:
	# 				parsimony_informative += 1

	# 			for char in column:
	# 				aminoacid_frequency[char] += 1


	def species_histogram (self, table, highlight_taxa=None):
		""" Creates a simple histogram with the frequency of species in a data set """

		species_values = dict ((taxon, 0) for taxon in self.taxaNames)

		for alignment in self.alignment_list:
			for taxon in species_values:
				species_values[taxon] += list(alignment.keys()).count(taxon)

		species_values = [(sp, val) for sp,val in species_values.items()]
		species_values = sorted(species_values, key=lambda x: x[1],reverse=True)

		if table != True:
			output_handle = open("Report_results/"+table+"_speciesFrequency.csv","w")
			output_handle.write("Species; Frequency\n")

			for sp, val in species_values:
				output_handle.write("%s; %s\n" % (sp, val))

			output_handle.close()

			if highlight_taxa != None:
				output_handle = open("Report_results/"+table+"_speciesFrequencyHighlighted_taxa.csv","w")
				output_handle.write("Species; Frequency\n")

				for sp, val in species_values:
					if sp in highlight_taxa:
						output_handle.write("%s; %s\n" % (sp, val))

				output_handle.close()

		species_chart = pygal.Bar(x_label_rotation=45,width=1200, legend_at_bottom=True, height=800,label_font_size=8,legend_font_size=20,margin=50,major_label_font_size=10)
		species_chart.title = "Species frequency in the data set"
		species_chart.x_labels = [element[0] for element in species_values]

		if highlight_taxa != None:
			major_taxa_subset = []
			for taxon in highlight_taxa:
				if taxon in species_chart.x_labels:
					major_taxa_subset.append(taxon)

			species_chart.x_labels_major = major_taxa_subset

		species_chart.add("Species frequency", [element[1] for element in species_values])

		return species_chart

	def gene_histogram (self,major_taxa=None):
		""" Creates a simple histrogram with the frequency of species contained in genes """

		def getMajorTaxaCount (alignment, major_taxa):
			""" Simply counts the number of major taxa in a given alignment """
			taxa = 0
			for species in alignment:
				if species in major_taxa:
					taxa+=1
			return taxa

		gene_chart = pygal.Bar(width=1200, legend_at_bottom=True, height=800,label_font_size=12,legend_font_size=20,margin=50,major_label_font_size=10)
		gene_chart.title = "Distribution of the species occurence in the data set genes"

		gene_values = []

		for alignment in self.alignment_list:
			if major_taxa == None:
				gene_values.append(len(alignment))
				gene_chart.x_labels = [str(int(round(float(x)/float(len(self.taxaNames)),2)*100))+"%" for x in range(max(gene_values), min(gene_values), -1)]
			else:
				gene_values.append(getMajorTaxaCount(alignment, major_taxa))
				gene_chart.x_labels = [str(int(round(float(x)/float(len(major_taxa)),2)*100))+"%" for x in range(max(gene_values), min(gene_values), -1)]

		gene_values = sorted(gene_values,reverse=True)

		histogram_values = [gene_values.count(x) for x in range(max(gene_values), min(gene_values), -1)]

		gene_chart.add("Gene species content", histogram_values)

		return gene_chart
		
	def avg_gene_length (self,table):

		species_dic = dict((sp,[]) for sp in self.taxaNames)

		for alignment in self.alignment_list:
			for species, sequence in alignment.items():
				seq = sequence.upper().replace("X","").replace("-","")
				species_dic[species].append(len(seq))

		# Get averages
		for species, values in species_dic.items():
			species_dic[species] = np.mean(values)

		output_handle = open("Report_results/"+table+"_averageGeneLength.csv","w")

		for species, value in species_dic.items():
			output_handle.write("%s;%s\n" % (species, value))

		output_handle.close()

		return 0


	def writePlot (self, plotObj, plotName):
		""" Creates a new folder and writes the results in it """
		subprocess.Popen(["mkdir Results"],shell=True).wait()

		plotObj.render_to_file("Results/"+plotName+".svg")