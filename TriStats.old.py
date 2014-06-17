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
import pygal
from base import html_creator
import subprocess
import stats.report_utils as ru

parser = argparse.ArgumentParser(description="Filters alignment files and creates statistics and graphics for alignments")

parser.add_argument("-in",dest="infile",nargs="+",required=True,help="Provide the input files")
parser.add_argument("-plot", dest="plot",action="store_const",const=True,help="Generate html plots")
parser.add_argument("-table", dest="table", action="store_const", const=True, help="Export tables in csv format")
parser.add_argument("-o", dest="project_name", required=True, help="Provide a name for the project")

modes = parser.add_argument_group("Report options")
modes.add_argument("-full", dest="full_report",action="store_const",const=True, help="Generate full report")
modes.add_argument("-missing", dest="missing_report",action="store_const",const=True, help="Generate missind data related report")
modes.add_argument("-sizes", dest="sizes_report", action="store_const", const=True, help="Generate report on sequence sizes distribution")
modes.add_argument("-sim", dest="similarity_report", action="store_const", const=True, help="Generate report on sequence similarity distribution")
modes.add_argument("-count",dest="count_report", action="store_const", const=True, help="Generate a report on gene and species frequency on the data set")

misc = parser.add_argument_group("Miscelaneous")
misc.add_argument("-taxa", dest="taxa_set", nargs="+", help="Provide a taxa set to be highlighted in the analysis")

arg = parser.parse_args()

def parse_alignments (file_list):
	""" Main alignment parsing function """

	alignment_list = []

	# Automatically find file format (from the three available: Fasta, Phylip, Nexus)
	input_format, character_code = ep.autofinder(file_list[0]) # Character code refers to DNA/Protein

	#Initializing parsing instance
	alignment_instance = ep.SeqUtils()

	#Parsing input alignments
	for input_alignment in file_list:
		alignment_list.append(alignment_instance.read_alignment(input_alignment, input_format)[0])

	return alignment_list, character_code

class MainReport ():

	def __init__ (self, alignment_dic, project_name, plot=True, table=True):
		self.alignment_dic = alignment_dic
		self.project_name = project_name
		self.plot = plot
		self.table = table
		# Initialization of report instance
		self.alignmentReportInstance = ru.Report(self.alignment_dic)
		# Initialization of html instance
		self.htmlInstance = html_creator.HTML_template()
		# Creation of results folder
		subprocess.Popen(["mkdir Report_results"],shell=True).wait()

	def define_taxa_subset (self, new_taxa_set=None):
		""" Use this function to define taxa set to be highlighted in the plots """
		EST_taxa = ["Armillaria_tabescens", "Puccinia_graminis_f._sp._tritici", "Suillus_luteus", "Uromyces_viciaefabae", "Pisolithus_microcarpus", "Pisolithus_tinctorius", "Melampsora_medusae_f._sp._tremuloidis", "Auricularia_auricula-judae", "Agrocybe_cylindracea", "Amanita_muscaria", "Amylostereum_areolatum", "Asterotremella_humicola", "Athelia_rolfsii", "Cronartium_quercuum_f._sp._fusiforme", "Cryptococcus_laurentii", "Cryptococcus_vishniacii", "Flammulina_velutipes", "Fomitopsis_palustris", "Ganoderma_lucidum", "Hemileia_vastatrix", "Hericium_sp", "Lactarius_quietus", "Lentinula_edodes", "Leucosporidium_scottii", "Melampsora_medusae_f._sp._deltoidis", "Melampsora_occidentalis", "Microbotryum_violaceum", "Moniliophthora_perniciosa","Moniliophthora_roreri","Phakopsora_pachyrhizi", "Phanerochaete_chrysosporium_EST", "Phellinidium_sulphurascens", "Pholiota_nameko", "Pisolithus_albus", "Pleurotus_sp._Florida","Puccinia_coronata_var._lolii", "Puccinia_striiformis", "Puccinia_striiformis_f._sp._tritici", "Rhizoctonia_solani", "Taiwanofungus_camphoratus", "Thanatephorus_cucumeris", "Uromyces_appendiculatus", "Volvariella_volvacea"]
		Pucciniomycotina_taxa = ["Melampsora_laricis-populina", "Puccinia_graminis", "Puccinia_triticina", "Sporobolomyces_roseus", "Cronartium_quercuum_f._sp._fusiforme", "Hemileia_vastatrix", "Leucosporidium_scottii", "Melampsora_medusae_f._sp._deltoidis", "Melampsora_medusae_f._sp._tremuloidis", "Melampsora_occidentalis", "Microbotryum_violaceum", "Phakopsora_pachyrhizi", "Puccinia coronata_var._lolii", "Puccinia_graminis_f._sp._tritici", "Puccinia_striiformis", "Puccinia_striiformis_f._sp._tritici", "Uromyces_appendiculatus, Uromyces_viciae-fabae"]
		if new_taxa_set == ["EST"]:
			self.taxaSet = EST_taxa
		elif new_taxa_set == ["Pucciniomycotina"]:
			self.taxaSet = Pucciniomycotina_taxa
		else:
			self.taxaSet = new_taxa_set

	def missing_report (self, missing_symbol, major_taxa=None):
		""" Main report generator for missing data. Requires a report instance as input """

		if self.table == True:
			returnGenerator = self.alignmentReportInstance.getSpeciesMissingData(plot=self.plot, missing_data_symbol=missing_symbol, table=self.project_name+"_missingDataTable")
		else:
			returnGenerator = self.alignmentReportInstance.getSpeciesMissingData(plot=self.plot, missing_data_symbol=missing_symbol)

		if returnGenerator != 0:
			charPlot, genePlot = returnGenerator

			if major_taxa != None:
				self.define_taxa_subset (major_taxa)
				major_taxa_subset = []

				for taxon in self.taxaSet:
					if taxon in charPlot.x_labels:
						major_taxa_subset.append(taxon)

				charPlot.x_labels_major = major_taxa_subset
				genePlot.x_labels_major = major_taxa_subset

			self.writePlot (charPlot, self.project_name+"_characterMissingPlot")
			self.writePlot (genePlot , self.project_name+"_geneMissingPlot")

			self.htmlInstance.addSinglePlot("Character missing data ", self.project_name+"_characterMissingPlot.svg", heading_level=3)
			self.htmlInstance.addSinglePlot("Gene missing data", self.project_name+"_geneMissingPlot.svg", heading_level=3)

	def sizes_report (self):
		""" Main report generator for sequence sizes distribution """

		if self.table == True:
			returnGenerator = self.alignmentReportInstance.getSizeDistribution(plot=self.plot, table=self.project_name+"_SeqSizeDistribution")
		else:
			returnGenerator = self.alignmentReportInstance.getSizeDistribution(plot=self.plot)

		if returnGenerator != 0:
			sizesPlot = returnGenerator

			self.writePlot (sizesPlot, self.project_name+"_sizeDistribution")

			self.htmlInstance.addSinglePlot("Sequence size distribution", self.project_name+"_sizeDistribution.svg", heading_level=3)

		return 0

	def similarity_report (self):
		""" Main report generator for sequence similarity distribution """

		if self.table == True:
			returnGenerator = self.alignmentReportInstance.getSimilarityDistribution(plot=self.plot, table=self.project_name+"_SeqSimilarityDistribution")
		else:
			returnGenerator = self.alignmentReportInstance.getSimilarityDistribution(plot=self.plot)

		if returnGenerator != 0:
			similarityPlot = returnGenerator

			self.writePlot(similarityPlot, self.project_name+"_similarityDistribution")

			self.htmlInstance.addSinglePlot("Sequence similarity distribution", self.project_name+"_similarityDistribution.svg", heading_level=3)

		return 0

	def count_report (self, major_taxa=None):
		""" Main report generator for counting methods """

		if major_taxa != None:
			self.define_taxa_subset(major_taxa)

		self.alignmentReportInstance.avg_gene_length(self.project_name)

		SpcountPlot = self.alignmentReportInstance.species_histogram(self.project_name, highlight_taxa=self.taxaSet)
		GncountPlot = self.alignmentReportInstance.gene_histogram(self.taxaSet)

		self.writePlot(SpcountPlot, self.project_name+"_speciesHistogram")
		self.writePlot(GncountPlot, self.project_name+"_geneHistogram")

		self.htmlInstance.addSinglePlot("Species frequency on the data set", self.project_name+"_speciesHistogram.svg", heading_level=3)
		self.htmlInstance.addSinglePlot("Gene species cotente", self.project_name+"_geneHistogram.svg",heading_level=3)

		return 0

	def writePlot (self,pygalObj, plotName):
		""" Simple function that renders pygalObj """
		pygalObj.render_to_file("Report_results/"+plotName+".svg")

		return 0

	def write2HTML (self):
		""" Main function that writes report contents to html """

		self.htmlInstance.addTitle(self.project_name)
		
		self.htmlInstance.write_file("Report_results/"+self.project_name+"_report")

		return 0

def main():

	### Arguments
	input_file_list = arg.infile
	project_name = arg.project_name
	taxa_subset = arg.taxa_set
	# Report options = [Missing]
	if arg.full_report == True:
		report_options = [True]*4
	else:
		report_options = [arg.missing_report, arg.sizes_report, arg.similarity_report, arg.count_report]

	# Check if there are any options specified
	if list(set(report_options)) == [False] or list(set(report_options)) == [None]:
		print ("No report options have been specified. Exiting...\n")
		raise SystemExit

	### Alignment parsing
	alignments = parse_alignments(input_file_list)
	alignment_dic = alignments[0]
	missing_code = alignments[1][1]

	### Report options
	# Initializing report instance
	alignmentReportInstance = MainReport(alignment_dic, project_name)

	# Report_execution
	if report_options[0] == True:
		alignmentReportInstance.missing_report(missing_code, taxa_subset)
	if report_options[1] == True:
		alignmentReportInstance.sizes_report()
	if report_options[2] == True:
		alignmentReportInstance.similarity_report()
	if report_options[3] == True:
		alignmentReportInstance.count_report(taxa_subset)

	alignmentReportInstance.write2HTML()
	
	return 0 

if __name__ == '__main__':
	main()

