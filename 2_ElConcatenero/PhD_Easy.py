#!/usr/bin/python3

# ElConcatenero v4.0.0.2
# Author: Diogo N Silva
# Last update: 03/12/2013
# ElConcatenero is tool to convert and concatenate several commonly used data format types. Currently, supported input formats include Nexus, FastA and Phylip. Output may be in Nexus, Phylip (wiht part file for RaXML) or FastA. Please type "ElConcatenero -h" or read the README.md file for information on usage.

#  Copyright 2012 Diogo N Silva <diogo@arch>
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
#  along with this program; if not,  If not, see <http://www.gnu.org/licenses/>.

import argparse
#import ElParsito3 as ep
from wingman import Alignment,Data
from wingman.ErrorHandling import *


##### ARGUMENT LIST ######

parser = argparse.ArgumentParser(description="Concatenates DNA data matrices")

# Main execution
main_exec = parser.add_argument_group("Main execution")
main_exec.add_argument("-in",dest="infile",nargs="+",help="Provide the input file name. If multiple files are provided, plase separated the names with spaces")
main_exec.add_argument("-if",dest="input_format",default="guess",choices=["fasta","nexus","phylip","guess"],help="Format of the input file(s). The default is 'guess' in which the program tries to guess the input format and genetic code automatically")
main_exec.add_argument("-of",dest="output_format",nargs="+",default="nexus",choices=["nexus","phylip","fasta","mcmctree"],help="Format of the ouput file(s). You may select multiple output formats simultaneously (default is '%(default)s')")
main_exec.add_argument("-o",dest="outfile",help="Name of the output file")

# Alternative modes
alternative = parser.add_argument_group("Alternative execution modes")
alternative.add_argument("-c",dest="conversion",action="store_const",const=True,help="Used for convertion of the input files passed as arguments with the -in option. This flag precludes the usage of the -o option, as the output file name is automatically generated based on the input file name")
alternative.add_argument("-r",dest="reverse",help="Reverse a concatenated file into its original single locus alignments. A partition file similar to the one read by RAxML must be provided")
alternative.add_argument("-z","--zorro-suffix",dest="zorro",type=str, help="Use this option if you wish to concatenate auxiliary Zorro files associated with each alignment. Provide the sufix for the concatenated zorro file")
alternative.add_argument("-p","--partition-file", dest="partition_file", type=str, help="Using this option and providing the partition file will convert it between a RAxML or Nexus format")
alternative.add_argument("-collapse", dest="collapse",action="store_const",const=True, default=False, help="Use this flag if you would like to collapse the input alignment(s) into unique haplotypes")
alternative.add_argument("-gcoder",dest="gcoder", action="store_const", const=True, default=False, help="Use this flag to code the gaps of the alignment into a binary state matrix that is appended to the end of the alignment")
alternative.add_argument("-filter", dest="filter", nargs=2, help="Use this option if you wish to filter the alignment's missing data. Along with this option provide the threshold percentages for gap and missing data, respectively (e.g. -filter 50 75 - filters alignments columns with more than 50%% of gap+missing data and columns with more than 75%% of true missing data)")

# Formatting options
formatting = parser.add_argument_group("Formatting options")
formatting.add_argument("-model",dest="model_phy",default="LG",choices=["DAYHOFF","DCMUT","JTT","MTREV","WAG","RTREV","CPREV","VT","BLOSUM62","MTMAM","LG"],help="This option only applies for the concatenation of protein data into phylip format. Specify the model for all partitions defined in the partition file (default is '%(default)s')")
formatting.add_argument("-interleave",dest="interleave",action="store_const",const="interleave",help="Specificy this option to write output files in interleave format (currently only supported for nexus files")
#formatting.add_argument("-g",dest="gap",default="-",help="Symbol for gap (default is '%(default)s')")
#formatting.add_argument("-m",dest="missing",default="n",help="Symbol for missing data (default is '%(default)s')")

# Data manipulation
manipulation = parser.add_argument_group("Data manipultation")
manipulation.add_argument("-rm",dest="remove",nargs="*",help="Removes the specified taxa from the final alignment. Unwanted taxa my be provided in a csv file containing 1 column with a species name in each line or they may be specified in the command line and separated by whitespace")
manipulation.add_argument("-outgroup", dest="outgroup_taxa", nargs="*", help="Provide taxon names/number for the outgroup (This option is only supported for NEXUS output format files)")

miscellaneous = parser.add_argument_group("Miscellaneous")
miscellaneous.add_argument("-quiet", dest="quiet", action="store_const", const=True,default=False, help="Removes all terminal output")

arg = parser.parse_args()

##### MAIN FUNCTIONS ######

def main_parser(alignment_list):
	""" Function with the main operations of ElConcatenero """
	
	# Defining main variables
	#gap = arg.gap
	#missing_sym = arg.missing
	conversion = arg.conversion
	input_format = arg.input_format
	output_format = arg.output_format
	outfile = arg.outfile
	interleave = arg.interleave
	model_phy = arg.model_phy
	outgroup_taxa = arg.outgroup_taxa

	# Setting leave/interleave format
	if interleave == None:
		sequence_format = "leave"
	else:
		sequence_format = "interleave"

	# Defining output file name
	if arg.conversion == None and arg.outfile != None:
		outfile = "".join(arg.outfile)
	elif arg.conversion != None and arg.outfile != None:
		outfile = "".join(arg.outfile)
	elif arg.conversion != None and arg.outfile == None:
		outfile = "".join(alignment_list).split(".")[0]

	# The input file at this stage is not necessary
	# If just converting the partition file format do this and exit
	if arg.partition_file != None:
		# Initializing Partitions instance
		partition = Data.Partitions(arg.partition_file)
		if partition.partition_format == "nexus":
			partition.write_to_file("raxml", outfile, model_phy)
		else:
			partition.write_to_file("nexus", outfile)
		return 0


	# From here, the input file is mandatory
	if len(alignment_list) == 1:

		# In case only one alignment
		alignment = Alignment.Alignment("".join(alignment_list))

		# Check if input format is the same as output format. If so, and no output file name has been provided, update the default output file name
		if alignment.input_format in output_format and output_format == None:
			outfile = "".join(alignment_list).split(".")[0]+"_conv"

		# If only to reverse a concatenated alignment into individual loci do this and exit
		if arg.reverse != None:
			partition = Data.Partitions(arg.reverse)
			reverse_alignments = alignment.reverse_concatenate(partition)
			reverse_alignments.write_to_file(output_format,form=sequence_format, outgroup_list=outgroup_taxa)
			return 0

	else:

		# With many alignments
		alignments = Alignment.AlignmentList(alignment_list)

		if arg.conversion != None:

			# In case multiple files are to be converted and an alignment filter is to be carried out
			if arg.filter != None:

				alignments.filter_missing_data64(arg.filter[0], arg.filter[1], verbose=True)

			# In case taxa are to be removed while converting
			if arg.remove != None:
				if arg.quiet is False:
					alignments.remove_taxa(arg.remove, verbose=True)
				else:
					alignments.remove_taxa(arg.remove)

			alignments.write_to_file(output_format, form=sequence_format, outgroup_list=outgroup_taxa)
			return 0

		else:

			alignment = alignments.concatenate()

			# If zorro weigth files are provided, concatenate them as well
			if arg.zorro != None:
				zorro = Data.Zorro(alignment_list, arg.zorro)

	# Removing taxa
	if arg.remove != None:
		if arg.quiet is False: print ("\rRemoving taxa", end="")
		alignment.remove_taxa(arg.remove)

	# Collapsing the alignment
	if arg.collapse != False:
		if arg.quiet is False: print ("\rCollapsing alignment", end="")
		alignment.collapse(haplotypes_file=outfile)

	# Codes gaps into binary states
	if arg.gcoder != False:
		if arg.quiet is False: print ("\rCoding gaps", end="")
		if output_format != ["nexus"]:
			raise OutputFormatError("Alignments with gaps coded can only be written in Nexus format")
		alignment.code_gaps()

	if arg.filter != None:
		alignment.filter_missing_data(arg.filter[0], arg.filter[1])

	## Writing files
	if arg.quiet is False: print ("\rWritting output file(s)",end="")
	alignment.write_to_file (output_format, outfile, form=sequence_format, outgroup_list=outgroup_taxa)

	# In case zorro weigth files are provide, write the concatenated file 
	if arg.zorro != None:
		zorro.write_to_file(outfile)

		
def main_check ():
	if arg.partition_file != None and arg.outfile == None:
		raise ArgumentError("An output file must be provided with option '-o'")
		
	if arg.partition_file != None:
		return 0
		
	if arg.conversion == None and arg.outfile == None and arg.reverse == None:
		raise ArgumentError("If you wish to concatenate provide the output file name using the '-o' option. If you wish to convert a file, specify it using the '-c' option")
		
	if len(arg.infile) == 1 and arg.conversion == None and arg.reverse == None and arg.collapse == None:
		raise ArgumentError ("Cannot perform concatenation of a single file. Please provide additional files to concatenate, or specify the conversion '-c' option")
		
	if arg.zorro != None and len(arg.infile) == 1:
		raise ArgumentError ("The '-z' option cannot be invoked when only a single input file is provided. This option is reserved for concatenation of multiple alignment files")

	else:
		return 0
				
def main():
	main_check()
	main_parser(arg.infile)

	if arg.quiet is False: 
		print ("\rProgram done!", end="")

##### EXECUTION ######

main()
