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

class Partitions ():
	""" The Partitions class is used to parse and work with partition files. Currently the supported input formats are RAxML partition files and Nexus chartset blocks. To initialize a class, only a file partition file is required. """
	def __init__ (self, partition_file, model_nexus="LG"):

		self.model_nexus = model_nexus

		self.partitions = self._get_partitions(partition_file)

	def _get_format (self,partition_file):
		""" Tries to guess the format of the partition file (Whether it is Nexus of RAxML's) """
		file_handle = open(partition_file)
		
		# Skips first empty lines, if any 
		header = file_handle.readline()
		while header.startswith("\n"):
			header = next(file_handle)
			
		fields = header.split()
		if fields[0].lower() == "charset":
			partition_format = "nexus"
		else:
			partition_format = "raxml"
		
		return partition_format

	def _get_partitions (self, partitions_file):
		""" This function parses a file containing partitions. Supports partitions files similar to RAxML's and NEXUS charset blocks. The NEXUS file, however, must only contain the charset block. The model_nexus argument provides a namespace for the model variable in the nexus format, since this information is not present in the file. However, it assures consistency on the Partition object """
		
		# Get the format of the partition file
		self.partition_format = self._get_format (partitions_file)

		part_file = open(partitions_file)
		partition_storage = []
		
		if self.partition_format == "raxml":
			for line in part_file:
				fields = line.split(",")
				model = fields[0]
				partition_name = fields[1].split("=")[0]
				partition_range_temp = fields[1].split("=")[1]
				partition_range = partition_range_temp.strip().split("-")
				partition_storage.append((model, partition_name, partition_range)) # Format of partition storage: ["str","str",["str","str"]]
				
		elif self.partition_format == "nexus":
			for line in part_file:
				fields = line.split("=")
				partition_name = fields[0].split()[1]
				partition_range = fields[1].replace(";","").strip().split("-")
				partition_storage.append((self.model_nexus, partition_name, partition_range)) # Format of partition storage: ["str", str", ["str","str"]]
			
		return partition_storage

	def write_to_file (self, output_format, output_file, model="LG"):
		""" Writes the Partitions object into an output file according to the output_format. The supported output formats are RAxML and Nexus. The model option is for the RAxML format """

		if output_format == "raxml":
			outfile_handle = open(output_file+".part.File","w")
			for part in self.partitions:
				partition_name = part[0]
				partition_range = "-".join([x for x in part[1]])
				outfile_handle.write("%s, %s = %s\n" % (model, part[0], partition_range))

		elif output_format == "nexus":
			outfile_handle = open(output_file+".charset","w")
			for part in self.partitions:
				outfile_handle.write("charset %s = %s;\n" % (part[1],"-".join(part[2])))

		outfile_handle.close()
		return 0

class Zorro ():

	def __init__ (self, alignment_list, suffix="_zorro.out"):

		self.suffix = suffix
		self.weigth_values = self._zorro2rax(alignment_list)

		def _zorro2rax (self, alignment_list):
			""" Function that converts the floating point numbers contained in the original zorro output files into intergers that can be interpreted by RAxML. If multiple alignment files are provided, it also concatenates them in the same order """
			weigths_storage = []
			for alignment_file in alignment_list:
				zorro_file = alignment_file.split(".")[0]+self.suffix # This assumes that the prefix of the alignment file is shared with the corresponding zorro file
				zorro_handle = open(zorro_file)
				weigths_storage += [round(float(weigth.strip())) for weigth in zorro_handle]
			return weigths_storage

		def write_to_file (self, output_file):
			""" Creates a concatenated file with the zorro weigths for the corresponding alignment files """
			outfile = output_file+"_zorro.out"
			outfile_handle = open(outfile,"w")
			for weigth in zorro_weigths:
				outfile_handle.write("%s\n" % weigth)
			outfile_handle.close()