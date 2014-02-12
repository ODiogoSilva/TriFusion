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

import sys
import re
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Pipeline for the OrthoMCL software")
parser.add_argument("-in",dest="infile",nargs="+",help="Input proteome files in fasta format")
parser.add_argument("-a",action="store_const",const=True,dest="adjust",help="Run only the AdjustFasta program")
parser.add_argument("-na",action="store_const",const=True,dest="no_adjust",help="Run only the AdjustFasta program")
parser.add_argument("-c",action="store_const",const=True,dest="code",help="The proteome file names are already in code (e.g. Homo_sapiens.fas -> HoSap.fas). Not advisable because it will overwrite the originals!")
parser.add_argument("-p",action="store_const",const=True,dest="check",help="Checks for duplicates and other potential errors")
parser.add_argument("-n",action="store_const",const=True,dest="normal",help="Normal run of the pipeline")

arg = parser.parse_args()

## PARAMETER DEFINITIONS ##

# Configuration file for orthocml
config_file = "orthomcl.config"

# For adjustFast
Unique_ID = 1 # Specify the field of the sequence header where the unique protein ID is. If the field is different depending on the proteomes, perhaps it will be better to group proteomes with unique ID in the same field number and run the script with the -a flag
Name_separator = "_" # Specify the name separator in the input files (e.g., the separator is "_" if the file name is Homo_sapiens.fasta). This parameter only applies if the file names are not already in code format (e.g., for the "Homo_sapiens" a code format could be "HoSap").

# For filterFasta
min_length = 10 # Minimum allowed length of proteins.  (suggested: 10)
max_percent_stop = 20 # Maximum percent stop codons.  (suggested 20)

# For allvsallBLAST
database_name = "goodProteins_db"
blast_out_name = "AllVsAll_BLAST.out"
evalue_cutoff = "0.00001" # 1E-5 - Recommended
CPUs = "15" # Number of CPU's for multiprocessing

# For mcl
inflation = "1.5"

# For mclGroups
prefix = "Basidiomycota" # Arbitrary string for the name of the groups
start_ID = "1000" # Startin number for the groups
groups_file = "groups.txt"

def loading (current_state,size,prefix,width,proteome):
	""" Function that prints the loading progress of the script """
	percentage = int((current_state/size)*100)
	complete = int(width*percentage*0.01)
	sys.stdout.write("\r%s [%s%s] %s%%  %s" % (prefix,"#"*complete,"."*(width-complete),percentage,proteome))
	sys.stdout.flush()

def lineCounter (infile):
	file_read = open(infile,"r")
	x = -1
	for line in file_read:
		x += 1
	return x

def ID_duplicateCheck(proteome):
	#print ("\t Checking duplicates for "+proteome)
	previous_IDs,duplicate_IDs = [], {}
	subprocess.Popen(["mv "+proteome+" "+proteome+".old"],shell=True).wait()
	with open(proteome+".old","r") as proteome_read, open(proteome,"w") as proteome_out, open("ID_check.log","a") as ID_check_log:
		x = 1
		ID_check_log.write("## ID duplicate check for "+proteome+"\n\n")
		for line in proteome_read:
			if line.startswith(">") and line not in previous_IDs:
				previous_IDs.append(line)
				proteome_out.write(line)
			elif line.startswith(">") and line in previous_IDs:
				duplicate_IDs[line] = line[:-1]+"_"+str(x)+"\n"
				previous_IDs.append(line[:-1]+"_"+str(x)+"\n")
				proteome_out.write(line[:-1]+"_"+str(x)+"\n")
				x += 1
			elif "(" in line or ")" in line:
				pass
			else:
				proteome_out.write(line)
		if duplicate_IDs == {}:
			ID_check_log.write("No duplicates found!\n")
		else:
			for k,v in duplicate_IDs.items():
				ID_check_log.write("Found duplicate entry: "+k[:-1]+"; replaced by "+v+"\n")
	subprocess.Popen(["rm "+proteome+".old"],shell=True).wait()
	
def prepBlastDB (proteome_file):
	print ("\t Preparing file for BLAST database")
	subprocess.Popen(["mv "+proteome_file+" "+proteome_file+".old"],shell=True).wait()
	with open (proteome_file+".old","r") as file_in, open(proteome_file,"w") as file_out:
		for line in file_in:
			if line.startswith(">"):
				file_out.write(">gnl|"+line[1:])
			else:
				file_out.write(line)
	subprocess.Popen(["rm "+proteome_file+".old"],shell=True).wait()
	
def adjustFasta(Proteome_list):
	print ("Running orthomclAdjustFasta")
	Proteome_code_list =  []
	for proteome in Proteome_list:
		if arg.code:
			code_name = protome.split(".")[0]
		else:
			code_temp = proteome.split(Name_separator)
			code_name = code_temp[0][:2]+code_temp[1][:3]
		Proteome_code_list.append(code_name)
		print ("\tProcessing proteome "+proteome)
		subprocess.Popen(["orthomclAdjustFasta "+code_name+" "+proteome+" "+str(Unique_ID)],shell=True).wait()
		ID_duplicateCheck(code_name+".fasta")
		prepBlastDB(code_name+".fasta")
	subprocess.Popen(["mkdir compliantFasta/ originalFasta/"],shell=True).wait()
	for code,file_i in zip(Proteome_code_list,Proteome_list):
		subprocess.Popen(["mv "+code+".fasta compliantFasta/"],shell=True).wait()
		subprocess.Popen(["mv "+file_i+" originalFasta/"],shell=True).wait()
	return Proteome_code_list

def checkFasta(Proteome_list):
	print("Check fasta files for duplicates and errors")
	for proteome in Proteome_list:
		loading (Proteome_list.index(proteome),len(Proteome_list),"Checking proteomes: ",50,proteome)
		ID_duplicateCheck(proteome)

def filterFasta():
	print ("Filtering proteome fasta files")
	subprocess.Popen(["\torthomclFilterFasta compliantFasta/ "+str(min_length)+" "+str(max_percent_stop)],shell=True).wait()

def makeBlastDB(goodproteins):
	print("Making BLAST database")
	subprocess.Popen(["makeblastdb -in "+goodproteins+" -dbtype prot -parse_seqids -input_type fasta -out "+database_name],shell=True).wait()
	
def allvsallBLAST(goodproteins,goodproteins_db):
	print ("BLASTing all the way (may take a while...)")
	subprocess.Popen(["blastp -query "+goodproteins+" -db "+database_name+" -out "+blast_out_name+" -evalue "+evalue_cutoff+" -outfmt 6 -num_threads "+CPUs],shell=True).wait()
	
def blastParser(blast_out):
	print ("Parsing BLAST output")
	subprocess.Popen(["orthomclBlastParser "+blast_out_name+" compliantFasta/ >> similarSequences.txt"],shell=True).wait()
	
def loadBlast(config_file):
	print ("Loading BLAST output into orthoMCL database")
	subprocess.Popen(["orthomclLoadBlast "+config_file+" similarSequences.txt"],shell=True).wait()
	
def pairs(config_file):
	print ("Finding pairs for orthoMCL")
	subprocess.Popen(["orthomclPairs "+config_file+" pairs.log cleanup=yes"],shell=True).wait()
	
def dumpPairs(config_file):
	print ("Dump files from the database produced by the orthomclPairs program")
	subprocess.Popen(["orthomclDumpPairsFiles "+config_file],shell=True).wait()

def mcl(inflation):
	print ("Running mcl algorithm")
	subprocess.Popen(["mcl mclInput --abc -I "+inflation+" -o mclOutput"],shell=True).wait()
	
def mclGroups(prefix,start_ID,groups_file):
	print ("Dumping groups")
	subprocess.Popen(["orthomclMclToGroups "+prefix+" "+start_ID+" < mclOutput > "+groups_file],shell=True).wait()


if arg.adjust:
	adjustFasta(arg.infile)
	
elif arg.no_adjust:
	filterFasta()
	makeBlastDB("goodProteins.fasta")
	allvsallBLAST("goodProteins.fasta","goodProteins.fasta_db")
	blastParser(blast_out_name)
	loadBlast(config_file)
	pairs(config_file)
	dumpPairs(config_file)
	mcl(inflation)
	mclGroups(prefix,start_ID,groups_file)
	
elif arg.check:
	checkFasta(arg.infile)
elif arg.normal:
#	adjustFasta(arg.infile)
#	filterFasta()
#	makeBlastDB("goodProteins.fasta")
	allvsallBLAST("goodProteins.fasta","goodProteins.fasta_db")
	blastParser(blast_out_name)
	loadBlast(config_file)
	pairs(config_file)
	dumpPairs(config_file)
	mcl(inflation)
	mclGroups(prefix,start_ID,groups_file)

