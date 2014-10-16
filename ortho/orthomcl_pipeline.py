#!/usr/bin/python3

import os
import sys
import codecs
import subprocess

#For systems without argparse installed, provide the path to the module
#sys.path.append("/home/diogo/Python/Modules")

import argparse

parser = argparse.ArgumentParser(description="Pipeline for the OrthoMCL software")
parser.add_argument("-in", dest="infile", type=str, help="Provide the path to the directory containing the proteome "
					"files")
parser.add_argument("-a", action="store_const", const=True, dest="adjust", help="Run only the adjust_fasta program")
parser.add_argument("-na", action="store_const", const=True, dest="no_adjust", help="Do not run only the adjust_fasta "
					"program")
parser.add_argument("-c", action="store_const", const=True, dest="code", help="The proteome file names are already in "
					"code (e.g. Homo_sapiens.fas -> HoSap.fas). Not advisable because it will overwrite the originals!")
parser.add_argument("-p", action="store_const", const=True, dest="check", help="Checks for duplicates and other "
					"potential errors")
parser.add_argument("-n", action="store_const", const=True, dest="normal", help="Normal run of the pipeline")

arg = parser.parse_args()

## PARAMETER DEFINITIONS ##

# Configuration file for orthocml
config_file = "orthomcl.config"

# For adjustFast
#Unique_ID = 1  	# Specify the field of the sequence header where the unique protein ID is. If the field is different
				# depending on the proteomes, perhaps it will be better to group proteomes with unique ID in the same
				# field number and run the script with the -a flag
Name_separator = "_"  	# Specify the name separator in the input files (e.g., the separator is "_" if the file name is
						# Homo_sapiens.fasta). This parameter only applies if the file names are not already in code
						# format (e.g., for the "Homo_sapiens" a code format could be "HoSap").

# For filter_fasta
min_length = 10  # Minimum allowed length of proteins.  (suggested: 10)
max_percent_stop = 20  # Maximum percent stop codons.  (suggested 20)

# For allvsall BLAST
database_name = "goodProteins_db"
blast_out_name = "AllVsAll_BLAST.out"
evalue_cutoff = "0.00001"  # 1E-5 - Recommended
CPUs = "15"  # Number of CPU's for multiprocessing

# For mcl
inflation = ["1.5", "2", "3", "4", "5"]

# For mcl_groups
prefix = "Basidiomycota"  # Arbitrary string for the name of the groups
start_ID = "1000"  # Starting number for the groups
groups_file = "groups.txt"


def loading(current_state, size, prefix, width, proteome):
	""" Function that prints the loading progress of the script """
	percentage = int((current_state / size) * 100)
	complete = int(width * percentage * 0.01)
	sys.stdout.write("\r%s [%s%s] %s%%  %s" % (prefix, "#" * complete, "." * (width - complete), percentage, proteome))
	sys.stdout.flush()


def install_schema():
	""" Install the schema for the mySQL database """
	print("Installing mySQL schema")
	subprocess.Popen(["orthomclInstallSchema " + config_file], shell=True).wait()


def id_duplicate_check(proteome):
	print("\t Checking duplicates for " + proteome)
	previous_ids, duplicate_ids = [], {}
	subprocess.Popen(["mv " + proteome + " " + proteome + ".old"], shell=True).wait()
	with codecs.open(proteome + ".old", "r", "cp1252") as proteome_read, open(proteome, "w") as proteome_out, \
			open("ID_check.log",
		"a") as ID_check_log:
		x = 1
		ID_check_log.write("## ID duplicate check for " + proteome + "\n\n")
		for line in proteome_read:
			if line.startswith(">") and line not in previous_ids:
				previous_ids.append(line)
				proteome_out.write(line)
			elif line.startswith(">") and line in previous_ids:
				duplicate_ids[line] = line[:-1] + "_" + str(x) + "\n"
				previous_ids.append(line[:-1] + "_" + str(x) + "\n")
				proteome_out.write(line[:-1] + "_" + str(x) + "\n")
				x += 1
			elif "(" in line or ")" in line:
				pass
			else:
				proteome_out.write(line)
		if duplicate_ids == {}:
			ID_check_log.write("No duplicates found!\n")
		else:
			for k, v in duplicate_ids.items():
				ID_check_log.write("Found duplicate entry: " + k[:-1] + "; replaced by " + v + "\n")
	subprocess.Popen(["rm " + proteome + ".old"], shell=True).wait()


def check_unique_field(proteome_file):
	"""  Checks the original proteome file for a field in the fasta header that is unique to all sequences"""

	# Some files may have uf8 encoding problems so I used codecs here
	file_handle = codecs.open(proteome_file, "r", "cp1252")
	header_list = []

	for line in file_handle:
		if line.startswith(">"):
			header = line[1:].strip()
			# Store header in list format
			header_list.append(header.split("|"))
	else:
		# Get the size of the header fields
		header_field_size = len(header.split("|"))

	for i in range(header_field_size):
		temp_list = []
		for header in header_list:
			temp_list.append(header[i])

		if len(temp_list) == len(set(temp_list)) and len(set(temp_list)) == len(header_list):
			# The orthoMCL program uses an index starting from 1, so the +1 is a necessary adjustment
			return i + 1

	else:
		return None


def prep_blastdb(proteome_file):
	print("\t Preparing file for BLAST database")
	subprocess.Popen(["mv " + proteome_file + " " + proteome_file + ".old"], shell=True).wait()
	with open(proteome_file + ".old", "r") as file_in, open(proteome_file, "w") as file_out:
		for line in file_in:
			if line.startswith(">"):
				file_out.write(">gnl|" + line[1:])
			else:
				file_out.write(line)
	subprocess.Popen(["rm " + proteome_file + ".old"], shell=True).wait()


def adjust_fasta(proteome_dir):
	print("Running orthomcladjust_fasta")
	proteome_file_list = os.listdir(proteome_dir)
	proteome_code_list = []

	for proteome in proteome_file_list:
		if arg.code:
			code_name = proteome.split(".")[0]
		else:
			code_temp = proteome.split(Name_separator)
			code_name = str(code_temp[0][:2] + code_temp[1][:3]).lower()

		proteome_code_list.append(code_name)

		unique_id = check_unique_field(proteome_dir + "/" + proteome)

		subprocess.Popen(["orthomclAdjustFasta " + code_name + " " + proteome_dir + "/" + proteome + " " +
						str(unique_id)], shell=True).wait()

		id_duplicate_check(code_name + ".fasta")
		prep_blastdb(code_name + ".fasta")

	subprocess.Popen(["mkdir compliantFasta/"], shell=True).wait()

	for code in proteome_code_list:
		subprocess.Popen(["mv " + code + ".fasta compliantFasta/"], shell=True).wait()

	return proteome_code_list


def check_fasta(proteome_list):
	print("Check fasta files for duplicates and errors")
	for proteome in proteome_list:
		loading(proteome_list.index(proteome), len(proteome_list), "Checking proteomes: ", 50, proteome)
		id_duplicate_check(proteome)


def filter_fasta():
	print("Filtering proteome fasta files")
	subprocess.Popen(["\torthomclFilterFasta compliantFasta/ " + str(min_length) + " " + str(max_percent_stop)],
					shell=True).wait()


def make_blastdb(goodproteins):
	print("Making BLAST database")
	subprocess.Popen(["makeblastdb -in " + goodproteins + " -dbtype prot -parse_seqids -input_type fasta -out "
					"" + database_name], shell=True).wait()


def allvsall_blast(goodproteins):
	print("BLASTing all the way (may take a while...)")
	subprocess.Popen(["blastp -query " + goodproteins + " -db " + database_name + " -out " + blast_out_name + " "
					"-evalue " + evalue_cutoff + " -outfmt 6 -num_threads " + CPUs], shell=True).wait()


def blast_parser():
	print("Parsing BLAST output")
	subprocess.Popen(["orthomclBlastParser " + blast_out_name + " compliantFasta/ >> similarSequences.txt"],
					shell=True).wait()


def load_blast():
	print("Loading BLAST output into orthoMCL database")
	subprocess.Popen(["orthomclLoadBlast " + config_file + " similarSequences.txt"], shell=True).wait()


def pairs():
	print("Finding pairs for orthoMCL")
	subprocess.Popen(["orthomclPairs " + config_file + " pairs.log cleanup=yes"], shell=True).wait()


def dump_pairs():
	print("Dump files from the database produced by the orthomclPairs program")
	subprocess.Popen(["orthomclDumpPairsFiles " + config_file], shell=True).wait()


def mcl():
	print("Running mcl algorithm")
	for val in inflation:
		subprocess.Popen(["mcl mclInput --abc -I " + val + " -o mclOutput_" + val.replace(".", "")], shell=True).wait()


def mcl_groups(prefix, start_id, groups_file):
	print("Dumping groups")
	for val in inflation:
		subprocess.Popen(["orthomclMclToGroups " + prefix + " " + start_id + " < mclOutput_" + val.replace(".", "") +
						  " > " + groups_file], shell=True).wait()


if arg.adjust:
	adjust_fasta(arg.infile)
	
elif arg.no_adjust:
	install_schema()
	filter_fasta()
	make_blastdb("goodProteins.fasta")
	allvsall_blast("goodProteins.fasta")
	blast_parser()
	load_blast()
	pairs()
	dump_pairs()
	mcl()
	mcl_groups(prefix, start_ID, groups_file)
	
elif arg.check:
	check_fasta(arg.infile)
elif arg.normal:
	install_schema()
	adjust_fasta(arg.infile)
	filter_fasta()
	make_blastdb("goodProteins.fasta")
	allvsall_blast("goodProteins.fasta")
	blast_parser()
	load_blast()
	pairs()
	dump_pairs()
	mcl()
	mcl_groups(prefix, start_ID, groups_file)
