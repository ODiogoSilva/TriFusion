#!/usr/bin/python3

import os
import sys
import codecs
import subprocess
import shutil

#For systems without argparse installed, provide the path to the module
#sys.path.append("/home/diogo/Python/Modules")

import argparse
from ortho import OrthomclToolbox as OT

parser = argparse.ArgumentParser(description="Pipeline for the OrthoMCL "
                                 "software")
parser.add_argument("-in", dest="infile", type=str, help="Provide the path to "
                    "the directory containing the proteome files")
parser.add_argument("-a", action="store_const", const=True, dest="adjust",
                    help="Run only the adjust_fasta program")
parser.add_argument("-na", action="store_const", const=True, dest="no_adjust",
                    help="Do not run only the adjust_fasta program")
parser.add_argument("-c", action="store_const", const=True, dest="code",
                    help="The proteome file names are already in "
                    "code (e.g. Homo_sapiens.fas -> HoSap.fas). Not advisable"
                    " because it will overwrite the originals!")
parser.add_argument("-p", action="store_const", const=True, dest="check",
                    help="Checks for duplicates and other potential errors")
parser.add_argument("-n", action="store_const", const=True, dest="normal",
                    help="Normal run of the pipeline")

arg = parser.parse_args()

## PARAMETER DEFINITIONS ##

# Configuration file for orthocml
config_file = "orthomcl.config"

# Output directory
output_dir = "./"

name_separator = "_"  # Specify the name separator in the input files (e.g.,
                      # the separator is "_" if the file name is
                      # Homo_sapiens.fasta). This parameter only applies if
                      # the file names are not already in code format (e.g.,
                      # for the "Homo_sapiens" a code format could be "HoSap").

# For filter_fasta
min_length = 10  # Minimum allowed length of proteins.  (suggested: 10)
max_percent_stop = 20  # Maximum percent stop codons.  (suggested 20)

# For allvsall BLAST or USEARCH (recommended)
database_name = "goodProteins_db"
usearch_out_name = "AllVsAll.out"
evalue_cutoff = "0.00001"  # 1E-5
CPUs = "4"  # Number of CPU's for multiprocessing

# For mcl
inflation = ["1.5", "2", "3", "4", "5"]

# For mcl_groups
prefix = "Test"  # Arbitrary string for the name of the groups
start_ID = "1000"  # Starting number for the groups
groups_file = "groups"


def loading(current_state, size, prefix_txt, width, proteome):
    """ Function that prints the loading progress of the script """
    percentage = int((current_state / size) * 100)
    complete = int(width * percentage * 0.01)
    sys.stdout.write("\r%s [%s%s] %s%%  %s" % (prefix_txt, "#" * complete,
                                               "." * (width - complete),
                                               percentage, proteome))
    sys.stdout.flush()


def install_schema(cfg_file, verbose=False):
    """ Install the schema for the mySQL database """

    if verbose:
        print("Installing mySQL schema")

    x = subprocess.Popen(["orthomclInstallSchema", cfg_file]).wait()


def check_unique_field(proteome_file):
    """  Checks the original proteome file for a field in the fasta header
     that is unique to all sequences"""

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

        if len(temp_list) == len(set(temp_list)) and len(set(temp_list)) ==\
                len(header_list):

            # The orthoMCL program uses an index starting from 1, so the +1 is
            #  a necessary adjustment
            return i

    else:
        return None


def prep_fasta(proteome_file, code, unique_id, verbose=False):

    if verbose:
        print("\t Preparing file for USEARCH")

    # Storing header list to check for duplicates
    header_list = []

    # Will prevent writing
    lock = True

    # File handles
    file_in = open(proteome_file)
    file_out = open(proteome_file.split(".")[0] + "_mod.fas", "w")

    for line in file_in:
        if line.startswith(">"):
            if line not in header_list:
                header_list.append(line)
                fields = line.split("|")
                file_out.write(">%s|%s\n" % (code, fields[unique_id]))
                lock = True
            else:
                lock = False
        elif lock:
            file_out.write(line)

    # Close file handles:
    file_in.close()
    file_out.close()


def adjust_fasta(proteome_files, code=True, verbose=False):

    if verbose:
        print("Running orthomcladjust_fasta")

    # Create compliant fasta directory
    if not os.path.exists(os.path.join(output_dir, "compliantFasta")):
        os.makedirs(os.path.join(output_dir, "compliantFasta"))

    for proteome in proteome_files:
        # Get code for proteome
        code_name = proteome.split(os.path.sep)[-1].split(".")[0]

        # Check the unique ID field
        unique_id = check_unique_field(proteome)

        # Adjust fasta
        prep_fasta(proteome, code_name, unique_id)

        protome_file_name = proteome.split(os.path.sep)[-1].split(".")[0] + \
                            ".fasta"

        shutil.move(proteome.split(".")[0] + "_mod.fas",
                    os.path.join(output_dir, "compliantFasta",
                                 protome_file_name))


def check_fasta(proteome_list):
    print("Check fasta files for duplicates and errors")
    for proteome in proteome_list:
        loading(proteome_list.index(proteome), len(proteome_list),
                "Checking proteomes: ", 50, proteome)


def filter_fasta(min_len, max_stop, verbose=False):

    if verbose:
        print("Filtering proteome fasta files")

    x = subprocess.Popen(["orthomclFilterFasta", "compliantFasta", str(min_len),
                      str(max_stop)]).wait()


def allvsall_usearch(goodproteins, eval, cpus, usearch_outfile, verbose=False):

    if verbose:
        print("Perfoming USEARCH allvsall")

    x = subprocess.Popen(["usearch", "-ublast", goodproteins, "-db",
                          goodproteins, "-blast6out", usearch_outfile,
                          "-evalue", str(eval), "--maxaccepts", "0",
                          "-threads", str(cpus)]).wait()


def blast_parser(usearch_ouput, verbose=False):

    if verbose:
        print("Parsing BLAST output")

    x = subprocess.Popen(["orthomclBlastParser " + usearch_ouput +
                      " compliantFasta/ >> similarSequences.txt"],
                      shell=True).wait()


def remove_duplicate_entries(verbose=False):

    if verbose:
        print("Removing possible dupplicate entries")

    shutil.move("similarSequences.txt", "similarSequences.txt.old")

    file_handle = open("similarSequences.txt.old")
    output_handle = open("similarSequences.txt", "w")

    storage = {}

    for line in file_handle:
        fields = line.split()
        id1 = fields[0]
        id2 = fields[1]

        if id1 not in storage:
            storage[id1] = [id2]
        else:
            if id2 not in storage[id1]:
                output_handle.write(line)
                storage[id1].append(id2)

    file_handle.close()
    output_handle.close()
    os.remove("similarSequences.txt.old")


def load_blast(cfg_file, verbose=False):

    if verbose:
        print("Loading BLAST output into orthoMCL database")

    x = subprocess.Popen(["orthomclLoadBlast", cfg_file,
                      "similarSequences.txt"]).wait()


def pairs(cfg_file, verbose=False):

    if verbose:
        print("Finding pairs for orthoMCL")

    x = subprocess.Popen(["orthomclPairs", cfg_file,
                      "pairs.log", "cleanup=yes"]).wait()


def dump_pairs(cfg_file, verbose=False):

    if verbose:
        print("Dump files from the database produced by the orthomclPairs "
              "program")

    x = subprocess.Popen(["orthomclDumpPairsFiles", cfg_file]).wait()


def mcl(inflation_list, verbose=False):

    if verbose:
        print("Running mcl algorithm")

    for val in inflation_list:
        x = subprocess.Popen(["mcl mclInput --abc -I " + val + " -o mclOutput_" +
                          val.replace(".", "")], shell=True).wait()


def mcl_groups(inflation_list, mcl_prefix, start_id, group_file, verbose=False):

    if verbose:
        print("Dumping groups")

    for val in inflation_list:
        x = subprocess.Popen(["orthomclMclToGroups " + mcl_prefix + " " +
                          start_id + " < mclOutput_" + val.replace(".", "")
                          + " > " + group_file + "_" + str(val) + ".txt"],
                         shell=True).wait()


def export_filtered_groups(inflation_list, group_prefix, gene_t, sp_t, db):

    stats_storage = {}

    for val in inflation_list:
        group_obj = OT.Group(group_prefix + "_%s.txt" % val, gene_t, sp_t)
        stats = group_obj.export_filtered_group(
            output_file_name="Filtered_groups_%s.txt" % val)
        group_obj.retrieve_fasta(db)
        stats_storage[val] = stats

    return stats


if __name__ == '__main__':

    # Get proteome files
    proteome_files = os.listdir(arg.infile)

    # Change working directory
    os.chdir(output_dir)

    if arg.adjust:
        adjust_fasta(proteome_files, verbose=True)

    elif arg.no_adjust:
        install_schema(config_file, verbose=True)
        filter_fasta(output_dir, min_length)
        allvsall_usearch("goodProteins.fasta")
        blast_parser()
        remove_duplicate_entries()
        load_blast()
        pairs()
        dump_pairs()
        mcl()
        mcl_groups(prefix, start_ID, groups_file)

    elif arg.check:
        check_fasta(arg.infile)
    elif arg.normal:
        install_schema(config_file, verbose=True)
        adjust_fasta(proteome_files, verbose=True)
        filter_fasta(min_length, max_percent_stop)
        allvsall_usearch("goodProteins.fasta")
        blast_parser()
        remove_duplicate_entries()
        load_blast(config_file)
        pairs(config_file)
        dump_pairs(config_file)
        mcl()
        mcl_groups(prefix, start_ID, groups_file)
