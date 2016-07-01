#!/usr/bin/env python2
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

import sys
from process.error_handling import *
from collections import OrderedDict

dna_chars = ["a", "t", "g", "c"]

aminoacid_table = OrderedDict({"a": ["Alanine", "nonpolar", "neutral"],
                               "r": ["Arginine", "polar", "positive"],
                               "n": ["Asparagine", "polar", "neutral"],
                               "d": ["Aspartate", "polar", "negative"],
                               "c": ["Cysteine", "nonpolar", "neutral"],
                               "e": ["Glutamate", "polar", "negative"],
                               "q": ["Glutamine", "polar", "neutral"],
                               "g": ["Glycine", "nonpolar", "neutral"],
                               "h": ["Histidine", "polar", "neutral"],
                               "i": ["Isoleucine", "nonpolar", "neutral"],
                               "l": ["Leucine", "nonpolar", "neutral"],
                               "k": ["Lysine", "polar", "positive"],
                               "m": ["Methionine", "nonpolar", "neutral"],
                               "f": ["Phenylalanine", "nonpolar", "neutral"],
                               "p": ["Proline", "nonpolar", "neutral"],
                               "s": ["Serine", "polar", "neutral"],
                               "t": ["Threonine", "polar", "neutral"],
                               "w": ["Tryptophan", "nonpolar", "neutral"],
                               "y": ["Tyrosine", "polar", "neutral"],
                               "v": ["Valine", "nonpolar", "neutral"],
                               "u": ["Selenocysteine", "", ""],
                               "o": ["Pyrrolysine", "", ""],
                               "x": ["Missing", "", ""]})

iupac = {"ag": "r", "ct": "y", "cg": "s", "at": "w", "gt": "k", "ac": "m",
         "cgt": "b", "agt": "d", "act": "h", "acg": "v", "acgt": "n"}

iupac_rev = {"v": "acg", "r": "ag", "m": "ac", "s": "cg", "d": "agt",
             "b": "cgt", "n": "acgt", "h": "act", "y": "ct", "w": "at",
             "k": "gt"}

iupac_conv = {"v": "acg", "r": "ag", "m": "ac", "s": "cg", "d": "agt",
              "b": "cgt", "n": "acgt", "h": "act", "y": "ct", "w": "at",
              "k": "gt", "a": "aa", "t": "tt", "c": "cc", "g": "gg"}


def merger(ranges):
    """
    Generator that merges ranges in a list of tuples. For example,
    if ranges is [(1, 234), (235, 456), (560, 607), (607,789)]
    this generator will yield [(1, 456), (560, 789)]
    """
    previous = 0
    last_start = 0
    for st, en in ranges:
        if not previous:
            last_start = st
            previous = en
        elif st - 1 == previous:
            previous = en
        else:
            yield last_start, previous
            previous = en
            last_start = st
    else:
        yield last_start, en


def has_colours(stream):
    if not hasattr(stream, "isatty"):
        return False
    if not stream.isatty():
        return False # auto color only on TTYs
    try:
        import curses
        curses.setupterm()
        return curses.tigetnum("colors") > 2
    except:
        # guess false in case of error
        return False

# Support for terminal colors
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
has_colours = has_colours(sys.stdout)


def print_col(text, color, i=0):
    p = ["TriSeq", "OrthoMCl Pipeline", "TriStats", "TriOrtho"]
    suf = {GREEN: "[%s] " % p[i], YELLOW: "[%s-Warning] " % p[i],
           RED: "[%s-Error] " % p[i]}
    if has_colours:
        seq = "\x1b[1;%dm" % (30 + color) + suf[color] + "\x1b[0m" + text
        print(seq)
    else:
        print(text)

    if color is RED:
        raise SystemExit


class Base:

    def autofinder(self, reference_file):
        """ Autodetect the type of file to be parsed. Based on headers """
        sequence = ""
        file_handle = open(reference_file, "r")

        # If input file is not a simple text file, which means it"s invalid,
        # handle this exception
        try:
            header = file_handle.readline()
        except UnicodeDecodeError:
            return InputError("Invalid input file.")

        # Skips first empty lines, if any
        while header.startswith("\n"):
            header = next(file_handle)

        # Recognition of NEXUS files is based on the existence of the string
        # "#NEXUS" in the first non-empty line
        if header.upper().strip().startswith("#NEXUS"):
            autofind = "nexus"
            while True:
                line = file_handle.readline()
                if line.strip().lower() == "matrix":
                    next_line = file_handle.readline()
                    sequence = "".join(next_line.split()[1:]).strip()
                    break

        # Recognition of Stockhold files is based on the existence of the
        # string "# stockholm" in the first non-empty line (case insensitive)
        elif header.upper().strip().startswith("# STOCKHOLM") or \
                header.upper().strip().startswith("#STOCKHOLM"):
            autofind = "stockholm"
            while True:
                line = file_handle.readline()
                if not line.startswith("#") and line.strip() != "":
                    sequence = line.split()[1]
                    break

        # Recognition of FASTA or .loci files is based on the existence of a ">"
        # character as the first character of a non-empty line
        elif header.strip().startswith(">"):
            next_line = next(file_handle)
            if next_line.strip().startswith(">"):
                autofind = "loci"
                sequence = header.split()[-1].strip()
            else:
                autofind = "fasta"
                sequence = next_line.strip()
                for line in file_handle:
                    if line.strip() != "" and line.strip()[0] != ">":
                        sequence += line.strip()
                    elif line.strip() != "" and line.strip()[0] == ">":
                        break

        # Recognition of Phylip files is based on the existence of two
        # integers separated by whitespace on the first non-empy line
        elif len(header.strip().split()) == 2 and header.strip().split()[0]\
                .isdigit() and header.strip().split()[1].isdigit():

            autofind = "phylip"
            sequence = "".join(file_handle.readline().split()[1:]).strip()

        # Check if there is any sequence. If not, the alignment file has no
        # sequence
        else:
            return InputError("Unknown input file format.")
        if sequence.replace("-", "") == "":
            print("\nAlignment file %s has no sequence or the first sequence "
                  "is empty. Please check the file." % reference_file)
            raise SystemExit

        # Guessing the genetic code
        code = self.guess_code(sequence)

        return autofind, code

    def get_loci_taxa(self, loci_file):
        """
        Gets a taxa list from a loci file. This is required prior to parsing
        the alignment in order to correctly add missing data when certain
        taxa are not present in a locus
        :param loci_file: string, path to loci file
        """

        file_handle = open(loci_file)
        taxa_list = []

        for line in file_handle:
            if line.strip().startswith(">"):
                taxon = line.strip().split()[0][1:]
                if taxon not in taxa_list:
                    taxa_list.append(taxon)

        return taxa_list

    def partition_format(self, partition_file):
        """ Tries to guess the format of the partition file (Whether it is
        Nexus of RAxML"s) """
        file_handle = open(partition_file)

        # Skips first empty lines, if any
        header = file_handle.readline()
        while header.startswith("\n"):
            header = next(file_handle)

        fields = header.split()
        if fields[0].lower() == "charset":
            p_format = "nexus"
        else:
            p_format = "phylip"

        return p_format

    def guess_code(self, sequence):
        """ Function that guesses the code of the molecular sequences (i.e.,
        DNA or Protein) based on the first sequence of a reference file """

        # Removes gaps from the sequence so that the frequencies are not biased
        sequence = sequence.upper().replace("-", "")

        dna_count = sequence.count("A") + sequence.count("T") + \
                    sequence.count("G") + sequence.count("C") + \
                    sequence.count("N")
        dna_proportion = float(dna_count) / float(len(sequence))
        if dna_proportion > 0.9:  # The 0.9 cut-off has been effective so far
            code = ("DNA", "n")
        else:
            code = ("Protein", "x")
        return code

    def rm_illegal(self, string):
        """ Function that removes illegal characters from taxa names """

        # Additional illegal characters are added here
        illegal_chars = [":", ",", ")", "(", ";", "[", "]", """, """]

        clean_name = "".join([char for char in string if char not in
                              illegal_chars])

        return clean_name

    def duplicate_taxa(self, taxa_list):
        """ Function that identifies duplicated taxa """
        import collections
        duplicated_taxa = [x for x, y in collections.Counter(taxa_list).items()
                           if y > 1]
        return duplicated_taxa

    def check_format(self, input_alignment, alignment_format):
        """ This function performs some very basic checks to see if the format
         of the input file is in accordance to the input file format
         specified when the script is executed """
        input_handle = open(input_alignment)
        line = input_handle.readline()
        while line.strip() == "":
            line = next(input_handle)

        if alignment_format == "fasta":
            if line.strip()[0] != ">":
                print("File not in Fasta format. First non-empty line of the"
                      " input file %s does not start with ">". Please verify "
                      "the file, or the input format settings\nExiting..." %
                      input_alignment)
                raise SystemExit

        elif alignment_format == "nexus":
            if line.strip().lower() != "#nexus":
                print("File not in Nexus format. First non-empty line of the"
                      " input file %s does not start with "#NEXUS". Please "
                      "verify the file, or the input format settings\n"
                      "Exiting..." % input_alignment)
                raise SystemExit

        elif alignment_format == "phylip":
            try:
                header = line.strip().split()
                int(header[0])
                int(header[1])
            except:
                print("File not in correct Phylip format. First non-empty "
                      "line of the input file %s does not start with two "
                      "integers separated by whitespace. Please verify the "
                      "file, or the input format settings\nExiting..." %
                      input_alignment)
                raise SystemExit

    def check_sizes(self, alignment_dic, current_file):
        """ This will make two sanity checks of the alignment contained in
        the alignment_dic object: First, it will check if none of the
        sequences is empty; If True, it will raise an error informing which
        taxa have empty sequences. If False, this will also test whether all
        sequences are of the same size and, if not, which are different """

        # Checking for taxa with empty sequences
        empty_taxa = []
        seq_list = []
        for taxa, f in alignment_dic.items():

            with open(f) as fh:
                seq = "".join(fh.readlines())

            if seq == "":
                empty_taxa.append(taxa)
            else:
                seq_list.append(len(seq))

        if empty_taxa is []:

            print("\nInputError: The following taxa contain empty sequences "
                  "in the file %s: %s\nPlease verify and re-run the program. "
                  "Exiting...\n" % (current_file, " ".join(empty_taxa)))
            raise SystemExit

        # Checking sequence lengths
        # Determine the most common length
        # commonseq = max(set([v for v in alignment_dic.values()]),
        #                 key=[v for v in alignment_dic.values()].count)
        # # Creates a dictionary with the sequences, and respective length,
        # # of different length
        # diflength = dict((key, len(value)) for key, value in alignment_dic.items()
        #                  if len(commonseq) != len(value))

        if len(set(seq_list)) > 1:
            return False
        else:
            return True

    def read_basic_csv(self, file_handle):
        """ This will parse a simples csv file with only one column and one
        or more lines. It returns a list containing the contents of each line
        as an element. This can be used by any class/function of the process
        submodule granted that a previous check for the presence of the file
        is made """

        storage = []

        for line in file_handle:

            storage.append(line.strip())

        return storage


class Progression():

    def record(self, name, obj_size, window_size=50):

        self.name = name
        self.size = obj_size
        self.width = window_size

    def progress_bar(self, position):
        """ this function requires the record method to be previously defined,
        as it will need its attributes. It will print a progress bar with a
        specified weight according to the position on the current data
        structure """

        # If there is a previous message in the output, erase it
        try:
            self.msg
            sys.stdout.write("\r" + " " * len(self.msg))
        except AttributeError:
            pass

        # The progress bar
        position_proportion = int((position / self.size) * self.width)

        msg = "\r%s [%s%s] %s%%" % (self.name, "#" * position_proportion,
                                    "-" * (self.width - position_proportion),
                                    int((position_proportion / self.width) *
                                        100))

        # Erase the last message
        if int((position_proportion / self.width) * 100) == 100:
            sys.stdout.write("\r" + " " * len(msg))

    def write(self, msg):
        """ This will simply write a provided string to the terminal """

        self.msg = msg


__author__ = "Diogo N. Silva"
