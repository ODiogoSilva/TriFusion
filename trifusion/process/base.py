#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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

try:
    from process.error_handling import InputError, EmptyAlignment
except ImportError:
    from trifusion.process.error_handling import InputError, EmptyAlignment

import os
import sys
import time
import shutil
import traceback
from collections import OrderedDict, Counter

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


class CleanUp(object):

    def __init__(self, func):
        self.func = func
        self.temp_dir = ".trifusion-temp"
        self.idx = 0 if self.func.__name__ == "main_parser" else 2

    def __call__(self, *args):

        try:
            start_time = time.time()
            self.func(*args)
            if not args[0].quiet:
                print_col("Program execution successfully completed in %s "
                          "seconds" %
                          (round(time.time() - start_time, 2)), GREEN,
                          self.idx)
        except KeyboardInterrupt:
            # Removing temporary directory, if any
            if os.path.exists(self.temp_dir):
                shutil.rmtree(self.temp_dir)
            print_col("Interrupting, by your command", RED, self.idx)
        # The broad exception handling is used to remove the temporary
        # directory under any circumstances
        except Exception as e:
            print(e)
            traceback.print_exc()
            # Removing temporary directory, if any
            if os.path.exists(self.temp_dir):
                shutil.rmtree(self.temp_dir)

            if not args[0].quiet:
                print_col("Program exited with errors!", RED, self.idx)

        # Removing temporary directory, if any
        if os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)


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


def print_col(text, color, i=0, quiet=False):
    if not quiet:
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


class Base(object):
    """
    Collection of basic methods that can be inherited by Alignment-like
    objects
    """

    def autofinder(self, reference_file):
        """
        Autodetects format, missing data symbol and sequence type
        
        Parameters
        ----------
        reference_file : str
            Path to sequence file
        
        Returns
        -------
        fmt : str
            File format of `reference_file`
        code : tuple
            The sequence type as string in first element and missing data
            symbol as string in the second element
        
        """

        file_handle = open(reference_file, "r")

        # Set to True when the format has been detected
        format_found = False

        # If input file is not a simple text file, which means it's invalid,
        # handle this exception
        try:
            header = file_handle.readline()
        except UnicodeDecodeError:
            return InputError("Invalid input file.")

        try:
            # Skips first empty lines, if any
            while header.startswith("\n"):
                header = next(file_handle)

            # Recognition of NEXUS files is based on the existence of the
            # string "#NEXUS" in the first non-empty line
            if header.upper().strip().startswith("#NEXUS"):
                fmt = "nexus"
                format_found = True
                while True:
                        line = next(file_handle)
                        if line.strip().lower() == "matrix":
                            next_line = next(file_handle)
                            while next_line.startswith("\n"):
                                next_line = next(file_handle)
                            sequence = "".join(next_line.split()[1:]).strip()
                            break

            # Recognition of Stockhold files is based on the existence of the
            # string "# stockholm" in the first non-empty line
            # (case insensitive)
            elif header.upper().strip().startswith("# STOCKHOLM") or \
                    header.upper().strip().startswith("#STOCKHOLM"):
                fmt = "stockholm"
                format_found = True
                while True:
                    line = file_handle.readline()
                    if not line.startswith("#") and line.strip() != "":
                        sequence = line.split()[1]
                        break

            # Recognition of FASTA or .loci files is based on the existence
            # of a ">" character as the first character of a non-empty line
            elif header.strip().startswith(">"):
                next_line = next(file_handle)
                if next_line.strip().startswith(">"):
                    fmt = "loci"
                    format_found = True
                    sequence = header.split()[-1].strip()
                else:
                    fmt = "fasta"
                    format_found = True
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

                fmt = "phylip"
                format_found = True
                sequence = "".join(file_handle.readline().split()[1:]).strip()

            # Recognition of ipyrad loci file, which does not start with a ">"
            # character
            if not format_found:
                while not header.startswith("//"):
                    try:
                        header = next(file_handle)
                    except StopIteration:
                        break
                    if header.startswith("//"):
                        if header.count("|") == 2:
                            fmt = "loci"
                            format_found = True
                            sequence = next(file_handle).split()[1]

            # Check if there is any sequence. If not, the alignment file has no
            # sequence
            if not format_found:
                return InputError("Unknown input file format.")
            if sequence.replace("-", "") == "":
                return EmptyAlignment("Alignment is empty")

            # Guessing the genetic code
            code = self.guess_code(sequence)

        except StopIteration:
            return EmptyAlignment("Alignment is empty")

        return fmt, code

    @staticmethod
    def get_loci_taxa(loci_file):
        """
        Get the list of taxa from a .loci file
        
        This is required prior to parsing the alignment in order to 
        correctly add missing data when certain taxa are not present in 
        a locus
        
        Parameters
        ----------
        loci_file : str
            Path to the .loci file.
        
        Returns
        -------
        taxa_list : list
            List of the taxon names as strings.
        
        
        """

        file_handle = open(loci_file)
        taxa_list = []

        for line in file_handle:
            # Skip empty lines and lines starting with //
            # The .lstrip(">") supports the new format of pyRAD's .loci
            # where taxon names start with a ">" character
            if not line.strip().startswith("//") and line.strip() != "":
                taxon = line.strip().split()[0].lstrip(">")
                if taxon not in taxa_list:
                    taxa_list.append(taxon)

        return taxa_list

    @staticmethod
    def guess_code(sequence):
        """
        Guess the sequence type, i.e. protein
        
        Guesses the code of the provided `sequence`, that is, if it is a
        DNA or Protein sequence based on the first sequence of 
        the reference file. The second item of the returned `code` variable 
        is a placeholder for the missing data symbol (can be either "?", "n"
        or "x"). This symbol will be evaluated during the Alignment parsing.

        Parameters
        ----------
        sequence : string
            Sequence string that will be used to guess the type
        
        Returns
        -------
        code : list
            Provides the (<sequence type>, <missing symbol character>).
            sequence type can be either "DNA" or "Protein". The
            missing symbol character may be "?", "n", "x" or None (See Notes).
        
        See also
        --------
        autofinder
        process.sequence.Alignment.read_alignment
        
        Notes
        -----
        The missing symbol character that is provided in the returned `code`
        tuple may be None if the missing character cannot be determined
        from the first reference sequence that is used to guess the 
        sequence type. If none of the potential missing data symbol is 
        present in this reference sequence, it will remain undetermined and
        will be evaluated during the alignment parsing. 
        
        """

        # Removes gaps from the sequence so that the frequencies are not biased
        sequence = sequence.upper().replace("-", "")

        # Try to evaluate whether missing data is in ? or N/X format. Since
        # the '?' character should only occur in sequences as missing data
        # (unlike the N character, which can be an amino acid residue),
        # we assume that if it exists in the sequence, that's the missing
        # data character
        if sequence.count("?"):
            missing = "?"
        else:
            missing = None

        dna_count = sequence.count("A") + sequence.count("T") + \
            sequence.count("G") + sequence.count("C") + \
            sequence.count("N")

        dna_proportion = float(dna_count) / float(len(sequence))

        # The 0.9 cut-off has been effective so far
        if dna_proportion > 0.9:
            if sequence.count("N"):
                missing = "n"
            code = ["DNA", missing]
        else:
            if sequence.count("X"):
                missing = "x"
            code = ["Protein", missing]
        return code

    @staticmethod
    def rm_illegal(taxon_string):
        """
        Removes illegal characters from taxon name
        
        The illegal characters are defined in the `illegal_chars`
        variable
        
        Parameters
        ----------
        taxon_string : str
            Taxon name
        
        Returns
        -------
        clean_name : str
            Taxon name without illegal characters
        
        """

        # Additional illegal characters are added here
        illegal_chars = [":", ",", ")", "(", ";", "[", "]", """, """]

        clean_name = "".join([char for char in taxon_string if char not in
                              illegal_chars])

        return clean_name

    @staticmethod
    def duplicate_taxa(taxa_list):
        """
        Identified duplicate items in a list
        
        Used to detect, for instance, duplicated taxa in the Alignment object
        
        Parameters
        ----------
        taxa_list : list
            List with taxon names
        
        Returns
        -------
        duplicated_taxa : list
            List with duplicated taxa from `taxa_list`
        """

        duplicated_taxa = [x for x, y in Counter(taxa_list).items()
                           if y > 1]

        return duplicated_taxa

    @staticmethod
    def read_basic_csv(file_handle):
        """
        Reads a basic CSV into a list
        
        Parses a simples CSV file with only one column and one
        or more lines while stripping whitespace.
        
        Parameters
        ----------
        file_handle : file object
            File object of the CSV file
        
        Returns
        -------
        result : list
            Contents of the CSV file with each line in a different entry
        
        """

        result = []

        for line in file_handle:

            result.append(line.strip())

        return result


class Progression(object):

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
