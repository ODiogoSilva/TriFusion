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
#  Version: 0.1.1
#  Last update: 11/02/14

# TriFusion imports
import process
from process.base import *
from process.data import Partitions

# Other imports
import numpy as np
from collections import OrderedDict, Counter
import itertools
import re
import os
import shutil
from os import sep
from os.path import join, basename, splitext
from itertools import compress
import fileinput
import multiprocessing
import cPickle as pickle
# import pickle
# TODO: Create a SequenceSet class for sets of sequences that do not conform
# to an alignment, i.e. unequal length.This would eliminate the problems of
# applying methods designed for alignments to sets of sequences with unequal
# length would allows these sets of sequences to have methods of their own.
# However, the __init__ of the current Alignment object should apply to both
# SequenceSet and Alignment classes. So, I'll have to re-structure the code
# somehow.
# TODO After creating the SequenceSet class, an additional class should be
# used to make the triage of files to either the Alignment or SequenceSet
# classes


class AlignmentException(Exception):
    pass


class AlignmentUnequalLength(Exception):
    pass


def read_alns(l):
    """
    Worker function that reads alignments from a list of file names.
    :param l: list, [0] - List of file names
                    [1] - dest variable, where the temporary sequences files
                    will be saved
                    [2] - int, work number
                    [3] - Namespace object, to share counter and processed files
    """

    # Open handle for this worker, where individual Alignment objects will be
    #  stored
    fh = open(join(l[1], "Worker{}.pc".format(l[2])), "wb")

    for i in l[0]:
        # Update the counter for the progress bar. The index is used instead
        # of a simple counter, because race conditioning of multiple workers
        # does not produce the expected behaviour with a int counter.
        if l[3]:
            count = (l[0].index(i) + 1) * multiprocessing.cpu_count()
            l[3].progress = count if count > l[3].progress else l[3].progress
            # Update the file being processed
            l[3].m = "Processing file %s" % basename(i)
        # Pickle the Alignment object
        pickle.dump(Alignment(i, dest=l[1]), fh)

    fh.close()


class Alignment(Base):

    def __init__(self, input_alignment, input_format=None, alignment_name=None,
                 partitions=None, dest=None, locus_length=None):
        """
        The basic Alignment instance requires only an alignment file or an
        OrderedDict object. In case the class is initialized with a dictionary
        object information on the partitions must be provide using the
        partitions argument.

        :param input_alignment: string or OrderedDict. If string, it must be the
        input file name; if OrderedDict, it must contain the alignment in an
        ordered dictionary format where keys are taxon names and values are
        sequence strings

        :param input_format: string. Input format of the Alignment object. If
        input_alignment is a file name, the input format is automatically
        detected. If it is an OrderedDict, it must be specified

        :param alignment_name: string. It sets the self.name attribute only when
        input_alignment is an OrderedDict. Otherwise, the input file name will
        be used.

        :param partitions: Partitions object. If provided, it will overwrite
        self.partitions.

        :param dest: string. Path where the temporary sequence files will be
        stored

        :param locus_length: Manually sets the length of the locus. Usually
        provided when input_alignment is a dict object and the result of
        concatenation.
        """

        self.log_progression = Progression()

        """
        Initializing a Partitions instance for the current alignment. By
        default, the Partitions instance will assume that the whole Alignment
        object is one single partition. However, if the current alignment is the
        result of a concatenation, or if a partitions file is provided, the
        partitions argument can be used. Partitions can be later changed using
        the set_partitions method. Substitution models objects are associated
        with each partition and they default to None
        """
        if isinstance(partitions, Partitions):
            self.partitions = partitions
        else:
            self.partitions = Partitions()

        """
        The length of the alignment object. Even if the current alignment object
        is partitioned, this will return the length of the entire alignment
        """
        if not locus_length:
            self.locus_length = 0
        else:
            self.locus_length = locus_length

        """
        This option is only relevant when gaps are coded. This will store a
        string with the range of the restriction-type data that will encode gaps
        and will only be used when nexus is in the output format
        """
        self.restriction_range = None

        """
        Attribute informing whether the current object is an actual alignment
        (defined as a sequence set with sequences of the same length), in which
        case it is set to True, or a sequence set (with sequences of unequal
        length), in which case it is set to False. This is automatically set
        in the read_alignment method
        """
        self.is_alignment = None

        """
        The alignment attribute is an ordered dictionary where the keys are
        taxon names and the values are file names containing the path to the
        temporary sequence files.
        """
        self.alignment = {}

        """
        The e attribute will store any exceptions that occur during the
        parsing of the alignment object. It remains None unless something
        wrong happens.
        """
        self.e = None

        # """
        # This attribute is derived from self.alignment. To prevent permanent
        # modifications to the original alignment, this attribute will be
        # populated when the Process execution is triggered. In this event,
        # the original temporary .seq files will be duplicated, their file paths
        # added to self.action_alignment and all modifications will be made in
        # this attribute. After the execution is done, the duplicated files
        # will be removed and this attribute will return to its empty state.
        # """
        # self.action_alignment = OrderedDict()

        self.path = None

        # If the object is initialized with a string
        if isinstance(input_alignment, str):

            """
            Sets the alignment object name based on the input alignment file
            name. The name attribute will remove the extension and a preceding
            path, if it exists. The name attribute will retain the
            extension
            """
            self.path = input_alignment
            # Short name - No extension
            self.sname = basename(splitext(input_alignment)[0])
            # Full name - with extension
            self.name = basename(input_alignment)

            if dest:
                self.dest = dest
            else:
                self.dest = ""

            # Create temporary directory to store sequences for the
            # current alignment object
            if not os.path.exists(join(self.dest, self.sname)):
                os.makedirs(join(self.dest, self.sname))

            # Get alignment format and code. Sequence code is a tuple of
            # (DNA, N) or (Protein, X)
            finder_content = self.autofinder(input_alignment)
            # Handles the case where the input format is invalid and
            # finder_content is an Exception
            if isinstance(finder_content, Exception) is False:
                self.input_format, self.sequence_code = self.autofinder(
                    input_alignment)

                # In case the input format is specified, overwrite the attribute
                if input_format:
                    self.input_format = input_format

                # parsing the alignment
                self.read_alignment(input_alignment, self.input_format)
            else:
                # If the input file is invalid, self.alignment will be an
                # Exception instance instead of an OrderedDict()
                self.alignment = None
                # Setting the sequence code attribute for seq type checking
                # in AlignmentList
                self.sequence_code = None
                self.e = finder_content

        # In case the class is initialized with a dictionary object
        elif isinstance(input_alignment, OrderedDict):

            # The name of the alignment (str)
            self.name = alignment_name
            # Gets several attributes from the dictionary alignment
            self._init_dicobj(input_alignment)
            # The input format of the alignment (str)
            self.input_format = input_format
            # Short name - No extension
            try:
                self.sname = basename(splitext(alignment_name)[0])
            except AttributeError:
                pass
            self.dest = dest

    def __iter__(self):
        """
        Iterate over Alignment objects
        """

        for tx, f in self.alignment.items():
            with open(f) as fh:
                yield tx, "".join(fh.readlines())

    def _set_format(self, input_format):
        """
        Manually sets the input format associated with the Alignment object

        :param input_format: string. Input format. It can be one out of
        "fasta", "nexus" and "phylip"
        """

        self.input_format = input_format

    def _set_alignment(self, alignment_dict):
        """
        Sets a new alignment dictionary to the Alignment object. This may be
        useful when only the alignment dict of the object has to be modified
        through other objects/functions

        :param alignment_dict: OrderedDict, containing taxa names as keys
        and sequences as values
        """

        if isinstance(alignment_dict, OrderedDict):
            self.alignment = alignment_dict
        else:
            raise AlignmentException("Alignments manually added to the "
                                     "Alignment object must be OrderedDict")

    def _init_dicobj(self, dictionary_obj):
        """
        Internal method to set the alignment and other attributed when the
        Alignment object is instantiated with an OrderedDict

        :param dictionary_obj: OrderedDict, containing the taxa names as keys
        and corresponding sequences as values
        """

        with open(list(dictionary_obj.values())[0]) as fh:
            seq = "".join(fh.readlines())

        self.sequence_code = self.guess_code(seq)
        self.alignment = dictionary_obj

    def start_action_alignment(self):
        """
        Creates temporary action files and directs self.alignment values to them
        """

        # Copies temporary files
        for f in self.alignment.values():
            shutil.copyfile(f, splitext(f)[0] + ".temp")

        # Populates action_alignment
        self.alignment = OrderedDict(
            [(tx, splitext(f)[0] + ".temp")
             for tx, f in self.alignment.items()])

    def stop_action_alignment(self):
        """
        Removes temporary seq files for modification and restores the original
        alignment files
        """

        for f in self.alignment.values():
            os.remove(f)

        self.alignment = OrderedDict(
            [(tx, splitext(f)[0] + ".seq")
             for tx, f in self.alignment.items()])

    def sequences(self):
        """
        Generator for sequence data of the alignment object. Akin to
        values() method of a dictionary
        """

        for f in self.alignment.values():
            with open(f) as fh:
                yield "".join(fh.readlines())

    def get_sequence(self, taxon):
        """
        Convenince method that returns the sequences for a given taxon
        :param taxa: string, taxon name. Must be in self.alignments
        """

        if taxon in self.alignment:
            with open(self.alignment[taxon]) as fh:
                return "".join(fh.readlines())
        else:
            raise KeyError

    def remove_alignment(self, remove_active=False):
        """
        Removes all temporary sequence files for the alignment object
        :param remove_active: Boolean. If False, it will remove all temporary
        files (both active and stopped). If True, it will only remove the active
        """

        if remove_active:
            for f in self.alignment.values():
                os.remove(f)
        else:
            shutil.rmtree(join(self.dest, self.sname))

    def read_alignment(self, input_alignment, alignment_format,
                       size_check=True):
        """
        The read_alignment method is run when the class is initialized to
        parse an alignment and set all the basic attributes of the class.

        :param input_alignment: string. File name containing the input alignment

        :param alignment_format: string. Format of the input file. It can be
        one of three: "fasta", "nexus", "phylip"

        :param size_check: Boolean. If True it will check the size consistency
        of the sequences in the alignment
        """

        file_handle = open(input_alignment)

        # ======================================================================
        # PARSING PHYLIP FORMAT
        # ======================================================================

        if alignment_format == "phylip":
            # Get the number of taxa and sequence length from the file header
            header = file_handle.readline().split()
            self.locus_length = int(header[1])
            self.partitions.set_length(self.locus_length)
            for line in file_handle:
                try:
                    taxa = line.split()[0].replace(" ", "")
                    taxa = self.rm_illegal(taxa)

                    # Create temporary directory to store sequences for the
                    # current alignment object
                    if not os.path.exists(join(self.dest, self.sname)):
                        os.makedirs(join(self.dest, self.sname))

                    self.alignment[taxa] = join(self.dest, self.sname,
                                                taxa + ".seq")
                    try:
                        sequence = line.split()[1].strip().lower()
                    except IndexError:
                        sequence = ""
                    with open(self.alignment[taxa], "w") as fh:
                        fh.write(sequence)
                except IndexError:
                    pass

                    # TODO: Read phylip interleave

            # Updating partitions object
            self.partitions.add_partition(self.name, self.locus_length,
                                          file_name=self.path)

        # ======================================================================
        # PARSING FASTA FORMAT
        # ======================================================================
        elif alignment_format == "fasta":
            for line in file_handle:
                if line.strip().startswith(">"):
                    taxa = line[1:].strip()
                    taxa = self.rm_illegal(taxa)

                    self.alignment[taxa] = join(self.dest, self.sname,
                                                taxa + ".seq")

                elif line.strip() != "" and taxa:
                    with open(self.alignment[taxa], "a") as fh:
                        fh.write(line.strip().lower().
                                 replace(" ", "").replace("*", ""))

            with open(self.alignment[taxa]) as fh:
                self.locus_length = len("".join(fh.readlines()))

            self.partitions.set_length(self.locus_length)

            # Updating partitions object
            self.partitions.add_partition(self.name, self.locus_length,
                                          file_name=self.path)

        # ======================================================================
        # PARSING LOCI FORMAT
        # ======================================================================
        elif alignment_format == "loci":
            taxa_list = self.get_loci_taxa(self.path)

            # Create empty dict
            self.alignment = dict([(x, open(join(self.dest, self.sname,
                 x + ".seq"), "a")) for x in taxa_list])

            # Add a counter to name each locus
            locus_c = 1
            present_taxa = []
            for line in file_handle:
                if line.strip().startswith(">"):
                    fields = line.strip().split()
                    taxon = fields[0][1:]
                    present_taxa.append(taxon)
                    self.alignment[taxon].write(fields[1].lower())

                elif line.strip().startswith("//"):

                    locus_len = len(fields[1])
                    self.locus_length += locus_len

                    # Adding missing data
                    for tx in taxa_list:
                        if tx not in present_taxa:
                            self.alignment[tx].write(self.sequence_code[1] *
                                                     locus_len)

                    present_taxa = []

                    self.partitions.add_partition("locus_{}".format(locus_c),
                                                  locus_len,
                                                  file_name=self.path)
                    locus_c += 1

            for tx, fh in self.alignment.items():
                fh.close()
                self.alignment[tx] = join(self.dest, self.sname, tx + ".seq")

            self.partitions.set_length(self.locus_length)

        # ======================================================================
        # PARSING NEXUS FORMAT
        # ======================================================================
        elif alignment_format == "nexus":
            counter = 0
            for line in file_handle:
                # Skips the nexus header
                if line.strip().lower() == "matrix" and counter == 0:
                    counter = 1
                # Stop sequence parser here
                elif line.strip() == ";" and counter == 1:
                    counter = 2

                    # Close file handles
                    for tx, fh in self.alignment.items():
                        fh.close()
                        self.alignment[tx] = join(self.dest, self.sname,
                                                  tx + ".seq")

                    with open(self.alignment[tx]) as fh:
                        self.locus_length = len("".join(fh.readlines()))

                    self.partitions.set_length(self.locus_length)
                # Start parsing here
                elif line.strip() != "" and counter == 1:
                    taxa = line.strip().split()[0].replace(" ", "")
                    taxa = self.rm_illegal(taxa)
                    # This accommodates for the interleave format
                    if taxa in self.alignment:
                        self.alignment[taxa].write("".join(
                            line.strip().split()[1:]).lower())
                    else:
                        self.alignment[taxa] = open(join(self.dest, self.sname,
                                                taxa + ".seq"), "a")
                        self.alignment[taxa].write("".join(line.strip().
                                                           split()[1:]).lower())

                # If partitions are specified using the charset command, this
                # section will parse the partitions
                elif line.strip().startswith("charset"):
                    self.partitions.read_from_nexus_string(line,
                                                    file_name=self.path)

                # If substitution models are specified using the lset or prset
                # commands, this will parse the model parameters
                if line.lower().strip().startswith("lset") or \
                        line.lower().strip().startswith("prset"):
                    self.partitions.parse_nexus_model(line)

            # If no partitions have been added during the parsing of the nexus
            # file, set a single partition
            if self.partitions.partitions == OrderedDict():
                self.partitions.add_partition(self.name, self.locus_length,
                                              file_name=self.path)

        file_handle.close()

        # Checks the size consistency of the alignment
        if size_check is True:
            self.is_alignment = self.check_sizes(self.alignment,
                                                 input_alignment)
            if not self.is_alignment:
                self.e = AlignmentUnequalLength()
                return

        # Checks for duplicate taxa
        if len(list(self.alignment)) != len(set(list(self.alignment))):
            taxa = self.duplicate_taxa(self.alignment.keys())
            self.log_progression.write("WARNING: Duplicated taxa have been "
                                       "found in file %s (%s). Please correct "
                                       "this problem and re-run the program\n"
                                       % (input_alignment, ", ".join(taxa)))
            raise SystemExit

    def iter_taxa(self):
        """
        Generator for taxa names
        """

        for sp in self.alignment:
            yield sp

    def iter_sequences(self):
        """
        Generator for sequences
        """

        for seq in self.alignment.values():
            yield seq

    def remove_taxa(self, taxa_list_file, mode="remove"):
        """
        Removes specified taxa from the alignment. As taxa_list, this
        method supports a python list or an input csv file with a single
        column containing the unwanted species in separate lines. It
        currently supports two modes:
            ..:remove: removes the specified taxa
            ..:inverse: removes all but the specified taxa

        :param taxa_list_file: list/string. A list of taxa names or a csv file
        with taxa names in each line

        :param mode: string. Mode of execution. It can be either "remove" or
        "inverse
        """

        new_alignment = OrderedDict()

        def remove(list_taxa):
            for taxa, seq in self.alignment.items():
                if taxa not in list_taxa:
                    new_alignment[taxa] = seq
            self.alignment = new_alignment

        def inverse(list_taxa):
            for taxa, seq in self.alignment.items():
                if taxa in list_taxa:
                    new_alignment[taxa] = seq
            self.alignment = new_alignment

        # Checking if taxa_list is an input csv file:
        try:
            file_handle = open(taxa_list_file[0])

            taxa_list = self.read_basic_csv(file_handle)

        # If not, then the method's argument is already the final list
        except (IOError, IndexError):
            taxa_list = taxa_list_file

        if mode == "remove":
            remove(taxa_list)
        if mode == "inverse":
            inverse(taxa_list)

    def change_taxon_name(self, old_name, new_name):

        self.alignment = OrderedDict([(new_name, seq) if taxon == old_name
                                      else (taxon, seq) for taxon, seq in
                                      list(self.alignment.items())])

    def collapse(self, write_haplotypes=True, haplotypes_file=None,
                 haplotype_name="Hap", dest="./"):
        """
        Collapses equal sequences into haplotypes. This method changes
        the alignment variable and only returns a dictionary with the
        correspondence between the haplotypes and the original taxa names

        :param write_haplotypes: Boolean, If true, a haplotype list
        mapping the haplotype names file will be created for each individual
        input alignment.
        :param haplotypes_file: String, Name of the haplotype list mapping file
        referenced in write_haplotypes
        :param haplotype_name: String, Custom name of the haplotypes
        """

        collapsed_dic, correspondence_dic = OrderedDict(), OrderedDict()
        counter = 1

        for taxa, seq in self.alignment.items():
            if seq in collapsed_dic:
                collapsed_dic[seq].append(taxa)
            else:
                collapsed_dic[seq] = [taxa]

        self.alignment = OrderedDict()
        for seq, taxa_list in collapsed_dic.items():
            haplotype = "%s_%s" % (haplotype_name, counter)
            self.alignment[haplotype] = seq
            correspondence_dic[haplotype] = taxa_list
            counter += 1

        if write_haplotypes is True:
            # If no output file for the haplotype correspondence is provided,
            # use the input alignment name as reference
            if haplotypes_file is None:
                haplotypes_file = self.name.split(".")[0]
            self.write_loci_correspondence(correspondence_dic, haplotypes_file,
                                           dest)

    def consensus(self, consensus_type):
        """
        Converts the current Alignment object dictionary into a single consensus
         sequence. The consensus_type argument determines how variation in the
        original alignment is handled for the generation of the consensus
        sequence. The options are:
            ..:iupac: Converts variable sites according to the corresponding
            IUPAC symbols
            ..:soft mask: Converts variable sites into missing data
            ..:remove: Removes variable sites
            ..:first sequence: Uses the first sequence in the dictionary
        :param consensus_type: string, from the list above.
        """

        from process.base import iupac

        # If sequence type is first sequence
        if consensus_type == "First sequence":
            self.alignment = OrderedDict([("consensus",
                                           list(self.alignment.values())[0])])
            return

        # Empty consensus sequence
        consensus_seq = []

        for column in zip(*[self.get_sequence(x) for x in self.alignment]):

            column = list(set(column))

            # If invariable, include in consensus and skip. Increases speed
            # when there are no missing or gap characters
            if len(column) == 1:
                consensus_seq.append(column[0])
                continue

            # Remove missing data and gaps
            column = sorted([x for x in column if
                             x != self.sequence_code[1] and x != "-"])

            # In case there is only missing/gap characters
            if not column:
                consensus_seq.append(self.sequence_code[1])
                continue

            # Check for invariable without missing data
            if len(column) == 1:
                consensus_seq.append(column[0])
                continue

            if consensus_type == "IUPAC":
                # Convert variation into IUPAC code
                try:
                    consensus_seq.append(iupac["".join(column)])
                except KeyError:
                    iupac_code = sorted(set(list(itertools.chain.from_iterable(
                        [x if x in dna_chars else iupac_rev[x]
                         for x in column]))))
                    consensus_seq.append(iupac["".join(iupac_code)])

            elif consensus_type == "Soft mask":
                consensus_seq.append(self.sequence_code[1])
                continue

            elif consensus_type == "Remove":
                continue

        # Replace alignment attribute with consensus sequence
        self.remove_alignment(remove_active=True)
        self.alignment = OrderedDict([("consensus",
                                       join(self.dest, self.sname,
                                            "consensus.seq"))])
        with open(self.alignment["consensus"], "w") as fh:
            fh.write("".join(consensus_seq))

    @staticmethod
    def write_loci_correspondence(dic_obj, output_file, dest="./"):
        """
        This function supports the collapse method by writing the
        correspondence between the unique haplotypes and the loci into a
        new file
        """

        output_handle = open(join(dest, output_file + ".haplotypes"), "w")

        for haplotype, taxa_list in dic_obj.items():
            output_handle.write("%s: %s\n" % (haplotype, "; ".join(taxa_list)))

        output_handle.close()

    def _check_partitions(self, partition_obj):
        """
        Internal. Makes a consistency check for the self.partitions attribute
        """

        # Checks if total lenght of partitions matches the lenght of the
        # current alignment

        if partition_obj.counter != self.locus_length:
            return process.data.InvalidPartitionFile("Partitions in partition"
                   "file are inconsistency with current alignment")

    def set_partitions(self, partitions):
        """
        Updates the Partitions object of the current alignment.

        :param partitions: Partitions object. Use one of the Partition parsers
        to retrieve partitions information from files or python data structures.
        See process.data.Partitions documentation
        """

        # Checks partition's consistency
        er = self._check_partitions(partitions)

        if isinstance(er, process.data.InvalidPartitionFile):
            return er
        else:
            self.partitions = partitions

    def reverse_concatenate(self):
        """
        This function divides a concatenated file according to the
        partitions set in self.partitions and returns an AlignmentList object
        """

        concatenated_aln = AlignmentList([], dest=self.dest)
        alns = []

        for name, part_range in self.partitions:

            current_dic = OrderedDict()
            for taxon, f in self.alignment.items():
                seq = self.get_sequence(taxon)
                sub_seq = seq[part_range[0][0]:part_range[0][1] + 1]

                # If sub_seq is not empty (only gaps or missing data)
                if sub_seq.replace(self.sequence_code[1], "") != "":
                    if not os.path.exists(join(self.dest, name)):
                        os.makedirs(join(self.dest, name))

                    current_dic[taxon] = join(self.dest, name, taxon + ".seq")

                    with open(current_dic[taxon], "w") as fh:
                        fh.write(sub_seq)

            current_partition = Partitions()
            current_partition.add_partition(name, part_range[0][1] -
                                            part_range[0][0])

            # Check if current alignment object is not empty. This may occur
            # when reverse concatenating an alignment with a taxa subset
            # selected.
            if current_dic == OrderedDict():
                continue

            current_aln = Alignment(current_dic, input_format=self.input_format,
                                    partitions=current_partition,
                                    alignment_name=name, dest=self.dest)

            alns.append(current_aln)

        concatenated_aln.add_alignments(alns, ignore_paths=True)

        return concatenated_aln

    def code_gaps(self):
        """
        This method codes gaps present in the alignment in binary format,
        according to the method of Simmons and Ochoterena (2000), to be read
        by phylogenetic programs such as MrBayes. The resultant alignment,
        however, can only be output in the Nexus format
        """

        def gap_listing(sequence, gap_symbol="-"):
            """ Function that parses a sequence string and returns the
            position of indel events. The returned list is composed of
            tuples with the span of each indel """
            gap = "%s+" % gap_symbol
            span_regex = ""
            gap_list, seq_start = [], 0
            while span_regex is not None:
                span_regex = re.search(gap, sequence)
                if span_regex is not None and seq_start == 0:
                    gap_list.append(span_regex.span())
                    sequence = sequence[span_regex.span()[1] + 1:]
                    seq_start = span_regex.span()[1] + 1
                elif span_regex is not None and seq_start != 0:
                    gap_list.append((span_regex.span()[0] + seq_start,
                                     span_regex.span()[1] + seq_start))
                    sequence = sequence[span_regex.span()[1] + 1:]
                    seq_start += span_regex.span()[1] + 1
            return gap_list

        def gap_binary_generator(sequence, gap_list):
            """ This function contains the algorithm to construct the binary
             state block for the indel events """
            for cur_gap in gap_list:
                cur_gap_start, cur_gap_end = cur_gap
                if sequence[cur_gap_start:cur_gap_end] == "-" * \
                        (cur_gap_end - cur_gap_start) and \
                        sequence[cur_gap_start - 1] != "-" and \
                        sequence[cur_gap_end] != "-":
                    sequence += "1"

                elif sequence[cur_gap_start:cur_gap_end] == "-" * \
                     (cur_gap_end - cur_gap_start):

                    if sequence[cur_gap_start - 1] == "-" or \
                    sequence[cur_gap_end] == "-":
                        sequence += "-"

                elif sequence[cur_gap_start:cur_gap_end] != "-" * \
                        (cur_gap_end - cur_gap_start):
                    sequence += "0"
            return sequence

        complete_gap_list = []

        # Get the complete list of unique gap positions in the alignment
        for taxa, seq in self.alignment.items():

            current_list = gap_listing(seq)
            complete_gap_list += [gap for gap in current_list if gap not in
                                  complete_gap_list]

        # This will add the binary matrix of the unique gaps listed at the
        # end of each alignment sequence
        for taxa, seq in self.alignment.items():
            self.alignment[taxa] = gap_binary_generator(seq, complete_gap_list)

        self.restriction_range = "%s-%s" % (int(self.locus_length),
                                            len(complete_gap_list) +
                                            int(self.locus_length) - 1)

    def _filter_terminals(self):
        """
        This will replace the gaps in the extremities of the alignment with
        missing data symbols
        :return:
        """

        for taxa, f in self.alignment.items():

            with open(f) as fh:
                seq = "".join(fh.readlines())

            # Condition where the sequence only has gaps
            if not seq.strip("-"):
                with open(f, "w") as fh:
                    fh.write(str(self.sequence_code[1]) * len(seq))
                    continue

            trim_seq = list(seq)
            counter, reverse_counter = 0, -1

            while trim_seq[counter] == "-":
                trim_seq[counter] = self.sequence_code[1]
                counter += 1

            while trim_seq[reverse_counter] == "-":
                trim_seq[reverse_counter] = self.sequence_code[1]
                reverse_counter -= 1

            with open(f, "w") as fh:
                fh.write("".join(trim_seq))

    def _filter_columns(self, gap_threshold, missing_threshold):
        """ Here several missing data metrics are calculated, and based on
         some user defined thresholds, columns with inappropriate missing
         data are removed """

        taxa_number = float(len(self.alignment))
        self.old_locus_length = len(list(self.alignment.values())[0])

        filtered_cols = []

        fhs = fileinput.input(files=[x for x in self.alignment.values()])

        # Creating the column list variable
        for column in zip(*fhs):

            cadd = column.count

            # Calculating metrics
            gap_proportion = (float(cadd("-")) /
                              taxa_number) * float(100)
            missing_proportion = (float(cadd(self.sequence_code[1])) /
                                  taxa_number) * float(100)
            total_missing_proportion = gap_proportion + missing_proportion

            if total_missing_proportion <= gap_threshold or \
                    missing_proportion <= missing_threshold:

                filtered_cols.append(1)

            else:
                filtered_cols.append(0)

        fhs.close()

        # Open file handles in write mode
        for f in self.alignment.values():

            with open(f) as fh:
                seq = "".join(fh.readlines())

            with open(f, "w") as fh:
                fh.write("".join(compress(seq, filtered_cols)))

        self.locus_length = len(seq)

    def filter_missing_data(self, gap_threshold, missing_threshold,
                            fork=False):
        """
        Filters gaps and true missing data from the alignment using tolerance
        thresholds for each type of missing data. Both thresholds are maximum
        percentages of sites in an alignment column containing the type of
        missing data. If gap_threshold=50, for example, alignment columns with
        more than 50% of sites with gaps are removed.

        This method will create a new temporary sequence file for each taxa,
        with the modified sequence. The modified file will have the _filtered
        suffix. All temporary files should only be created when creating the
        output files and removed as soon as possible

        :param gap_threshold: int ranging from 0 to 100.
        :param missing_threshold: int ranging from 0 to 100.
        :param fork: boolean. If True, in addition to performing the
        filtering in the action_alignment it will also perform it on the
        original alignment and write the output file.
        """

        self._filter_terminals()
        self._filter_columns(gap_threshold, missing_threshold)

    def write_to_file(self, output_format, output_file, new_alignment=None,
                      seq_space_nex=40, seq_space_phy=30, seq_space_ima2=10,
                      cut_space_nex=50, cut_space_phy=258, cut_space_ima2=8,
                      interleave=False, gap="-", model_phylip=None,
                      outgroup_list=None, ima2_params=None, use_charset=True,
                      partition_file=True, output_dir=None,
                      phy_truncate_names=False, ld_hat=None):
        """ Writes the alignment object into a specified output file,
        automatically adding the extension, according to the output format
        This function supports the writing of both converted (no partitions)
        and concatenated (partitioned files). The choice between these modes
        is determined by the Partitions object associated with the Alignment
        object. If it contains multiple partitions, it will produce a
        concatenated alignment and the auxiliary partition files where
        necessary. Otherwise it will treat the alignment as a single partition.

        :param output_format: string. Format of the output file. It can be one
        of five: "fasta", "nexus", "phylip", "mcmctree" and "ima2"

        :param output_file: string. Name of the output file. It will overwrite
        files with the same name.

        :param new_alignment: OrderedDict. An option to provide an alternative
        alignment to write. It is set to None by default, in which case it
        uses self.alignment.

        :param interleave: Boolean. Determines whether the output alignment
        will be in leave (False) or interleave (True) format. Not all
        output formats support this option.

        :param gap: string. Symbol for gap data.

        :param model_phylip. string. Substitution model for the auxiliary
        partition file of phylip format, compliant with RAxML.

        :param outgroup_list. list. The outgroup_list argument is used only for
        Nexus output format and consists in writing a line defining the
        outgroup. This may be useful for analyses with MrBayes or other
        software that may require outgroups

        :param ima2_params: The ima2_params argument is used to provide
        information for the ima2 output format. If the argument is used,
        it should be in a list format and contain the following information:
          [[str, file_name containing the species and populations],
          [str, the population tree in newick format, e.g. (0,1):2],
          [mut_model:[str, mutational model for all alignments],
          [str, inheritance scalar]]

        :param use_charset: Boolean. If true, partitions from the Partitions
        object will be written in the nexus output format

        :param partition_file: Boolean. If true, the auxiliary partitions file
        will be writen.

        :param output_dir: String. If provided, the output file will be written
        on the specified path

        :param phy_truncate_names: Boolean. Whether names in phylip output
        format should be truncated to 10 characters or not.

        :param ld_hat: Boolean: If not None, the Fasta output format will
        include a first line compliant with the format of LD Hat and will
        truncate sequence names and sequence lenght per line accordingly.
        """

        # If this function is called in the AlignmentList class, there may
        # be a need to specify a new alignment dictionary, such as a
        # concatenated one
        if new_alignment is not None:
            alignment = new_alignment
        else:
            alignment = self.alignment

        # This will determine the default model value. GTR for nucleotides
        # and LG for proteins
        if not model_phylip:
            if self.sequence_code[0] == "DNA":
                model_phylip = "GTR"
            else:
                model_phylip = "LG"

        # If a specific output directory is provided, the output file will be
        # written there
        if output_dir:
            output_file = join(output_dir, output_file)

        # Checks if there is any other format besides Nexus if the
        # alignment's gap have been coded
        if self.restriction_range is not None:
            if output_format != ["nexus"]:
                self.log_progression.write("OutputFormatError: Alignments "
                                           "with gaps coded can only be written"
                                           " in Nexus format")
                return 0
        else:
            pass

        # Writes file in IMa2 format
        if "ima2" in output_format:

            population_file = ima2_params[0]
            population_tree = ima2_params[1]
            mutational_model = ima2_params[2]
            inheritance_scalar = ima2_params[3]

            # Get information on which species belong to each population from
            #  the populations file
            population_handle = open(population_file)
            population_storage = OrderedDict()
            for line in population_handle:
                taxon, population = line.strip().split()
                try:
                    population_storage[population.strip()].append(taxon)
                except KeyError:
                    population_storage[population.strip()] = [taxon]

            # Write the general header of the IMa2 input file
            out_file = open(output_file + ".txt", "w")
            # First line with general description
            out_file.write("Input file for IMa2 using %s alignments\n"
                        "%s\n"  # Line with number of loci
                        "%s\n"  # Line with name of populations
                        "%s\n"  # Line with population string
                        % (len(self.partitions.partitions),
                           len(population_storage),
                           " ".join(population_storage.keys()),
                           population_tree))

            if not self.partitions.is_single():
                # Write each locus
                for partition, lrange in self.partitions:

                    # Retrieving taxon names and sequence data. This step is
                    # the first because it will enable the removal of species
                    # containing only missing data.
                    new_alignment = []

                    # This temporary ordered dictionary is created so that
                    # the number of taxa per populations is corrected in
                    # each locus
                    current_locus_populations = OrderedDict(
                        (x, []) for x in population_storage)

                    for population, taxa_list in population_storage.items():
                        for taxon in taxa_list:
                            # This try statement catches common errors, such as
                            #  providing a species in the mapping file that
                            # does not exist in the alignment
                            try:
                                seq = self.get_sequence(taxon)[
                                    lrange[0][0]:lrange[0][1]].upper()
                            except KeyError:
                                print("Taxon %s provided in auxiliary "
                                      "population mapping file is not found "
                                      "in the alignment")
                                raise SystemExit

                            if seq.replace("N", "") != "":
                                new_alignment.append((taxon[:cut_space_ima2]
                                                      .ljust(seq_space_ima2),
                                                      seq))

                                current_locus_populations[population]\
                                    .append(taxon)

                    # Write the header of each partition
                    out_file.write("%s %s %s %s %s\n" % (
                        partition,
                        " ".join([str(len(x)) for x in
                                  list(current_locus_populations.values())]),
                        (lrange[0][1]) - lrange[0][0],
                        mutational_model,
                        inheritance_scalar))

                    # Write sequence data according to the order of the
                    # population mapping file
                    for taxon, seq in new_alignment:
                        out_file.write("%s%s\n" % (taxon, seq))

            if self.partitions.is_single():
                # Write the header for the single
                out_file.write("%s %s %s %s %s\n" % (
                               partition,
                               " ".join(population_storage.values()),
                               self.locus_length,
                               mutational_model,
                               inheritance_scalar))

                # Write sequence data
                for population, taxa_list in population_storage.items():
                    for taxon in taxa_list:
                        seq = self.get_sequence(taxon).upper()
                        out_file.write("%s%s\n" %
                            (taxon[:cut_space_ima2].ljust(seq_space_ima2), seq))

        # Writes file in phylip format
        if "phylip" in output_format:

            # Change taxa space if phy_truncate_names option is set to True
            if phy_truncate_names:
                cut_space_phy = 10

            out_file = open(output_file + ".phy", "w")
            out_file.write("%s %s\n" % (len(alignment), self.locus_length))
            if interleave:
                counter = 0
                for i in range(90, self.locus_length, 90):
                    for key, f in alignment.items():

                        with open(f) as fh:
                            seq = "".join(fh.readlines())[counter:i]

                        out_file.write("%s %s\n" % (
                            key[:cut_space_phy].ljust(seq_space_phy),
                            seq.upper()))
                    else:
                        out_file.write("\n")
                        counter = i
            else:
                for key, f in alignment.items():
                        with open(f) as fh:
                            seq = "".join(fh.readlines())
                        out_file.write("%s %s\n" % (
                            key[:cut_space_phy].ljust(seq_space_phy),
                            seq.upper()))

            # In case there is a concatenated alignment being written
            if not self.partitions.is_single() and partition_file:
                partition_file = open(output_file + "_part.File", "w")
                for name, lrange in self.partitions:
                    # Get model from app if it exists and there are no codon
                    # positions
                    model = model_phylip if self.partitions.models[name] == \
                        [None] or len(self.partitions.models[name][1]) > 1 \
                        else self.partitions.models[name][1][0]
                    partition_file.write("%s, %s = %s\n" % (
                                         model if model else
                                         self.sequence_code[0], name,
                                         "-".join([str(x + 1) for x in
                                                   lrange[0]])))
                partition_file.close()

            out_file.close()

        if "mcmctree" in output_format:

            out_file = open(output_file + "_mcmctree.phy", "w")
            taxa_number = len(self.alignment)

            if self.partitions.is_single() is False:
                for lrange in self.partitions.partitions.values():
                    lrange = lrange[0]
                    out_file.write("%s %s\n" % (taxa_number,
                                                (lrange[1] - (lrange[0]))))

                    for taxon, f in alignment.items():

                        with open(f) as fh:
                            seq = "".join(fh.readlines())

                        out_file.write("%s  %s\n" % (
                                       taxon[:cut_space_phy].ljust(
                                         seq_space_phy),
                                       seq[lrange[0]:lrange[1]].upper()))
            else:
                out_file.write("%s %s\n" % (taxa_number, self.locus_length))
                for taxon, f in alignment.items():

                    with open(f) as fh:
                        seq = "".join(fh.readlines())

                    out_file.write("%s  %s\n" % (
                                   taxon[:cut_space_phy].ljust(seq_space_phy),
                                   seq.upper()))

            out_file.close()

        # Writes file in nexus format
        if "nexus" in output_format:

            out_file = open(output_file + ".nex", "w")

            # This writes the output in interleave format
            if interleave:
                if self.restriction_range is not None:
                    out_file.write("#NEXUS\n\nBegin data;\n\tdimensions "
                                   "ntax=%s nchar=%s ;\n\tformat datatype="
                                   "mixed(%s:1-%s, restriction:%s) interleave="
                                   "yes gap=%s missing=%s ;\n\tmatrix\n" %
                                   (len(alignment),
                                    self.locus_length,
                                    self.sequence_code[0].upper(),
                                    self.locus_length - 1,
                                    self.restriction_range,
                                    gap,
                                    self.sequence_code[1]))
                else:
                    out_file.write("#NEXUS\n\nBegin data;\n\tdimensions "
                                   "ntax=%s nchar=%s ;\n\tformat datatype=%s "
                                   "interleave=yes gap=%s missing=%s ;\n\t"
                                   "matrix\n" %
                                   (len(alignment),
                                    self.locus_length,
                                   self.sequence_code[0].upper(),
                                    gap,
                                    self.sequence_code[1]))
                counter = 0
                for i in range(90, self.locus_length, 90):
                    for key, f in alignment.items():

                        with open(f) as fh:
                            seq = "".join(fh.readlines())

                        out_file.write("%s %s\n" % (
                                       key[:cut_space_nex].ljust(
                                         seq_space_nex),
                                       seq[counter:i].upper()))
                    else:
                        out_file.write("\n")
                        counter = i
                else:
                    for key, f in alignment.items():

                        with open(f) as fh:
                            seq = "".join(fh.readlines())

                        out_file.write("%s %s\n" % (
                                       key[:cut_space_nex].ljust(
                                         seq_space_nex),
                                       seq[i:self.locus_length].upper()))
                    else:
                        out_file.write("\n")
                out_file.write(";\n\tend;")

            # This writes the output in leave format (default)
            else:
                if self.restriction_range is not None:
                    out_file.write("#NEXUS\n\nBegin data;\n\tdimensions "
                                   "ntax=%s nchar=%s ;\n\tformat datatype=mixed"
                                   "(%s:1-%s, restriction:%s) interleave=yes "
                                   "gap=%s missing=%s ;\n\tmatrix\n" %
                                   (len(alignment),
                                    self.locus_length,
                                    self.sequence_code[0],
                                    self.locus_length - 1,
                                    self.restriction_range,
                                    gap,
                                    self.sequence_code[1]))
                else:
                    out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s"
                                   " nchar=%s ;\n\tformat datatype=%s "
                                   "interleave=no gap=%s missing=%s ;\n\t"
                                   "matrix\n" % (
                                    len(alignment),
                                    self.locus_length,
                                    self.sequence_code[0],
                                    gap, self.sequence_code[1]))

                for key, f in alignment.items():

                    with open(f) as fh:
                        seq = "".join(fh.readlines())

                    out_file.write("%s %s\n" % (key[:cut_space_nex].ljust(
                        seq_space_nex), seq))
                out_file.write(";\n\tend;")

            if use_charset:
                # Writing partitions, if any
                if not self.partitions.is_single():
                    out_file.write("\nbegin mrbayes;\n")
                    # Full gene partitions
                    for name, lrange in self.partitions:
                        # If there are codon partitions, write those
                        if lrange[1]:
                            for i in lrange[1]:
                                out_file.write("\tcharset %s_%s = %s-%s\\3;\n" %
                                       (name, i + 1, i + 1, lrange[0][1] + 1))
                        else:
                            out_file.write("\tcharset %s = %s-%s;\n" %
                                       (name, lrange[0][0] + 1,
                                        lrange[0][1] + 1))

                    out_file.write("\tpartition part = %s: %s;\n\tset "
                                   "partition=part;\nend;\n" %
                                   (len(self.partitions.partitions),
                                    ", ".join([name for name in
                                    self.partitions.get_partition_names()])))

            # In case outgroup taxa are specified
            if outgroup_list is not None:

                # This assures that only the outgroups present in the current
                #  file are written
                compliant_outgroups = [taxon for taxon in outgroup_list
                                       if taxon in self.iter_sequences()]
                if compliant_outgroups is not []:
                    out_file.write("\nbegin mrbayes;\n\toutgroup %s\nend;\n" %
                                   (" ".join(compliant_outgroups)))

            out_file.close()

        # Writes file in fasta format
        if "fasta" in output_format:
            out_file = open(output_file + ".fas", "w")

            # If LD HAT sub format has been specificed, write the first line
            # containing the number of sequences, sites and genotype phase
            if ld_hat:
                out_file.write("{} {} {}\n".format(len(self.alignment),
                                                 self.locus_length,
                                                 "2"))

            for key, f in alignment.items():

                with open(f) as fh:
                    seq = "".join(fh.readlines())

                if ld_hat:
                    # Truncate sequence name to 30 characters
                    out_file.write(">%s\n" % (key[:30]))
                    # Limit each sequence line to 2000 characters
                    if len(seq) > 2000:
                        for i in range(0, len(seq), 2000):
                            out_file.write("%s\n" % (seq[i:i + 2000].upper()))
                else:
                    out_file.write(">%s\n%s\n" % (key, seq.upper()))

            out_file.close()


class AlignmentList(Base):
    """
    At the most basic instance, this class contains a list of Alignment
    objects upon which several methods can be applied. It only requires either
    a list of alignment files or. It inherits methods from Base and
    Alignment classes for the write_to_file methods.
    """

    def __init__(self, alignment_list, dest=None, shared_namespace=None):
        """
        :param alignment_list: List of Alignment objects
        :param dest: String. Path to temporary directory that will store
        the sequence data of each alignment object
        :param shared_namespace: Namespace object, used to share information
        between subprocesses
        """

        self.log_progression = Progression()

        """
        Stores the "active" Alignment objects for the current AlignmentList.
        Keys will be the Alignment.path for quick lookup of Alignment object
        values
        """
        self.alignments = OrderedDict()

        """
        Stores the "inactive" or "shelved" Alignment objects. All AlignmentList
        methods will operate only on the alignments attribute, unless explicitly
        stated otherwise. Key-value is the same as self.alignments
        """
        self.shelve_alignments = OrderedDict()

        """
        Attribute list that stores the Alignment.name attribute of badly
        formatted alignments
        """
        self.bad_alignments = []

        """
        Attribute list that stores duplicate Alignment.name.
        """
        self.duplicate_alignments = []

        """
        Attribute list that stores sequence sets of unequal lenght
        """
        self.non_alignments = []

        """
        List with the name of the taxa included in the AlignmentList object
        """
        self.taxa_names = []

        """
        Lists the Alignment.name attributes of the current AlignmentList object
        """
        self.filename_list = []
        self.path_list = []

        """
        Tuple with the AlignmentList sequence code. Either ("DNA", "n") or
        ("Protein", "x")
        """
        self.sequence_code = None
        self.gap_symbol = "-"

        # Set partitions object
        self.partitions = Partitions()

        c = 0
        # if type(alignment_list[0]) is str:
        if alignment_list:
            for alignment in alignment_list:

                if shared_namespace:
                    shared_namespace.progress = c
                    shared_namespace.m = alignment.split(sep)[-1]
                    c += 1

                alignment_object = Alignment(alignment, dest=dest)

                # Check for badly formatted alignments
                if isinstance(alignment_object.e, InputError):
                    self.bad_alignments.append(alignment_object.path)
                elif isinstance(alignment_object.e,
                                AlignmentUnequalLength):
                    self.non_alignments.append(alignment_object.path)
                    print("Warning: Sequences of unequal length detected"
                          " in file {}".format(alignment_object.name))

                # Check for duplicate alignments
                elif alignment_object.path in [x.path for x in
                                             self.alignments.values()]:
                    self.duplicate_alignments.append(alignment_object.name)
                else:
                    # Get seq code
                    if not self.sequence_code:
                        self.sequence_code = alignment_object.sequence_code

                    self.alignments[alignment_object.name] = alignment_object
                    self.set_partition(alignment_object)
                    self.filename_list.append(alignment_object.name)

        self.taxa_names = self._get_taxa_list()

    def __iter__(self):
        """
        Iterate over Alignment objects
        """
        return iter(self.alignments.values())

    def clear_alignments(self):
        """
        Clears the current AlignmentList object
        :return:
        """

        # Clear temporary sequence files
        for aln in self.alignments.values():
            aln.remove_alignment()

        self.alignments = {}
        self.shelve_alignments = {}
        self.bad_alignments = []
        self.duplicate_alignments = []
        self.partitions = Partitions()
        self.filename_list = []
        self.taxa_names = []
        self.non_alignments = []
        self.sequence_code = None
        self.path_list = []

    def update_active_alignments(self, aln_list):
        """
        Updates the self.alignments and self.shelve_alignments attributes.
        The Alignment.name's provided by the argument will populate
        self.alignments and the remaining will be
        """

        for aln_name in self.filename_list:
            if aln_name in aln_list:
                if aln_name in self.shelve_alignments:
                    self.alignments[aln_name] = self.shelve_alignments[aln_name]
                    del self.shelve_alignments[aln_name]
            else:
                if aln_name in self.alignments:
                    self.shelve_alignments[aln_name] = self.alignments[aln_name]
                    del self.alignments[aln_name]

    def update_active_alignment(self, aln_name, direction):
        """
        Same as update_active_alignments but for a single aln_name, so that
        the whole list does not need to be iterated
        :param aln_name: string, name of the alignment to move
        :param direction: string, can be either 'shelve' or 'active'
        """

        if direction == "shelve":
            self.shelve_alignments[aln_name] = self.alignments[aln_name]
            del self.alignments[aln_name]

        else:
            self.alignments[aln_name] = self.shelve_alignments[aln_name]
            del self.shelve_alignments[aln_name]

    def format_list(self):
        """
        :return: List with the unique sequence types of the Alignment objects
        """

        return list(set([x.sequence_code[0] for x in
                         self.alignments.values() if x]))

    def _get_taxa_list(self):
        """
        Gets the full taxa list of all alignment objects
        :return full_taxa. List of taxa names in the AlignmentList
        """

        full_taxa = []

        for alignment in self.alignments.values():
            diff = set(alignment.iter_taxa()) - set(full_taxa)
            if diff != set():
                full_taxa.extend(diff)

        return full_taxa

    def _get_filename_list(self):
        """
        Returns a list with the input file names
        """
        return (alignment.name for alignment in self.alignments.values())

    def start_action_alignment(self):
        """
        Wrapper of the start_action_alignment of Alignment object
        """

        for aln_obj in self.alignments.values():
            aln_obj.start_action_alignment()

    def stop_action_alignment(self):
        """
        Wrapper of the stop_action_alignment of the Alignment object
        """

        for aln_obj in self.alignments.values():
            aln_obj.stop_action_alignment()

    def set_partition(self, alignment_obj):
        """
        Updates the partition object with the provided alignment_obj
        :param alignment_obj: Alignment object
        :return:
        """

        # Update partitions object
        if not alignment_obj.partitions.is_single():
            for k, v in alignment_obj.partitions:
                self.partitions.add_partition(k, locus_range=v[0], codon=v[1],
                            use_counter=False, file_name=alignment_obj.path,
                            model_cls=alignment_obj.partitions.models[k])
        else:
            self.partitions.add_partition(alignment_obj.name,
                                use_counter=True,
                                file_name=alignment_obj.path,
                                length=alignment_obj.locus_length,
                                model_cls=alignment_obj.partitions.models[
                                    alignment_obj.name])

    def add_alignments(self, alignment_obj_list, ignore_paths=False):
        """
        Adds a new Alignment object
        :param alignment_obj_list: list with Alignment objects
        """

        for alignment_obj in alignment_obj_list:

            # Check for badly formatted alignments
            if isinstance(alignment_obj.e, InputError):
                self.bad_alignments.append(alignment_obj.path)
            elif isinstance(alignment_obj.e,
                            AlignmentUnequalLength):
                self.non_alignments.append(alignment_obj.path)

            # if not ignore_paths:
            if not ignore_paths:
                if alignment_obj.path in [x.path for x in
                                          self.alignments.values()]:
                    self.duplicate_alignments.append(alignment_obj.name)
                else:
                    # Get seq code
                    if not self.sequence_code:
                        self.sequence_code = alignment_obj.sequence_code

                    self.alignments[alignment_obj.name] = alignment_obj
                    self.set_partition(alignment_obj)
                    self.filename_list.append(alignment_obj.name)
            else:
                # Get seq code
                if not self.sequence_code:
                    self.sequence_code = alignment_obj.sequence_code

                self.alignments[alignment_obj.name] = alignment_obj
                self.set_partition(alignment_obj)
                self.filename_list.append(alignment_obj.name)

        self.taxa_names = self._get_taxa_list()

    def add_alignment_files(self, file_name_list, dest=None,
                            shared_namespace=None):
        """
        Adds a new alignment based on a file name
        :param file_name_list: list, with the path to the alignment files
        :param dest: string. Path to the temporary directory that will store
        the sequence data of the Alignment object
        """

        # Check for duplicates
        for i in set(self.path_list).intersection(set(file_name_list)):
            self.duplicate_alignments.append(i)
            file_name_list.remove(i)

        if shared_namespace:
            shared_namespace.progress = 0

        njobs = len(file_name_list) if len(file_name_list) <= \
            multiprocessing.cpu_count() else multiprocessing.cpu_count()

        if shared_namespace:
            jobs = [[x.tolist(), dest, y, shared_namespace] for y, x in
                enumerate(np.array_split(np.array(file_name_list), njobs))]
        else:
            jobs = [[x.tolist(), dest, y, None] for y, x in enumerate(
                np.array_split(np.array(file_name_list), njobs))]
        # Execute alignment reading in parallel
        multiprocessing.Pool(njobs).map(
            read_alns, jobs)

        # Read the pickle files with the saved Alignment objects
        for i in range(njobs):
            fh = open(join(dest, "Worker{}.pc".format(i)), "rb")
            while 1:
                try:
                    aln = pickle.load(fh)

                    if isinstance(aln.e, InputError):
                        self.bad_alignments.append(aln.path)
                    elif isinstance(aln.e, AlignmentUnequalLength):
                        self.non_alignments.append(aln.path)

                    else:
                        # Get seq code
                        if not self.sequence_code:
                            self.sequence_code = aln.sequence_code

                        self.alignments[aln.name] = aln
                        self.set_partition(aln)
                        self.filename_list.append(aln.name)
                        self.path_list.append(aln.path)

                except EOFError:
                    fh.close()
                    os.remove(join(dest, "Worker{}.pc".format(i)))
                    break

        self.taxa_names = self._get_taxa_list()

        if shared_namespace:
            shared_namespace.m = "Updating App structures"

    def retrieve_alignment(self, name):
        """
        :param name: string. Name of the input alignment
        :return: Returns an Alignment object with a given name attribute
        """

        if name in self.alignments:
            return self.alignments[name]
        elif name in self.shelve_alignments:
            return self.shelve_alignments[name]
        else:
            return None

    def iter_alignment_dic(self):
        """
        :return: List of the dictionary alignments
        """

        return iter(alignment.alignment for alignment in
                    self.alignments.values())

    def write_taxa_to_file(self):
        """
        Compiles the taxa names of all alignments and writes them in a single
        column .csv file
        """

        output_handle = open("Taxa_list.csv", "w")

        for taxon in self.taxa_names:
            output_handle.write(taxon + "\n")

        output_handle.close()

    def concatenate(self, alignment_name=None, dest=None, remove_temp=True):
        """
        Concatenates multiple sequence alignments creating a single alignment
        object and the auxiliary Partitions object defining the partitions
        of the concatenated alignment
        :param alignment_name: string. Optional. Name of the new concatenated
        alignment object. This should be used when collapsing the alignment
        afterwards.
        :param dest: Path to temporary directory where sequence data will be
        stored
        :param remove_temp: boolean. If True, it will remove active temporary
        sequence files.
        :return concatenated_alignment: Alignment object
        """

        # If dest is provided, this will ensure that the directory exists.
        # Take caution to provide a unique path, so that other temp
        # files are not overwritten. If dest is not provided, temp files will
        # be created on the current working directory
        if dest and alignment_name:
            if not os.path.exists(join(dest, alignment_name)):
                os.makedirs(join(dest, alignment_name))
        else:
            dest = "./"

        # Variable that will store the lenght of the concatenated alignment
        # and provided it when initializing the Alignment object
        locus_length = None

        # Initializing alignment dict to store the alignment information
        if alignment_name:
            concatenation = OrderedDict(
                [(key, join(dest, alignment_name, key + ".seq")) for key in
                 self.taxa_names])
        else:
            concatenation = OrderedDict(
                [(key, join(dest, key + ".seq")) for key in self.taxa_names])

        # Concatenation is performed for each taxon at a time
        for taxon in self.taxa_names:

            # This variable will remain False unless some data is provided for
            # the current taxa. If there is only missing data, this taxon
            # will be removed from the concatenation dict
            has_data = False

            with open(concatenation[taxon], "w") as fh:

                for aln in self.alignments.values():

                    # Setting missing data symbol
                    missing = aln.sequence_code[1]

                    if taxon in aln.alignment:
                        fh.write(aln.get_sequence(taxon))
                        has_data = True
                    else:
                        fh.write(missing * aln.locus_length)

            # Remove taxa that do not contain any data
            if not has_data:
                os.remove(concatenation[taxon])
                del concatenation[taxon]
            # Retrieve locus length
            else:
                if not locus_length:
                    with open(concatenation[taxon]) as fh:
                        locus_length = len(fh.read().strip())

        # Removes partitions that are currently in the shelve
        for aln_obj in self.shelve_alignments.values():
            self.partitions.remove_partition(file_name=aln_obj.path)

        # Remove previous temp sequence files
        if remove_temp:
            for aln in self.alignments.values():
                aln.remove_alignment(remove_active=True)

        # Create the concatenated file in an Alignment object
        concatenated_alignment = Alignment(concatenation,
                                           dest=dest,
                                           partitions=self.partitions,
                                           alignment_name=alignment_name,
                                           locus_length=locus_length)

        return concatenated_alignment

    def filter_min_taxa(self, min_taxa):
        """
        Filters Alignment objects based on a minimum taxa representation
        threshold. Alignments with less that the specified minimum taxa
        percentage will be moved to the filtered_alignments attribute.

        NOTE: Since this filtering is meant to be performed when executing
        the process operations it will permanently change the AlignmentList
        object, which means both self.alignments and self.partitions. Not doing
        so and removing/adding the partitions would create a great deal of
        conflicts that can be easily avoided by simply copying the
        AlignmentList object and modifying this object for the process execution

        :param min_taxa: integer, percentage of minimum taxa below which
        alignments are moved to the filtered_alignments attribute
        """

        for k, alignment_obj in list(self.alignments.items()):
            if len(alignment_obj.alignment) < \
                    (min_taxa / 100) * len(self.taxa_names):
                del self.alignments[k]
                self.partitions.remove_partition(file_name=alignment_obj.path)

    def filter_by_taxa(self, filter_mode, taxa_list):
        """
        Filters the alignments attribute by taxa list. The filtering may be done
        to exclude or include a particular set of taxa
        :param filter_mode: string, determines the filtering mode. Can be either
        'Contain' or 'Exclude'
        :param taxa_list: list, contains the list of taxa to be used for
        filtering
        """

        for k, alignment_obj in list(self.alignments.items()):

            # Filter alignments that do not contain at least all taxa in
            # taxa_list
            if filter_mode == "Contain":
                if set(taxa_list) - set(list(alignment_obj.alignment)) != set():
                    del self.alignments[k]

            # Filter alignments that contain the taxa in taxa list
            if filter_mode == "Exclude":
                if any((x for x in taxa_list
                        if x in list(alignment_obj.alignment))):
                    del self.alignments[k]

        # If the resulting alignment is empty, raise an Exception
        if self.alignments == {}:
            raise EmptyAlignment("Alignment is empty after taxa filter")

    def filter_codon_positions(self, position_list):
        """
        Filter codon positions from DNA alignments.
        :param position_list: list containing a boolean value for each codon
        position. Ex. [True, True, True] will save all positions while
        [True, True, False] will exclude the third codon position
        """

        def index(length, pos):
            """
            index generator
            """
            for _ in range(0, length, 3):
                for j in pos:
                    if j:
                        yield 1
                    else:
                        yield 0

        # Reset partitions
        self.partitions = Partitions()

        for alignment_obj in list(self.alignments.values()):

            for taxon, seq in alignment_obj:
                filtered_seq = "".join(list(itertools.compress(seq,
                                            index(alignment_obj.locus_length,
                                                  position_list))))
                alignment_obj.alignment[taxon] = filtered_seq

            alignment_obj.locus_length = len(filtered_seq)

            self.set_partition(alignment_obj)

    def filter_missing_data(self, gap_threshold, missing_threshold):
        """
        Wrapper of the filter_missing_data method of the Alignment object.
        See the method's documentation.
        :param gap_threshold: integer, percentage of gap symbols below which
        the alignment column should be filtered
        :param missing_threshold: integer, percentage of missing data (gaps +
        true missing data) below which the alignment column should be fitered
        """

        for alignment_obj in list(self.alignments.values()):

            alignment_obj.filter_missing_data(gap_threshold=gap_threshold,
                                        missing_threshold=missing_threshold)

    def remove_taxa(self, taxa_list, mode="remove"):
        """
        Wrapper of the remove_taxa method of the Alignment object for
        multiple alignments. It current supports two modes:

            ..:remove: removes specified taxa
            ..:inverse: removes all but the specified taxa
        """

        for alignment_obj in list(self.alignments.values()):
            alignment_obj.remove_taxa(taxa_list, mode=mode)

        # Updates taxa names
        for tx in taxa_list:
            try:
                self.taxa_names.remove(tx)
            except ValueError:
                # TODO: log a warning
                pass

    def change_taxon_name(self, old_name, new_name):
        """
        Changes the name of a taxon
        """

        for alignment_obj in list(self.alignments.values()):
            alignment_obj.change_taxon_name(old_name, new_name)

        # update taxa names
        self.taxa_names = [new_name if x == old_name else x
                           for x in self.taxa_names]

    def remove_file(self, filename_list):
        """
        Removes alignment objects based on their name attribute
        :param filename_list: list with the names of the alignment objects to
        be removed
        """

        for nm_path in filename_list:
            nm = nm_path.split(sep)[-1]
            if nm in self.alignments:
                self.alignments[nm].remove_alignment()
                del self.alignments[nm]
            elif nm in self.shelve_alignments:
                self.alignments[nm].remove_alignment()
                del self.shelve_alignments[nm]
            self.partitions.remove_partition(file_name=nm_path)

        # Updates taxa names
        self.taxa_names = self._get_taxa_list()

    def shelve_file(self, filename_list):
        """
        Instead of completely removing the Alignment object, these are moved
        to the shelve_alignments list.
        :param filename_list: list with the names of the alignment objects to
        be removed
        """

        for nm in filename_list:
            nm_wext = basename(nm)
            nm = basename(nm).split(".")[0]
            if nm in self.alignments:
                self.shelve_alignments[nm] = self.alignments[nm]
                del self.alignments[nm]
            self.partitions.remove_partition(file_name=nm_wext)

        # Updates taxa names
        self.taxa_names = self._get_taxa_list()

    def select_by_taxa(self, taxa_list, mode="strict"):
        """
        This method is used to selected gene alignments according to a list
        of taxa.

        :param taxa_list. List of taxa names

        :param mode. String. Modes can be the following:
            ..:strict: The taxa of the alignment must be exactly the same as the
        specified taxa.
            ..:inclusive: The taxa of the alignment must contain all specified
        taxa.
            ..:relaxed: At least on of the specified taxa must be in the taxa of
        the alignment.
        """

        selected_alignments = []

        # taxa_list may be a file name (string) or a list containing the name
        # of the taxa. If taxa_list is a file name this code will parse the
        # csv file and return a list of the taxa. Otherwise, the taxa_list
        # variable remains the same.
        try:
            file_handle = open("".join(taxa_list))
            taxa_list = self.read_basic_csv(file_handle)
        except EnvironmentError:
            pass

        for alignment_obj in self.alignments.values():

            alignment_taxa = list(alignment_obj.alignment)

            # Selected only the alignments with the exact same taxa
            if mode == "strict":
                if set(taxa_list) == set(alignment_taxa):
                    selected_alignments.append(alignment_obj)

            # Selected alignments that include the specified taxa
            if mode == "inclusive":
                if set(taxa_list) - set(alignment_taxa) == set():
                    selected_alignments.append(alignment_obj)

            if mode == "relaxed":
                for taxon in taxa_list:
                    if taxon in alignment_taxa:
                        selected_alignments.append(alignment_obj)
                        continue

        return selected_alignments

    def code_gaps(self):
        """
        Wrapper for the code_gaps method of the Alignment object.
        """

        for alignment_obj in self.alignments.values():
            alignment_obj.code_gaps()

    def collapse(self, write_haplotypes=True, haplotypes_file="",
                 haplotype_name="Hap", dest="./"):
        """
        Wrapper for the collapse method of the Alignment object. If
        write_haplotypes is True, the haplotypes file name will be based on the
        individual input file

        :param write_haplotypes: Boolean, if True, a haplotype list
        mapping the haplotype names file will be created for each individual
        input alignment.

        :param haplotype_name: String, Custom name of the haplotypes

        :param dest: string. Path to write the .haplotypes file
        """

        for alignment_obj in self.alignments.values():
            if write_haplotypes:
                # Set name for haplotypes file
                output_file = alignment_obj.name.split(".")[0] + haplotypes_file
                alignment_obj.collapse(haplotypes_file=output_file,
                                       haplotype_name=haplotype_name, dest=dest)
            else:
                alignment_obj.collapse(write_haplotypes=False,
                                       haplotype_name=haplotype_name, dest=dest)

    def consensus(self, consensus_type, single_file=False):
        """
        Wrapper for the consensus method of the Alignment object.
        :param consensus_type: string. How variation is handled when creating
        the consensus. See the Alignment method for more information
        """

        self.taxa_names = ["consensus"]

        if single_file:
            single_consensus = OrderedDict()

        for alignment_obj in self.alignments.values():
            alignment_obj.consensus(consensus_type)

            if single_file:
                single_consensus[alignment_obj.name] = \
                    alignment_obj.alignment["consensus"]

        if single_file:
            consensus_aln = Alignment(single_consensus)
            return consensus_aln

    def reverse_concatenate(self, dest=None):
        """
        Internal function to reverse concatenate an alignment according to
        defined partitions in a Partitions object

        This will only work if alignment_object_list has one alignment, as it
        is intended to be a wrapper of sorts for the Alignment object method

        :param partition_obj: Partitions object, containing the partitions of
        the input file
        :return: AlignmentList object with individual alignments
        """

        concatenated_aln = self.concatenate(alignment_name="concatenated",
                                            dest=dest)

        reverted_alns = concatenated_aln.reverse_concatenate()

        return reverted_alns

    def write_to_file(self, output_format, output_suffix="", interleave=False,
                      outgroup_list=None, partition_file=True, output_dir=None,
                      use_charset=True, phy_truncate_names=False, ld_hat=None,
                      ima2_params=None):
        """
        Wrapper of the write_to_file method of the Alignment object for multiple
        alignments.

        :param output_format: string, format of the output file
        :param output_suffix: string, optional suffix that is added at the end
        of the original file name
        :param interleave: boolean, Whether the output alignment will be in
        leave (False) or interleave (True) format. Not all output formats
        respect this option.
        :param outgroup_list: list, containing the taxa names of the outgroup.
        (Nexus output format only)
        :param partition_file: boolean, If true, the auxiliary partitions file
        will be writen.
        :param output_dir: string, if provided, the output file will be written
        on the specified path
        :param use_charset: boolean, if true, partitions from the Partitions
        object will be written in the nexus output format
        :param phy_truncate_names: Boolean. Whether names in phylip output
        format should be truncated to 10 characters or not.
        """

        for alignment_obj in self.alignments.values():

            if alignment_obj.input_format in output_format:
                output_file_name = alignment_obj.name.split(".")[0] \
                                   + output_suffix + "_conv"
            else:
                output_file_name = alignment_obj.name.split(".")[0] + \
                                   output_suffix

            alignment_obj.write_to_file(output_format,
                                        output_file=output_file_name,
                                        interleave=interleave,
                                        outgroup_list=outgroup_list,
                                        partition_file=partition_file,
                                        output_dir=output_dir,
                                        use_charset=use_charset,
                                        phy_truncate_names=phy_truncate_names,
                                        ld_hat=ld_hat,
                                        ima2_params=ima2_params)

    # Stats methods
    def gene_occupancy(self):
        """
        Creates data for an interpolation plot to visualize the amount of
        missing genes in a phylogenomics data set
        """

        data = []

        for alignment in self.alignments.values():
            data.append([1 if x in alignment.alignment.keys() else 0
                         for x in self.taxa_names])

        data = np.transpose(data)

        return {"data": data}

    def missing_data_distribution(self):
        """
        Creates data for an overall distribution of missing data
        """

        legend = ["Gaps", "Missing data", "Data"]

        data_storage = OrderedDict((x, []) for x in legend)

        for aln in self.alignments.values():
            gaps_g, missing_g, data_g = [], [], []
            for tx in self.taxa_names:
                if tx in aln.alignment:
                    seq = aln.get_sequence(tx)
                    # Get gaps
                    gaps = float(seq.count("-"))
                    gaps_g.append(gaps / float(aln.locus_length))
                    # Get missing data
                    missing = float(seq.count(aln.sequence_code[1]))
                    missing_g.append(missing / float(aln.locus_length))
                    # Get actual data
                    actual_data = (float(aln.locus_length) - gaps - missing) / \
                        float(aln.locus_length)
                    data_g.append(actual_data)
                else:
                    gaps_g.append(0)
                    missing_g.append(1)
                    data_g.append(0)

            data_storage["Gaps"].append(np.mean(gaps_g))
            data_storage["Missing data"].append(np.mean(missing_g))
            data_storage["Data"].append(np.mean(data_g))

        data = [x for x in data_storage.values()]

        return {"data": data,
                "legend": legend,
                "ax_names": ["Proportion", "Number of genes"],
                "table_header": ["Bin"] + legend}

    def missing_data_per_species(self):
        """
        Creates data for a distribution of missing data per species
        """

        # Data for a stacked bar plot. First element for gaps, second for
        # missing, third for actual data
        data_storage = OrderedDict((taxon, [0, 0, 0]) for taxon in
                                   self.taxa_names)
        total_len = 0

        legend = ["Gaps", "Missing", "Data"]

        for aln in self.alignments.values():
            total_len += aln.locus_length
            for key in data_storage:
                if key in aln.alignment:
                    # Get gaps
                    seq = aln.get_sequence(key)
                    gaps = seq.count("-")
                    data_storage[key][0] += gaps
                    # Get missing
                    missing = seq.count(aln.sequence_code[1])
                    data_storage[key][1] += missing
                    # Get actual data
                    actual_data = aln.locus_length - gaps - missing
                    data_storage[key][2] += actual_data
                else:
                    data_storage[key][1] += aln.locus_length

        data_storage = OrderedDict(sorted(data_storage.items(),
                                          key=lambda x: x[1][1] + x[1][0],
                                          reverse=True))

        data = np.array([[float(x[0]) for x in
                          data_storage.values()],
                         [float(x[1]) for x in
                          data_storage.values()],
                         [float(x[2]) for x in
                          data_storage.values()]])

        return {"data": data,
                "labels": list(data_storage.keys()),
                "legend": legend,
                "table_header": ["Taxon", "Gaps", "%", "Missing", "%", "Data",
                                 "%"],
                "normalize": True,
                "normalize_factor": total_len}

    def missing_genes_per_species(self):
        """
        Creates data for the distribution of missing genes per species
        :return: dictionary with arguments for plotting functions
        """

        data_storage = OrderedDict((taxon, 0) for taxon in self.taxa_names)

        for aln in self.alignments.values():
            for key in data_storage:
                if key not in aln.alignment:
                    data_storage[key] += 1

        # Sort data in descending order of missing genes
        data_storage = OrderedDict(sorted(data_storage.items(), reverse=True,
                                          key=lambda t: t[1]))

        return {"data": [list(data_storage.values())],
                "labels": list(data_storage.keys()),
                "title": "Distribution of missing genes per species",
                "ax_names": [None, "Frequency"],
                "table_header": ["Taxon", "Missing genes"]
                }

    def missing_genes_average(self):
        """
        Creates histogram data for average mssing genes
        """

        data = []

        for aln in self.alignments.values():
            data.append(len(set(self.taxa_names) - set(aln.alignment.keys())))

        return {"data": data,
                "title": "Distribution of missing genes",
                "ax_names": ["Number of missing genes", "Frequency"],
                "table_header": ["Number of missing genes", "Frequency"]}

    def average_seqsize_per_species(self):
        """
        Creates data for the average sequence size for each taxa
        :return: dictionary with arguments for plotting functions
        """

        data_storage = OrderedDict((taxon, []) for taxon in self.taxa_names)

        for aln in self.alignments.values():
            for sp, seq in aln:
                data_storage[sp].append(len(seq.replace("-", "").
                                        replace(aln.sequence_code[1], "")))

        # Adapt y-axis label according to sequence code
        seq_code = aln.sequence_code[0]
        ax_ylabel = "Size (bp)" if seq_code == "DNA" else "Size (residues)"

        data_storage = OrderedDict(sorted(data_storage.items(), reverse=True,
                                   key=lambda t: np.mean(t[1])))

        return {"data": list(data_storage.values()),
                "labels": list(data_storage.keys()),
                "title": "Sequence size distribution per species",
                "ax_names": [None, ax_ylabel]}

    def average_seqsize(self):
        """
        Creates data for the average sequence size for the entire data set
        :return:
        """

        data_storage = []

        for aln in self.alignments.values():
            data_storage.append(aln.locus_length)

        # Adapt y-axis label according to sequence code
        seq_code = aln.sequence_code[0]
        ax_xlabel = "Size (bp)" if seq_code == "DNA" else "Size (residues)"

        return {"data": data_storage,
                "title": "Average sequence size distribution",
                "ax_names": [ax_xlabel, "Frequency"],
                "table_header": [ax_xlabel, "Frequency"]}

    def characters_proportion(self):
        """
        Creates data for the proportion of nucleotides/residues for the data set
        """

        data_storage = Counter()

        for aln in self.alignments.values():
            for seq in aln.sequences():
                data_storage += Counter(seq.replace("-", "").
                                        replace(self.sequence_code[1], ""))

        # Determine total number of characters
        chars = float(sum(data_storage.values()))

        # Valid characters list
        valid_chars = dna_chars if self.sequence_code[0] == "DNA" else \
            list(aminoacid_table.keys())

        data, xlabels = zip(*[(float(x) / chars, y.upper()) for y, x in
                              data_storage.items() if y in valid_chars])

        title = "Nucleotide proportions" if self.sequence_code[0] == "DNA" \
            else "Amino acid proportions"
        ax_xlabel = "Nucleotide" if self.sequence_code[0] == "DNA" \
            else "Amino acid"

        return {"data": [data],
                "labels": xlabels,
                "title": title,
                "ax_names": [ax_xlabel, "Proportion"],
                "table_header": [ax_xlabel, "Proportion"]}

    def characters_proportion_per_species(self):
        """
        Creates data for the proportion of nucleotides/residures per species
        """

        data_storage = OrderedDict((x, Counter()) for x in self.taxa_names)

        for aln in self.alignments.values():
            for sp, seq in aln:
                data_storage[sp] += Counter(seq.replace("-", "").
                                            replace(self.sequence_code[1], ""))

        legend = dna_chars if self.sequence_code[0] == "DNA" else \
            list(aminoacid_table.keys())

        data = [[] for _ in legend]

        for p, char in enumerate(legend):
            for c in data_storage.values():
                chars = float(sum([x for y, x in c.items() if y in legend]))
                data[p].append(float(c[char]) / chars)

        data = np.array(data)

        ax_ylabel = "Nucleotide" if self.sequence_code[0] == "DNA" \
            else "Amino acid"

        return {"data": data,
                "labels": list(data_storage.keys()),
                "legend": legend,
                "ax_names": ["Taxa", ax_ylabel],
                "table_header": ["Taxon"] + legend}

    def _get_similarity(self, seq1, seq2):
        """
        Gets the similarity between two sequences in proportion. The
        proportion is calculated here so that sites with only missing
        data/gaps are not accounted in the total length.
        :param seq1: string
        :param seq2: string
        """

        similarity = 0.0
        total_len = 0.0

        missing = [self.sequence_code[1], self.gap_symbol]

        for c1, c2 in zip(*[seq1, seq2]):
            # Ignore comparisons with ONLY missing data / gaps
            if c1 in missing and c2 in missing:
                continue
            elif c1 == c2:
                similarity += 1.0
                total_len += 1.0
            else:
                total_len += 1.0

        return similarity / total_len

    def _get_differences(self, seq1, seq2):
        """
        Returns the number of differences between two sequences
        """

        s = 0

        for c1, c2 in zip(*[seq1, seq2]):
            # Ignore comparisons with gaps or missing data
            if self.sequence_code[1] in [c1, c2] or self.gap_symbol in [c1, c2]:
                continue
            if c1 != c2:
                s += 1

        return s

    def _get_differences_prop(self, seq1, seq2):
        """
        Returns the proportion of differences between two sequences,
        accounting for sites with only missing data / gaps
        :param seq1: string, sequence 1
        :param seq2: string sequence 2
        :return: s: float, proportion of segregating sites
        """

        s = 0.0
        total_len = 0.0

        for c1, c2 in zip(*[seq1, seq2]):
            # Ignore comparisons with ONLY missing data / gaps
            if c1 + c2.replace(self.sequence_code[1], "").replace(
                    self.gap_symbol, "") == "":
                continue
            if c1 != c2:
                s += 1.0
                total_len += 1.0
            else:
                total_len += 1.0

        return s / total_len

    def _get_informative_sites(self, aln):
        """
        Determines the number of informative sites in an Alignment object
        :param aln: Alignment object
        :return:
        """

        informative_sites = 0

        for column in zip(*aln.sequences()):

            column = Counter([x for x in column if x != aln.sequence_code[1] and
                              x != self.gap_symbol])

            if len(column) > 1:
                # Remove most common and check the length of the remaining
                del column[column.most_common()[0][0]]
                if sum(column.values()) > 1:
                    informative_sites += 1

        return informative_sites

    def sequence_similarity(self):
        """
        Creates average sequence similarity data
        """

        data = []

        for aln in self.alignments.values():

            aln_similarities = []

            for seq1, seq2 in itertools.combinations(aln.sequences(), 2):

                x = self._get_similarity(seq1, seq2)

                if x:
                    aln_similarities.append(x)

            if aln_similarities:
                data.append(np.mean(aln_similarities) * 100)

        return {"data": data,
                "ax_names": ["Similarity (%)", "Frequency"]}

    def sequence_similarity_per_species(self):
        """
        Creates data for a triangular matrix of sequence similarity for pairs
        of taxa
        """

        # Create matrix for parwise comparisons
        data = [np.empty((len(self.taxa_names), 0)).tolist() for _ in
                range(len(self.taxa_names))]

        taxa_pos = OrderedDict((x, y) for y, x in enumerate(self.taxa_names))

        for aln in self.alignments.values():

            for tx1, tx2 in itertools.combinations(taxa_pos.keys(), 2):

                try:
                    seq1, seq2 = aln.get_sequence(tx1), aln.get_sequence(tx2)
                except KeyError:
                    continue

                similarity = self._get_similarity(seq1, seq2)
                data[taxa_pos[tx1]][taxa_pos[tx2]].append(similarity)

        data = np.array([[np.mean(y) if y else 0. for y in x] for x in data])
        mask = np.tri(data.shape[0], k=0)
        data = np.ma.array(data, mask=mask)

        return {"data": data,
                "labels": list(taxa_pos)}

    def sequence_similarity_gene(self, gene_name, window_size):

        aln_obj = self.alignments[gene_name]

        data = []

        for i in range(0, aln_obj.locus_length, window_size):

            window_similarities = []

            seqs = np.array([[y for y in x[i:i + window_size]] for x in
                             aln_obj.sequences()])

            for seq1, seq2 in itertools.combinations(seqs, 2):
                window_similarities.append(self._get_similarity(seq1, seq2)
                                           * 100)

            if window_similarities:
                data.append(np.mean(window_similarities))

        return {"data": data,
                "window_size": window_size,
                "ax_names": ["Sequence (bp)", "Similarity (%)"],
                "table_header": ["Sequence (bp)", "Similarity (%)"]}

    def sequence_segregation(self, proportions=False):
        """
        Generates data for distribution of segregating sites

        :param proportions: Boolean. If True, use proportions instead of
        absolute values
        """

        data = []

        for aln in self.alignments.values():

            segregating_sites = 0

            for column in zip(*aln.sequences()):

                # Remove gaps and missing characters
                column = set([x for x in column if x != aln.sequence_code[1]
                              and x != self.gap_symbol])

                if len(column) > 1:
                    segregating_sites += 1

            if proportions:
                data.append(
                    float(segregating_sites / float(aln.locus_length)))
                ax_names = ["Segregating sites", "Percentage"]
                real_bin = False
            else:
                data.append(segregating_sites)
                ax_names = ["Segregating sites", "Frequency"]
                real_bin = True

        return {"data": data,
                "ax_names": ax_names,
                "title": "Distribution of segregating sites",
                "table_header": ax_names,
                "real_bin_num": real_bin}

    def sequence_segregation_per_species(self):
        """
        Creates a data for a triangular matrix of sequence segregation for
        pairs of taxa
        """

        # Create matrix for parwise comparisons
        data = [np.empty((len(self.taxa_names), 0)).tolist() for _ in
                range(len(self.taxa_names))]

        taxa_pos = OrderedDict((x, y) for y, x in enumerate(self.taxa_names))

        for aln in self.alignments.values():

            for tx1, tx2 in itertools.combinations(taxa_pos.keys(), 2):

                try:
                    seq1, seq2 = aln.get_sequence(tx1), aln.get_sequence(tx2)
                except KeyError:
                    continue

                aln_diff = self._get_differences(seq1, seq2)
                data[taxa_pos[tx1]][taxa_pos[tx2]].append(aln_diff)

        data = np.array([[np.mean(y) if y else 0. for y in x] for x in data])
        mask = np.tri(data.shape[0], k=0)
        data = np.ma.array(data, mask=mask)

        return {"data": data,
                "labels": list(taxa_pos),
                "color_label": "Segregating sites"}

    def sequence_segregation_gene(self, gene_name, window_size):
        """
        Generates data for a sliding window analysis of segregating sites
        :param gene_name: string, name of gene in self.alignments
        :param window_size: size of the sliding window
        """

        aln_obj = self.alignments[gene_name]

        data = []

        for i in range(0, aln_obj.locus_length, window_size):

            segregating_sites = 0

            seqs = np.array([[y for y in x[i:i + window_size]] for x in
                             aln_obj.sequences()])

            for column in zip(*seqs):
                column = set([x for x in column if x != aln_obj.sequence_code[1]
                              and x != "-"])

                if len(column) > 1:
                    segregating_sites += 1

            data.append(segregating_sites)

        return {"data": data,
                "window_size": window_size,
                "ax_names": ["Sequence (bp)", "Segregating sites"],
                "table_header": ["Sequence (bp)", "Segregating sites"]}

    def length_polymorphism_correlation(self):
        """
        Generates data for a scatter plot and correlation analysis between
        alignment length and informative sites (polymorphic sites in at least
        two taxa)
        """

        data_length = []
        data_inf = []

        for aln in self.alignments.values():

            # Get informative sites for alignment
            inf_sites = self._get_informative_sites(aln)
            data_inf.append(inf_sites)

            # Get size
            data_length.append(aln.locus_length)

        return {"data": [data_length, data_inf],
                "ax_names": ["Alignment length", "Informative sites"],
                "table_header": ["Alignment length", "Informative sites"],
                "correlation": True}

    def taxa_distribution(self):
        """
        Generates data for a distribution of taxa frequency across alignments
        """

        data = []

        for aln in self.alignments.values():

            # Get number of taxa
            data.append(len(aln.alignment))

        return {"data": data,
                "ax_names": ["Number of taxa", "Frequency"],
                "title": "Distribution of taxa frequency",
                "table_header": ["Number of taxa", "Frequency"],
                "real_bin_num": True}

    @staticmethod
    def _mad_based_outlier(p, threshold=3.5):
        """
        An outlier detection method based on median absolute deviation. This
        code was adapted from http://stackoverflow.com/a/22357811/1990165.
        The usage of media is much less biased that the mean and is robust to
        smaller data sets.
        :param p: Num of observations
        :param threshold: modified Z-score to use as a threshold.
        """

        if len(p.shape) == 1:
            p = p[:, None]
        # Get sample median.
        median = np.median(p, axis=0)
        diff = np.sum((p - median) ** 2, axis=-1)
        diff = np.sqrt(diff)
        med_abs_deviation = np.median(diff)

        # The 0.6745is a constant that makes values roughly equivalent in
        # units to standard deviations.
        z_score = 0.6745 * diff / med_abs_deviation

        return z_score > threshold

    def outlier_missing_data(self):
        """
        Get data for outlier detection of genes based on the distribution of
        average missing data per gene. Data points will be based on the
        proportion of missing data symbols out of the possible total. For
        example, in an alignment with three taxa, each with 100 sites,
        the total possible missing data is 300 (100 * 3). Here, missing data
        will be gathered from all taxa and a proportion will be calculated
        based n the total possible
        """

        data_labels = []
        data_points = []

        for gn, aln in self.alignments.items():

            total_len = aln.locus_length * len(aln.alignment)
            gn_data = 0

            for seq in aln.sequences():
                gn_data += seq.count(self.sequence_code[1]) + \
                          seq.count(self.gap_symbol)

            m_data = float(gn_data) / float(total_len)
            data_points.append(m_data)
            data_labels.append(gn)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier labels
        outliers_labels = list(data_labels[self._mad_based_outlier(
            data_points)])

        return {"data": data_points,
                "outliers": outliers_points,
                "outliers_labels": outliers_labels,
                "ax_names": ["Proportion of missing data", "Frequency"]}

    def outlier_missing_data_sp(self):
        """
        Gets data for outlier detection of species based on missing data. For
        this analysis, genes for which a taxa is completely absent will be
        ignored for the calculations of that taxa. The reason for this,
        is that including genes where the taxon is absent would bias the
        outlier detection towards taxa that have low prevalence in the data
        set, even if they have low missing data in the alignments where they
        are present.
        """

        data = dict((tx, []) for tx in self.taxa_names)

        for aln in self.alignments.values():

            total_len = aln.locus_length

            # Get missing data for every taxon
            for tx in data:
                if tx in aln.alignment:
                    seq = aln.get_sequence(tx)
                    m_data = float(seq.count(self.sequence_code[1]) +
                                   seq.count(self.gap_symbol)) / \
                        float(total_len)
                    data[tx].append(m_data)

        # Get average for each taxon
        for tx, vals in data.items():
            data[tx] = np.mean(vals)

        # Prepara data for plotting
        data_points = []
        data_labels = []
        for tx, vals in data.items():
            data_points.append(vals)
            data_labels.append(tx)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier taxa
        outlier_labels = list(data_labels[self._mad_based_outlier(data_points)])

        return {"data": data_points,
                "outliers": outliers_points,
                "outliers_labels": outlier_labels,
                "ax_names": ["Proportion of missing data", "Frequency"]}

    def outlier_segregating(self):
        """
        Generates data for the outlier detection of genes based on
        segregating sites. The data will be based on the number of alignments
        columns with a variable number of sites, excluding gaps and missing
        data
        """

        data_points = []
        data_labels = []

        for aln in self.alignments.values():

            segregating_sites = 0

            for column in zip(*aln.sequences()):

                # Remove gaps and missing characters
                column = set([x for x in column if x != aln.sequence_code[1] and
                              x != self.gap_symbol])

                if len(column) > 1:
                    segregating_sites += 1

            # Get proportion of segregating sites for current alignment
            data_points.append(float(segregating_sites) /
                               float(aln.locus_length))
            data_labels.append(aln.name)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier taxa
        outlier_labels = list(data_labels[self._mad_based_outlier(data_points)])

        return {"data": data_points,
                "outliers": outliers_points,
                "outliers_labels": outlier_labels,
                "ax_names": ["Proportion of segregating sites", "Frequency"]}

    def outlier_segregating_sp(self):
        """
        Generates data for the outlier detection of species based on their
        average pair-wise proportion of segregating sites. Comparisons
        with gaps or missing data are ignored
        """

        data = OrderedDict((tx, []) for tx in self.taxa_names)

        for aln in self.alignments.values():

            for tx1, tx2 in itertools.combinations(data.keys(), 2):
                try:
                    seq1, seq2 = aln.get_sequence(tx1), aln.get_sequence(tx2)
                except KeyError:
                    continue

                s_data = self._get_differences_prop(seq1, seq2)

                data[tx1].append(s_data)
                data[tx2].append(s_data)

        for tx, vals in data.items():
            data[tx] = np.mean(vals)

        # Prepara data for plotting
        data_points = []
        data_labels = []
        for tx, vals in data.items():
            data_points.append(vals)
            data_labels.append(tx)

        data_points = np.asarray(data_points)
        data_labels = np.asarray(data_labels)

        # Get outliers
        outliers_points = data_points[self._mad_based_outlier(data_points)]
        # Get outlier taxa
        outlier_labels = list(data_labels[self._mad_based_outlier(data_points)])

        return {"data": data_points,
                "outliers": outliers_points,
                "outliers_labels": outlier_labels,
                "ax_names": ["Proportion of segregating sites", "Frequency"]}

__author__ = "Diogo N. Silva"
__copyright__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva", "Tiago F. Jesus"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Diogo N. Silva"
__email__ = "o.diogosilva@gmail.com"
__status__ = "Prototype"
