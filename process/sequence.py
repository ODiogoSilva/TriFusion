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
from process.base import *
from process.missing_filter import MissingFilter
from process.data import Partitions

# Other imports
from collections import OrderedDict
import re
from os import sep
from os.path import join

## TODO: Create a SequenceSet class for sets of sequences that do not conform
# to an alignment, i.e. unequal length.This would eliminate the problems of
# applying methods designed for alignments to sets of sequences with unequal
# length would allows these sets of sequences to have methods of their own.
# However, the __init__ of the current Alignment object should apply to both
# SequenceSet and Alignment classes. So, I'll have to re-structure the code
# somehow.
## TODO After creating the SequenceSet class, an additional class should be
# used to make the triage of files to either the Alignment or SequenceSet
# classes


class AlignmentException(Exception):
    pass


class Alignment (Base):

    def __init__(self, input_alignment, input_format=None, alignment_name=None,
                 partitions=None):
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
        self.locus_length = 0

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

        # If the object is initialized with a string
        if isinstance(input_alignment, str):

            """
            Sets the alignment object name based on the input alignment file
            name. This will remove the extension and a preceding path, if
            it exists
            """
            self.name = input_alignment.split(sep)[-1].split(".")[0]

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
                self.alignment = finder_content
                # Setting the sequence code attribute for seq type checking
                # in AlignmentList
                self.sequence_code = None

        # In case the class is initialized with a dictionary object
        elif isinstance(input_alignment, OrderedDict):

            # The name of the alignment (str)
            self.name = alignment_name
            # Gets several attributes from the dictionary alignment
            self._init_dicobj(input_alignment)
            #The input format of the alignment (str)
            self.input_format = input_format

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

        self.sequence_code = self.guess_code(list(dictionary_obj.values())[0])
        self.alignment = dictionary_obj
        self.locus_length = len(list(dictionary_obj.values())[0])

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

        # Storage taxa names and corresponding sequences in an ordered
        # Dictionary
        self.alignment = OrderedDict()

        file_handle = open(input_alignment)

        #=======================================================================
        # PARSING PHYLIP FORMAT
        #=======================================================================
        if alignment_format == "phylip":
            # Get the number of taxa and sequence length from the file header
            header = file_handle.readline().split()
            self.locus_length = int(header[1])
            self.partitions.set_length(self.locus_length)
            for line in file_handle:
                try:
                    taxa = line.split()[0].replace(" ", "")
                    taxa = self.rm_illegal(taxa)
                    try:
                        sequence = line.split()[1].strip().lower()
                    except IndexError:
                        sequence = ""

                    self.alignment[taxa] = sequence
                except IndexError:
                    pass

                    ## TO DO: Read phylip interleave

            # Updating partitions object
            self.partitions.add_partition(self.name, self.locus_length)

        #=======================================================================
        # PARSING FASTA FORMAT
        #=======================================================================
        elif alignment_format == "fasta":
            for line in file_handle:
                if line.strip().startswith(">"):
                    taxa = line[1:].strip().replace(" ", "_")
                    taxa = self.rm_illegal(taxa)
                    self.alignment[taxa] = ""
                elif line.strip() != "" and taxa is not None:
                    self.alignment[taxa] += line.strip().lower().\
                        replace(" ", "")
            self.locus_length = len(list(self.alignment.values())[0])
            self.partitions.set_length(self.locus_length)

            # Updating partitions object
            self.partitions.add_partition(self.name, self.locus_length)

        #=======================================================================
        # PARSING NEXUS FORMAT
        #=======================================================================
        elif alignment_format == "nexus":
            counter = 0
            for line in file_handle:
                # Skips the nexus header
                if line.strip().lower() == "matrix" and counter == 0:
                    counter = 1
                # Stop sequence parser here
                elif line.strip() == ";" and counter == 1:
                    counter = 2
                    self.locus_length = len(list(self.alignment.values())[0])
                    self.partitions.set_length(self.locus_length)
                # Start parsing here
                elif line.strip() != "" and counter == 1:
                    taxa = line.strip().split()[0].replace(" ", "")
                    taxa = self.rm_illegal(taxa)
                    # This accommodates for the interleave format
                    if taxa in self.alignment:
                        self.alignment[taxa] += "".join(
                            line.strip().split()[1:]).lower()
                    else:
                        self.alignment[taxa] = "".join(
                            line.strip().split()[1:]).lower()

                # If partitions are specified using the charset command, this
                # section will parse the partitions
                elif line.strip().startswith("charset"):
                    self.partitions.read_from_nexus_string(line)

            # If no partitions have been added during the parsing of the nexus
            # file, set a single partition
            if self.partitions.partitions == OrderedDict():
                self.partitions.add_partition(self.name,
                                                     self.locus_length)

        # Checks the size consistency of the alignment
        if size_check is True:
            self.is_alignment = self.check_sizes(self.alignment,
                                                 input_alignment)

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
        Returns a list with the taxa contained in the alignment
        """

        taxa = [sp for sp in self.alignment]

        return taxa

    def iter_sequences(self):
        """
        Returns a list with the sequences contained in the alignment
        """

        sequences = [seq for seq in self.alignment]

        return sequences

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
        except (FileNotFoundError, IndexError):
            taxa_list = taxa_list_file

        if mode == "remove":
            remove(taxa_list)
        if mode == "inverse":
            inverse(taxa_list)

    def collapse(self, write_haplotypes=True, haplotypes_file=None,
                 haplotype_name="Hap"):
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
            #  use the input alignment name as reference
            if haplotypes_file is None:
                haplotypes_file = self.name.split(".")[0]
            self.write_loci_correspondence(correspondence_dic, haplotypes_file)

    @staticmethod
    def write_loci_correspondence(dic_obj, output_file):
        """
        This function supports the collapse method by writing the
        correspondence between the unique haplotypes and the loci into a
        new file
        """

        output_handle = open(output_file + ".haplotypes", "w")

        for haplotype, taxa_list in dic_obj.items():
            output_handle.write("%s: %s\n" % (haplotype, "; ".join(taxa_list)))

        output_handle.close()

    def set_partitions(self, partitions):
        """
        Updates the Partitions object of the current alignment.

        :param partitions: Partitions object. Use one of the Partition parsers
        to retrieve partitions information from files or python data structures.
        See process.data.Partitions documentation
        """

        self.partitions = partitions

    def reverse_concatenate(self):
        """
        This function divides a concatenated file according to the
        partitions set in self.partitions and returns an AlignmentList object
        """

        concatenated_aln = AlignmentList([])

        for name, part_range in self.partitions:

            current_dic = OrderedDict()
            for taxon, seq in self.alignment.items():
                sub_seq = seq[part_range[0]:part_range[1]]

                # If sub_seq is not empty (only gaps or missing data)
                if sub_seq.replace(self.sequence_code[1], "") != "":
                    current_dic[taxon] = sub_seq

            current_partition = Partitions()
            current_partition.add_partition(name, part_range[1] -
                                                   part_range[0])
            current_aln = Alignment(current_dic, input_format=self.input_format,
                                    partitions=current_partition,
                                    alignment_name=name)

            concatenated_aln.add_alignment(current_aln)

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

    def filter_missing_data(self, gap_threshold, missing_threshold):
        """
        Filters gaps and true missing data from the alignment using tolerance
        thresholds for each type of missing data. Both thresholds are maximum
        percentages of sites in an alignment column containing the type of
        missing data. If gap_threshold=50, for example, alignment columns with
        more than 50% of sites with gaps are removed.

        :param gap_threshold: int ranging from 0 to 100.
        :param missing_threshold: int ranging from 0 to 100.
        """

        # When the class is initialized, it performs the basic filtering
        # operations based on the provided thresholds
        alignment_filter = MissingFilter(self.alignment,
                                         gap_threshold=gap_threshold,
                                         missing_threshold=missing_threshold,
                                         gap_symbol="-",
                                         missing_symbol=self.sequence_code[1])

        # Replace the old alignment by the filtered one
        self.alignment = alignment_filter.alignment
        self.locus_length = alignment_filter.locus_length

    def write_to_file(self, output_format, output_file, new_alignment=None,
                      seq_space_nex=40, seq_space_phy=30, seq_space_ima2=10,
                      cut_space_nex=50, cut_space_phy=50, cut_space_ima2=8,
                      interleave=False, gap="-", model_phylip="LG",
                      outgroup_list=None, ima2_params=None, use_charset=True,
                      partition_file=True, output_dir=None):
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
        """

        # If this function is called in the AlignmentList class, there may
        # be a need to specify a new alignment dictionary, such as a
        # concatenated one
        if new_alignment is not None:
            alignment = new_alignment
        else:
            alignment = self.alignment

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
                taxon, population = re.split(r'[\t;,]', line)
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

            if self.partitions.is_single() is False:
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
                                seq = self.alignment[taxon][
                                      lrange[0]:lrange[1]].upper()
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
                        (lrange[1]) - lrange[0],
                        mutational_model,
                        inheritance_scalar))

                    # Write sequence data according to the order of the
                    # population mapping file
                    for taxon, seq in new_alignment:
                        out_file.write("%s%s\n" % (taxon, seq))

            if self.partitions.is_single():
                #Write the header for the single
                out_file.write("%s %s %s %s %s\n" % (
                               partition,
                               " ".join(population_storage.values()),
                               self.locus_length,
                               mutational_model,
                               inheritance_scalar))

                #Write sequence data
                for population, taxa_list in population_storage.items():
                    for taxon in taxa_list:
                        seq = self.alignment[taxon].upper()
                        out_file.write("%s%s\n" %
                                      (taxon[:cut_space_ima2].ljust(
                                         seq_space_ima2),
                                       seq))

        # Writes file in phylip format
        if "phylip" in output_format:

            out_file = open(output_file + ".phy", "w")
            out_file.write("%s %s\n" % (len(alignment), self.locus_length))
            for key, seq in alignment.items():
                    out_file.write("%s %s\n" % (
                                   key[:cut_space_phy].ljust(seq_space_phy),
                                   seq.upper()))

            # In case there is a concatenated alignment being written
            if self.partitions.is_single() is False and partition_file:
                partition_file = open(output_file + "_part.File", "w")
                for name, lrange in self.partitions:
                    partition_file.write("%s, %s = %s\n" % (
                                         model_phylip,
                                         name,
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

                    for taxon, seq in self.alignment.items():
                        out_file.write("%s  %s\n" % (
                                       taxon[:cut_space_phy].ljust(
                                         seq_space_phy),
                                       seq[lrange[0]:lrange[1]].upper()))
            else:
                out_file.write("%s %s\n" % (taxa_number, self.locus_length))
                for taxon, seq in self.alignment.items():
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
                                    self.sequence_code[0],
                                    self.locus_length - 1,
                                    self.restriction_range,
                                    gap,
                                    self.sequence_code[1]))
                else:
                    out_file.write("#NEXUS\n\nBegin data;\n\tdimensions "
                                   "ntax=%s nchar=%s ;\n\tformat datatype=%s "
                                   "interleave=yes gap=%s missing=%s ;\n\t"
                                   "matrix\n" % (
                                   len(alignment),
                                   self.locus_length,
                                   self.sequence_code[0], gap,
                                   self.sequence_code[1]))
                counter = 0
                for i in range(90, self.locus_length, 90):
                    for key, seq in alignment.items():
                        out_file.write("%s %s\n" % (
                                       key[:cut_space_nex].ljust(
                                         seq_space_nex),
                                       seq[counter:i].upper()))
                    else:
                        out_file.write("\n")
                        counter = i
                else:
                    for key, seq in alignment.items():
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

                for key, seq in alignment.items():
                    out_file.write("%s %s\n" % (key[:cut_space_nex].ljust(
                        seq_space_nex), seq))
                out_file.write(";\n\tend;")

            if use_charset:
                # Writing partitions, if any
                if self.partitions.is_single() is False:
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

                    out_file.write("\tpartition part = %s: %s;\n\tset partition="
                                   "part;\nend;\n" %
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
            for key, seq in self.alignment.items():
                out_file.write(">%s\n%s\n" % (key, seq.upper()))

            out_file.close()


class AlignmentList (Base):
    """
    At the most basic instance, this class contains a list of Alignment
    objects upon which several methods can be applied. It only requires either
    a list of alignment files or. It inherits methods from Base and
    Alignment classes for the write_to_file methods.
    """

    def __init__(self, alignment_list, verbose=True):
        """
        :param alignment_list: List of Alignment objects
        :param verbose: Boolean. If True, it prints a progression bar to the
        terminal when loading data
        """

        self.log_progression = Progression()

        self.alignment_object_list = []
        # Stores badly formatted/invalid input alignments
        self.bad_alignments = []

        # if type(alignment_list[0]) is str:
        if alignment_list:
            self.log_progression.record("Parsing file", len(alignment_list))
            for alignment in alignment_list:

                if verbose is True:
                    self.log_progression.progress_bar(
                        alignment_list.index(alignment) + 1)

                alignment_object = Alignment(alignment)
                if isinstance(alignment_object.alignment, Exception):
                    self.bad_alignments.append(alignment_object)
                else:
                    self.alignment_object_list.append(alignment_object)

        # elif type(alignment_list[0]) is OrderedDict:
        #
        #     # Setting class flag so that the object is aware that is a result
        #     #  of reverse concatenation and not an actual list of alignments
        #     self.reverse_concatenation = True
        #
        #     for alignment, model, name in zip(alignment_list, model_list,
        #                                       name_list):
        #
        #         alignment_object = Alignment(alignment, model_list=[model],
        #                                      alignment_name=name)
        #         # In some cases, like reverse concatenating, the input format
        #         #  must be manually set
        #         if input_format is not None:
        #             alignment_object._set_format(input_format)
        #         self.alignment_object_list.append(alignment_object)

        #### SETTING GENERAL ATTRIBUTES

        # list of file names, complete with path
        self.filename_list = self._get_filename_list()

        # list of taxon names
        self.taxa_names = self._get_taxa_list()

    def __iter__(self):
        """
        Iterate over Alignment objects
        """
        return iter(self.alignment_object_list)

    def format_list(self):
        """
        :return: List with the unique sequence types of the Alignment objects
        """

        return list(set([x.sequence_code[0] for x in
                         self.alignment_object_list if x]))

    def _get_format(self):
        """
        Gets the input format of the first alignment in the list
        """
        # TODO: To support different input format, provide information on how
        # many alignment are in a certain input format

        return self.alignment_object_list[0].input_format

    def _get_taxa_list(self):
        """
        Gets the full taxa list of all alignment objects
        :return full_taxa. List of taxa names in the AlignmentList
        """

        full_taxa = []

        for alignment in self.alignment_object_list:
            diff = set(alignment.iter_taxa()) - set(full_taxa)
            if diff != set():
                full_taxa.extend(diff)

        return full_taxa

    def _get_filename_list(self):
        """
        Returns a list with the input file names
        """
        return [alignment.name for alignment in self.alignment_object_list]

    def add_alignment(self, alignment_obj):
        """
        Adds a new Alignment object
        :param alignment_obj: Alignment object
        """

        self.alignment_object_list.append(alignment_obj)

        # Update taxa names with the new alignment
        self.taxa_names = self._get_taxa_list()

    def retrieve_alignment(self, name):
        """
        :param name: string. Name of the input alignment
        :return: Returns an Alignment object with a given name attribute
        """

        name = name.split(sep)[-1].split(".")[0]

        alignment_obj = [x for x in self.alignment_object_list if x.name ==
                         name][0]

        if alignment_obj:
            return alignment_obj
        else:
            return None

    def iter_alignment_dic(self):
        """
        :return: List of the dictionary alignments
        """

        return [alignment.alignment for alignment in self.alignment_object_list]

    def write_taxa_to_file(self):
        """
        Compiles the taxa names of all alignments and writes them in a single
        column .csv file
        """

        output_handle = open("Taxa_list.csv", "w")

        for taxon in self._get_taxa_list():
            output_handle.write(taxon + "\n")

        output_handle.close()

    def concatenate(self, progress_stat=True):
        """
        Concatenates multiple sequence alignments creating a single alignment
        object and the auxiliary Partitions object defining the partitions
        of the concatenated alignment
        :param progress_stat: Boolean. If True, it will print a progression bar
        to the terminal
        :return concatenated_alignment: Alignment object
        """

        self.log_progression.record("Concatenating file", len(
                                    self.alignment_object_list))

        taxa_list = self._get_taxa_list()

        # Initializing alignment dict to store the alignment information
        concatenation = OrderedDict([(key, []) for key in taxa_list])
        # Instanciating a Partitions object to store information on the
        # partitions and substitution models
        partitions = Partitions()

        for alignment_object in self.alignment_object_list:

            # When set to True, this statement produces a progress status on
            # the terminal
            if progress_stat is True:
                self.log_progression.progress_bar(self.alignment_object_list.
                                                  index(alignment_object) + 1)

            # Setting the missing data symbol
            missing = alignment_object.sequence_code[1]

            # Appending sequence data from current alignment to concatenated
            # alignment. Appending each sequence data in a list and in the end
            # converting that list into a string is much much faster than
            # concatenating strings at this point
            for taxa in taxa_list:
                try:
                    sequence = alignment_object.alignment[taxa]
                    concatenation[taxa].append(sequence)
                except KeyError:
                    concatenation[taxa].append(missing
                                             * alignment_object.locus_length)

            # UPDATING PARTITIONS OBJECT
            # This is intended to support the concatenation of both individual
            # and concatenated files. If one of the alignment objects is already
            # a concatenated alignment, this will add each partition to the
            # new partitions object
            if not alignment_object.partitions.is_single():
                for k, v in \
                        alignment_object.partitions.partitions.items():
                            partitions.add_partition(k, locus_range=v[0],
                                                     codon=v[1],
                                                     use_counter=True)
            else:
                partitions.add_partition(alignment_object.name,
                                         length=alignment_object.locus_length,
                                         use_counter=True)

        # Each taxa is a list of strings for each alignment object, so here
        # the list is being converted into a string
        for taxa, seq in concatenation.items():
            concatenation[taxa] = "".join(seq)

        # Create the concatenated file in an Alignment object
        concatenated_alignment = Alignment(concatenation,
                                           input_format=self._get_format(),
                                           partitions=partitions)

        return concatenated_alignment

    def filter_missing_data(self, gap_threshold, missing_threshold,
                            verbose=True):
        """
        Wrapper of the filter_missing_data method of the Alignment object.
        See the method's documentation.
        """

        self.log_progression.record("Filtering file",
                                    len(self.alignment_object_list))
        for alignment_obj in self.alignment_object_list:

            if verbose is True:
                self.log_progression.progress_bar(
                    self.alignment_object_list.index(alignment_obj))

            alignment_obj.filter_missing_data(gap_threshold=gap_threshold,
                          missing_threshold=missing_threshold)

    def remove_taxa(self, taxa_list, verbose=False, mode="remove"):
        """
        Wrapper of the remove_taxa method of the Alignment object for
        multiple alignments. It current supports two modes:

            ..:remove: removes specified taxa
            ..:inverse: removes all but the specified taxa
        """

        if verbose is True:
            self.log_progression.write("Removing taxa")

        for alignment_obj in self.alignment_object_list:
            alignment_obj.remove_taxa(taxa_list, mode=mode)

        # Updates taxa names
        self.taxa_names = self._get_taxa_list()

    def remove_file(self, filename_list):
        """
        Removes alignment objects based on their name attribute
        :param filename_list: list with the names of the alignment objects to
        be removed
        """

        for nm in filename_list:
            nm = nm.split(sep)[-1].split(".")[0]
            self.alignment_object_list = [x for x in self.alignment_object_list
                                          if nm != x.name]

        # Updates taxa names
        self.taxa_names = self._get_taxa_list()

    def clear_files(self):
        """
        Removes all Alignment objects from the AlignmentList
        """

        self.alignment_object_list = []

    def select_by_taxa(self, taxa_list, mode="strict", verbose=True):
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

        :param verbose. Boolean. If True, a progression bar will be printed
        on the terminal
        """

        selected_alignments = []

        if verbose is True:
            self.log_progression.write("Selecting alignments")

        # taxa_list may be a file name (string) or a list containing the name
        #  of the taxa. If taxa_list is a file name this code will parse the
        # csv file and return a list of the taxa. Otherwise, the taxa_list
        # variable remains the same.
        try:
            file_handle = open("".join(taxa_list))
            taxa_list = self.read_basic_csv(file_handle)
        except FileNotFoundError:
            pass

        for alignment_obj in self.alignment_object_list:

            alignment_taxa = alignment_obj.iter_taxa()

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

        for alignment_obj in self.alignment_object_list:
            alignment_obj.code_gaps()

    def collapse(self, write_haplotypes=True, haplotypes_file="",
                 haplotype_name="Hap"):
        """
        Wrapper for the collapse method of the Alignment object. If
        write_haplotypes is True, the haplotypes file name will be based on the
        individual input file

        :param write_haplotypes: Boolean, if True, a haplotype list
        mapping the haplotype names file will be created for each individual
        input alignment.

        :param haplotype_name: String, Custom name of the haplotypes
        """

        for alignment_obj in self.alignment_object_list:
            if write_haplotypes:
                # Set name for haplotypes file
                output_file = alignment_obj.name.split(".")[0] + haplotypes_file
                alignment_obj.collapse(haplotypes_file=output_file,
                                       haplotype_name=haplotype_name)
            else:
                alignment_obj.collapse(write_haplotypes=False,
                                       haplotype_name=haplotype_name)

    def _reverse_concatenate(self, partition_obj):
        """
        Internal function to reverse concatenate an alignment according to
        defined partitions in a Partitions object

        This will only work if alignment_object_list has one alignment, as it
         is intended to be a wrapper of sorts for the Alignment object method

        :param partition_obj: Partitions object, containing the partitions of
        the input file
        :return: AlignmentList object with individual alignments
        """

        if len(self.alignment_object_list) == 1:
            aln_obj = self.alignment_object_list[0]
            return aln_obj.reverse_concatenate(partition_obj)

    def write_to_file(self, output_format, output_suffix="", interleave=False,
                      outgroup_list=None, partition_file=True, output_dir=None,
                      use_charset=True):
        """ This method writes a list of alignment objects or a concatenated
         alignment into a file """

        for alignment_obj in self.alignment_object_list:
            if alignment_obj.input_format in output_format:
                output_file_name = alignment_obj.name.split(".")[0]\
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
                                        use_charset=use_charset)

__author__ = "Diogo N. Silva"
__copyright__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva", "Tiago F. Jesus"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Diogo N. Silva"
__email__ = "o.diogosilva@gmail.com"
__status__ = "Prototype"