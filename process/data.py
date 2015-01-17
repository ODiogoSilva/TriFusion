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

from collections import OrderedDict


class PartitionException(Exception):
    pass


class Partitions():
    """
    The Partitions class is used to define partitions for Alignment objects and
    associate substitution models for each partition. Partitions may be set
    in two ways:

    ..: Partition files: Being Nexus charset blocks and RAxML partition files
    currently supported
    ..: Tuple-like objects: Containing the ranges and names of the partitions

    A SubstitutionModels object will be associated to each partition, and by
    default there will be no substitution model selected.
    """

    def __init__(self):
        """
        Setting the self._partitions private attribute. This will contain an
        ordered dictionary with the partition names as keys and information on
        their range and substitution model object as values. The ranges will be
        in tuple format with the initial position as the first element and
        final position as the second element

        e.g. self.partitions["GeneA"] = [(0, 953), SubstitutionModels]
        Defines the partition GeneA whose sequence spans from 0 to the 953rd
        character
        """
        self.partitions = OrderedDict()

        """
        partitions_nice is similar to partitions but the locus ranges are
        already in an appropriate format to be used in output files. This is
        necessary because the sequence index for python has an offset of 0,
        while all subsequent programs that use partition ranges use an index
        of 1.
        """
        self.partitions_nice = OrderedDict()

        """
        The private self._models attribute will contain the same key list as
        self._partitions and will associate the substitution models to each
        partitions
        """
        self._models = OrderedDict()

        """
        The counter attribute will be used as an indication of where the last
        partition ends when one or more partitions are added
        """
        self.counter = 0

    def __iter__(self):
        """
        The class iterator will iterate over a list containing the partition
        names and a modified version of their ranges that is compatible with
        other software (unlike the 0 offset of python)
        :return:
        """

        return iter(self.partitions.items())

    #===========================================================================
    # Parsers
    #===========================================================================

    @staticmethod
    def _get_file_format(partition_file):
        """ Tries to guess the format of the partition file (Whether it is
         Nexus of RAxML's) """
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

    def read_from_file(self, partitions_file):
        """
        This function parses a file containing partitions. Supports
        partitions files similar to RAxML's and NEXUS charset blocks. The
        NEXUS file, however, must only contain the charset block. The
        model_nexus argument provides a namespace for the model variable in
        the nexus format, since this information is not present in the file.
        However, it assures consistency on the Partition object
        :param partitions_file: string, file name of the file containing the
        partitions
        """

        # Get the format of the partition file
        partition_format = self._get_file_format(partitions_file)

        part_file = open(partitions_file)

        if partition_format == "raxml":
            for line in part_file:
                fields = line.split(",")
                # Get model name as string
                model_name = fields[0]
                # Get partition name as string
                partition_name = fields[1].split("=")[0].strip()
                # Get partition range as list of int
                partition_range_temp = fields[1].split("=")[1]
                partition_range = [int(x) for x in
                                   partition_range_temp.strip().split("-")]
                # Add information to partitions storage
                self.add_partition(partition_name, locus_range=partition_range,
                                   model=model_name)

        elif partition_format == "nexus":
            for line in part_file:
                self.read_from_nexus_string(line)

    def read_from_nexus_string(self, string):
        """
        Parses the partition defined in a charset command
        :param string: string with the charset command.
        :return:
        """

        fields = string.split("=")
        partition_name = fields[0].split()[1].strip()

        # In case the charset command is setting a codon specific partition,
        # such as 'charset EF1a_2nd = 1245-1611\3;' no partition will be added
        try:
            partition_range = [int(x) for x in
                               fields[1].replace(";", "").strip().split("-")]
        except ValueError:
            return 0

        self.add_partition(partition_name, locus_range=partition_range)

    def read_from_dict(self, dict_obj):
        """
        Parses partitions defined and stored in a special OrderedDict. The
        values of dict_obj should be the partition names and their corresponding
        values should contain the loci range and substitution model, if any

        Example

        dict_obj = OrderedDict(("GeneA", [(0,234), "GTR"]), ("GeneB", [(235,
                                865), "JC"))
        :param dict_obj: And OrderedDict object
        """

        for k, v in dict_obj:
            # Determining if value contains only the range or the substitution
            # model as well
            if len(v) > 1:
                self.add_partition(k, locus_range=v[0], model=v[1])
            else:
                self.add_partition(k, locus_range=v[0])

    def is_single(self):
        """
        :return: Boolean. Returns True is there is only a single partition
        defined, and False if there are multiple partitions
        """

        if len(self.partitions) == 1:
            return True
        else:
            return False

    def set_single_partition(self, name, length):
        """
        Sets the self.partitions attribute consisting of a single partition.
        Usually used by single alignments.

        Adding or setting partitions always associates them with a substitution
        model object

        :param length: int. Lenght of the alignment
        :param name: string. Name of the alignment
        """

        # Add partition
        self.partitions[name] = (0, length)
        self.partitions_nice[name] = (1, length)
        # Associate substitution model
        self._models[name] = SubstitutionModels()

        #Update counter
        self.counter += length

    def add_partition(self, name, length=None, locus_range=None, model=None):
        """
        Adds a new partition providing the length or the range of current
        alignment. If both are provided, the length takes precedence.
        :param name: string. Name of the alignment
        :param length: int. Length of the alignment
        :param locus_range: list/tuple. Range of the alignment
        :param model: string. [optional] Name of the substitution model
        """

        if name in self.partitions:
            raise PartitionException("Partition name %s is already in partition"
                                     "table" % name)

        if length:
            self.partitions[name] = (self.counter, self.counter + length)
            self.partitions_nice[name] = (self.counter + 1, self.counter +
                                          length)
            self.counter += length
        elif locus_range:
            self.partitions[name] = (locus_range[0], locus_range[1])
            self.partitions_nice[name] = (locus_range[0] + 1, locus_range[1])
            self.counter = locus_range[1]

        self._models[name] = SubstitutionModels()

    # def write_to_file(self, output_format, output_file, model="LG"):
    #     """ Writes the Partitions object into an output file according to the
    #      output_format. The supported output formats are RAxML and Nexus.
    #      9The model option is for the RAxML format """
    #
    #     if output_format == "raxml":
    #         outfile_handle = open(output_file + ".part.File", "w")
    #         for part in self.partitions:
    #             partition_name = part[0]
    #             partition_range = "-".join([x for x in part[1]])
    #             outfile_handle.write("%s, %s = %s\n" % (model,
    #                                                     partition_name,
    #                                                     partition_range))
    #
    #         outfile_handle.close()
    #
    #     elif output_format == "nexus":
    #         outfile_handle = open(output_file + ".charset", "w")
    #         for part in self.partitions:
    #             outfile_handle.write("charset %s = %s;\n" % (
    #                                  part[1],
    #                                  "-".join(part[2])))
    #
    #         outfile_handle.close()
    #
    #     return 0


class Zorro ():

    def __init__(self, alignment_list, suffix="_zorro.out"):

        def zorro2rax(alignment_list):
            """ Function that converts the floating point numbers contained
            in the original zorro output files into integers that can be
            interpreted by RAxML. If multiple alignment files are provided,
            it also concatenates them in the same order """
            weigths_storage = []
            for alignment_file in alignment_list:
                # This assumes that the prefix of the
                zorro_file = alignment_file.split(".")[0] + self.suffix
                # alignment file is shared with the corresponding zorro file
                zorro_handle = open(zorro_file)
                weigths_storage += [round(float(weigth.strip())) for
                                    weigth in zorro_handle]
            return weigths_storage

        self.suffix = suffix
        self.weigth_values = zorro2rax(alignment_list)

    def write_to_file(self, output_file):
        """ Creates a concatenated file with the zorro weights for the
        corresponding alignment files """
        outfile = output_file + "_zorro.out"
        outfile_handle = open(outfile, "w")
        for weigth in self.weigth_values:
            outfile_handle.write("%s\n" % weigth)
        outfile_handle.close()


class SubstitutionModels():
    """
    This class handles the storage, parsing and retrieval of substitution models
    in several formats.
    """

    _models = {"mrbayes": {}}

    #===========================================================================
    #   MrBayes models
    #===========================================================================
    """
    MrBayes substitution models are stored in the dictionary _models["mrbayes"].
    The keys of the dictionary are the name of the substitution models (usually
    in capital letters) and the values will contain the instructions to
    specific such model in a list. Each element of the list corresponds to one
    line
    """

    # GTR
    _models["GTR"] = ["lset nst=6"]

    # SYM
    _models["SYM"] = ["lset nst=6", "prset statefreqpr=fixed(equal)"]

    # HKY
    _models["HKY"] = ["lset nst=2"]

    # K2P
    _models["K2P"] = ["lset nst=2", "prset statefreqpr=fixed(equal)"]

    # F81
    _models["F81"] = ["lset nst=1"]

    # JC
    _models["JC"] = ["lset nst=1", "prset statefreqpr=fixed(equal)"]


__author__ = "Diogo N. Silva"
__copyright__ = "Diogo N. Silva"
__credits__ = ["Diogo N. Silva"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Diogo N. Silva"
__email__ = "o.diogosilva@gmail.com"
__status__ = "Prototype"