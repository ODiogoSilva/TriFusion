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

import re
from os.path import basename, splitext, join
from os import sep
from collections import OrderedDict


class PartitionException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class InvalidPartitionFile(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Partitions(object):
    """Alignment partitions interface for `Alignment` and `AlignmentList`.

    The Partitions class is used to define partitions for `Alignment`
    and `AlignmentList` objects and associate substitution models to
    each partition. After instantiating, partitions may be set in two ways:

      - Partition files: Being Nexus charset blocks and RAxML partition files
        currently supported
      - Tuple-like objects: Containing the ranges and names of the partitions

    Attributes
    ----------
    partition_length : int
        Length of the total partitions.
    partitions : OrderedDict
        Storage of partition names (key) and their range (values).
    partitions_index : list
        The index (starting point) for each partition, including codon
        partitions.
    partitions_alignments : OrderedDict
        Storage of the partition names (key) and their corresponding
        alignment files (values).
    alignments_range : OrderedDict
        Storage of the alignment names (key) and their range (values).
    models : OrderedDict
        Storage of partition names (key) and their models (values).
    merged_files : dict
        Storage of the original range (values) of every alignment file (key).
    counter : int
        Indicator of where the last partition ended.
    partition_format : str
        Format of the original partition file, if any.
    """

    _models = {"mrbayes": {}}

    # =========================================================================
    #   MrBayes models
    # =========================================================================
    """
    MrBayes substitution models are stored in the dictionary 
    _models["mrbayes"]. The keys of the dictionary are the name of the 
    substitution models (usually in capital letters) and the values will 
    contain the instructions to specific such model in a list. Each element 
    of the list corresponds to one line
    """

    # GTR
    _models["mrbayes"]["GTR"] = ["nst=6", "statefreqpr=dirichlet(1,1,1,1)"]

    # SYM
    _models["mrbayes"]["SYM"] = ["nst=6", "statefreqpr=fixed(equal)"]

    # HKY
    _models["mrbayes"]["HKY"] = ["nst=2", "statefreqpr=dirichlet(1,1,1,1)"]

    # K2P
    _models["mrbayes"]["K2P"] = ["nst=2", "statefreqpr=fixed(equal)"]

    # F81
    _models["mrbayes"]["F81"] = ["nst=1", "statefreqpr=dirichlet(1,1,1,1)"]

    # JC
    _models["mrbayes"]["JC"] = ["nst=1", "statefreqpr=fixed(equal)"]

    def __init__(self):

        self.partition_length = 0
        """
        The length of the locus may be necessary when partitions are defined
        in the input files using the "." notation, meaning the entire locus.
        Therefore, to convert this notation into workable integers, the size
        of the locus must be provided using the set_length method.
        """

        self.partitions = OrderedDict()
        """
        partitions will contain the name and range of the partitions for a given
        alignment object. Both gene and codon partitions will be stored in this
        attribute, but gene partitions are the main entries. An example of
        different stored partitions is::

            partitions = {"partitionA": ((0, 856), False),
                          "partitionB": ((857, 1450), [857,858,859] }

        "partitionA" is a simple gene partition ranging from 0 to 856, while
        "partitionB" is an assembly of codon partitions. The third element of
        the tuple is destined to codon partitions. If there are none, it should
        be False. If there are codon partitions, a list should be provided with
        the desired initial codons. In the example above, "partitionB" has
        actually 3 partitions starting at the first, second and third sequence
        nucleotide of the main partition.
        """

        self.partitions_index = []
        """
        partitions_index will remember the index of all added partitions. This
        attribute was created because codon models are added to the same parent
        partitions, thus losing their actual index. This is important for
        Nexus files, where models are applied to the index of the partition.
        This will simply store the partition names, which can be accessed using
        their index, or searched to return their index. To better support codon
        partitions, each entry in the partitions_index will consist in a list,
        in which the first element is the partition name, and the second element
        is the index of the subpartition. An example would be::

            partitions_index = [["partA", 0], ["partA", 1], ["partA", 2],
                                ["partB", 0]]

        in which, partA has 3 codon partitions, and partB has only one partition

        """

        self.partitions_alignments = OrderedDict()
        """
        The partitions_alignments attribute will associate the partition with
        the corresponding alignment files. For single alignment partitions,
        this will provide information on the file name. For multiple 
        alignments, besides the information of the file names, it will 
        associate which alignments are contained in a given partition and 
        support multi alignment partitions. An example would be::

            partitions_alignments = {"PartitionA": ["FileA.fas"],
                                     "PartitionB": ["FileB.fas", "FileC.fas"]}

        """

        self.alignments_range = OrderedDict()
        """
        """

        self.models = OrderedDict()
        """
        The self.models attribute will contain the same key list as
        self.partitions and will associate the substitution models to each
        partitions. For each partition, the format should be as follows::

            models["partA"] = [[[..model_params..]],[..model_names..],
                               ["12", "3"]]

        The first element is a list that may contain the substitution model
        parameters for up to three subpartitions, the second element is also
        a list with the corresponding names of the substitution models and
        the third list will store any links between models
        """

        self.merged_files = {}
        """
        This attribute will keep a record of the original ranges of every file
        that was merged. This is useful to split partitions according to files
        or to undo any changes. Each entry should be::

            {"alignment_file1": (0, 1234), "alignment_file2": (3444, 6291)}
        """

        self.counter = 0
        """
        The counter attribute will be used as an indication of where the last
        partition ends when one or more partitions are added
        """

        self.partition_format = None

    def __iter__(self):
        """Iterator behavior for `Partitions`.

        The class iterator will iterate over a list containing the partition
        names and a modified version of their ranges that is compatible with
        other software (unlike the 0 offset of python)

        Returns
        _ : iter
            Iterator of `partitions.items()`.
        """

        return iter(self.partitions.items())

    def reset(self, keep_alignments_range=False):
        """Clears partitions and attributes

        Clears partitions and resets object to __init__ state. The original
        alignment range can be retained by setting the `keep_alignments_range`
        argument to True.

        Parameters
        ----------
        keep_alignments_range : bool
            If True, the `alignments_range` attribute will not be reset.
        """

        self.partitions = OrderedDict()
        self.partitions_index = []
        self.partitions_alignments = OrderedDict()
        self.models = OrderedDict()
        self.counter = 0
        if not keep_alignments_range:
            self.alignments_range = OrderedDict()

    def iter_files(self):
        """Iterates over `partitions_alignments.items()`.

        Returns
        -------
        _ : iter
            Iterator of `partitions_alignments.items()`.
        """

        return iter(self.partitions_alignments.items())

    def set_length(self, length):
        """Set total length of current locus (over all partitions).

        Sets the length of the locus. This may be important to convert certain
        partition defining nomenclature, such as using the "." to indicate
        whole length of the alignment

        Parameters
        ----------
        length : int
            Integer that will be set as `partition_length`.
        """

        self.partition_length = length

    #===========================================================================
    # Parsers
    #===========================================================================

    @staticmethod
    def _get_file_format(partition_file):
        """Guesses the format of the partition file (Nexus or RAxML's).

        Returns
        -------
        partition_format : str
            Format of the partition file ("nexus" or "raxml").
        """

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
        """Parses partitions from file

        This method parses a file containing partitions. It supports
        partitions files similar to RAxML's and NEXUS charset blocks. The
        NEXUS file, however, must only contain the charset block. The
        model_nexus argument provides a namespace for the model variable in
        the nexus format, since this information is not present in the file.
        However, it assures consistency on the Partition object.

        Parameters
        ----------
        partitions_file : str
            Path to partitions file.

        Raises
        ------
        PartitionException
            When one partition definition cannot be parsed.
        """

        # Resets previous partitions (except alignments_range)
        self.reset(keep_alignments_range=True)

        # Get the format of the partition file
        self.partition_format = self._get_file_format(partitions_file)

        part_file = open(partitions_file)

        # In order to support unsorted partition ranges, the complete
        # partition set will be stored temporary in memory. Even very large
        # partition files should result in relatively small data structures.
        # Once this variable is populated, it will be sorted according to the
        # first element of the range.
        temp_ranges = []

        # TODO: Add support for codon partitions in raxml format
        if self.partition_format == "raxml":
            for p, line in enumerate(part_file):

                # Ignore empty lines
                if line.strip() == "":
                    continue

                # A wrongly formatted raxml partition file may be provided, in
                # which case an IndexError exception will be raised. This will
                # handle that exception
                try:
                    fields = line.split(",")
                    # Get partition name as string
                    partition_name = fields[1].split("=")[0].strip()
                    # Get partition range as list of int
                    partition_range_temp = fields[1].split("=")[1]
                    try:
                        partition_range = [
                            int(x) - 1 for x in
                            partition_range_temp.strip().split("-")]

                    except ValueError as e:
                        # A ValueError may be raise when there is a "."
                        # notation in the partition range. If so, convert
                        # the "." to the sequence lenght. If no sequence lenght
                        # has been provided raise another exception
                        pr = partition_range_temp.strip().split("-")
                        if pr[1] == ".":
                            if self.partition_length:
                                partition_range = [int(pr[0]) - 1,
                                                   self.partition_length - 1]
                            else:
                                raise PartitionException(
                                    "The length of the locus must be "
                                    "provided when partitions are "
                                    "defined using '.' notation to "
                                    "mean full length")
                        else:
                            raise e

                    # Check which alignment file contains the current partition
                    if self.alignments_range:
                        try:
                            file_name = \
                                [x for x, y in self.alignments_range.items() if
                                 y[0] <= partition_range[0] < y[1]][0]
                        except IndexError:
                            file_name = None
                    else:
                        file_name = None

                    temp_ranges.append([partition_name, file_name,
                                        partition_range])

                except (IndexError, ValueError):
                    return InvalidPartitionFile(
                        "Badly formatted partitions file in line {} "
                        "with:\n\n{}".format(p + 1, line))

        elif self.partition_format == "nexus":
            for line in part_file:
                # Ignore empty lines
                if line.strip() != "":
                    try:
                        res = self.read_from_nexus_string(line,
                                                          return_res=True)
                    except PartitionException as e:
                        return e
                    if res:
                        temp_ranges.append(res)

        # Sort partition ranges according to the first element of the range
        temp_ranges.sort(key=lambda part: part[2][0])

        for name, file_name, part_range in temp_ranges:
            # Add information to partitions storage
            try:
                self.add_partition(name,
                                   locus_range=part_range,
                                   file_name=file_name)
            except InvalidPartitionFile as e:
                return e

    def read_from_nexus_string(self, nx_string, file_name=None,
                               return_res=False):
        """Parses a single nexus string with partition definition.

        Parameters
        ----------
        nx_string : str
            String with partition definition
        file_name : str, optional
            String with name of the file corresponding to the partition.
        return_res : bool
            If True, it will only return the parsed partition information.
            If False, it will add the parsed partition to the `Partitions`
            object.
        """

        try:
            fields = nx_string.split("=")
            partition_name = fields[0].split()[1].strip()

            # If this list has 2 elements, it should be a simple gene partition
            # If it has 3 elements, it should be a codon partition
            partition_full = re.split(r"[-\\]", fields[1].strip().
                                      replace(";", "").replace("/", "\\"))

            # If partition is defined using "." notation to mean full length
            if partition_full[1] == ".":
                if self.partition_length:
                    partition_range = [int(partition_full[0]) - 1,
                                       self.partition_length - 1]
                else:
                    raise PartitionException("The length of the locus must be "
                                             "provided when partitions are "
                                             "defined using '.' notation to "
                                             "mean full length")
            else:
                partition_range = [int(partition_full[0]) - 1,
                                   int(partition_full[1]) - 1]

            # Check which alignment file contains the current partition
            if self.alignments_range:
                try:
                    file_name = [x for x, y in self.alignments_range.items() if
                                 partition_range[0] in xrange(*y) and
                                 partition_range[1] in xrange(*y)][0]
                except IndexError:
                    pass

            if return_res:
                return [partition_name, file_name, partition_range]
            else:
                self.add_partition(partition_name, locus_range=partition_range,
                                   file_name=file_name)
        # If, for some reason, the current line cannot be interpreted as a
        # charset line, ignore it.
        except IndexError:
            if return_res:
                return None
            else:
                pass

    def get_partition_names(self):
        """Returns a list with the name of the partitions

        Returns
        -------
        names : list
            List with names of the partitions. When a parent
            partition has multiple codon partitions, it returns a partition
            name for every codon starting position present.
        """

        names = []

        for part, vals in self.partitions.items():
            if vals[1]:
                names.extend([part + "_%s" % (x + 1) for x in vals[1]])
            else:
                names.append(part)

        return names

    def read_from_dict(self, dict_obj):
        """Reads partition information from a dict object

        Parses partitions defined and stored in a special OrderedDict. The
        values of dict_obj should be the partition names and their
        corresponding values should contain the loci range and substitution
        model, if any.

        Parameters
        ----------
        dict_obj : OrderedDict
            Ordered dictionary with the definition of the partitions

        Examples
        --------
        Here is an example of a `dict_obj`::

            dict_obj = OrderedDict(("GeneA", [(0,234), "GTR"]),
                                   ("GeneB", [(235, 865), "JC"))
        """

        for k, v in dict_obj:
            # Determining if value contains only the range or the substitution
            # model as well
            if len(v) > 1:
                self.add_partition(k, locus_range=v[0])
            else:
                self.add_partition(k, locus_range=v[0])

    def is_single(self):
        """Returns whether the current `Partitions` has single or multiple
        partitions.

        Returns
        -------
            _ : bool
            Returns True is there is only a single partition defined,
            and False if there are multiple partitions.
        """

        if len(self.partitions) == 1:
            if not [x for x in self.partitions.values()][0][1]:
                return True
            else:
                return False
        else:
            return False

    def _find_parent(self, max_range):
        """Finds the parent partition from a specified range.

        Finds a parent partition of a codon partition.

        Parameters
        ----------
        max_range : int
            The maximum range of the codon partition.

        Returns
        -------
        part : str
            The name of the parent partition, from the `partitions` attribute.
        """

        for part, vals in self.partitions.items():
            lrange = vals[0]
            if lrange[1] == max_range:
                return part

    def add_partition(self, name, length=None, locus_range=None, codon=False,
                      use_counter=False, file_name=None, model_cls=None,
                      auto_correct_name=True):
        """Adds a new partition.

        Adds a new partition providing the length or the range of current
        alignment. If both are provided, the length takes precedence.The range
        of the partition should be in python index, that is, the first position
        should be 0 and not 1.

        Parameters
        ----------
        name : str
            Name of the partition.
        length : int, optional
            Length of the alignment.
        locus_range : list or tuple, optional
            Range of the partition.
        codon : list
            If the codon partitions are already defined, provide the
            starting points in list format, e.g: [1,2,3].
        use_counter : bool
            If True, `locus_range` will be updated according to the `counter`
            attribute.
        file_name : str
            Name of the alignment file.
        model_cls :
            Specified the substitution model that will be set in `models`.
        auto_correct_name : bool
            If set to True, when a partition name already exist, add a counter
            to the end of the name.

        Notes
        -----
        IMPORTANT NOTE on self.model: The self.model attribute was designed
        in a way that allows the storage of different substitution models
        inside the same partition name. This is useful for codon partitions that
        share the same parent partition name. So, for example, a parent
        partition named "PartA" with 3 codon partitions can have a different
        model for each one like this::

            self.models["PartA"] = [[[..model1_params..], [..model2_params..],
                [..model3_params..]], [GTR, GTR, GTR], ["1", "2", "3"]]

        """

        # Check for duplicate names in partitions
        if name in self.partitions:
            if auto_correct_name:
                c = 1
                while "{}_{}".format(name, c) in self.partitions:
                    c += 1
            else:
                raise PartitionException("Partition name %s is already in "
                                         "partition table" % name)

        # When length is provided
        if length:
            # Add to or update alignments_range attribute. This will store the
            # original range of the alignment
            if file_name:
                if file_name in self.alignments_range:
                    current_range = [self.counter, self.counter + (length - 1)]
                    # If start position is earlier than before, update
                    if current_range[0] < self.alignments_range[file_name][0]:
                        self.alignments_range[file_name][0] = current_range[0]
                    # If stop position if later than before, update
                    if current_range[1] > self.alignments_range[file_name][1]:
                        self.alignments_range[file_name][1] = current_range[1]
                else:
                    self.alignments_range[file_name] = [
                        self.counter, self.counter + (length - 1)]

            # Add partition to index list
            self.partitions_index.append([name, 0])
            # Add partition to alignment list
            try:
                self.partitions_alignments[name].append(file_name if file_name
                                                        else name)
            except KeyError:
                self.partitions_alignments[name] = [file_name if file_name else
                                                    name]
            # Create empty model attribute for a single partition
            self.models[name] = [[[]], [None], []]

            self.partitions[name] = [(self.counter,
                                      self.counter + (length - 1)), codon]
            self.counter += length
            self.partition_length += length

        # When a list/tuple range is provided
        elif locus_range:

            if use_counter:
                locus_range = (self.counter,
                               self.counter + locus_range[1] - locus_range[0])

            # Add to or update alignments_range attribute. This will store the
            # original range of the alignment
            if file_name:
                if file_name in self.alignments_range:
                    if locus_range[0] < self.alignments_range[file_name][0]:
                        self.alignments_range[file_name][0] = locus_range[0]
                    if locus_range[1] > self.alignments_range[file_name][1]:
                        self.alignments_range[file_name][1] = locus_range[1]
                else:
                    self.alignments_range[file_name] = list(locus_range)

            # If the maximum range of the current partition is already included
            # in some other partition, and no codon partitions were provided
            # using the "codon" argument, then it should be an undefined codon
            # partition and should be added to an existing partition
            if locus_range[1] <= self.counter and not codon:

                # Find the parent partition
                parent_partition = self._find_parent(locus_range[1])

                if not parent_partition:
                    raise InvalidPartitionFile(
                        "Could not find parent partition of {}. Check the"
                        " ranges of your partitions to ensure no range "
                        "overlaps".format(name))

                # If no codon partition is present in the parent partition,
                # create one
                if not self.partitions[parent_partition][1]:
                    # Add partition to index list
                    self.partitions_index.append([parent_partition, 1])
                    # Create empty model attribute for two partitions
                    self.models[parent_partition] = [[[], []], [None, None], []]

                    parent_start = self.partitions[parent_partition][0][0]
                    self.partitions[parent_partition][1] = [parent_start,
                                                            locus_range[0]]
                else:
                    # Create empty model attribute for additional partitions
                    self.models[parent_partition][0].append([])
                    self.models[parent_partition][1].append(None)

                    # Add partition to index list
                    self.partitions_index.append([parent_partition, 2])

                    self.partitions[parent_partition][1].append(locus_range[0])

            # If the start of the current partition is already within the range
            # of a previous partitions, raise an exception
            elif locus_range[0] < self.counter:
                raise InvalidPartitionFile(
                    "Badly formatted partition with range [{}-{}] starts "
                    "inside the range of a previous partitions ({})".format(
                        locus_range[0], locus_range[1], self.counter))

            # Else, create the new partition. If codon is provided, the codon
            # information is automatically added
            else:
                if model_cls:
                    self.models[name] = model_cls
                else:
                    # Create empty model attribute for a single partition
                    self.models[name] = [[[]], [None], []]
                if codon:
                    self.partitions_index = [[name, x] for x in codon]
                else:
                    # Add partition to index list
                    self.partitions_index.append([name, 0])
                try:
                    self.partitions_alignments[name].append(file_name if
                                                            file_name else name)
                except KeyError:
                    self.partitions_alignments[name] = [file_name if file_name
                                                        else name]

                self.partitions[name] = [(locus_range[0], locus_range[1]),
                                         codon]

                self.counter = locus_range[1] + 1
                self.partition_length = locus_range[1] + 1

    def remove_partition(self, partition_name=None, file_name=None):
        """Removes partitions.

        Removes a partitions by a given partition or file name. This will
        handle any necessary changes on the remaining partitions. The changes
        will be straightforward for most attributes, such as partitions_index,
        partitions_alignments and models, but it will require a re-structuring
        of partitions because the ranges of the subsequent partitions will
        have to be adjusted.

        Parameters
        ----------
        partition_name : str
            Name of the partition.
        file_name : str
            Name of the alignment file.
        """

        def rm_part(nm):
            """
            Remove a partition from self.partitions and update the ranges of
            the remaining partitions
            """

            del self.partitions[nm]

            new_dic = OrderedDict()

            counter = 0
            for nm, vals in self.partitions.items():
                # Check if the starting position of the next partition is the
                # same as the counter. If so, add the vals to the new dict.
                # Else, correct the ranges based on the counter
                if vals[0][0] == counter:
                    new_dic[nm] = vals
                    counter = vals[0][1] + 1
                else:
                    # Get lenght of the partition
                    part_len = vals[0][1] - vals[0][0]
                    # Create corrected range
                    part_range = (counter, counter + part_len)
                    # Correct codon position start if any
                    if vals[1]:
                        codon = [counter, counter + 1, counter + 2]
                    else:
                        codon = False
                    new_dic[nm] = [part_range, codon]
                    counter = counter + part_len + 1

            return new_dic

        def remove_routine(part_name):
            """
            Routine that removes a partition based on its name. It ca be used
            when calling the remove_partition method with the partition_name
            argument, or with the file_name argument when the partition only
            contains that file name
            """
            # Remove partition from partition_index
            self.partitions_index = [x for x in self.partitions_index if x[0] !=
                                     part_name]

            # Remove partitions_alignments
            del self.partitions_alignments[part_name]

            # Remove models
            del self.models[part_name]

            # Remove from partitions
            self.partitions = rm_part(part_name)

        if partition_name:
            # Raise exception if partition name does not exist
            if partition_name not in self.partitions:
                raise PartitionException("%s is not a partition name" %
                                         partition_name)

            remove_routine(partition_name)

        if file_name:

            # Set file_found to True, when there is a match. If no match is
            # found, raise a PartitionException at the end of the loop.
            file_found = False
            for part, file_list in self.partitions_alignments.items():
                if file_name in file_list:
                    file_found = True
                    # If the partitions consists only of the provided file,
                    # Remove the entire partition
                    if len(file_list) == 1:
                        remove_routine(part)
                    # If the partition contains other files, then only remove
                    # the current file from the partition
                    else:
                        self.partitions_alignments[part].remove(file_name)

            if not file_found:
                raise PartitionException("%s file does not belong to any"
                                         "partition" % file_name)

    def change_name(self, old_name, new_name):
        """Changes name of a partition.

        Parameters
        ----------
        old_name : str
            Original partition name.
        new_name : str
            New partition name.
        """

        self.partitions[new_name] = self.partitions.pop(old_name)
        self.partitions_alignments[new_name] = \
            self.partitions_alignments.pop(old_name)
        self.models[new_name] = self.models.pop(old_name)

    def merge_partitions(self, partition_list, name):
        """Merges multiple partitions into a single one.

        Parameters
        ----------
        partition_list : list
            List with partition names to be merged.
        name : str
            Name of new partition
        """

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

        def flatter(s):
            """
            Creates a flat iterator of tuples. If s is [[(1,2), (2,3)], (4,5)]
            this will yield ((1,2), (2,3), (4,5))
            """
            for i in s:
                if isinstance(i, tuple):
                    yield i
                else:
                    for j in i:
                        yield j

        # Get new range
        new_range = [x for x in merger(flatter((y[0] for x, y in
                                               self.partitions.items()
                                               if x in partition_list)))]

        # Add entries for new partition
        self.partitions[name] = [new_range[0] if len(new_range) == 1 else
            new_range, False]
        self.partitions_alignments[name] = list(set([i for x, y in
                                            self.partitions_alignments.items()
                                            if x in partition_list for i in y]))
        self.models[name] = [[[]], [None], []]

        # Delete previous partitions and update merged dict
        for p in partition_list:
            if len(self.partitions_alignments[p]) == 1:
                self.merged_files[self.partitions_alignments[p][0]] = \
                    self.partitions[p][0]
            del self.partitions[p]
            del self.partitions_alignments[p]
            del self.models[p]

    def split_partition(self, name, new_range=None, new_names=None):
        """Splits one partition into two.

        Splits a partitions with `name` into two with the tuple list provided
        by `new_range`. If new_range is None, This will split the partition
        by its alignment files instead.

        Parameters
        ----------
        name : str
            Name of the partition to be split.
        new_range : list or tuple, optional
            List of two tuples, containing the ranges of the new partitions.
        new_names : list, optional
            The names of the new partitions.
        """

        if new_range:

            # Add new partitions
            for n, r in zip(new_names, new_range):
                self.partitions[n] = [r, False]
                # Create new partitions_alignments. Keep the original alignment
                # file for both
                self.partitions_alignments[n] = self.partitions_alignments[name]
                self.models[n] = [[[]], [None], []]

        else:

            for aln in self.partitions_alignments[name]:
                #  Get original range of alignment file
                new_range = self.merged_files[aln]
                # Add new partitions
                aln_name = basename(aln)
                self.partitions[aln_name] = [new_range, False]
                self.partitions_alignments[aln_name] = [aln]
                self.models[aln_name] = [[[]], [None], []]

        # Delete original partition
        del self.partitions[name]
        del self.partitions_alignments[name]
        del self.models[name]

    # ==========================================================================
    # Model handling
    # ==========================================================================

    def parse_nexus_model(self, string):
        """Parses a substitution model defined in a prset and/or lset command.

        Parameters
        ----------
        string : str
            String with the prset or lset command.
        """

        string = string.lower()

        # Find out which partitions the current parameters apply to. If
        # detected, it should be something like "applyto=(1,2)"
        applyto = re.findall(r"applyto=\(.*\)", string)
        # Find parameters
        nst = re.findall(r"nst=[0-9]", string)
        statefreqpr = re.findall(r"statefreqpr=.*\)", string)

        # Collect params
        params = [x[0] for x in [nst, statefreqpr] if x]

        if applyto:
            if applyto == ["applyto=(all)"]:
                for partition in self.partitions:
                    self.models[partition][0] += params
            else:
                # Get target partitions
                part_index = [int(x) for x in
                              re.split("[()]", applyto[0])[1].split(",")]
                for i in part_index:
                    part = self.partitions_index[i - 1]
                    # Get partition name
                    part_name = part[0]
                    # Get subpartition index. 0 if single partition, other if
                    # multiple subpartition
                    part_subpart = part[1]
                    self.models[part_name][0][part_subpart] += params

    def get_model_name(self, params):
        """Given a list of parameters, return the name of the model

        Parameters
        ----------
        p : list
            List of prset/lset parameters

        Returns
        -------
        model : str or None
            Returns the name of the model if it finds. Else, returns None.
        """

        for model, p in self._models["mrbayes"].items():
            if params == p:
                return model
        else:
            return None

    def set_model(self, partition, models, links=None, apply_all=False):
        """Sets substitution model for a given partition.

        Parameters
        ----------
        partition : str
            Partition name.
        models : list
            Model names for each of the three codon partitions. If there
            are no codon partitions, provide only a single element to the list.
        links : list
            Provide potential links between codon models. For
            example, if codon 1 and 2 are to be linked, it should be:
            links=["12", "3"]
        apply_all : bool
            If True, the current model will be applied to all partitions.
        """

        # Get list with partitions to be changed
        if apply_all:
            plist = [x for x in self.partitions]
        else:
            plist = [partition]

        # Replace "No model" string with None
        models = [None if x == "No model" else x for x in models]

        # Set model to the whole partition
        if len(models) == 1:
            # If the current partition was previously defined as having codon
            # partitions, revert it
            for p in plist:
                if self.partitions[p][1]:
                    self.partitions[p][1] = False
                self.models[p][1] = models

        # Set codon models
        else:
            for p in plist:
                # Change the partition in self.partitions to have codon
                # partitions
                st_idx = self.partitions[p][0][0]
                self.partitions[p][1] = [st_idx + x for x in range(3)]
                self.models[p][1] = models
                self.models[p][2] = links

    def write_to_file(self, output_format, output_file, model="LG"):
        """Writes partitions to a file.

        Writes the Partitions object into an output file according to the
        output_format. The supported output formats are RAxML and Nexus. The
        `model` option is for the RAxML format only.

        Parameters
        ----------
        output_format : str
            Output format of partitions file. Can be either "nexus" or
            "raxml".
        output_file : str
            Path to output file.
        model : str
            Name of the model for the partitions. "raxml" format only.
        """

        if output_format == "raxml":
            outfile_handle = open(output_file + ".part.File", "w")
            for part, rge in self.partitions.items():
                partition_range = "-".join([str(x + 1) for x in rge[0]])
                outfile_handle.write("%s, %s = %s\n" % (model,
                                                        part,
                                                        partition_range))

            outfile_handle.close()

        elif output_format == "nexus":
            outfile_handle = open(output_file + ".charset", "w")
            for part, rge in self.partitions.items():
                outfile_handle.write("charset %s = %s;\n" % (
                                     part,
                                     "-".join([str(x + 1) for x in rge[0]])))

            outfile_handle.close()

        return 0


class Zorro(object):
    """
    Class that handles the concatenation of zorro weights.

    Parameters
    ----------
    alignment_list : trifusion.process.sequence.AlignmentList
        AlignmentList object.
    suffix : str
        Suffix of the zorro weight files, based on the corresponding
        input alignments.
    zorro_dir : str
        Path to directory where zorro weight files are stored.
    """

    def __init__(self, alignment_list, suffix="_zorro.out", zorro_dir=None):

        self.weigth_values = []
        self.suffix = suffix

        for file_path in [x.path for x in alignment_list.alignments.values()]:
            # If zorro_dir is provided, use the specified path
            if zorro_dir:
                zorro_file = splitext(file_path.split(sep)[-1])[0]
                zorro_file = "{}{}.txt".format(join(zorro_dir, zorro_file), suffix)
            # If zorro_dir is not provided, use the same path as the input alignment
            else:
                zorro_file = file_path.split(".")[0] + self.suffix + ".txt"
            # alignment file is shared with the corresponding zorro file
            zorro_handle = open(zorro_file)
            self.weigth_values += [int(round(float(weigth.strip()))) for
                                   weigth in zorro_handle]

    def write_to_file(self, output_file):
        """ Creates a concatenated file with the zorro weights for the
        corresponding alignment files."""
        outfile = output_file + "_zorro.out"
        outfile_handle = open(outfile, "w")
        for weigth in self.weigth_values:
            outfile_handle.write("%s\n" % weigth)
        outfile_handle.close()


__author__ = "Diogo N. Silva"
