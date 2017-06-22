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
    from process.sequence import Alignment
    from base.plotter import bar_plot, multi_bar_plot
    from process.error_handling import KillByUser
except ImportError:
    from trifusion.process.sequence import Alignment
    from trifusion.base.plotter import bar_plot, multi_bar_plot
    from trifusion.process.error_handling import KillByUser

from collections import OrderedDict, Counter
import pickle
import os
import sqlite3
from os.path import join
import random
import string
import copy


class Cluster(object):
    """ Object for clusters of the OrthoMCL groups file. It is useful to set a
     number of attributes that will make subsequent filtration and
     processing much easier """

    def __init__(self, line_string):
        """
        To initialize a Cluster object, only a string compliant with the
        format of a cluster in an OrthoMCL groups file has to be provided.
        This line should contain the name of the group, a colon, and the
        sequences belonging to that group separated by whitespace
        :param line_string: String of a cluster
        """

        # Initializing attributes for parse_string
        self.name = None
        self.sequences = None
        self.species_frequency = {}

        # Initializing attributes for apply filter
        # If the value is different than None, this will inform downstream
        # objects of whether this cluster is compliant with the specified
        # gene_threshold
        self.gene_compliant = None
        # If the value is different than None, this will inform downstream
        # objects of whether this cluster is compliant with the specified
        # species_threshold
        self.species_compliant = None

        self.parse_string(line_string)

    def parse_string(self, cluster_string):
        """
        Parses the string and sets the group name and sequence list attributes
        """

        fields = cluster_string.split(":")
        # Setting the name and sequence list of the clusters
        self.name = fields[0].strip()
        self.sequences = fields[1].strip().split()

        # Setting the gene frequency for each species in the cluster
        self.species_frequency = Counter([field.split("|")[0] for field in
                                          self.sequences])

    def remove_taxa(self, taxa_list):
        """
        Removes the taxa contained in taxa_list from self.sequences and
        self.species_frequency
        :param taxa_list: list, each element should be a taxon name
        """

        self.sequences = [x for x in self.sequences if x.split("|")[0]
                          not in taxa_list]

        self.species_frequency = dict((taxon, val) for taxon, val in
                                      self.species_frequency.items()
                                      if taxon not in taxa_list)

    def apply_filter(self, gene_threshold, species_threshold):
        """
        This method will update two Cluster attributes, self.gene_flag and
        self.species_flag, which will inform downstream objects if this
        cluster respects the gene and species threshold
        :param gene_threshold: Integer for the maximum number of gene copies
        per species
        :param species_threshold: Integer for the minimum number of species
        present
        """

        # Check whether cluster is compliant with species_threshold
        if len(self.species_frequency) >= species_threshold and \
                species_threshold:
            self.species_compliant = True
        else:
            self.species_compliant = False

        # Check whether cluster is compliant with gene_threshold
        if max(self.species_frequency.values()) <= gene_threshold and \
                gene_threshold:
            self.gene_compliant = True
        else:
            self.gene_compliant = False


class OrthoGroupException(Exception):
    pass


class GroupLight(object):
    """
    Analogous to Group object but with several changes to reduce memory usage
    """

    def __init__(self, groups_file, gene_threshold=None,
                 species_threshold=None, ns=None):

        self.gene_threshold = gene_threshold if gene_threshold else None
        self.species_threshold = species_threshold if species_threshold \
            else None

        # Attribute containing the list of included species
        self.species_list = []
        # Attribute that will contain taxa to be excluded from analyses
        self.excluded_taxa = []
        self.species_frequency = []

        # Attributes that will store the number (int) of cluster after gene and
        # species filter
        self.all_clusters = 0
        self.num_gene_compliant = 0
        self.num_species_compliant = 0
        self.all_compliant = 0

        # Attribute containing the total number of sequences
        self.total_seqs = 0
        # Attribute containing the maximum number of extra copies found in the
        # clusters
        self.max_extra_copy = 0

        # Attribute with name of the group file, which will be an ID
        self.name = os.path.abspath(groups_file)
        self.table = groups_file.split(os.sep)[-1].split(".")[0]

        # Initialize atribute containing the groups filtered using the gene and
        # species threshold. This attribute can be updated at any time using
        # the update_filtered_group method
        self.filtered_groups = []

        self._parse_groups(ns)

        if type(self.species_threshold) is float:
            self._get_sp_proportion()

    def groups(self):
        """
        Generator for group file. This replaces the self.groups attribute of
        the original Group Object. Instead of loading the whole file into
        memory, a generator is created to iterate over its contents. It may
        run a bit slower but its a lot more memory efficient.
        :return:
        """

        file_handle = open(self.name)
        for line in file_handle:
            if line.strip() != "":
                yield line.strip()

    def iter_species_frequency(self):
        """
        In order to prevent permanent changes to the species_frequency
        attribute due to the filtering of taxa, this iterable should be used
        instead of the said variable. This creates a temporary deepcopy of
        species_frequency which will be iterated over and eventually modified.
        """

        # Since the items of species_frequency are mutable, this ensures
        # that even those objects are correctly cloned
        sp_freq = copy.deepcopy(self.species_frequency)

        for cl in sp_freq:
            yield cl

    def _remove_tx(self, line):
        """
        Given a group line, remove all references to the excluded taxa
        :param line: raw group file line
        """

        new_line = "{}:".format(line.split(":")[0])

        tx_str = "\t".join([x for x in line.split(":")[1].split() if
                           x.split("|")[0] not in self.excluded_taxa])

        return new_line + tx_str

    def _apply_filter(self, cl):
        """
        Sets or updates the basic group statistics, such as the number of
        orthologs compliant with the gene copy and minimum taxa filters.

        :param cl: dictionary. Contains the number of occurrences for each
         taxon present in the ortholog cluster
         (e.g. {"taxonA": 2, "taxonB": 1).
        """

        # First, remove excluded taxa from the cl object since this will
        # impact all other filters
        for tx in self.excluded_taxa:
            cl.pop(tx, None)

        if cl:

            self.all_clusters += 1

            extra_copies = max(cl.values())
            if extra_copies > self.max_extra_copy:
                self.max_extra_copy = extra_copies

            if extra_copies <= self.gene_threshold and self.gene_threshold and\
                len(cl) >= self.species_threshold and  \
                    self.species_threshold:
                self.num_gene_compliant += 1
                self.num_species_compliant += 1
                self.all_compliant += 1

            elif (extra_copies <= self.gene_threshold and
                      self.gene_threshold) or self.gene_threshold == 0:
                self.num_gene_compliant += 1

            elif len(cl) >= self.species_threshold and \
                    self.species_threshold:
                self.num_species_compliant += 1

    def _get_compliance(self, cl):
        """
        Determines whether an ortholog cluster is compliant with the specified
        ortholog filters.
        :param ccl: dictionary. Contains the number of occurrences for each
         taxon present in the ortholog cluster
         (e.g. {"taxonA": 2, "taxonB": 1).

        :return: tuple. The first element refers to the gene copy filter
        while the second refers to the minimum taxa filter. Values of 1
        indicate that the ortholg cluster is compliant.
        """

        for tx in self.excluded_taxa:
            cl.pop(tx, None)

        if cl:

            cp = max(cl.values())

            if not self.gene_threshold and not self.species_threshold:
                return 1, 1

            if cp <= self.gene_threshold and self.gene_threshold and\
                len(cl) >= self.species_threshold and  \
                    self.species_threshold:
                return 1, 1

            elif (cp <= self.gene_threshold and self.gene_threshold) or \
                    not self.gene_threshold:
                return 1, 0

            elif (len(cl) >= self.species_threshold and
                    self.species_threshold) or not self.species_threshold:
                return 0, 1

            else:
                return 0, 0

    def _reset_counter(self):

        self.all_clusters = 0
        self.num_gene_compliant = 0
        self.num_species_compliant = 0
        self.all_compliant = 0

    def _parse_groups(self, ns=None):

        for cl in self.groups():

            if ns:
                if ns.stop:
                    raise KillByUser("")

            # Retrieve the field containing the ortholog sequences
            sequence_field = cl.split(":")[1]

            # Update species frequency list
            sp_freq = Counter((x.split("|")[0] for x in
                              sequence_field.split()))

            self.species_frequency.append(sp_freq)

            # Update number of sequences
            self.total_seqs += len(sequence_field)

            # Update max number of extra copies
            extra_copies = max(sp_freq.values())
            if extra_copies > self.max_extra_copy:
                self.max_extra_copy = max(sp_freq.values())

            self.species_list.extend([x for x in sp_freq if x not in
                                      self.species_list])

            # Apply filters, if any
            # gene filter
            if self.species_threshold and self.gene_threshold:
                self._apply_filter(sp_freq)

    def exclude_taxa(self, taxa_list, update_stats=False):
        """
        Updates the excluded_taxa attribute and updates group statistics if
        update_stats is True. This does not change the Group object data
        permanently, only sets an attribute that will be taken into account
        when plotting and exporting data.
        :param taxa_list: list. List of taxa that should be excluded from
        downstream operations
        :param update_stats: boolean. If True, it will update the group
        statistics
        """

        # IF the taxa_list is the same as the excluded_taxa attribute,
        # there is nothing to do
        if sorted(taxa_list) == sorted(self.excluded_taxa):
            return

        self.species_list = [x for x in self.species_list + self.excluded_taxa
                             if x not in taxa_list]

        self.excluded_taxa = taxa_list

        if update_stats:
            self._reset_counter()
            for cl in self.iter_species_frequency():
                self._apply_filter(cl)

    def basic_group_statistics(self, update_stats=True):

        if update_stats:
            self._reset_counter()
            for cl in self.iter_species_frequency():
                self._apply_filter(cl)

        return len(self.species_frequency), self.total_seqs, \
            self.num_gene_compliant, self.num_species_compliant, \
            self.all_compliant

    def _get_sp_proportion(self):
        """
        When the species filter is a float value between 0 and 1, convert
        this proportion into absolute values (rounded up), since filters were
        already designed for absolutes.
        """

        self.species_threshold = int(self.species_threshold *
                                     len(self.species_list))

    def update_filters(self, gn_filter, sp_filter, update_stats=False):
        """
        Updates the group filter attributes and group summary stats if
        update_stats is True. This method does not change the
        data of the Group object, only sets attributes that will be taken into
        account when plotting or exporting data
        :param gn_filter: integer. Maximum number of gene copies allowed in an
        ortholog cluster
        :param sp_filter: integer/float. Minimum number/proportion of taxa
        representation
        :param update_stats: boolean. If True it will update the group summary
        statistics
        """

        # If the provided filters are the same as the current group attributes
        # there is nothing to do
        if (gn_filter, sp_filter) == (self.gene_threshold,
                                      self.species_threshold):
            return

        self.gene_threshold = gn_filter
        self.species_threshold = sp_filter

        if type(self.species_threshold) is float:
            self._get_sp_proportion()

        if update_stats:
            self._reset_counter()
            for cl in self.iter_species_frequency():
                self._apply_filter(cl)

    def retrieve_sequences(self, sqldb, protein_db, dest="./",
                             shared_namespace=None, outfile=None):
        """
        :param sqldb: srting. Path to sqlite database file
        :param protein_db: string. Path to protein database file
        :param dest: string. Directory where sequences will be exported
        :param shared_namespace: Namespace object to communicate with
        TriFusion's main process
        :param outfile: If set, all sequeces will be instead saved in a
        single output file. This is used for the nucleotide sequence export
        :return:
        """

        if not os.path.exists(dest) and not outfile:
            os.makedirs(dest)

        if shared_namespace:
            shared_namespace.act = shared_namespace.msg = "Creating database"
            # Stores sequences that could not be retrieved
            shared_namespace.missed = shared_namespace.counter = 0
            shared_namespace.progress = 0
            # Get number of lines of protein database
            p = 0
            with open(protein_db) as fh:
                for p, _ in enumerate(fh):
                    pass
            shared_namespace.max_pb = shared_namespace.total = p + 1

        # Connect to database
        conn = sqlite3.connect(sqldb, timeout=100000)
        c = conn.cursor()
        table_name = "".join([x for x in protein_db if x.isalnum()]).encode(
            "utf8")

        # Create table if it does not exist
        if not c.execute("SELECT name FROM sqlite_master WHERE type='table' "
                         "AND name='{}'".format(table_name)).fetchall():

            c.execute("CREATE TABLE {} (seq_id text PRIMARY KEY, seq text)".
                      format(table_name))

            # Populate database
            with open(protein_db) as ph:
                seq = ""
                for line in ph:

                    # Kill switch
                    if shared_namespace:
                        if shared_namespace.stop:
                            conn.close()
                            raise KillByUser("")
                        shared_namespace.progress += 1
                        shared_namespace.counter += 1

                    if line.startswith(">"):
                        seq_id = line.strip()[1:]
                        if seq != "":
                            c.execute("INSERT INTO {} VALUES (?, ?)".
                                      format(table_name), (seq_id, seq))
                        seq = ""
                    else:
                        seq += line.strip()

            conn.commit()

        if shared_namespace:
            shared_namespace.act = shared_namespace.msg = "Fetching sequences"
            shared_namespace.good = shared_namespace.counter = 0
            shared_namespace.progress = 0
            shared_namespace.max_pb = shared_namespace.total = \
                self.all_compliant

        # Set single output file, if option is set
        if outfile:
            output_handle = open(join(dest, outfile), "w")

        # Fetching sequences
        for line, cl in zip(self.groups(), self.iter_species_frequency()):

            # Kill switch
            if shared_namespace:
                if shared_namespace.stop:
                    conn.close()
                    raise KillByUser("")

            # Filter sequences
            if self._get_compliance(cl) == (1, 1):

                if shared_namespace:
                    shared_namespace.good += 1
                    shared_namespace.progress += 1
                    shared_namespace.counter += 1

                # Retrieve sequences from current cluster
                if self.excluded_taxa:
                    line = self._remove_tx(line)
                fields = line.split(":")

                # Open file
                if not outfile:
                    cl_name = fields[0]
                    output_handle = open(os.path.join(dest, cl_name) + ".fas",
                                         "w")

                seqs = fields[-1].split()
                for i in seqs:
                    # Query database
                    c.execute("SELECT * FROM {} WHERE seq_id = ?".
                              format(table_name), (i,))
                    vals = c.fetchone()
                    # Handles cases where the sequence could not be retrieved
                    # If outfile is set, output_handle will be a single file
                    # for all groups. If not, it will represent an individual
                    # group file
                    try:
                        output_handle.write(">{}\n{}\n".format(vals[0],
                                                               vals[1]))
                    except TypeError:
                        pass

                if not outfile:
                    output_handle.close()

        if outfile:
            output_handle.close()

        conn.close()

    def export_filtered_group(self, output_file_name="filtered_groups",
                                 dest="./", shared_namespace=None):

        if shared_namespace:
            shared_namespace.act = "Exporting filtered orthologs"
            shared_namespace.missed = 0
            shared_namespace.good = 0

        output_handle = open(os.path.join(dest, output_file_name), "w")

        for p, (line, cl) in enumerate(zip(self.groups(),
                                           self.iter_species_frequency())):

            if shared_namespace:
                if shared_namespace.stop:
                    raise KillByUser("")

            if shared_namespace:
                shared_namespace.progress = p

            if self._get_compliance(cl) == (1, 1):
                if shared_namespace:
                    shared_namespace.good += 1
                if self.excluded_taxa:
                    l = self._remove_tx(line)
                else:
                    l = line
                output_handle.write("{}\n".format(l))

        output_handle.close()

    def bar_species_distribution(self, filt=False):

        if filt:
            data = Counter((len(cl) for cl in self.iter_species_frequency() if
                           self._get_compliance(cl) == (1, 1)))
        else:
            data = Counter((len(cl) for cl in self.species_frequency))

        x_labels = [x for x in list(data)]
        data = list(data.values())

        # When data is empty, return an exception
        if not data:
            return {"data": None}

        # Sort lists
        x_labels = [list(x) for x in zip(*sorted(zip(x_labels, data)))][0]

        # Convert label to strings
        x_labels = [str(x) for x in x_labels]

        title = "Taxa frequency distribution"
        ax_names = ["Number of taxa", "Ortholog frequency"]

        return {"data": [data],
                "title": title,
                "ax_names": ax_names,
                "labels": x_labels,
                "table_header": ["Number of species",
                                      "Ortholog frequency"]}

    def bar_genecopy_distribution(self, filt=False):
        """
        Creates a bar plot with the distribution of gene copies across
        clusters
        :param filt: Boolean, whether or not to use the filtered groups.
        """

        if filt:
            data = Counter((max(cl.values()) for cl in
                            self.iter_species_frequency() if
                            self._get_compliance(cl) == (1, 1)))
        else:
            data = Counter((max(cl.values()) for cl in self.species_frequency
                           if cl))

        x_labels = [x for x in list(data)]
        data = list(data.values())

        # When data is empty, return an exception
        if not data:
            return {"data": None}

        x_labels, data = (list(x) for x in zip(*sorted(zip(x_labels, data))))

        # Convert label to strings
        x_labels = [str(x) for x in x_labels]

        title = "Gene copy distribution"
        ax_names = ["Number of gene copies", "Ortholog frequency"]

        return {"data": [data],
                "labels": x_labels,
                "title": title,
                "ax_names": ax_names,
                "table_header": ["Number of gene copies",
                                 "Ortholog frequency"]}

    def bar_species_coverage(self, filt=False):
        """
        Creates a stacked bar plot with the proportion of
        :return:
        """

        data = Counter(dict((x, 0) for x in self.species_list))

        self._reset_counter()

        for cl in self.iter_species_frequency():

            self._apply_filter(cl)
            if filt:
                data += Counter(dict((x, 1) for x, y in cl.items() if y > 0 and
                           self._get_compliance(cl) == (1, 1)))
            else:
                data += Counter(dict((x, 1) for x, y in cl.items() if y > 0))

        data = data.most_common()

        # When data is empty, return an exception
        if not data:
            return {"data": None}

        x_labels = [str(x[0]) for x in data]
        data = [[x[1] for x in data], [self.all_clusters - x[1] if not
                                      filt else self.all_compliant - x[1]
                                      for x in data]]

        lgd_list = ["Available data", "Missing data"]
        ax_names = [None, "Ortholog frequency"]

        return {"data": data,
                "labels": x_labels,
                "lgd_list": lgd_list,
                "ax_names": ax_names}

    def bar_genecopy_per_species(self, filt=False):

        data = Counter(dict((x, 0) for x in self.species_list))

        self._reset_counter()

        for cl in self.iter_species_frequency():

            self._apply_filter(cl)
            if filt:
                data += Counter(dict((x, y) for x, y in cl.items() if y > 1 and
                           self._get_compliance(cl) == (1, 1)))
            else:
                data += Counter(dict((x, y) for x, y in cl.items() if y > 1))

        data = data.most_common()

        # When data is empty, return an exception
        if not data:
            return {"data": None}

        x_labels = [str(x[0]) for x in data]
        data = [[x[1] for x in data]]
        ax_names = [None, "Gene copies"]

        return {"data": data,
                "labels": x_labels,
                "ax_names": ax_names}


class Group(object):
    """ This represents the main object of the orthomcl toolbox module. It is
     initialized with a file name of a orthomcl groups file and provides
     several methods that act on that group file. To process multiple Group
     objects, see MultiGroups object """

    def __init__(self, groups_file, gene_threshold=None,
                 species_threshold=None, project_prefix="MyGroups"):

        # Initializing thresholds. These may be set from the start, or using
        #  some method that uses them as arguments
        self.gene_threshold = gene_threshold
        self.species_threshold = species_threshold

        # Attribute containing the list of included species
        self.species_list = []
        # Attribute that will contain taxa to be excluded from analyses
        self.excluded_taxa = []

        # Attributes that will store the number (int) of cluster after gene and
        # species filter
        self.all_compliant = 0
        self.num_gene_compliant = 0
        self.num_species_compliant = 0

        # Attribute containing the total number of sequences
        self.total_seqs = 0
        # Attribute containing the maximum number of extra copies found in the
        # clusters
        self.max_extra_copy = 0

        # Attribute with name of the group file, which will be an ID
        self.group_name = groups_file
        # Initialize the project prefix for possible output files
        self.prefix = project_prefix
        # Initialize attribute containing the original groups
        self.groups = []
        # Initialize atribute containing the groups filtered using the gene and
        # species threshold. This attribute can be updated at any time using
        # the update_filtered_group method
        self.filtered_groups = []
        self.name = None
        # Parse groups file and populate groups attribute
        self.__parse_groups(groups_file)

    def __parse_groups(self, groups_file):
        """
        Parses the ortholog clusters in the groups file and populates the
         self.groups list with Cluster objects for each line in the groups file.
        :param groups_file: File name for the orthomcl groups file
        :return: populates the groups attribute
        """

        self.name = groups_file
        self.species_list = []
        groups_file_handle = open(groups_file)

        for line in groups_file_handle:
            cluster_object = Cluster(line)

            # Add cluster to general group list
            self.groups.append(cluster_object)

            # Update total sequence counter
            self.total_seqs += len(cluster_object.sequences)

            # Update maximum number of extra copies, if needed
            if max(cluster_object.species_frequency.values()) > \
                    self.max_extra_copy:
                self.max_extra_copy = \
                    max(cluster_object.species_frequency.values())

            # Update species_list attribute
            self.species_list = list(set(self.species_list).union(
                set(cluster_object.species_frequency.keys())))

            # If thresholds have been specified, update self.filtered_groups
            # attribute
            if self.species_threshold and self.gene_threshold:
                cluster_object.apply_filter(self.gene_threshold,
                                            self.species_threshold)
                if cluster_object.species_compliant and \
                        cluster_object.gene_compliant:
                    # Add cluster to the filtered group list
                    self.filtered_groups.append(cluster_object)
                    self.all_compliant += 1

                # Update num_species_compliant attribute
                if cluster_object.species_compliant:
                    self.num_species_compliant += 1
                # Update num_gene_compliant attribute
                if cluster_object.gene_compliant:
                    self.num_gene_compliant += 1

    def exclude_taxa(self, taxa_list):
        """
        Adds a taxon_name to the excluded_taxa list and updates the
        filtered_groups list
        """

        self.excluded_taxa.extend(taxa_list)

        # Storage variable for new filtered groups
        filtered_groups = []

        # Reset max_extra_copy attribute
        self.max_extra_copy = 0

        for cl in self.groups:
            cl.remove_taxa(taxa_list)
            if cl.iter_sequences and cl.species_frequency:
                filtered_groups.append(cl)

                # Update maximum number of extra copies, if needed
                if max(cl.species_frequency.values()) > self.max_extra_copy:
                    self.max_extra_copy = max(cl.species_frequency.values())

        # Update species_list
        self.species_list = sorted(list(set(self.species_list) -
                                        set(taxa_list)))

        self.filtered_groups = self.groups = filtered_groups

    def get_filters(self):
        """
        Returns a tuple with the thresholds for max gene copies and min species
        """

        return self.gene_threshold, self.species_threshold

    def basic_group_statistics(self):
        """
        This method creates a basic table in list format containing basic
        information of the groups file (total number of clusters, total number
        of sequences, number of clusters below the gene threshold, number of
        clusters below the species threshold and number of clusters below the
        gene AND species threshold)
        :return: List containing number of

        [total clusters,
         total sequences,
         clusters above gene threshold,
         clusters above species threshold,
         clusters above gene and species threshold]
        """

        # Total number of clusters
        total_cluster_num = len(self.groups)

        # Total number of sequenes
        total_sequence_num = self.total_seqs

        # Gene compliant clusters
        clusters_gene_threshold = self.num_gene_compliant

        # Species compliant clusters
        clusters_species_threshold = self.num_species_compliant

        clusters_all_threshold = len(self.filtered_groups)

        statistics = [total_cluster_num, total_sequence_num,
                      clusters_gene_threshold, clusters_species_threshold,
                      clusters_all_threshold]

        return statistics

    def paralog_per_species_statistic(self, output_file_name=
                                      "Paralog_per_species.csv", filt=True):
        """
        This method creates a CSV table with information on the number of
        paralog clusters per species
        :param output_file_name: string. Name of the output csv file
        :param filt: Boolean. Whether to use the filtered groups (True) or
        total groups (False)
        """

        # Setting which clusters to use
        if filt:
            groups = self.filtered_groups
        else:
            groups = self.groups

        paralog_count = dict((species, 0) for species in self.species_list)

        for cluster in groups:
            for species in paralog_count:
                if cluster.species_frequency[species] > 1:
                    paralog_count[species] += 1

        # Writing table
        output_handle = open(output_file_name, "w")
        output_handle.write("Species; Clusters with paralogs\n")

        for species, val in paralog_count.items():
            output_handle.write("%s; %s\n" % (species, val))

        output_handle.close()

    def export_filtered_group(self, output_file_name="filtered_groups",
                              dest="./", get_stats=False,
                              shared_namespace=None):
        """
        Export the filtered groups into a new file.
        :param output_file_name: string, name of the filtered groups file
        :param dest: string, path to directory where the filtered groups file
        will be created
        :param get_stats: Boolean, whether to return the basic count stats or
        not
        :param shared_namespace: Namespace object, for communicating with
        main process.
        """

        if self.filtered_groups:

            if shared_namespace:
                shared_namespace.act = "Exporting filtered orthologs"

            output_handle = open(os.path.join(dest, output_file_name), "w")

            if get_stats:
                all_orthologs = len(self.groups)
                sp_compliant = 0
                gene_compliant = 0
                final_orthologs = 0

            for cluster in self.filtered_groups:

                if shared_namespace:
                    shared_namespace.progress = \
                        self.filtered_groups.index(cluster)

                if cluster.species_compliant and cluster.gene_compliant:
                    output_handle.write("%s: %s\n" % (
                                    cluster.name, " ".join(cluster.iter_sequences)))
                    if get_stats:
                        final_orthologs += 1
                if get_stats:
                    if cluster.species_compliant:
                        sp_compliant += 1
                    if cluster.gene_compliant:
                        gene_compliant += 1

            output_handle.close()

            if get_stats:
                return all_orthologs, sp_compliant, gene_compliant,\
                           final_orthologs

        else:
            raise OrthoGroupException("The groups object must be filtered "
                                      "before using the export_filtered_group"
                                      "method")

    def update_filters(self, gn_filter, sp_filter):
        """
        Sets new values for the self.species_threshold and self.gene_threshold
        and updates the filtered_group
        :param gn_filter: int. Maximum value for gene copies in cluster
        :param sp_filter:  int. Minimum value for species in cluster
        """

        self.species_threshold = int(sp_filter)
        self.gene_threshold = int(gn_filter)

        self.update_filtered_group()

    def update_filtered_group(self):
        """
        This method creates a new filtered group variable, like
        export_filtered_group, but instead of writing into a new file, it
        replaces the self.filtered_groups variable
        """

        self.filtered_groups = []

        # Reset gene and species compliant counters
        self.num_gene_compliant = 0
        self.num_species_compliant = 0

        for cluster in self.groups:
            cluster.apply_filter(self.gene_threshold, self.species_threshold)
            if cluster.species_compliant and cluster.gene_compliant:
                self.filtered_groups.append(cluster)

            # Update num_species_compliant attribute
            if cluster.species_compliant:
                self.num_species_compliant += 1
            # Update num_gene_compliant attribute
            if cluster.gene_compliant:
                self.num_gene_compliant += 1

    def retrieve_sequences(self, database, dest="./", mode="fasta",
                             filt=True, shared_namespace=None):
        """
        When provided with a database in Fasta format, this will use the
        Alignment object to retrieve sequences
        :param database: String. Fasta file
        :param dest: directory where files will be save
        :param mode: string, whether to retrieve sequences to a file ('fasta'),
        or a dictionary ('dict')
        :param filt: Boolean. Whether to use the filtered groups (True) or
        total groups (False)
        :param shared_namespace: Namespace object. This argument is meant for
        when fast are retrieved in a background process, where there is a need
        to update the main process of the changes in this method
        :param dest: string. Path to directory where the retrieved sequences
        will be created.
        """

        if mode == "dict":
            seq_storage = {}

        if filt:
            groups = self.filtered_groups
        else:
            groups = self.groups

        if not os.path.exists("Orthologs"):
            os.makedirs("Orthologs")

        # Update method progress
        if shared_namespace:
            shared_namespace.act = "Creating database"
            shared_namespace.progress = 0

        print("Creating db")
        # Check what type of database was provided
        #TODO: Add exception handling if file is not parsed with Aligment
        if isinstance(database, str):
            try:
                db_aln = pickle.load(open(database, "rb"))
            except (EnvironmentError, pickle.UnpicklingError):
                db_aln = Alignment(database)
                db_aln = db_aln.alignment
        elif isinstance(database, dict):
            db_aln = database
        else:
            raise OrthoGroupException("The input database is neither a string"
                                      "nor a dictionary object")

        print("Retrieving seqs")
        # Update method progress
        if shared_namespace:
            shared_namespace.act = "Retrieving sequences"
        for cluster in groups:

            if shared_namespace:
                shared_namespace.progress += 1

            if mode == "dict":
                seq_storage[cluster.name] = []

            output_handle = open(join(dest, cluster.name + ".fas"), "w")
            for sequence_id in cluster.iter_sequences:
                seq = db_aln[sequence_id]
                if mode == "fasta":
                    output_handle.write(">%s\n%s\n" % (sequence_id, seq))
                elif mode == "dict":
                    seq_storage[cluster.name].append([sequence_id.split("|")[0],
                                                      seq])
            output_handle.close()

        if mode == "dict":
            return seq_storage

    def bar_species_distribution(self, dest="./", filt=False, ns=None,
                                 output_file_name="Species_distribution"):
        """
        Creates a bar plot with the distribution of species numbers across
        clusters
        :param dest: string, destination directory
        :param filt: Boolean, whether or not to use the filtered groups.
        :param output_file_name: string, name of the output file
        """

        data = []

        # Determine which groups to use
        if filt:
            groups = self.filtered_groups
        else:
            groups = self.groups

        for i in groups:
            if ns:
                if ns.stop:
                    raise KillByUser("")
            data.append(len([x for x, y in i.species_frequency.items()
                             if y > 0]))

        # Transform data into histogram-like
        transform_data = Counter(data)
        x_labels = [x for x in list(transform_data)]
        y_vals = list(transform_data.values())

        # Sort lists
        x_labels, y_vals = (list(x) for x in zip(*sorted(zip(x_labels,
                                                             y_vals))))
        # Convert label to strings
        x_labels = [str(x) for x in x_labels]

        if ns:
            if ns.stop:
                raise KillByUser("")

        # Create plot
        b_plt, lgd, _ = bar_plot([y_vals], x_labels,
                                 title="Taxa frequency distribution",
                                 ax_names=["Number of taxa", "Ortholog frequency"])
        b_plt.savefig(os.path.join(dest, output_file_name), bbox_inches="tight",
                      dpi=400)

        # Create table
        table_list = [["Number of species", "Ortholog frequency"]]
        for x, y in zip(x_labels, y_vals):
            table_list.append([x, y])

        return b_plt, lgd, table_list

    def bar_genecopy_distribution(self, dest="./", filt=False,
                                output_file_name="Gene_copy_distribution.png"):
        """
        Creates a bar plot with the distribution of gene copies across
        clusters
        :param dest: string, destination directory
        :param filt: Boolean, whether or not to use the filtered groups.
        :param output_file_name: string, name of the output file
        """

        data = []

        # Determin which groups to use
        if filt:
            groups = self.filtered_groups
        else:
            groups = self.groups

        for cl in groups:
            # Get max number of copies
            max_copies = max(cl.species_frequency.values())

            data.append(max_copies)

        # Transform data into histogram-like
        transform_data = Counter(data)
        x_labels = [x for x in list(transform_data)]
        y_vals = list(transform_data.values())

        # Sort lists
        x_labels, y_vals = (list(x) for x in zip(*sorted(zip(x_labels,
                                                             y_vals))))
        # Convert label to strings
        x_labels = [str(x) for x in x_labels]

        # Create plot
        b_plt, lgd, _ = bar_plot([y_vals], x_labels,
                    title="Gene copy distribution",
                    ax_names=["Number of gene copies", "Ortholog frequency"],
                    reverse_x=False)
        b_plt.savefig(os.path.join(dest, output_file_name), bbox_inches="tight",
                      figsize=(8 * len(x_labels) / 4, 6), dpi=200)

        # Create table
        table_list = [["Number of gene copies", "Ortholog frequency"]]
        for x, y in zip(x_labels, y_vals):
            table_list.append([x, y])

        return b_plt, lgd, table_list

    def bar_species_coverage(self, dest="./", filt=False, ns=None,
                             output_file_name="Species_coverage"):
        """
        Creates a stacked bar plot with the proportion of
        :return:
        """

        # Determine which groups to use
        if filt:
            groups = self.filtered_groups
        else:
            groups = self.groups

        data = Counter(dict((x, 0) for x in self.species_list))

        for cl in groups:
            if ns:
                if ns.stop:
                    raise KillByUser("")
            data += Counter(dict((x, 1) for x, y in cl.species_frequency.items()
                            if y > 0))

        xlabels = [str(x) for x in list(data.keys())]
        data = [list(data.values()), [len(groups) - x for x in
                                      data.values()]]

        lgd_list = ["Available data", "Missing data"]

        if ns:
            if ns.stop:
                raise KillByUser("")

        b_plt, lgd, _ = bar_plot(data, xlabels, lgd_list=lgd_list,
                              ax_names=[None, "Ortholog frequency"])
        b_plt.savefig(os.path.join(dest, output_file_name), bbox_inches="tight",
                      dpi=200)

        return b_plt, lgd, ""


class MultiGroups(object):
    """ Creates an object composed of multiple Group objects """

    def __init__(self, groups_files=None, gene_threshold=None,
                 species_threshold=None, project_prefix="MyGroups"):
        """
        :param groups_files: A list containing the file names of the multiple
        group files
        :return: Populates the self.multiple_groups attribute
        """

        # If a MultiGroups is initialized with duplicate Group objects, these
        # will be stored in a list. If all Group objects are unique, the list
        # will remain empty
        self.duplicate_groups = []

        # Initializing thresholds. These may be set from the start, or using
        # some method that uses them as arguments
        self.gene_threshold = gene_threshold
        self.species_threshold = species_threshold

        self.prefix = project_prefix

        self.multiple_groups = {}
        self.filters = {}

        if groups_files:
            for group_file in groups_files:

                # If group_file is already a Group object, just add it
                if not isinstance(group_file, Group):
                    # Check for duplicate group files
                    group_object = Group(group_file, self.gene_threshold,
                                         self.species_threshold)
                else:
                    group_object = group_file

                if group_object.name in self.multiple_groups:
                    self.duplicate_groups.append(group_object.name)
                else:
                    self.multiple_groups[group_object.name] = group_object
                    self.filters[group_object.name] = (1,
                                                len(group_object.species_list))

    def __iter__(self):

        return iter(self.multiple_groups)

    def iter_gnames(self):

        return (x.name for x in self.multiple_groups)

    def get_gnames(self):

        return [x.name for x in self.multiple_groups]

    def add_group(self, group_obj):
        """
        Adds a group object
        :param group_obj: Group object
        """

        # Check for duplicate groups
        if group_obj.name in self.multiple_groups:
            self.duplicate_groups.append(group_obj.name)
        else:
            self.multiple_groups[group_obj.name] = group_obj

    def remove_group(self, group_id):
        """
        Removes a group object according to its name
        :param group_id: string, name matching a Group object name attribute
        """

        if group_id in self.multiple_groups:
            del self.multiple_groups[group_id]

    def get_group(self, group_id):
        """
        Returns a group object based on its name. If the name does not match
        any group object, returns None
        :param group_id: string. Name of group object
        """

        try:
            return self.multiple_groups[group_id]
        except KeyError:
            return

    def add_multigroups(self, multigroup_obj):
        """
        Merges a MultiGroup object
        :param multigroup_obj: MultiGroup object
        """

        for group_obj in multigroup_obj:
            self.add_group(group_obj)

    def update_filters(self, gn_filter, sp_filter, group_names=None,
                       default=False):
        """
        This will not change the Group object themselves, only the filter
        mapping. The filter is only applied when the Group object is retrieved
        to reduce computations
        :param gn_filter: int, filter for max gene copies
        :param sp_filter: int, filter for min species
        :param group_names: list, with names of group objects
        """

        if group_names:
            for group_name in group_names:
                # Get group object
                group_obj = self.multiple_groups[group_name]
                # Define filters
                gn_filter = gn_filter if not default else 1
                sp_filter = sp_filter if not default else \
                    len(group_obj.species_list)
                # Update Group object with new filters
                group_obj.update_filters(gn_filter, sp_filter)
                # Update filter map
                self.filters[group_name] = (gn_filter, sp_filter)

            for group_name, group_obj in self.multiple_groups.items():
                # Define filters
                gn_filter = gn_filter if not default else 1
                sp_filter = sp_filter if not default else \
                    len(group_obj.species_list)
                # Update Group object with new filters
                group_obj.update_filters(gn_filter, sp_filter)
                # Update filter map
                self.filters[group_name] = (gn_filter, sp_filter)

    def basic_multigroup_statistics(self, output_file_name=
                                    "multigroup_base_statistics.csv"):
        """
        :param output_file_name:
        :return:
        """

        # Creates the storage for the statistics of the several files
        statistics_storage = OrderedDict()

        for group in self.multiple_groups:
            group_statistics = group.basic_group_statistics()
            statistics_storage[group.name] = group_statistics

        output_handle = open(self.prefix + "." + output_file_name, "w")
        output_handle.write("Group file; Total clusters; Total sequences; "
                            "Clusters below gene threshold; Clusters above "
                            "species threshold; Clusters below gene and above"
                            " species thresholds\n")

        for group, vals in statistics_storage.items():
            output_handle.write("%s; %s\n" % (group, ";".join([str(x) for x
                                                               in vals])))

        output_handle.close()

    def bar_orthologs(self, output_file_name="Final_orthologs",
                             dest="./", stats="total"):
        """
        Creates a bar plot with the final ortholog values for each group file
        :param output_file_name: string. Name of output file
        :param dest: string. output directory
        :param stats: string. The statistics that should be used to generate
        the bar plot. Options are:
            ..: "1": Total orthologs
            ..: "2": Species compliant orthologs
            ..: "3": Gene compliant orthologs
            ..: "4": Final orthologs
            ..: "all": All of the above
            Multiple combinations can be provided, for instance: "123" will
            display bars for total, species compliant and gene compliant stats
        """

        # Stores the x-axis labels
        x_labels = []
        # Stores final ortholog values for all 4 possible data sets
        vals = [[], [], [], []]
        lgd = ["Total orthologs", "After species filter", "After gene filter",
               "Final orthologs"]

        # Get final ortholog values
        for g_obj in self.multiple_groups:

            x_labels.append(g_obj.name.split(os.sep)[-1])
            # Populate total orthologs
            if "1" in stats or stats == "all":
                vals[0].append(len(g_obj.groups))
            # Populate species compliant orthologs
            if "2" in stats or stats == "all":
                vals[1].append(g_obj.num_species_compliant)
            # Populate gene compliant orthologs
            if "3" in stats or stats == "all":
                vals[2].append(g_obj.num_gene_compliant)
            # Populate final orthologs
            if "4" in stats or stats == "all":
                vals[3].append(len(g_obj.filtered_groups))

        # Filter valid data sets
        lgd_list = [x for x in lgd if vals[lgd.index(x)]]
        vals = [l for l in vals if l]

        # Create plot
        b_plt, lgd = multi_bar_plot(vals, x_labels, lgd_list=lgd_list)
        b_plt.savefig(os.path.join(dest, output_file_name),
                      bbox_extra_artists=(lgd,), bbox_inches="tight")

        # Create table list object
        table_list = []
        # Create header
        table_list.append([""] + x_labels)
        # Create content
        for i in range(len(vals)):
            table_list += [x for x in [[lgd_list[i]] + vals[i]]]

        return b_plt, lgd, table_list

    def group_overlap(self):
        """
        This will find the overlap of orthologs between two group files.
        THIS METHOD IS TEMPORARY AND EXPERIMENTAL
        """

        def parse_groups(group_obj):
            """
            Returns a list with the sorted ortholog clusters
            """

            storage = []

            for cluster in group_obj.groups:
                storage.append(set(cluster.iter_sequences))

            return storage

        if len(self.multiple_groups) != 2:
            raise SystemExit("This method can only be used with two group "
                             "files")

        group1 = self.multiple_groups[0]
        group2 = self.multiple_groups[1]

        group1_list = parse_groups(group1)
        group2_list = parse_groups(group2)

        counter = 0
        for i in group1_list:
            if i in group2_list:
                counter += 1

        print(counter)


class MultiGroupsLight(object):
    """
    Creates an object composed of multiple Group objects like MultiGroups.
    However, instead of storing the groups in memory, these are shelved in
    the disk
    """

    # The report calls available
    calls = ['bar_genecopy_distribution',
             'bar_species_distribution',
             'bar_species_coverage',
             'bar_genecopy_per_species']

    def __init__(self, db_path, groups=None, gene_threshold=None,
                 species_threshold=None, project_prefix="MyGroups",
                 ns=None):
        """
        :param groups: A list containing the file names of the multiple
        group files
        :return: Populates the self.multiple_groups attribute
        """

        self.db_path = db_path

        # If a MultiGroups is initialized with duplicate Group objects, their
        # names will be stored in a list. If all Group objects are unique, the
        # list will remain empty
        self.duplicate_groups = []

        self.groups = {}

        self.groups_stats = {}

        # Attribute that will store the paths of badly formated group files
        self.bad_groups = []

        # Initializing thresholds. These may be set from the start, or using
        # some method that uses them as arguments
        self.gene_threshold = gene_threshold
        self.species_threshold = species_threshold

        # Initializing mapping of group filters to their names. Should be
        # something like {"groupA": (1, 10)}
        self.filters = {}

        self.taxa_list = {}
        self.excluded_taxa = {}

        # This attribute will contain a dictionary with the maximum extra copies
        # for each group object
        self.max_extra_copy = {}
        # This attribute will contain a list with the number of species for
        # each group object, excluding replicates. If a MultiGroupLight object
        # contains Group objects with different taxa numbers, this attribute
        # can be used to issue a warning
        self.species_number = []

        self.prefix = project_prefix

        if ns:
            ns.files = len(groups)

        if groups:
            for group_file in groups:
                # If group_file is already a Group object, just add it
                if not isinstance(group_file, GroupLight):
                    try:
                        if ns:
                            if ns.stop:
                                raise KillByUser("")
                            ns.counter += 1
                        group_object = GroupLight(group_file,
                                                  self.gene_threshold,
                                                  self.species_threshold,
                                                  ns=ns)
                    except Exception as e:
                        print(e.message)
                        self.bad_groups.append(group_file)
                        continue
                else:
                    group_object = group_file

                # Check for duplicate group files
                if group_object.name in self.groups:
                    self.duplicate_groups.append(group_file.name)
                else:
                    self.add_group(group_object)

    def __iter__(self):
        for k, val in self.groups.items():
            yield k, pickle.load(open(val, "rb"))

    def clear_groups(self):
        """
        Clears the current MultiGroupsLight object
        """

        for f in self.groups.values():
            os.remove(f)

        self.duplicate_groups = []
        self.groups = {}
        self.groups_stats = {}
        self.filters = {}
        self.max_extra_copy = {}
        self.species_number = []
        self.gene_threshold = self.species_threshold = 0

    def add_group(self, group_obj):
        """
        Adds a group object
        :param group_obj: Group object
        """

        # Check for duplicate groups
        if group_obj.name not in self.groups:
            gpath = os.path.join(self.db_path,
                    "".join(random.choice(string.ascii_uppercase) for _ in
                            range(15)))
            pickle.dump(group_obj, open(gpath, "wb"))
            self.groups[group_obj.name] = gpath
            self.filters[group_obj.name] = (1, len(group_obj.species_list), [])
            self.max_extra_copy[group_obj.name] = group_obj.max_extra_copy
            if len(group_obj.species_list) not in self.species_number:
                self.species_number.append(len(group_obj.species_list))
        else:
            self.duplicate_groups.append(group_obj.name)

    def remove_group(self, group_id):
        """
        Removes a group object according to its name
        :param group_id: string, name matching a Group object name attribute
        """

        if group_id in self.groups:
            os.remove(self.groups[group_id])
            del self.groups[group_id]

    def get_group(self, group_id):
        """
        Returns a group object based on its name. If the name does not match
        any group object, returns None
        :param group_id: string. Name of group object
        """

        try:
            return pickle.load(open(self.groups[group_id], "rb"))
        except KeyError:
            return

    def add_multigroups(self, multigroup_obj):
        """
        Merges a MultiGroup object
        :param multigroup_obj: MultiGroup object
        """

        for _, group_obj in multigroup_obj:
            self.add_group(group_obj)

    def update_filters(self, gn_filter, sp_filter, excluded_taxa,
                       group_names=None, default=False):
        """
        This will not change the Group object themselves, only the filter
        mapping. The filter is only applied when the Group object is retrieved
        to reduce computations

        :param gn_filter: int, filter for max gene copies
        :param sp_filter: int, filter for min species
        :param group_names: list, with names of group objects
        """

        # There are no groups to update
        if group_names == []:
            return

        if group_names:
            glist = group_names
        else:
            glist = self.groups

        for group_name in glist:
            # Get group object
            group_obj = pickle.load(open(self.groups[group_name], "rb"))

            # Define excluded taxa
            group_obj.exclude_taxa(excluded_taxa, True)

            # Define filters
            gn_filter = gn_filter if not default else 1
            sp_filter = sp_filter if not default else \
                len(group_obj.species_list)

            # Correct maximum filter values after excluding taxa
            gn_filter = gn_filter if gn_filter <= group_obj.max_extra_copy \
                else group_obj.max_extra_copy
            sp_filter = sp_filter if sp_filter <= len(group_obj.species_list) \
                else len(group_obj.species_list)

            # Update Group object with new filters
            group_obj.update_filters(gn_filter, sp_filter)
            # Update group stats

            self.get_multigroup_statistics(group_obj)
            pickle.dump(group_obj, open(self.groups[group_name], "wb"))
            # Update filter map

            self.filters[group_name] = (gn_filter, group_obj.species_threshold)
            self.taxa_list[group_name] = group_obj.species_list
            self.excluded_taxa[group_name] = group_obj.excluded_taxa

    def get_multigroup_statistics(self, group_obj):
        """
        :return:
        """

        stats = group_obj.basic_group_statistics()

        self.groups_stats[group_obj.name] = {"stats": stats,
                                        "species": group_obj.species_list,
                                        "max_copies": group_obj.max_extra_copy}

    def bar_orthologs(self, group_names=None, output_file_name="Final_orthologs",
                             dest="./", stats="all"):
        """
        Creates a bar plot with the final ortholog values for each group file
        :param group_names: list. If None, all groups in self.group_stats will
        be used to generate the plot. Else, only the groups with the names in
        the list will be plotted.
        :param output_file_name: string. Name of output file
        :param dest: string. output directory
        :param stats: string. The statistics that should be used to generate
        the bar plot. Options are:
            ..: "1": Total orthologs
            ..: "2": Species compliant orthologs
            ..: "3": Gene compliant orthologs
            ..: "4": Final orthologs
            ..: "all": All of the above
            Multiple combinations can be provided, for instance: "123" will
            display bars for total, species compliant and gene compliant stats
        """

        # Stores the x-axis labels
        x_labels = []
        # Stores final ortholog values for all 4 possible data sets
        vals = [[], [], [], []]
        lgd = ["Total orthologs", "After species filter", "After gene filter",
               "Final orthologs"]

        # Determine which groups will be plotted
        if group_names:
            groups_lst = group_names
        else:
            groups_lst = self.groups_stats.keys()

        for gname in groups_lst:

            gstats = self.groups_stats[gname]

            x_labels.append(gname.split(os.sep)[-1])
            # Populate total orthologs
            if "1" in stats or stats == "all":
                vals[0].append(gstats["stats"][0])
            # Populate species compliant orthologs
            if "2" in stats or stats == "all":
                vals[1].append(gstats["stats"][3])
            # Populate gene compliant orthologs
            if "3" in stats or stats == "all":
                vals[2].append(gstats["stats"][2])
            # Populate final orthologs
            if "4" in stats or stats == "all":
                vals[3].append(gstats["stats"][4])

        # Filter valid data sets
        lgd_list = [x for x in lgd if vals[lgd.index(x)]]
        vals = [l for l in vals if l]

        # Create plot
        b_plt, lgd = multi_bar_plot(vals, x_labels, lgd_list=lgd_list)
        b_plt.savefig(os.path.join(dest, output_file_name),
                      bbox_extra_artists=(lgd,), bbox_inches="tight", dpi=200)

        # Create table list object
        table_list = []
        # Create header
        table_list.append([""] + x_labels)
        # Create content
        for i in range(len(vals)):
            table_list += [x for x in [[lgd_list[i]] + vals[i]]]

        return b_plt, lgd, table_list


__author__ = "Diogo N. Silva"
