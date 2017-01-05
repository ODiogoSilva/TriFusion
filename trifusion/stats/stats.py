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

try:
    from process.sequence import Alignment
except ImportError:
    from trifusion.process.sequence import Alignment

from collections import OrderedDict
import numpy as np
import pygal


class UniReport(Alignment):
    """ This will create a report object for single alignments and will
    inherit from the Alignment class. As with its base class, it will be
    possible to create a UniReport object directly from an external file or
    from an ordered dictionary object
    """

    def base_statistics(self):
        """
        returns a list with the basic phylogenetic statistics, such as number
         of indels, overall percentage of indels variable sites,parsimony
         informative sites.
        """

        from collections import Counter

        parsimonious = 0
        variable = 0
        indel = 0

        for column_position in range(self.locus_length):

            column = [char[column_position].lower() for char in
                      self.alignment.values()]

            # Check for variable sites
            column_set = len(set(column) - {"n", self.sequence_code[1]})
            if column_set > 1:
                variable += 1

                # Check for parsimonious informative sites
                # Filters the Counter dictionary with key:val pairs with vals
                # that are not missing data or gap and higher than 1
                if len([val for char, val in Counter(column).items()
                        if val > 1 and char not in "-" and
                        char not in "n"]) > 1:
                    parsimonious += 1

            # Check for indels
            if self.sequence_code[1] in column:
                indel += 1

        cumulative_gap = 0
        for sequence in self.alignment.values():
            cumulative_gap += sequence.count(self.sequence_code[1]) +\
                              sequence.count("-")

        indel_perc = float(cumulative_gap) / (self.locus_length *
                                              len(self.alignment))

        return indel, indel_perc, variable, parsimonious

    def missing_data(self):
        """
        Returns a dictionary with the missing data information for each
        species. Each dictionary entry will be:
        "taxon":[missing_int, gap_int, gap_int+missing_int]
        """

        contents_missing_data = OrderedDict()
        # Creates the base dictionary with taxa names as keys and empty lists
        # as values
        # contents_missing_data = dict((taxa, []) for taxa in self.iter_taxa())

        for taxon, sequence in self.alignment.items():
            gap_value = sequence.count("-")
            missing_value = sequence.count(self.sequence_code[1])
            contents_missing_data[taxon] = [gap_value, missing_value,
                                            gap_value + missing_value]

        return contents_missing_data

    def species_gene_length(self):
        """
        :return: A dictionary with the gene length information for each
        species. The keys will be species names and the values floats with
        the gene length
        """

        species_gene_length = {}

        for taxon, sequence in self.alignment.items():

            strip_sequence = sequence.replace("-", "")
            species_gene_length[taxon] = float(len(strip_sequence))

        return species_gene_length


class MultiReport:
    """ This will create a report object for multiple alignments. Its most
     basic instance will be a list of UniReport objects
    """
    def __init__(self, alignment_list):

        self.report_list = []

        for report in alignment_list:
            report_object = UniReport(report)

            self.report_list.append(report_object)

    def get_species_set(self):
        """
        :return: A list containing the unique taxa names from all alignment
        files
        """

        species_set = set()

        for report in self.report_list:

            current_taxa = report.taxa_list
            if set(current_taxa) != species_set:
                species_set = set(list(species_set) + current_taxa)

        return list(species_set)

    def get_gene_set(self):
        """
        :return: A list containing the name of the gene alignments
        """

        return [report.input_alignment for report in self.report_list]

    def report_table(self, output_file):
        """
        :param output_file: output file name string
        Generates a csv table containing basic phylogenetic statistics for
        each alignment
        """

        output_handle = open(output_file, "w")
        output_handle.write("Gene; Number of sites; Number of indels; "
                            "Percentage of missing data; Variable sites; "
                            "Parsimonious sites\n")

        table_contents = OrderedDict()

        for report in self.report_list:
            name = report.input_alignment
            sites = report.locus_length
            indel, indel_perc, variable, parsimonious = report.base_statistics()

            # In case I'll need this in anything other than a table
            table_contents[name] = [sites, indel, indel_perc, variable,
                                    parsimonious]

        for gene, vals in table_contents.items():

            output_handle.write("%s; %s; %s; %s; %s; %s\n" % (gene, vals[0],
                                                              vals[1], vals[2],
                                                              vals[3], vals[4]))

        output_handle.close()

    def species_gene_length(self, output_file, table=False, plot=False):
        """
        :param output_file: String with the name of the output file
        :return: Creates a table and/or plot with information on the average
        gene length (and corresponding standard
        deviation) for each species
        """

        raw_data = dict((sp, []) for sp in self.get_species_set())

        for report in self.report_list:

            current_data = report.species_gene_length()

            for sp in raw_data:
                if sp not in current_data:
                    raw_data[sp].append(0)
                else:
                    raw_data[sp].append(current_data[sp])

        data = dict((sp, [np.mean(val), np.std(val)]) for sp, val in
                    raw_data.items())

        if table is not False:

            output_handle = open(output_file, "w")
            output_handle.write("Species; Average gene length; SD\n")

            for sp, val in data.items():
                output_handle.write("%s; %s; %s\n" % (sp, val[0], val[1]))

            output_handle.close()

        if plot is not False:

            gene_length_box_chart = pygal.Box(x_label_rotation=45, width=1200,
                                              legend_at_bottom=True, height=800,
                                              label_font_size=8,
                                              legend_font_size=20, margin=50,
                                              major_label_font_size=10)
            gene_length_box_chart.title = "Average gene length per species"

            for sp, vals in raw_data.items():

                gene_length_box_chart.add(sp, vals)

            return gene_length_box_chart

    def gene_variation_plot(self, output_file):
        """ Creates a bar plot with basic information on the variation and
        missing data for each gene. It is similar to the species_missing_data
         method, but the focus in on variation per gene instead of missing
         data per species """

        # UNFINISHED
        plot_contents = []

        for report in self.report_list:

            gene_statistics = report.base_statistics

            #Get absolute values
            parsimoniously_informative = float(gene_statistics[3])
            singletons = float(gene_statistics[2] - gene_statistics[3])
            missing_data = float(gene_statistics[0])
            total_characters = float(report.locus_length)
            gene_data_absolute = [parsimoniously_informative, singletons,
                                  missing_data]

            # Get proportions
            gene_data = [(val / total_characters) for val in gene_data_absolute]

            plot_contents.append((report.input_alignment, gene_data))

        variation_bar_chart = pygal.StackedBar(x_label_rotation=90, width=1200,
                                               legend_at_bottom=True,
                                               height=800, label_font_size=8,
                                               legend_font_size=20, margin=50,
                                                major_label_font_size=10,
                                                print_values=False,
                                                y_title='Proportion')

        variation_bar_chart.title = "Character missing data per species"
        variation_bar_chart.x_labels = [taxon[0] for taxon in
                                        self.get_gene_set()]

    def species_missing_data(self, output_file, table=False, plot=False):
        """
        :param table: Boolean. True will generate a csv table with information
        on the missing data for each species
        :param plot: Boolean. True will generate a stacked bar plot with
        information on the missing data for each species
        Depending on the table and plot argument values, this will generate a
        table and/or a plot with information on the missing data for each
        species
        """

        taxa_set = self.get_species_set()
        missing_data_contents = dict((taxa, [0, 0, 0]) for taxa in taxa_set)
        alignment_length = 0

        for report in self.report_list:
            report_missing_data = report.missing_data()
            alignment_length += report.locus_length

            for taxon, vals in report_missing_data.items():
                missing_data_contents[taxon] = [x + y for x, y in
                                                zip(missing_data_contents[
                                                        taxon], vals)]

        # Convert missing_data_contents to list and sort by missing data
        data_list = [(taxon, [vals[0], vals[1], vals[2]]) for taxon, vals in
                     missing_data_contents.items()]
        sorted_data_list = sorted(data_list, key=lambda x: x[1][2],
                                  reverse=True)

        # Table
        if table is True:

            try:
                output_handle = open(output_file, "w")
            except TypeError:
                print("TypeError: Please specify an output file name")
                raise SystemExit

            output_handle.write("Taxon; Gap chars; Missing chars; Total missing"
                                " chars; Effective chars; Total chars\n")

            for taxon, vals in sorted_data_list:
                output_handle.write("%s; %s; %s; %s; %s; %s\n" % (
                                    taxon,
                                    vals[0],
                                    vals[1],
                                    vals[2],
                                    alignment_length - vals[2],
                                    alignment_length))
            else:
                output_handle.close()

        # Plot
        if plot is True:

            missing_bar_chart = pygal.StackedBar(x_label_rotation=90,
                                                 width=1200,
                                                 legend_at_bottom=True,
                                                 height=800,
                                                label_font_size=8,
                                                legend_font_size=20, margin=50,
                                                major_label_font_size=10,
                                                print_values=False,
                                                y_title='Proportion')

            missing_bar_chart.title = "Character missing data per species"
            missing_bar_chart.x_labels = [taxon[0] for taxon in sorted_data_list]

            gap_list = []
            missing_list = []
            character_list = []

            for taxon, vals in sorted_data_list:
                # Adjusting missing data values to be in proportion
                gap_proportion = float(vals[0]) / float(alignment_length)
                missing_proportion = float(vals[1]) / float(alignment_length)
                character_proportion = (float(alignment_length) -
                                        float(vals[2])) / \
                                        float(alignment_length)

                # Creating chart variables
                gap_list.append({"value": gap_proportion,
                                 "label": str(vals[0])})
                missing_list.append({"value": missing_proportion,
                                     "label": str(vals[1])})
                character_list.append({"value": character_proportion,
                                       "label": str(alignment_length -
                                                    vals[2])})

            # Adding chart variables
            missing_bar_chart.add("Gaps", gap_list)
            missing_bar_chart.add("Missing data", missing_list)
            missing_bar_chart.add("Characters", character_list)

            return missing_bar_chart

__author__ = "Diogo N. Silva"
