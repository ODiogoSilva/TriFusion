#!/usr/bin/env python2
# -*- coding: utf-8 -*-
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


################################################################################

"""
This module contains the text content for the informative popups that can be
accessed by clicking info (?) buttons throughout the app. The text content
will be stored in a single dictionary variable, whose format should be the
following:
   ..: key: string id of the informative content
   ..: value: list containing the actual informative content
      ..: 1st element: string with the title
      ..: 2nd element: string with text body

Both list elements will be kivy Label objects that support markdown syntax
"""

informative_storage = {
    # Orthology
    "orthology_filters":
        ["[b]Orthology filters - Help[/b]",
         "These filters are applied at the end of the Orthology Search "
         "operation for each"
         " specified inflation value "
         "and determine which ortholog clusters are selected according to "
         "their number of gene copies and taxa. They can be later modified in"
         " the Orthology Explore tab to evaluate their effect on the final "
         "number of ortholog clusters."
         "\n\n- [b]Maximum number of genes copies[/b]: "
         "Sets the maximum number of gene copies in each ortholog cluster. To "
         "remove any upper limit on the number of gene copies, set this value "
         "to 0.\n\n- [b]Minimum number of taxa[/b]: Sets the minimum number of "
         "taxa that must be present in an ortholog cluster. Values can be "
         "either absolute or in proportion."],
    "orthology_output":
        ["[b]Orthology output directory - Help[/b]",
         "Sets the directory where all orthology search results will be "
         "stored. The directory structure will follow these general rules:"
         "\n\n- [b]'backstage_files'[/b]: Intermediate files that are "
         "necessary "
         "during the execution of orthology search. Can be useful for error "
         "checking. Also contains the protein data base file used for the "
         "all-vs-all search (goodProteins_db, by default)\n\n- ["
         "b]'Orthology_results'[/b]: "
         "Contains the final results of the orthology search, including the "
         "final group files. For each inflation parameter that was selected, "
         "a directory will be created with the results for the corresponding "
         "parameter value and containing the sequence files of the ortholog "
         "clusters after applying the orthology filters."],
    "database_name":
        ["[b]Database name - Help[/b]",
         "This option will set the name of the protein sequence database "
         "that will be used to performed all-vs-all searches with USEARCH. "
         "This database can also be subsequently used to convert OrthoMCL "
         "group files into protein sequences, by fetching the sequence "
         "references of the group clusters from the database file."],
    "inflation":
        ["[b]Inflation - Help[/b]",
         "Inflation is a parameter of the MCL algorithm that regulates "
         "the tightness of the ortholog clusters, with conceivable values "
         "ranging from 1.1 to 10, but rarely useful at the extremes. Using "
         "higher values (higher tightness) results in an increase of ortholog"
         " clusters but reduces the number of sequences included in each "
         "cluster. Conversely, lower values result in the inclusion of more "
         "sequences in fewer groups. According to the [u][ref=http://www.ncbi"
         ".nlm.nih.gov/pmc/articles/PMC403725/]OrthoMCL paper[/ref][/u], "
         "coarsed-grained clustering (1.5 value) provides sufficient "
         "tightness for identifying coherent ortholog clusters and appears to "
         "strike a balance between sensitivity and sensibility. However, "
         "they also note that this parameter may have little impact on the "
         "final number of clusters. A sensible approach would be to test "
         "several values and evaluate the sensibility of the data set to this "
         "parameter."],
    # Process
    "revert_concatenation":
        ["[b]Revert concatenation - Help[/b]",
         "The [b]revert concatenation[/b] option allows for the conversion of a"
         " single file with defined partitions into individual files for each "
         "partition. Partitions can be specified in two ways:\n\n- A "
         "partition file may be imported using the 'Select partition file' in "
         "the revert concatenation settings. Supported formats for the "
         "partition file include the Nexus charset block and a RAxML "
         "partition file;\n\n- Using defined partitions in "
         "the 'Partitions' tab in the side panel. Note that some input file "
         "formats, such as Nexus, may already contain defined partitions.\n\n"
         " This option takes precedence over all other Process' options, "
         "which means that all secondary operations (Collapse, Consensus, "
         "Filter, etc.) are applied to the individual files that result from "
         "this option."],
    "process_dataset":
        ["[b]Data set selection - Help[/b]",
         "The [b]Data set[/b] option offers a convenient way to perform "
         "operations on specific sets of taxa and/or files. User specified "
         "groups can be selected by toggling the taxon or file buttons in 'Menu"
         " > Open/View data' and selecting the 'Active taxa/files' option or "
         "by creating taxa groups in 'Menu > Dataset Groups' and selecting the "
         "name of the corresponding group."],
    "zorro":
        ["[b]ZORRO weights - Help[/b]",
         "(Concatenation only)\n\nThis option allows the user to"
         " concatenate alignment weight files that are generated by the "
         "ZORRO software [u][ref=http://journals.plos.org/plosone/article?id=10"
         ".1371/journal.pone.0030288](Wu M et al. 2012, PLoS ONE 7(1): "
         "e30288)[/ref][/u].\n\nIn order to maintain correspondence between "
         "the alignment and weight files, the weight file names should begin "
         "with the name of the corresponding alignment file and end with a "
         "suffix that must be provided. For example, if the files:\n\n "
         "-alignment_a.fas\n -alignment_b.fas\n\nare to be concatenated, "
         "then the corresponding zorro weight files should be:\n\n "
         "-alignment_a_zorro.txt\n -alignment_b_zorro.txt\n\nwhere '_zorro' "
         "is the suffix."],
    "collapse":
        ["[b]Collapse - Help",
         "The [b]Collapse[/b] option will merge all identical sequences into "
         "single haplotypes, so that the final alignment  contains only "
         "unique sequences. As a consequence, the taxa names will be "
         "replaced by haplotype names, which can be partially "
         "determined by the  user using the [b]Haplotype prefix[/b] option, "
         "and the correspondence between haplotype names and taxa names with "
         "the same sequence will be written into an *.haplotypes file."],
    "filter":
        ["[b]Gap/Missing Filter - Help",
         "The [b]Gap/Missing data[/b] option will perform two main tasks:\n\n "
         "[b]1. Within alignments:\n[/b]"
         "- Replacement of the gap symbols by missing data symbols at the "
         "beginning and end of the alignments created by the alignment "
         "software;\n"
         "- Removal of alignment columns containing a proportion of gaps "
         "and/or missing data higher than the thresholds specified using "
         "the [b]Gap threshold[/b] and [b]Missing data threshold[/b] sliders."
         "The [b]gap threshold[/b] refers only to gap characters found in the "
         "middle of the alignment. The [b]missing data threshold[/b] refers "
         "to the combination of gap and missing data characters found in any "
         "column of the alignment.\n\n"
         "[b]2. Among multiple alignments:[/b]\n"
         "- Alignments can be filtered based on the proportion of taxa "
         "represented. The total number of taxa is gathered from the full "
         "data set. If, for example, the minimum taxa representation value is "
         "set to 50%, only alignments with at least 50% of the full data "
         "set's taxa will be processed."],
    "taxa_filter":
        ["Taxa Filter - Help",
         "The [b]Taxa Filter[/b] option allows alignments to be filtered "
         "depending on whether they contain or exclude a given set of taxa. "
         "When choosing the filter mode  '[b]Contain[/b]', only those "
         "alignments containing at least the specified taxa in the taxa "
         "group will be saved. When choosing the '[b]Exclude[/b]' filter "
         "mode, alignments that contain [i]all[/i] the "
         "taxa specified in the taxa group will [i]not[/i] be saved."],
    "codon_filter":
        ["Codon Filter - Help",
         "The [b]Codon filter[/b] option is only available for nucleotide "
         "sequences, where it allows the filtering of codon positions in "
         "the alignments. By default, all positions are set to be saved but "
         "any combination can be specified to be filtered out."],
    "variation_filter":
        ["Variation Filter - Help",
         "The [b]variation filter[/b] option allows alignments to be "
         "filtered according their amount or type of sequence variation.\n\n["
         "b]Variable sites[/b] refer to any alignment column that contains at "
         "least one variante.\n\n[b]Informative sites[/b] refer to variable "
         "sites whose minimum frequency variant is present in at least "
         "in two taxa"],
    "consensus":
        ["Alignment consensus - Help",
         "The alignment consensus option will merge all sequences in an "
         "alignment into a single consensus sequence. Variation in alignment "
         "columns can be either converted into the corresponding IUPAC symbols,"
         " soft masked or removed from the final consensus"],
    "consensus_mode":
        ["Consensus variation handling - Help",
         "The consensus mode option determines how variation in an alignment "
         "will be handled to create a consensus sequence. The options are:"
         "\n\n[b]IUPAC:[/b] Converts variable sites using the IUPAC ambiguity"
         " code.\n\n[b]Soft mask:[/b] Replaces variable sites with missing "
         "data symbol.\n\n[b]Remove:[/b] Removes variable sites.\n\n"
         "[b]First sequence:[/b] Usesthe first sequence in each alignment"],
    "gcoder":
        ["[b]Gap coding - Help[/b]",
         "(Nexus output format only)\n\n"
         "This option allows for the codification of gap patterns as a binary "
         "matrix that is appended at the end of the sequences. The "
         "codification of the gaps follows the method selected in the '[b]Gap "
         "coding method[/b]' option. "],

    # Statistics
    "info_sequence_size":
        ["[b]Distribution of sequence size - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Per species; Average["
         "/b]\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Plots the "
         "distribution of average sequence size (excluding sites with gaps "
         "or missing data)\n\n"
         "[b]y-axis:[/b] Sequence size in nucleotides (DNA) or residues ("
         "Protein)\n"
         "[b]x-axis:[/b] Taxon name\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Plots the distribution of "
         "the average sequence size (excluding sites with gaps or missing "
         "data) for the entire active data set.\n\n"
         "[b]y-axis:[/b] Absolute frequency of alignments\n"
         "[b]x-axis:[/b] Sequence size in nucleotides (DNA) or residues ("
         "Protein)"],
    "info_nucleotide_proportion":
        ["[b]Nucleotide or residue proportion - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Per species; Average["
         "/b]\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Stacked bar "
         "plot with the proportions of each nucleotide (DNA) or residue ("
         "Protein) for each species in the active data set. A sorted color "
         "legend is provided above the plot.\n\n"
         "[b]y-axis:[/b] Proportion of nucleotide or residue\n"
         "[b]x-axis:[/b] Taxon name\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Bar plot "
         "with the averaged proportions for each nucleotide or residue across "
         "the active data set.\n\n"
         "[b]y-axis:[/b] Proportion\n"
         "[b]x-axis:[/b] Nucleotide or residue"],
    "info_taxa_frequency":
        ["[b]Distribution of taxa frequency - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Average[/b]\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Histogram with "
         "the frequency of alignments containing a given number of taxa. This "
         "plot can be useful in large data sets to evaluate the number of "
         "genes containing many or few taxa. Dense data sets (few missing "
         "taxa) tend to have distributions shifted to the right side, "
         "while sparse data sets tend to have distributions shifted to the "
         "left side.\n\n"
         "[b]y-axis:[/b] Absolute frequency of alignments\n"
         "[b]x-axis:[/b] Number of taxa present"],
    "pairwise_seq_similarity":
        ["[b]Pairwise sequence similarity - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Single gene; "
         "Per species; Average[/b]\n\n"
         "[color=37abc8ff][b]Single gene[/color][/b]: Sliding "
         "window line plot with the average sequence similarity across the "
         "length of the selected alignment.\n\n"
         "[b]y-axis:[/b] Sequence similarity, in percentage\n"
         "[b]x-axis:[/b] Alignment nucleotide/residue position\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Creates a triangular "
         "heat map matrix with the sequence similarity values for each "
         "pair-wise combination of taxa. Sequence similarity values are color"
         " coded according to the color bar at the right of the plot.\n\n"
         "[b]y-axis:[/b] Taxa names\n"
         "[b]x-axis:[/b] Taxa names\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Plots an histogram with the"
         " distribution of sequence similarity values across the "
         "active data set. Each entry in the histogram corresponds to "
         "average sequence similarity of one alignment.\n\n"
         "[b]y-axis[/b] Absolute frequency of alignments\n"
         "[b]x-axis[/b] Sequence similarity, in percentage"],
    "segregating_sites":
        ["[b]Segregating sites - Help",
         "[color=37abc8ff][b]Available options:[/color] Single gene; "
         "Per species; Average[/b]\n\n"
         "[color=37abc8ff][b]Single gene[/color][/b]: Sliding "
         "window line plot with the number of segregating sites across the "
         "length of the selected alignment.\n\n"
         "[b]y-axis:[/b] Number of segregating sites\n"
         "[b]x-axis:[/b] Alignment nucleotide/residue position\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Triangular "
         "heat map matrix with the number of segregating sites for each "
         "pair-wise combination of taxa. The number of segregating sites are "
         "color coded according to the color bar at the right of the plot.\n\n"
         "[b]y-axis:[/b] Taxa names\n"
         "[b]x-axis:[/b] Taxa names\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Plots an histogram with the"
         " distribution of the number of segregating sites, in absolute "
         "values or percentage, across the "
         "active data set. Each entry in the histogram corresponds to "
         "averaged number of segregating sites in one alignment.\n\n"
         "[b]y-axis:[/b] Absolute frequency of alignments\n"
         "[b]x-axis:[/b] Number of segregating sites, in absolute or "
         "percentage"],
    "len_pol_correlation":
        ["[b]Alignment length/polymorphism correlation[/b]",
         "[color=37abc8ff][b]Available options:[/color] Average[/b]\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Calculates the "
         "non-parametric spearman's rank correlation coefficient between the "
         "alignment length and number of informative sites for each "
         "alignment across the active data set. A "
         "scatter plot depicting the relationship between these two variables"
         " is generated along with the best fit line. The correlation "
         "coeficient and p-value are provided in the top right area of the "
         "plot.\n\n"
         "[b]x-axis:[/b] Number of informative sites\n"
         "[b]y-axis:[/b] Alignment length"],
    "afs":
        ["[b]Allele frequency spectrum - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Single gene, Average["
         "/b]\n\n"
         "[b](Nucleotide sequence only)[/b]\n\n"
         "[color=37abc8ff][b]Single gene[/color][/b]: Plots the "
         "distribution of allele frequencies for a single gene.\n\n"
         "[b]y-axis:[/b] Frequency of occurrence\n"
         "[b]x-axis:[/b] Derived allele frequency\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Creates the distribution of "
         "allele frequencies across the active data set. For each alignment, "
         "the distribution of the derived allele frequency is calculated and "
         "compiled into a single data matrix that depicts the distribution of "
         "allele frequencies for the entire data set.\n\n"
         "[b]y-axis:[/b] Frequency of occurrence\n"
         "[b]x-axis:[/b] Derived allele frequency\n"],
    "gene_occupancy":
        ["[b]Gene occupancy - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Average[/b]\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Plots an interpolation "
         "matrix that scores the presence/absence of each taxon for every "
         "alignment in the active data set. Each taxa is depicted as a line "
         "in the y-axis and each alignment as a column in the x-axis. "
         "Whenever a taxon is absent in a given alignment, the coordinates of "
         "this taxon-alignment combination receive a score of 0 and a blank "
         "space is represented. Conversely, if the taxa in present in a given "
         "alignment, the corresponding coordinate will receive a score of 1 "
         "and a black space is represent. Sparse data sets with taxa missing "
         "in several alignments will tend to have more white spaces. "
         "Conversely, dense data sets with very few missing taxa will have "
         "fewer white spaces.\n\n"
         "[b]y-axis:[/b] Taxa index, each line represents one taxon\n"
         "[b]x-axis:[/b] Alignment index, each column represents an alignment"],
    "missing_taxa":
        ["[b]Distribution of missing taxa - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Per species, Average["
         "/b]\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Plots the distribution "
         "of the number of genes missing for each taxon.\n\n"
         "[b]y-axis:[/b] Frequency, number of genes missing for a given taxon\n"
         "[b]x-axis:[/b] Taxa names.\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Plots the distribution of "
         "the number of taxa missing for each alignment across the entire "
         "active data set.\n\n"
         "[b]y-axis:[/b] Frequency, number of occurences\n"
         "[b]x-axis:[/b] Number of taxa missing"],
    "missing_data":
        ["[b]Distribution of missinf data - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Per species, Average["
         "/b]\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Stacked bar "
         "plot with the proportions of gaps, missing data and effective data "
         "for each taxon averaged over the active alignments. In cases "
         "where a given taxon is missing from an alignment, this will score "
         "100% of missing data for that alignment.\n\n"
         "[b]y-axis:[/b] Proportion\n"
         "[b]x-axis:[/b] Taxa names\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Smooth line "
         "histogram for three data categories: Gaps, Missing data "
         "and effective Data. For each alignment in the active data set, "
         "the proportion for each of these categories is calculated and then "
         "compiled into three data matrices. Each plot will then contain "
         "the distribution of the corresponding category across the entire "
         "active data set. As in the 'Per species' plot, taxa that are "
         "missing from a given alignment will score 100% of missing data for "
         "that alignment. Data sets "
         "with little missing data will tend to have 'Data' distributions "
         "abutting 100% values and 'Gaps' and 'Missing data' distributions "
         "abutting 0%. Increasing the amount of missing data will erase this "
         "pattern and tend to homogenise all three distribution or even to "
         "invert the trend.\n\n"
         "[b]y-axis: [/b]Frequency, number of genes\n"
         "[b]x-axis: [/b]Proportion, for each of the three categories"],
    "cumulative_missing":
        ["[b]Cumulative distribution of missing taxa - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Average[/b]\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Plots the number of genes "
         "available when applying a consecutive minimum taxa representation "
         "filters. This range spans from 0% (All values of missing data "
         "allowed) up to 100% (Alignments with no missing data) with steps "
         "of 5%. This plot can be useful when trying to decide the "
         "the best value for the minimum taxa representation filter that "
         "minimizes missing data and maximized effective data.\n\n"
         "[b]y-axis:[/b] Frequency, number of alignments\n"
         "[b]x-axis:[/b] Minimum taxa representation threshold, in percentage"],
    "outliers_missing":
        ["[b]Missing data outliers - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Per Species, Average["
         "/b]\n\n"
         "[b]All outlier detection methods employed in TriFusion use a "
         "slightly modified Median Absolute Deviation (MAD) method to find "
         "outliers. As the name implies, this method relies on the "
         "calculation of the median for the data set and then evaluate the "
         "distance from each data point to this value. The recommended "
         "threshold for this distance (above which a data point is considered "
         "an outlier) is 3.5, which is roughly equivalent to 3.5 standard "
         "deviations. The slight modification of the method is the inclusion "
         "of a consistency constant, which ensures that MAD provides a good "
         "estimate of the standard deviation independently of the sample "
         "size.[/b]\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Smooth line "
         "distribution of the amount of missing data for each taxon across the "
         "active data set. Data points will be based on the proportion of "
         "missing data symbols out of the possible total. For example, in an "
         "alignment with three taxa, each with 100 sites, the total "
         "possible missing data is 300 (100 * 3). Taxa that are absent from "
         "a given alignment are ignored, in order to avoid biasing outlier "
         "detection to missing taxa instead of missing data. Data points "
         "represent individual taxa and are depicted in a rug-like style near "
         "the x-axis line. Outliers are marked with red dots.\n\n"
         "[b]y-axis:[/b] Frequency, number of occurrences\n"
         "[b]x-axis:[/b] Proportion of missing data\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Smooth line "
         "distribution of the average missing data for each alignment across "
         "the active data set. Data points represent individual alignments "
         "and are depicted in a rug-like style near the x-axis line. Outliers "
         "are marked with red dots.\n\n"
         "[b]y-axis:[/b] Frequency, number of occurrences\n"
         "[b]x-axis:[/b] Proportion of missing data"],
    "outliers_segregating":
        ["[b]Segregating sites outliers - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Per Species, Average["
         "/b]\n\n"
         "[b]All outlier detection methods employed in TriFusion use a "
         "slightly modified Median Absolute Deviation (MAD) method to find "
         "outliers. As the name implies, this method relies on the "
         "calculation of the median for the data set and then evaluate the "
         "distance from each data point to this value. The recommended "
         "threshold for this distance (above which a data point is considered "
         "an outlier) is 3.5, which is roughly equivalent to 3.5 standard "
         "deviations. The slight modification of the method is the inclusion "
         "of a consistency constant, which ensures that MAD provides a good "
         "estimate of the standard deviation independently of the sample "
         "size.[/b]\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Smooth line "
         "distribution of the average pair-wise proportion of segregating "
         "sites of a single taxon when compared with every other taxon. "
         "Taxa outliers will be "
         "identified when a given taxon has an average excess or lack of "
         "segregating sites with every other taxon. Alignment columns "
         "containing gaps or missing data are ignored. Data points represent"
         " individual taxa and are depicted in a rug-like style near the "
         "x-axis line. Outliers are marked with red dots.\n\n"
         "[b]y-axis:[/b] Frequency, number of occurrences\n"
         "[b]x-axis:[/b] Proportion of segregating sites\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Smooth line "
         "distribution of the proportion of segregating sites for each "
         "alignment across the active data set. Alignment outliers will be"
         " identified when an alignment has an excess of lack of segregating "
         "sites compared to the remaining alignments. Alignment columns "
         "containing gaps or missing data are ignored. Data points represent "
         "individual alignments and are depicted in a rug-like style near "
         "the x-axis line. Outliers are marked with red dots.\n\n"
         "[b]y-axis:[/b] Frequency, number of occurrences\n"
         "[b]x-axis:[/b] Proportion of segregating sites"],
    "outliers_size":
        ["[b]Sequence size outliers - Help[/b]",
         "[color=37abc8ff][b]Available options:[/color] Per Species, Average["
         "/b]\n\n"
         "[b]All outlier detection methods employed in TriFusion use a "
         "slightly modified Median Absolute Deviation (MAD) method to find "
         "outliers. As the name implies, this method relies on the "
         "calculation of the median for the data set and then evaluate the "
         "distance from each data point to this value. The recommended "
         "threshold for this distance (above which a data point is considered "
         "an outlier) is 3.5, which is roughly equivalent to 3.5 standard "
         "deviations. The slight modification of the method is the inclusion "
         "of a consistency constant, which ensures that MAD provides a good "
         "estimate of the standard deviation independently of the sample "
         "size.[/b]\n\n"
         "[color=37abc8ff][b]Per species[/color][/b]: Smooth line "
         "distribution of the average sequence size for each taxon. The "
         "sequence size for each taxon excludes gaps and missing data. Taxa "
         "outliers will be identified when a given taxon has sequences much "
         "larger or smaller than the median of the active data set. Data "
         "points represent"
         " individual taxa and are depicted in a rug-like style near the "
         "x-axis line. Outliers are marked with red dots.\n\n"
         "[b]y-axis:[/b] Frequency, number of occurrences\n"
         "[b]x-axis:[/b] Sequence size in nucleotides/residues\n\n"
         "[color=37abc8ff][b]Average[/color][/b]: Smooth line "
         "distribution of the average sequence size for each alignment. The "
         "sequence size for each alignment excludes gaps and missing data. "
         "Outlier alignments will be identified when the average sequence "
         "size of an alignment is much larger or smaller than the median of "
         "the active data set. Data points represent "
         "individual alignments and are depicted in a rug-like style near "
         "the x-axis line. Outliers are marked with red dots.\n\n"
         "[b]y-axis:[/b] Frequency, number of occurrences\n"
         "[b]x-axis:[/b] Sequence size in nucleotides/residues"]
}

orthology_storage = {
    "Taxa distribution": "Distribution of the number of taxa across clusters",
    "Taxa coverage": "Shows available and missing ortholog data per taxon",
    "Gene copy distribution": "Distribution of gene copy numbers across "
        "clusters",
    "Taxa gene copies": "Cumulative number of multiple gene copies per taxon"
}

# Informative text for Orthology plots. Information for each plot is provided
#  as a tuple of three elements: (title, category, figure name, description)
orthology_plots = [
    ("Distribution of maximum gene copies",
     "Gene focused",
     "bar_genecopy_distribution.png",
     "Displays the number of orthologs in function of the maximum number of "
     "gene copies. Here, the number of gene copies in each ortholog cluster is "
     "obtained from the taxa with the highest number of gene copies. "
     "Therefore, if an ortholog cluster with 10 taxa has a single copy "
     "for 9 taxa and 10 copies for 1 taxon, that ortholog cluster scores a "
     "value of 10. "
     "Orthologs with a single gene copy effectively represent single "
     "copy genes. We recommend creating this plot without applying a maximum "
     "gene copy filter in order to get the full distribution of the data set."
     "<br/><br/>"
     "<b>x-axis:</b> Number of maximum gene copies<br/>"
     "<b>y-axis:</b> Frequency of ortholog groups"),
    ("Species distribution",
     "Species focused",
     "bar_species_distribution.png",
     "Displays the number of orthologs in function of the exact number of "
     "unique taxa represented in a given ortholog. Here, taxa represented by "
     "multiple copies are scored only once but the range of the number of "
     "taxa is determined by the minimum taxa representation filter. In the "
     "absence of a minimum taxa representation filter, it is entirely "
     "possible to have groups with a single taxon, which represent"
     "clusters of recent paralogs for that taxon.<br/><br/>"
     "<b>x-axis:</b> Number of taxa<br/>"
     "<b>y-axis:</b> Frequency of ortholog groups"),
    ("Data coverage per species",
     "Species focused",
     "bar_species_coverage.png",
     "Displays the amount of available and missing ortholog groups for each "
     "taxon. The number of available genes is displayed as dark blue bars, "
     "while the number of missing genes is displayed as light blue bars. Given "
     "the fact that we are dealing with non-aligned data, this will "
     "simply score the number of genes available for a given taxa in relation "
     "to the total number of ortholog cluster for the current filters. "
     "(Therefore, these amounts can also be viewed as proportions.) This "
     "analysis will NOT take into account the number of gaps that will be "
     "generated after the alignment procedure. The horizontal dashed line "
     "represents the mean available data for the data set.<br/><br/>"
     "<b>x-axis:</b> Taxon<br/>"
     "<b>y-axis:</b> Frequency of ortholog groups"),
    ("Distribution of gene copy per taxa",
     "Species focused",
     "bar_genecopy_per_species.png",
     "Displays the number of gene copies for each taxon. Here, only ortholog "
     "clusters with multiple copies for a given taxa are taken into account. "
     "If a taxon has only single copies in the entire data set, it scores 0."
     "<br/><br/>"
     "<b>x-axis:</b> Taxon<br/>"
     "<b>y-axis:</b> Number of gene copies (above single copy)")
]

__author__ = "Diogo N. Silva"
