__author__ = 'diogo'

from process import data
from process.error_handling import *
from process.sequence import AlignmentList, Alignment
from ortho import orthomcl_pipeline as ortho_pipe
from ortho import OrthomclToolbox as OrthoTool

from os.path import join, basename
from os import sep
from collections import OrderedDict
from copy import deepcopy
import logging
import cPickle as pickle


def load_proc(aln_list, file_list, nm, dest):
    try:
        if aln_list:
            aln_list.add_alignment_files(file_list, dest=dest,
                                         shared_namespace=nm)
            aln_obj = aln_list
        else:
            aln_obj = AlignmentList(file_list, dest=dest, shared_namespace=nm)

        # To speed sharing of the new AlignmentList object with the main
        # process, the class instance is being saved in disk using pickle,
        # and latter loaded into the add in the load_files method
        with open(join(dest, "alns.pc"), "wb") as fh:
            pickle.dump(aln_obj, fh)

    except:
        logging.exception("Unexpected error when loading input data")
        nm.exception = True


def background_process(f, ns, a):
    """
    Executes the func in the background and returns its value to the
    shared namespace
    :param f: bound method
    :param ns: shared namespace
    :param a: list arguments
    """
    try:
        if a:
            val = f(*a)
        else:
            val = f()
        ns.val = val
    except Exception:
        logging.exception("Unexpected exit in {}".format(f.__name__))
        ns.exception = True


def background_export_groups(f, nm, a):
    """
    Calls function f with a arguments
    :param f: bound method
    :param nm: shared namespace object to provide to bound method
    :param a: list, with arguments
    """
    try:
        f(*a, shared_namespace=nm)
    except:
        logging.exception("Unexpected error when exporting ortholog "
                          "groups")
        nm.exception = True


def orto_execution(nm, temp_dir, proteome_files, protein_min_len,
                   protein_max_stop, cur_dir, usearch_evalue,
                   usearch_threads, usearch_output, mcl_inflation,
                   ortholog_prefix, group_prefix, orto_max_gene,
                   orto_min_sp, sqldb, ortho_dir):
    """
    Executes all pipeline subprocesses sequentially and updates the
    Progess dialog label
    """

    try:
        if nm.k:
            nm.t = "Installing schema"
            nm.c = 1
            ortho_pipe.install_schema(temp_dir)
        if nm.k:
            nm.t = "Adjusting Fasta Files"
            nm.c = 2
            ortho_pipe.adjust_fasta(proteome_files)
        if nm.k:
            nm.t = "Filtering Fasta Files"
            nm.c = 3
            ortho_pipe.filter_fasta(protein_min_len, protein_max_stop,
                                    bin_path=join(cur_dir, "ortho",
                                                  "orthomclFilterFasta"))
        if nm.k:
            nm.t = "Running USearch. This may take a while..."
            nm.c = 4
            ortho_pipe.allvsall_usearch("goodProteins.fasta", usearch_evalue,
                                        usearch_threads, usearch_output,
                                        usearch_bin=join(cur_dir, "data",
                                                         "resources",
                                                         "usearch_linux"))
        if nm.k:
            nm.t = "Parsing USEARCH output"
            nm.c = 5
            ortho_pipe.blast_parser(usearch_output,
                                    db_dir=temp_dir)
        if nm.k:
            nm.t = "Obtaining Pairs"
            nm.c = 6
            ortho_pipe.pairs(temp_dir)
        if nm.k:
            ortho_pipe.dump_pairs(temp_dir)
        if nm.k:
            nm.t = "Running MCL"
            nm.c = 7
            ortho_pipe.mcl(mcl_inflation)
        if nm.k:
            nm.t = "Dumping groups"
            nm.c = 8
            ortho_pipe.mcl_groups(mcl_inflation, ortholog_prefix, "1000",
                                  group_prefix, bin_path=join(cur_dir, "ortho",
                                                "orthomclMclToGroups"))
        if nm.k:
            nm.t = "Filtering group files"
            nm.c = 9
            stats, groups_obj = ortho_pipe.export_filtered_groups(mcl_inflation,
                                                                  group_prefix,
                                                                  orto_max_gene,
                                                                  orto_min_sp,
                                                                  sqldb,
                                                                join(ortho_dir,
                                                         "backstage_files",
                                                         "goodProteins.fasta"),
                                                                  temp_dir)
            # stats is a dictionary containing the inflation value as
            #  key and a list with the orthologs as value
            nm.stats = stats
            nm.groups = groups_obj
    except:
        logging.exception("Unexpected exit in Orthology search")
        nm.exception = True


def update_active_fileset(aln_obj, set_name, file_list, file_groups):
    """
    This method is similar in purpose and functionality to the
    update_active_taxaset, but it updates the set of files. It should be
    used before the said method.
    :param aln_obj: The alignment object being used during execution of
    Process operations
    """
    # Determine the selected active taxa set from the dropdown menu
    if set_name == "All files":
        aln_obj.update_active_alignments([basename(x) for x in file_list])
        return aln_obj
    if set_name == "Active files":
        return aln_obj
    else:
        aln_obj.update_active_alignments(file_groups[set_name])
        return aln_obj


def update_active_taxaset(aln_obj, set_name, active_taxa_list, taxa_groups):
    """
    Do not use this method on the original self.alignment_list or
    self.active_alignment list, as it may cause unwanted permanent changes
    to the taxa set.
    This will take the complete taxa set from self.alignment_list.taxa_names
    and the currently active taxa set from self.active_taxa_list and remove
    the all taxa that are not present in the active taxa set from the
    alignment object passed as argument. If the lists are the same, no taxa
    will be removed
    """

    if set_name == "All taxa":
        return aln_obj
    if set_name == "Active taxa":
        tx_set = active_taxa_list
    else:
        tx_set = taxa_groups[set_name]

    # Remove taxa
    aln_obj.remove_taxa(list(set(aln_obj.taxa_names) - set(tx_set)))
    return aln_obj


def process_execution(aln_list, file_set_name, file_list, file_groups, 
                      taxa_set_name, active_taxa_list, ns, taxa_groups,
                      hap_prefix, secondary_operations, secondary_options,
                      missing_filter_settings, taxa_filter_settings,
                      codon_filter_settings, output_file, rev_infile, 
                      main_operations, zorro_suffix, partitions_file,
                      output_formats, create_partfile, use_nexus_partitions,
                      phylip_truncate_name, output_dir, use_app_partitions,
                      consensus_type, ld_hat, temp_dir, ima2_params):
    """
    Process execution function
    :param ns: Namespace object
    """

    def reverse_concatenation(aln):
        """
        Wrapper of the reverse concatenation operation
        :return: AlignmentList object
        """

        if not use_app_partitions:
            partition_obj = data.Partitions()
            # In case the partitions file is badly formatted or invalid, the
            # exception will be returned by the read_from_file method.
            er = partition_obj.read_from_file(partitions_file)
            aln = aln.retrieve_alignment(basename(rev_infile))

            aln.set_partitions(partition_obj)

        aln = aln.reverse_concatenate(dest=temp_dir)

        return aln

    def filter_aln(aln):
        """
        Wrapper for filtering operations, given an alignment object
        :param aln: AlignmentList object
        """

        # Check if a minimum taxa representation was specified
        if secondary_options["gap_filter"]:
            if missing_filter_settings[2]:
                aln.filter_min_taxa(missing_filter_settings[2])

        # Filter by taxa
        if secondary_options["taxa_filter"]:
            # Get taxa list from taxa groups
            taxa_list = taxa_groups[taxa_filter_settings[1]]
            aln.filter_by_taxa(taxa_filter_settings[0], taxa_list)

        # Filter codon positions
        if secondary_options["codon_filter"]:
            aln.filter_codon_positions(codon_filter_settings)

        # Filter missing data
        if secondary_options["gap_filter"]:
            aln.filter_missing_data(missing_filter_settings[0],
                                    missing_filter_settings[1])

        return aln

    def concatenation(aln, remove_temp=True):
        """
        Wrapper for concatenation operation
        :param aln: AlignmentList object
        """

        aln = aln.concatenate(alignment_name=basename(output_file),
                              dest=temp_dir, remove_temp=remove_temp)

        if secondary_options["zorro"]:
            ns.msg = "Concatenating ZORRO files"
            zorro_data = data.Zorro(aln, zorro_suffix)
            zorro_data.write_to_file(output_file)

        return aln

    def consensus(aln):
        """
        Wrapper of the consensus operation
        :param aln: AlignmentObject list
        """

        if secondary_options["consensus_single"]:
            aln = aln.consensus(consensus_type=consensus_type,
                                single_file=True)
        else:
            aln.consensus(consensus_type=consensus_type)

        return aln

    def writer(aln, filename=None, suffix_str=""):
        """
        Wrapper for the output writing operations
        :param aln: AlignmentList object
        :param filename: string. If provided, it will overwrite the output_file
        :param suffix_str: string. Provides the suffix for the AlignmentList
        write_to_file method
        argument
        """

        if filename:
            outfile = filename
        else:
            outfile = output_file

        # The output file(s) will only be written after all the required
        # operations have been concluded. The reason why there are two "if"
        # statement for "concatenation" is that the input alignments must be
        # concatenated before any other additional operations. If the
        # first if statement did not exist, then all additional options would
        # have to be manually written for both "conversion" and "concatenation".
        #  As it is, when "concatenation", the aln_obj is firstly converted
        # into the concatenated alignment, and then all additional
        # operations are conducted in the same aln_obj
        ns.msg = "Writing output"

        if isinstance(aln, Alignment):
            aln.write_to_file(output_formats,
                outfile if outfile else join(output_dir, "consensus"),
                interleave=secondary_options["interleave"],
                partition_file=create_partfile,
                use_charset=use_nexus_partitions,
                phy_truncate_names=phylip_truncate_name,
                ld_hat=ld_hat,
                ima2_params=ima2_params)
        elif isinstance(aln, AlignmentList):
            aln.write_to_file(
                output_formats,
                output_suffix=suffix_str,
                interleave=secondary_options["interleave"],
                partition_file=create_partfile,
                output_dir=output_dir,
                use_charset=use_nexus_partitions,
                phy_truncate_names=phylip_truncate_name,
                ld_hat=ld_hat,
                ima2_params=ima2_params)

    try:
        ns.msg = "Setting active data sets"
        # Setting the alignment to use.
        # Update active file set of the alignment object
        aln_object = update_active_fileset(aln_list,
                                           file_set_name,
                                           file_list,
                                           file_groups)

        # Update active taxa set of the alignment object
        aln_object = update_active_taxaset(aln_object, taxa_set_name,
                                           active_taxa_list,
                                           taxa_groups)

        ns.proc_files = len(aln_object.alignments)

        # The execution of the process module will begin with all the
        # operations on the main output alignment. Only after the main
        # output file has been created will the additional secondary output
        # files be processed. Since each output file requires the duplication
        # of the temporary sequence files for modification, this ensures that
        # there is only one set of "duplicate" temporary sequences files at
        # one time.

        #####
        # Perform operations on MAIN OUTPUT
        #####
        ns.msg = "Preparing data"
        main_aln = deepcopy(aln_object)
        main_aln.start_action_alignment()
        # Reverse concatenation
        if main_operations["reverse_concatenation"]:
            ns.msg = "Reverse concatenating"
            main_aln = reverse_concatenation(main_aln)
        # Filtering
        if secondary_operations["filter"] and not \
                secondary_options["filter_file"]:
            ns.msg = "Filtering alignment(s)"
            main_aln = filter_aln(main_aln)
        # Concatenation
        if main_operations["concatenation"]:
            ns.msg = "Concatenating"
            main_aln = concatenation(main_aln)
        # Collapsing
        if secondary_operations["collapse"] and not \
                secondary_options["collapse_file"]:
            ns.msg = "Collapsing alignment(s)"
            main_aln.collapse(haplotype_name=hap_prefix, dest=output_dir)
        # Gcoder
        if secondary_operations["gcoder"] and not \
                secondary_options["gcoder_file"]:
            ns.msg = "Coding gaps"
            main_aln.code_gaps()
        # Consensus
        if secondary_operations["consensus"] and not \
                secondary_options["consensus_file"]:
            ns.msg = "Creating consensus sequence(s)"
            main_aln = consensus(main_aln)

        # Writing main output
        writer(main_aln)

        main_aln.stop_action_alignment()

        #####
        # Perform operations on ADDITIONAL OUTPUTS
        #####

        # Stores operations that must be performed before concatenation,
        # if it was specified
        before_conc = ["filter_file", "consensus_file"]

        # Perform the filtering and consensus option separately, since these
        # must be done before concatenation

        for op in [x for x, y in secondary_options.items() if
                   x in before_conc and y]:

            main_aln = deepcopy(aln_object)
            main_aln.start_action_alignment()

            if op == "filter_file" and secondary_operations["filter"]:

                ns.msg = "Creating additional filtered alignments(s)"
                suffix = "_filtered"
                main_aln = filter_aln(main_aln)

            elif op == "consensus_file" and secondary_operations["consensus"]:
                suffix = "_consensus"
                main_aln = consensus(main_aln)
                if secondary_options["consensus_single"]:
                    filename = join(output_dir, "consensus")

            if main_operations["concatenation"] and \
                    isinstance(main_aln, AlignmentList):
                filename = output_file + suffix
                main_aln = concatenation(main_aln)
                writer(main_aln, filename=filename)
            else:
                writer(main_aln, suffix_str=suffix,
                       filename=filename if filename else None)

            main_aln.stop_action_alignment()

        if len([x for x, y in secondary_options.items() if
                x.endswith("_file") and y and x != "filter_file" and
                x != "consensus_file"]) >= 1:

            main_aln = deepcopy(aln_object)

            if main_operations["concatenation"]:
                main_aln = concatenation(main_aln, remove_temp=False)

            # Perform the remaining secondary operations
            for op in [x for x, y in secondary_options.items()
                       if x.endswith("_file") and y and x != "filter_file"]:

                main_aln.start_action_alignment()

                if op == "collapse_file":
                    if secondary_operations["collapse"]:
                        ns.msg = "Creating additional collapsed alignment(s)"
                        suffix = "_collapsed"

                        main_aln.collapse(haplotype_name=hap_prefix,
                                          dest=output_dir)

                elif op == "gcoder_file":
                    if secondary_operations["gcoder"]:
                        ns.msg = "Creating additional gap coded alignments(s)"
                        suffix = "_coded"
                        main_aln.code_gaps()

                if main_operations["concatenation"]:
                    filename = output_file + suffix
                    writer(main_aln, filename=filename)
                else:
                    writer(main_aln, suffix_str=suffix)

                main_aln.stop_action_alignment()

            else:
                if main_operations["concatenation"]:
                    main_aln.stop_action_alignment()

    except EmptyAlignment:
        logging.exception("Empty alignment")
        ns.exception = "EmptyAlignment"

    except:
        # Log traceback in case any unexpected error occurs. See
        # self.log_file for whereabouts of the traceback
        logging.exception("Unexpected exit in Process execution")
        ns.exception = "Unknown"


def load_group_files(group_files, temp_dir):
    og = OrthoTool.MultiGroupsLight(db_path=temp_dir,
                                    groups=group_files)
    return [og, og.filters]


def orto_update_filters(ortho_groups, gn_filter, sp_filter,
                        group_names=None, default=False):
    if group_names:
        ortho_groups.update_filters(gn_filter, sp_filter, group_names,
                                    default=default)
    else:
        ortho_groups.update_filters(gn_filter, sp_filter,
                                    default=default)
    return [ortho_groups]


def get_active_group(ortho_groups, old_active_group, active_group_name):

    if not old_active_group:
        active_group = ortho_groups.get_group(active_group_name)
    else:
        if active_group_name == old_active_group.name:
            return [None]
        else:
            active_group = ortho_groups.get_group(active_group_name)

    return [active_group]


def get_stats_data(aln_obj, stats_idx, active_file_set, active_taxa_set,
                   additional_args):
    """
    Given an aln_obj, this function will execute the according method to
    generate plot data

    :param aln_obj: AlignmentObject
    :param stats_idx: string, identifier that maps to an AlignmentObject method
    :return: data for plot production
    """

    # Update alignment object according to active file and taxa sets
    aln_obj.update_active_alignments(active_file_set)
    aln_obj.remove_taxa(list(set(aln_obj.taxa_names) - set(active_taxa_set)))

    # Check if active data sets are not empty. If so, raise an exception
    if aln_obj.alignments == OrderedDict() or not aln_obj.taxa_names:
        return [EmptyAlignment("Active alignment is empty")]

    # List of gene specific idx. These plots only have one gene for the footer
    gene_specific = ["Pairwise sequence similarity gn"]

    if stats_idx in gene_specific:
        footer = [1, len(aln_obj.taxa_names)]
    else:
        footer = [len(aln_obj.alignments), len(aln_obj.taxa_names)]

    methods = {"Gene occupancy": aln_obj.gene_occupancy,
               "Distribution of missing data": aln_obj.missing_data_per_species,
               "Distribution of missing orthologs":
                   aln_obj.missing_genes_per_species,
               "Distribution of missing orthologs avg":
                   aln_obj.missing_genes_average,
               "Distribution of sequence size":
                   aln_obj.average_seqsize_per_species,
               "Distribution of sequence size all": aln_obj.average_seqsize,
               "Proportion of nucleotides or residues":
                   aln_obj.characters_proportion,
               "Proportion of nucleotides or residues sp":
                   aln_obj.characters_proportion_per_species,
               "Pairwise sequence similarity": aln_obj.sequence_similarity,
               "Pairwise sequence similarity sp":
                   aln_obj.sequence_similarity_per_species,
               "Pairwise sequence similarity gn":
                   aln_obj.sequence_similarity_gene,
               "Segregating sites": aln_obj.sequence_segregation,
               "Segregating sites sp": aln_obj.sequence_segregation_per_species,
               "Segregating sites gn": aln_obj.sequence_segregation_gene}

    if additional_args:
        plot_data = methods[stats_idx](**additional_args)
    else:
        plot_data = methods[stats_idx]()

    return [plot_data, footer]

