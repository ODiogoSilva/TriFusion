try:
    from process import data
    from process.error_handling import *
    from ortho.error_handling import *
    from process.sequence import AlignmentList, Alignment
    import orthomcl_pipeline as ortho_pipe
    from ortho import OrthomclToolbox as OrthoTool
except ImportError:
    from trifusion.process import data
    from trifusion.process.error_handling import *
    from trifusion.ortho.error_handling import *
    from trifusion.process.sequence import AlignmentList, Alignment
    import trifusion.orthomcl_pipeline as ortho_pipe
    from trifusion.ortho import OrthomclToolbox as OrthoTool

from os.path import join, basename
from collections import OrderedDict
from copy import deepcopy
import logging
import shutil
import cPickle as pickle

import os
os.environ["KIVY_NO_ARGS"] = "1"

"""
Tasks that are executed in a background process are placed here and called
from the App class in TriFusion.py. The reason for this is that
multiprocessing in Windows does not rely on os.fork to spawn a new process
child with shared memory maps with the parent. Instead it relies on pickle to
transfer data. This reliance on pickle prevents application methods from being
sent to the background process because pickle cannot handle class methods.
The solution was to define the functions outside the App class and provide
them with the necessary arguments upon calling.
"""


def remove_tmp(temp_dir, sql_con):
    """
    Removes all temporary files in temp directory
    :param temp_dir: string, path to trifusion's temporary directory
    """

    sql_con.close()

    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)

    return 1


def load_proc(aln_list, file_list, nm, dest, queue):

    try:
        if aln_list:
            aln_list.add_alignment_files(file_list,
                                         shared_namespace=nm)
            aln_obj = aln_list
        else:
            aln_obj = AlignmentList(file_list, shared_namespace=nm)

        queue.put(aln_obj)

    except MultipleSequenceTypes:
        nm.exception = "multiple_type"

    except IOError:
        return

    except:
        logging.exception("Unexpected error when loading input data")
        nm.exception = True


def get_stats_summary(dest, aln_list, active_file_set, active_taxa_set,
                        ns):
    """
    Runs the get_summary_stats method in the background and writes the output
    in a pickle file
    :param aln_list: AlignmentList object
    :param dest: temporary file where stats will be written
    :param active_file_set: list, with files to be included in summary
    statistics
    :param active_taxa_set: list, with taxa to be included in summary statistics
    """

    try:

        # Update alignment object according to active file and taxa sets
        aln_list.update_active_alignments(active_file_set)
        aln_list.remove_taxa(list(set(aln_list.taxa_names) -
                                  set(active_taxa_set)))

        with open(join(dest, "stats.pc"), "wb") as fh_stats, \
                open(join(dest, "table.pc"), "wb") as fh_table:

            # Check if active data sets are not empty. If so, raise an exception
            if aln_list.alignments == OrderedDict() or not aln_list.taxa_names:
                for fh in [fh_stats, fh_table]:
                    pickle.dump({"exception": "Alignment is empty after file and "
                                              "taxa filters"}, fh)
                return

            stats = aln_list.get_summary_stats(ns=ns)
            table = aln_list.get_gene_table_stats()
            pickle.dump(stats, fh_stats)
            pickle.dump(table, fh_table)

    except KillByUser:
        pass


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
        logging.exception("Unexpected exit in %s" % f.__name__)
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
                   protein_max_stop, usearch_file, usearch_evalue,
                   usearch_threads, usearch_output, mcl_file, mcl_inflation,
                   ortholog_prefix, group_prefix, orto_max_gene,
                   orto_min_sp, sqldb, ortho_dir, usearch_db):
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
            ortho_pipe.adjust_fasta(proteome_files, ortho_dir)
        if nm.k:
            nm.t = "Filtering Fasta Files"
            nm.c = 3
            ortho_pipe.filter_fasta(protein_min_len, protein_max_stop,
                                    usearch_db, ortho_dir)
        if nm.k:
            nm.t = "Running USearch. This may take a while..."
            nm.c = 4
            ortho_pipe.allvsall_usearch(usearch_db, usearch_evalue, ortho_dir,
                                        usearch_threads, usearch_output,
                                        usearch_bin=usearch_file)
        if nm.k:
            nm.t = "Parsing USEARCH output"
            nm.c = 5
            ortho_pipe.blast_parser(usearch_output, ortho_dir,
                                    db_dir=temp_dir)
        if nm.k:
            nm.t = "Obtaining Pairs"
            nm.c = 6
            ortho_pipe.pairs(temp_dir)
        if nm.k:
            ortho_pipe.dump_pairs(temp_dir, ortho_dir)
        if nm.k:
            nm.t = "Running MCL"
            nm.c = 7
            ortho_pipe.mcl(mcl_inflation, ortho_dir, mcl_file=mcl_file)
        if nm.k:
            nm.t = "Dumping groups"
            nm.c = 8
            ortho_pipe.mcl_groups(mcl_inflation, ortholog_prefix, "1000",
                                  group_prefix, ortho_dir)
        if nm.k:
            nm.t = "Filtering group files"
            nm.c = 9
            stats, groups_obj = ortho_pipe.export_filtered_groups(
                mcl_inflation,
                group_prefix,
                orto_max_gene,
                orto_min_sp,
                sqldb,
                join(ortho_dir, "backstage_files", usearch_db),
                temp_dir,
                ortho_dir)
            # stats is a dictionary containing the inflation value as
            #  key and a list with the orthologs as value
            nm.stats = stats
            nm.groups = groups_obj

    except IOError as e:
        nm.exception = str(e)
        print(e)
        return

    except Exception as e:
        logging.exception("Unexpected exit in Orthology search")
        nm.exception = str(e)


def update_active_fileset(aln_obj, set_name, file_list, file_groups,
                          filename_map):
    """
    This method is similar in purpose and functionality to the
    update_active_taxaset, but it updates the set of files. It should be
    used before the said method.
    :param aln_obj: The alignment object being used during execution of
    Process operations
    """
    # Determine the selected active taxa set from the dropdown menu
    if set_name == "All files":
        aln_obj.update_active_alignments([x for x in file_list])
        return aln_obj
    if set_name == "Active files":
        return aln_obj
    else:
        aln_obj.update_active_alignments(
            [filename_map[x] for x in file_groups[set_name]])
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

    # Update active taxa
    aln_obj.update_taxa_names(tx_set)
    return aln_obj


def process_execution(aln_list, file_set_name, file_list, file_groups,
                      filename_map,
                      taxa_set_name, active_taxa_list, ns, taxa_groups,
                      hap_prefix, secondary_operations, secondary_options,
                      missing_filter_settings, taxa_filter_settings,
                      codon_filter_settings, variation_filter_settings,
                      output_file, rev_infile, main_operations, zorro_suffix,
                      partitions_file, output_formats, create_partfile,
                      use_nexus_partitions, use_nexus_models,
                      phylip_truncate_name, output_dir, use_app_partitions,
                      consensus_type, ld_hat, temp_dir, ima2_params,
                      conversion_suffix):
    """
    Process execution function
    :param ns: Namespace object
    """

    def reverse_concatenation(aln):
        """
        Wrapper of the reverse concatenation operation
        :return: AlignmentList object
        """

        con = aln.con

        if not use_app_partitions:
            partition_obj = data.Partitions()
            # In case the partitions file is badly formatted or invalid, the
            # exception will be returned by the read_from_file method.
            er = partition_obj.read_from_file(partitions_file)
            aln = aln.retrieve_alignment(rev_infile)

            aln.set_partitions(partition_obj)

        if isinstance(aln, AlignmentList):
            aln = aln.reverse_concatenate()
        else:
            aln = aln.reverse_concatenate(db_con=con)

        return aln

    def filter_aln(aln, table_out="_filter"):
        """
        Wrapper for filtering operations, given an alignment object
        :param aln: AlignmentList object
        :param table_out: string. Suffix of the sqlite table for operations
        that change alignments
        """

        # Check if a minimum taxa representation was specified
        if secondary_options["gap_filter"]:
            if missing_filter_settings[1][0]:
                aln.filter_min_taxa(missing_filter_settings[1][1])

        # Filter by taxa
        if secondary_options["taxa_filter"]:
            # Get taxa list from taxa groups
            taxa_list = taxa_groups[taxa_filter_settings[1]]
            aln.filter_by_taxa(taxa_filter_settings[0], taxa_list)

        # Filter codon positions
        if secondary_options["codon_filter"]:
            aln.filter_codon_positions(codon_filter_settings,
                                       table_out=table_out)

        # Filter missing data
        if secondary_options["gap_filter"]:
            if missing_filter_settings[0][0]:
                aln.filter_missing_data(missing_filter_settings[0][1],
                                        missing_filter_settings[0][2],
                                        table_out=table_out)

        # Filter variation
        if secondary_options["variation_filter"]:
            # Checks for variable site filter
            if variation_filter_settings[0] or variation_filter_settings[1]:
                aln.filter_segregating_sites(variation_filter_settings[0],
                                             variation_filter_settings[1],
                                             table_in=table_out)
            # Checks for informative site filter
            if variation_filter_settings[2] or variation_filter_settings[3]:
                aln.filter_informative_sites(variation_filter_settings[2],
                                             variation_filter_settings[3],
                                             table_in=table_out)

        # Pipe the information on the filtered alignments to the main process
        # only if it was applied a filter that changes the final alignments
        if set(aln.filtered_alignments.values()) != {None}:
            ns.filtered_alns = aln.filtered_alignments

        # Some filter configurations may result in empty final alignment
        # list. In such cases, return and issue warning
        if not aln.alignments:
            raise EmptyAlignment("Active alignment is empty")

        return aln

    def concatenation(aln, table_in=""):
        """
        Wrapper for concatenation operation
        :param aln: AlignmentList object
        """

        if secondary_options["zorro"]:
            ns.msg = "Concatenating ZORRO files"
            zorro_data = data.Zorro(aln, zorro_suffix)
            zorro_data.write_to_file(output_file)

        aln = aln.concatenate(alignment_name=basename(output_file),
                              table_in=table_in)

        return aln

    def consensus(aln, table_out):
        """
        Wrapper of the consensus operation
        :param aln: AlignmentObject list
        """

        if secondary_options["consensus_single"]:
            if isinstance(aln, AlignmentList):
                aln = aln.consensus(consensus_type=consensus_type,
                                    single_file=True,
                                    table_out=table_out)
            else:
                aln.consensus(consensus_type=consensus_type,
                              table_out=table_out)
        else:
            aln.consensus(consensus_type=consensus_type,
                          table_out=table_out)

        return aln

    def writer(aln, filename=None, suffix_str="", conv_suffix="",
               table_suffix=None, table_name=None):
        """
        Wrapper for the output writing operations
        :param aln: AlignmentList object
        :param filename: string. If provided, it will overwrite the output_file
        :param suffix_str: string. Provides the suffix for the AlignmentList
        write_to_file method argument
        :param conv_suffix: string. Provides the suffix provided for
        the conversion of files. This suffix will allway precede the
        suffix_str, which is meant to apply suffixes specific to secondary
        operations.
        :param table_suffix: string. Suffix of the table from where the
        sequence data will be retrieved
        :param table_name: string. Name of the table from where the
        sequence data will be retrived
        """

        try:
            if filename:
                outfile = filename
            else:
                outfile = output_file

            # The output file(s) will only be written after all the required
            # operations have been concluded. The reason why there are two "if"
            # statement for "concatenation" is that the input alignments must
            # be concatenated before any other additional operations. If the
            # first if statement did not exist, then all additional options
            # would have to be manually written for both "conversion" and
            # "concatenation". As it is, when "concatenation", the aln_obj is
            # firstly converted into the concatenated alignment, and then all
            # additional operations are conducted in the same aln_obj
            ns.msg = "Writing output"

            if isinstance(aln, Alignment):
                aln.write_to_file(output_formats,
                    outfile if outfile else join(output_dir, "consensus"),
                    interleave=secondary_options["interleave"],
                    partition_file=create_partfile,
                    use_charset=use_nexus_partitions,
                    phy_truncate_names=phylip_truncate_name,
                    ld_hat=ld_hat,
                    ima2_params=ima2_params,
                    use_nexus_models=use_nexus_models,
                    ns_pipe=ns,
                    table_suffix=table_suffix,
                    table_name=table_name)
            elif isinstance(aln, AlignmentList):
                aln.write_to_file(
                    output_formats,
                    output_suffix=suffix_str,
                    conversion_suffix=conv_suffix,
                    interleave=secondary_options["interleave"],
                    partition_file=create_partfile,
                    output_dir=output_dir,
                    use_charset=use_nexus_partitions,
                    phy_truncate_names=phylip_truncate_name,
                    ld_hat=ld_hat,
                    ima2_params=ima2_params,
                    use_nexus_models=use_nexus_models,
                    ns_pipe=ns,
                    table_suffix=table_suffix,
                    table_name=table_name)

        except IOError:
            pass

    try:

        aln_object = deepcopy(aln_list)
        # Restore database connections, since they are broken during the
        # deepcopy operation
        aln_object.set_database_connections(aln_list.cur, aln_list.con)

        ns.msg = "Setting active data sets"
        # Setting the alignment to use.
        # Update active file set of the alignment object
        aln_object = update_active_fileset(aln_object,
                                           file_set_name,
                                           file_list,
                                           file_groups,
                                           filename_map)

        # Update active taxa set of the alignment object
        main_aln = update_active_taxaset(aln_object, taxa_set_name,
                                           active_taxa_list,
                                           taxa_groups)

        ns.proc_files = len(aln_object.alignments)

        # Initialize attribute tha will store the number of filtered
        # alignments for reporting purposes
        ns.filtered_alns = None

        # The execution of the process module will begin with all the
        # operations on the main output alignment. Only after the main
        # output file has been created will the additional secondary output
        # files be processed.

        #####
        # Perform operations on MAIN OUTPUT
        #####

        # Set the suffix for the sqlite table harboring the main output
        # alignment if any of the secondary operations is specified
        if any(secondary_operations.values()):
            main_table = "main"
        else:
            main_table = ""

        # Reverse concatenation
        # Active table: Based on partition names
        if main_operations["reverse_concatenation"]:
            ns.msg = "Reverse concatenating"
            main_aln = reverse_concatenation(main_aln)
        # Filtering
        # Active table: * / *main
        if secondary_options["collapse_filter"] and not \
                secondary_options["collapse_file"]:
            ns.msg = "Filtering data for collapsing"
            # If the the collapse filter is active, perform this
            # filtering first. This is because the filter will allow 0% of
            # missing data, which will always be as stringent or more than any
            # missing data filter set.
            main_aln.filter_missing_data(0, 0, table_out=main_table)

        # Active table: * / *main
        if secondary_operations["filter"] and not \
                secondary_options["filter_file"]:
            ns.msg = "Filtering alignment(s)"
            main_aln = filter_aln(main_aln, table_out=main_table)
        # Concatenation
        # Active table: concatenation
        if main_operations["concatenation"]:
            ns.msg = "Concatenating"
            main_aln = concatenation(main_aln, table_in=main_table)
        # Collapsing
        # Active table: *main / concatenationmain
        if secondary_operations["collapse"] and not \
                secondary_options["collapse_file"]:
            ns.msg = "Collapsing alignment(s)"
            main_aln.collapse(haplotype_name=hap_prefix, dest=output_dir,
                                table_out=main_table)
        # Gcoder
        # Active table: *main / concatenationmain
        if secondary_operations["gcoder"] and not \
                secondary_options["gcoder_file"]:
            ns.msg = "Coding gaps"
            main_aln.code_gaps(table_out=main_table)
        # Consensus
        # Active table: *main / concatenationmain / consensus
        if secondary_operations["consensus"] and not \
                secondary_options["consensus_file"]:
            ns.msg = "Creating consensus sequence(s)"
            main_aln = consensus(main_aln, table_out=main_table)

        # ### Guide on possible final tables
        # --- Base types
        # 1: Conversion -> * (main -- checks)
        # 2: Concatenation -> concatenation (main -- checks)
        # 3: Reverse concatenation -> * (main --checks)
        # --- Conversion and Reverse concatenation + secondary ops
        # 4: Collapse/Filter/Gcoder -> *main (sec --checks)
        # 5: Consensus -> *main (sec -- checks)
        # 6: Consensus (single file) -> consensus (main -- fallback)
        # --- Concatenation + secondary ops
        # 7: Filter (only) -> concatenation (main -- fallback)
        # 8: Collapse/gcoder -> concatenationmain (sec -- checks)
        # 9: Consensus -> concatenationmain (sec -- checks)
        # 10: Consensus (single file) -> consensus (main -- fallback)

        # NOTE ON TABLE NAMES
        # Some combinations of operations actually create a table_suffix
        # that does not exist in the database. This happens in cases 6, 7 and
        # 10. However, in these cases the main table(s) of the Aligment
        # object should be used, so we let the sequence fetching methods
        # fail to find the suggested table and fallback to the main table.

        writer(main_aln, conv_suffix=conversion_suffix,
               table_suffix=main_table)

        #####
        # Perform operations on ADDITIONAL OUTPUTS
        #####

        # Stores operations that must be performed before concatenation,
        # if it was specified
        before_conc = ["filter"]

        # Perform the filtering and consensus option separately, since these
        # must be done before concatenation
        for op in [x for x, y in secondary_operations.items() if
                   x in before_conc and y and
                   secondary_options["%s_file" % x]]:

            if op == "filter":

                # Remove previous temporary tables
                aln_object.remove_tables(aln_object.get_tables())

                main_aln = deepcopy(aln_object)
                main_aln.set_database_connections(aln_list.cur, aln_list.con)

                ns.msg = "Creating additional filtered alignments(s)"
                suffix = "_filtered"
                main_aln = filter_aln(main_aln, table_out=suffix[1:])

            if main_operations["concatenation"] and \
                    isinstance(main_aln, AlignmentList):

                filename = output_file + suffix
                main_aln = concatenation(main_aln, table_in=suffix[1:])
                writer(main_aln, filename=filename, table_suffix=suffix[1:])
            else:
                writer(main_aln, suffix_str=suffix, table_suffix=suffix[1:],
                       conv_suffix=conversion_suffix)

        # Remove previous temporary tables
        aln_object.remove_tables(aln_object.get_tables())

        main_aln = deepcopy(aln_object)
        main_aln.set_database_connections(aln_list.cur, aln_list.con)

        concatenated = False

        for op in [x for x, y in secondary_operations.items() if
                   x not in before_conc and y and
                   secondary_options["%s_file" % x]]:

            # Filter data for collapsing
            if secondary_options["collapse_filter"] and op == "collapse":
                main_aln.filter_missing_data(0, 0)

            if main_operations["concatenation"]:
                # If suffix was specified, it means that the filter was
                # ON. In that case, use that suffix in the table input
                # for concatenation.
                if not concatenated:
                    try:
                        main_aln = concatenation(aln_object)
                    except NameError:
                        main_aln = concatenation(aln_object)

                    concatenated = True

            if op == "consensus":

                ns.msg = "Creating additional consensus alignment(s)"
                suffix = "_consensus"
                main_aln = consensus(main_aln, table_out=suffix[1:])

            elif op == "gcoder":
                ns.msg = "Creating additional gap coded alignments(s)"
                suffix = "_coded"
                main_aln.code_gaps(table_out=suffix[1:])

            elif op == "collapse":

                # In this case, filter the unconcatenated alignment
                # and concatenate again
                if secondary_options["collapse_filter"]:
                    ns.msg = "Filtering alignments before collapsing"
                    aln_object.remove_tables(aln_object.get_tables())
                    aln_object.filter_missing_data(0, 0)
                    main_aln = aln_object.concatenate(table_in="filter")

                ns.msg = "Creating additional collapsed alignment(s)"
                suffix = "_collapsed"

                main_aln.collapse(haplotype_name=hap_prefix,
                                  conversion_suffix=conversion_suffix,
                                  dest=output_dir,
                                  table_out=suffix[1:])

            if main_operations["concatenation"]:
                filename = output_file + suffix
                writer(main_aln, filename=filename, table_suffix=suffix[1:])
            else:
                writer(main_aln, suffix_str=suffix,
                       conv_suffix=conversion_suffix, table_suffix=suffix[1:])

        # Resets the taxa_names attribute of the aln_obj to include all taxa
        # aln_object.update_taxa_names(all_taxa=True)

    except IOError:
        # Resets the taxa_names attribute of the aln_obj to include all taxa
        # aln_object.update_taxa_names(all_taxa=True)
        return

    except EmptyAlignment:
        logging.exception("Empty alignment")
        # Resets the taxa_names attribute of the aln_obj to include all taxa
        # aln_object.update_taxa_names(all_taxa=True)
        # ns.exception = "EmptyAlignment"

    except Exception as e:
        print(e)
        # Log traceback in case any unexpected error occurs. See
        # self.log_file for whereabouts of the traceback
        logging.exception("Unexpected exit in Process execution")
        # Resets the taxa_names attribute of the aln_obj to include all taxa
        # aln_object.update_taxa_names(all_taxa=True)
        ns.exception = "Unknown"


def load_group_files(group_files, temp_dir):
    og = OrthoTool.MultiGroupsLight(db_path=temp_dir,
                                    groups=group_files)
    return [og, og.filters]


def orto_update_filters(ortho_groups, gn_filter, sp_filter,
                        group_names=None, default=False):
    try:
        if group_names:
            ortho_groups.update_filters(gn_filter, sp_filter, group_names,
                                        default=default)
        else:
            ortho_groups.update_filters(gn_filter, sp_filter,
                                        default=default)
    except EOFError:
        pass
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
               "Distribution of missing data sp":
                   aln_obj.missing_data_per_species,
               "Distribution of missing data":
                   aln_obj.missing_data_distribution,
               "Distribution of missing orthologs":
                   aln_obj.missing_genes_per_species,
               "Distribution of missing orthologs avg":
                   aln_obj.missing_genes_average,
               "Distribution of sequence size":
                   aln_obj.average_seqsize_per_species,
               "Distribution of sequence size all": aln_obj.average_seqsize,
               "Cumulative distribution of missing genes":
                   aln_obj.cumulative_missing_genes,
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
               "Segregating sites gn": aln_obj.sequence_segregation_gene,
               "Segregating sites prop": aln_obj.sequence_segregation,
               "Alignment length/Polymorphism correlation":
                   aln_obj.length_polymorphism_correlation,
               "Distribution of taxa frequency":
                   aln_obj.taxa_distribution,
               "Allele Frequency Spectrum":
                   aln_obj.allele_frequency_spectrum,
               "Allele Frequency Spectrum prop":
                   aln_obj.allele_frequency_spectrum,
               "Allele Frequency Spectrum gn":
                   aln_obj.allele_frequency_spectrum_gene,
               "Missing data outliers": aln_obj.outlier_missing_data,
               "Missing data outliers sp": aln_obj.outlier_missing_data_sp,
               "Segregating sites outliers": aln_obj.outlier_segregating,
               "Segregating sites outliers sp": aln_obj.outlier_segregating_sp,
               "Sequence size outliers sp": aln_obj.outlier_sequence_size_sp,
               "Sequence size outliers": aln_obj.outlier_sequence_size}

    if additional_args:
        plot_data = methods[stats_idx](**additional_args)
    else:
        plot_data = methods[stats_idx]()

    return [plot_data, footer]


__author__ = "Diogo N. Silva"
