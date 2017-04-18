try:
    from process import data
    from process.error_handling import KillByUser,\
        MultipleSequenceTypes, EmptyAlignment
    from process.sequence import AlignmentList, Alignment
    import orthomcl_pipeline as ortho_pipe
    from ortho import OrthomclToolbox as OrthoTool
except ImportError:
    from trifusion.process import data
    from trifusion.process.error_handling import KillByUser,\
        MultipleSequenceTypes, EmptyAlignment
    from trifusion.process.sequence import AlignmentList, Alignment
    import trifusion.orthomcl_pipeline as ortho_pipe
    from trifusion.ortho import OrthomclToolbox as OrthoTool

from os.path import join, basename
from collections import OrderedDict
from copy import deepcopy
import logging
import shutil
import cPickle as pickle
import time

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

    # Give some time to child threads to exit
    time.sleep(1)

    # Close database connection
    sql_con.close()

    # Remove temporary files
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)

    return 1


def load_proc(aln_list, file_list, nm, queue):
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

    except KillByUser:
        return

    except Exception as e:
        logging.exception("Unexpected error when loading input data")
        print(e)


def get_stats_summary(dest, aln_list, active_file_set, active_taxa_set,
                        ns):
    """
    Runs the get_summary_stats method in the background and writes the output
    in a pickle file
    :param aln_list: AlignmentList object
    :param dest: temporary file where stats will be written
    :param active_file_set: list, with files to be included in summary
    statistics
    :param active_taxa_set: list, with taxa to be included in summary
    statistics
    """

    try:
        # Creating deepcopy to perform changes  without impacting main
        # attribute
        main_aln = deepcopy(aln_list)
        main_aln.set_database_connections(aln_list.cur, aln_list.con)
        # Update alignment object according to active file and taxa sets
        main_aln.update_active_alignments(active_file_set)
        main_aln.update_taxa_names(active_taxa_set)

        with open(join(dest, "stats.pc"), "wb") as fh_stats, \
                open(join(dest, "table.pc"), "wb") as fh_table:

            # Check if active data sets are not empty. If so, raise an
            # exception
            if main_aln.alignments == OrderedDict() or not main_aln.taxa_names:
                for fh in [fh_stats, fh_table]:
                    pickle.dump({"exception": "Alignment is empty after file "
                                              "and taxa filters"}, fh)
                return

            stats = main_aln.get_summary_stats(ns=ns)
            table = main_aln.get_gene_table_stats()
            pickle.dump(stats, fh_stats)
            pickle.dump(table, fh_table)

    except KillByUser:
        return


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
            if "use_ns" in a:
                a.remove("use_ns")
                val = f(*a, ns=ns)
            else:
                val = f(*a)
        else:
            val = f()
        ns.val = val

    except KillByUser:
        return

    except IOError as e:
        print(e)

    except Exception:
        logging.exception("Unexpected exit in %s", f.__name__)
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

    except KillByUser:
        pass

    except IOError as e:
        print(e)

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
        nm.finished_tasks = []

        nm.task = "schema"
        ortho_pipe.install_schema(temp_dir)
        nm.finished_tasks = ["schema"]

        if nm.stop:
            raise KillByUser("")

        nm.task = "adjust"
        ortho_pipe.adjust_fasta(proteome_files, ortho_dir, nm)
        nm.finished_tasks = ["schema", "adjust"]

        if nm.stop:
            raise KillByUser("")

        nm.task = "filter"
        ortho_pipe.filter_fasta(protein_min_len, protein_max_stop,
                                usearch_db, ortho_dir, nm)
        nm.finished_tasks = ["schema", "adjust", "filter"]

        if nm.stop:
            raise KillByUser("")

        nm.task = "usearch"
        ortho_pipe.allvsall_usearch(usearch_db, usearch_evalue, ortho_dir,
                                    usearch_threads, usearch_output,
                                    usearch_bin=usearch_file, nm=nm)
        nm.finished_tasks = ["schema", "adjust", "filter", "usearch"]

        if nm.stop:
            raise KillByUser("")

        nm.task = "parse"
        ortho_pipe.blast_parser(usearch_output, ortho_dir,
                                db_dir=temp_dir, nm=nm)
        nm.finished_tasks = ["schema", "adjust", "filter", "usearch", "parse"]

        if nm.stop:
            raise KillByUser("")

        nm.task = "pairs"
        ortho_pipe.pairs(temp_dir, nm=nm)
        ortho_pipe.dump_pairs(temp_dir, ortho_dir, nm=nm)
        nm.finished_tasks = ["schema", "adjust", "filter", "usearch", "parse",
                             "pairs"]

        if nm.stop:
            raise KillByUser("")

        nm.task = "mcl"
        ortho_pipe.mcl(mcl_inflation, ortho_dir, mcl_file=mcl_file, nm=nm)
        nm.finished_tasks = ["schema", "adjust", "filter", "usearch", "parse",
                             "pairs", "mcl"]

        if nm.stop:
            raise KillByUser("")

        nm.task = "dump"
        ortho_pipe.mcl_groups(mcl_inflation, ortholog_prefix, "1000",
                              group_prefix, ortho_dir, nm=nm)
        nm.finished_tasks = ["schema", "adjust", "filter", "usearch", "parse",
                             "pairs", "mcl", "dump"]

        if nm.stop:
            raise KillByUser("")

        nm.task = "filter_groups"
        stats, groups_obj = ortho_pipe.export_filtered_groups(
            mcl_inflation,
            group_prefix,
            orto_max_gene,
            orto_min_sp,
            sqldb,
            join(ortho_dir, "backstage_files", usearch_db),
            temp_dir,
            ortho_dir, nm=nm)
        nm.finished_tasks = ["schema", "adjust", "filter", "usearch", "parse",
                             "pairs", "mcl", "dump", "filter_groups"]

        if nm.stop:
            raise KillByUser("")

        # stats is a dictionary containing the inflation value as
        #  key and a list with the orthologs as value
        nm.stats = stats
        nm.groups = groups_obj

    except KillByUser:
        return

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
                        consensus_type, ld_hat, ima2_params,
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
            # Retrieve the alignment object that will be reverted. This
            # is done first in order to retrieve the length of the locus,
            # which is provided to the Partitions object for checking and
            # conversion of "." notation
            aln = aln.retrieve_alignment(rev_infile)

            # Instanciate Partitions object and set its length attribute
            partition_obj = data.Partitions()
            partition_obj.set_length(aln.locus_length)

            # In case the partitions file is badly formatted or invalid, the
            # exception will be returned by the read_from_file method.
            e = partition_obj.read_from_file(partitions_file)
            if e:
                ns.exception = {"InvalidPartitionFile": e.value}
                raise data.InvalidPartitionFile("")

            # If there are no issues with the partitions file, set the new
            # partitions
            aln.set_partitions(partition_obj)

        if aln.__class__.__name__ == "AlignmentList":
            aln = aln.reverse_concatenate(ns=ns)
        else:
            aln = aln.reverse_concatenate(db_con=con, ns=ns)

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
                ns.main_msg = "Filter (minimum taxa)"
                aln.filter_min_taxa(missing_filter_settings[1][1], ns=ns)

        # Filter by taxa
        if secondary_options["taxa_filter"]:
            # Get taxa list from taxa groups
            ns.main_msg = "Filter (by taxa)"
            taxa_list = taxa_groups[taxa_filter_settings[1]]
            aln.filter_by_taxa(taxa_filter_settings[0], taxa_list, ns=ns)

        # Filter codon positions
        if secondary_options["codon_filter"]:
            ns.main_msg = "Filter (by codon)"
            aln.filter_codon_positions(codon_filter_settings,
                                       table_out=table_out, ns=ns)

        # Filter missing data
        if secondary_options["gap_filter"]:
            ns.main_msg = "Filter (by missing data)"
            if missing_filter_settings[0][0]:
                aln.filter_missing_data(missing_filter_settings[0][1],
                                        missing_filter_settings[0][2],
                                        table_out=table_out, ns=ns)

        # Filter variation
        if secondary_options["variation_filter"]:
            # Checks for variable site filter
            if variation_filter_settings[0] or variation_filter_settings[1]:
                ns.main_msg = "Filter (by variable sites)"
                aln.filter_segregating_sites(variation_filter_settings[0],
                                             variation_filter_settings[1],
                                             table_in=table_out, ns=ns)
            # Checks for informative site filter
            if variation_filter_settings[2] or variation_filter_settings[3]:
                ns.main_msg = "Filter (by informative sites)"
                aln.filter_informative_sites(variation_filter_settings[2],
                                             variation_filter_settings[3],
                                             table_in=table_out, ns=ns)

        # Pipe the information on the filtered alignments to the main process
        # only if it was applied a filter that changes the final alignments
        if set(aln.filtered_alignments.values()) != {None}:
            ns.filtered_alns = aln.filtered_alignments

        # Reset main label text
        ns.main_msg = None

        # Some filter configurations may result in empty final alignment
        # list. In such cases, return and issue warning
        if not aln.alignments:
            ns.exception = "EmptyAlignment"
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
                              table_in=table_in, ns=ns)

        # Sets the single alignment to True, for other method to be aware of
        # this
        ns.sa = True

        return aln

    def consensus(aln, table_out):
        """
        Wrapper of the consensus operation
        :param aln: AlignmentObject list
        """

        if secondary_options["consensus_single"]:
            if aln.__class__.__name__ == "AlignmentList":
                aln = aln.consensus(consensus_type=consensus_type,
                                    single_file=True,
                                    table_out=table_out,
                                    ns=ns)
                ns.sa = True
            else:
                aln.consensus(consensus_type=consensus_type,
                              table_out=table_out,
                              ns=ns)
        else:
            aln.consensus(consensus_type=consensus_type,
                          table_out=table_out,
                          ns=ns)

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

            if aln.__class__.__name__ == "Alignment":
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
            elif aln.__class__.__name__ == "AlignmentList":
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

        except IOError as e:
            logging.exception(e)

    try:

        aln_object = deepcopy(aln_list)
        # Restore database connections, since they are broken during the
        # deepcopy operation
        aln_object.set_database_connections(aln_list.cur, aln_list.con)

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
            ns.task = "reverse_concatenation"
            main_aln = reverse_concatenation(main_aln)
            ns.finished_tasks.append("reverse_concatenation")

        # Filtering
        # Active table: * / *main
        if secondary_options["collapse_filter"] and not \
                secondary_options["collapse_file"]:
            ns.task = "collapse"
            # If the the collapse filter is active, perform this
            # filtering first. This is because the filter will allow 0% of
            # missing data, which will always be as stringent or more than any
            # missing data filter set.
            main_aln.filter_missing_data(0, 0, table_out=main_table, ns=ns)
            ns.finished_tasks.append("collapse_filter")

        # Active table: * / *main
        if secondary_operations["filter"] and not \
                secondary_options["filter_file"]:
            ns.task = "filter"
            main_aln = filter_aln(main_aln, table_out=main_table)
            ns.finished_tasks.append("filter")
        # Concatenation
        # Active table: concatenation
        if main_operations["concatenation"]:
            ns.task = "concatenation"
            main_aln = concatenation(main_aln, table_in=main_table)
            ns.finished_tasks.append("concatenation")
        # Collapsing
        # Active table: *main / concatenationmain
        if secondary_operations["collapse"] and not \
                secondary_options["collapse_file"]:
            ns.task = "collapse"
            main_aln.collapse(haplotype_name=hap_prefix, dest=output_dir,
                                table_out=main_table, ns=ns)
            ns.finished_tasks.append("collapse")
        # Gcoder
        # Active table: *main / concatenationmain
        if secondary_operations["gcoder"] and not \
                secondary_options["gcoder_file"]:
            ns.task = "gcoder"
            main_aln.code_gaps(table_out=main_table, ns=ns)
            ns.finished_tasks.append("gcoder")
        # Consensus
        # Active table: *main / concatenationmain / consensus
        if secondary_operations["consensus"] and not \
                secondary_options["consensus_file"]:
            ns.task = "consensus"
            main_aln = consensus(main_aln, table_out=main_table)
            ns.finished_tasks.append("consensus")

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
        # 10. However, in these cases the main table(s) of the Alignment
        # object should be used, so we let the sequence fetching methods
        # fail to find the suggested table and fallback to the main table.

        ns.task = "write"
        writer(main_aln, conv_suffix=conversion_suffix,
               table_suffix=main_table)
        ns.finished_tasks.append("write")

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

                ns.task = "filter"

                # Remove previous temporary tables
                aln_object.remove_tables(aln_object.get_tables())

                main_aln = deepcopy(aln_object)
                main_aln.set_database_connections(aln_list.cur, aln_list.con)

                ns.msg = "Creating additional filtered alignments(s)"
                suffix = "_filtered"
                main_aln = filter_aln(main_aln, table_out=suffix[1:])

            if main_operations["concatenation"] and \
                    main_aln.__class__.__name__ == "AlignmentList":

                filename = output_file + suffix

                ns.main_msg = "Concatenating"
                main_aln = concatenation(main_aln, table_in=suffix[1:])
                ns.main_msg = "Writing output"
                writer(main_aln, filename=filename, table_suffix=suffix[1:])
            else:
                ns.main_msg = "Writing output"
                writer(main_aln, suffix_str=suffix, table_suffix=suffix[1:],
                       conv_suffix=conversion_suffix)

        # Remove previous temporary tables
        aln_object.remove_tables(aln_object.get_tables())

        main_aln = deepcopy(aln_object)
        main_aln.set_database_connections(aln_list.cur, aln_list.con)

        concatenated = False
        ns.sa = False

        for op in [x for x, y in secondary_operations.items() if
                   x not in before_conc and y and
                   secondary_options["%s_file" % x]]:

            ns.task = op

            # Filter data for collapsing
            # if secondary_options["collapse_filter"] and op == "collapse":
                # ns.task = "collapse"
                # ns.main_msg = "Filtering for collapse"
                # main_aln.filter_missing_data(0, 0, ns=ns)

            if main_operations["concatenation"]:
                # If suffix was specified, it means that the filter was
                # ON. In that case, use that suffix in the table input
                # for concatenation.
                if op == "collapse" and secondary_options["collapse_filter"]:
                    pass
                else:
                    if not concatenated:
                        try:
                            main_aln = concatenation(aln_object)
                        except NameError:
                            main_aln = concatenation(aln_object)

                        concatenated = True

            if op == "consensus":

                suffix = "_consensus"
                main_aln = consensus(main_aln, table_out=suffix[1:])

            elif op == "gcoder":
                suffix = "_coded"
                main_aln.code_gaps(table_out=suffix[1:], ns=ns)

            elif op == "collapse":

                suffix = "_collapsed"

                # In this case, filter the unconcatenated alignment
                # and concatenate again
                if secondary_options["collapse_filter"]:
                    ns.main_msg = "Filtering for collapse"
                    aln_object.remove_tables(aln_object.get_tables())
                    aln_object.filter_missing_data(0, 0, ns=ns,
                                                   table_out=suffix[1:])
                    if main_operations["concatenation"]:
                        ns.main_msg = "Concatenation"
                        main_aln = aln_object.concatenate(table_in=suffix[1:],
                                                          ns=ns)

                ns.main_msg = "Collapse"
                main_aln.collapse(haplotype_name=hap_prefix,
                                  conversion_suffix=conversion_suffix,
                                  dest=output_dir,
                                  table_out=suffix[1:], ns=ns)

            ns.main_msg = "Writing output"
            if main_operations["concatenation"]:
                filename = output_file + suffix
                writer(main_aln, filename=filename, table_suffix=suffix[1:])
            else:
                writer(main_aln, suffix_str=suffix,
                       conv_suffix=conversion_suffix, table_suffix=suffix[1:])

        # Resets the taxa_names attribute of the aln_obj to include all taxa
        # aln_object.update_taxa_names(all_taxa=True)

    except KillByUser:
        pass

    except IOError as e:
        logging.exception(e)
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
        if not hasattr(ns, "exception"):
            ns.exception = "Unknown"


def load_group_files(group_files, temp_dir, ns=None):
    og = OrthoTool.MultiGroupsLight(db_path=temp_dir,
                                    groups=group_files,
                                    ns=ns)
    return [og, og.filters]


def orto_update_filters(ortho_groups, gn_filter, sp_filter, excluded_taxa,
                        group_names=None, default=False):
    try:
        if group_names or group_names == []:
            ortho_groups.update_filters(gn_filter, sp_filter, excluded_taxa,
                                        group_names, default=default)
        else:
            ortho_groups.update_filters(gn_filter, sp_filter, excluded_taxa,
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


def get_orto_data(active_group, plt_idx, filt, exclude_taxa):
    """
    Given a GroupLight object, this function will execute the method that
    corresponds to plt_idx to generate data.
    :param active_group: GroupLight object. Active group that will be used
    to plot data
    :param plt_idx: string. Identifier of the plot type that must have a
    correspondence in the method dictionary below
    :return plot_data: dictionary. Has the required data for the
    plotting methods to generate the plot.
    """

    # Store the plot generation method in a dictionary where keys are
    # the text attributes of the plot spinner and the values are
    # bound methods
    methods = {
        "Taxa distribution": active_group.bar_species_distribution,
        "Taxa coverage": active_group.bar_species_coverage,
        "Gene copy distribution": active_group.bar_genecopy_distribution,
        "Taxa gene copies": active_group.bar_genecopy_per_species
    }

    # Check for excluded taxa. If any have been provided and are different from
    # the ones already set in the GroupLight object, then update the taxa
    # filter.
    if (exclude_taxa and exclude_taxa != active_group.excluded_taxa) or \
            (exclude_taxa == [] and
             exclude_taxa != active_group.excluded_taxa):
        # The update_stats flag of the exclude_taxa method is set to True
        # to update the group summary stats. This is important for eventual
        # corrections to the ortholog cluster filters.
        active_group.exclude_taxa(exclude_taxa, True)

    if filt:
        # If filters AND excluded taxa have been provided, the first thing
        # to do is check whether the provided filters are still correct.
        # Excluding taxa may lead to changes in the maximum values of the
        # ortholog filters, and this needs to be corrected here
        if exclude_taxa and exclude_taxa != active_group.excluded_taxa:
            # Correct gene copy filter maximum
            gn_filt = filt[0] if filt[0] <= active_group.max_extra_copy \
                else active_group.max_extra_copy
            # Correct minimum taxa representation maximum
            sp_filt = filt[1] if \
                filt[1] <= len(active_group.species_list) \
                else len(active_group.species_list)
        else:
            # No taxa have been excluded this time, so we'll keep the provided
            # filters
            gn_filt, sp_filt = filt

        active_group.update_filters(gn_filt, sp_filt, True)

    plot_data = methods[plt_idx](filt=filt)

    return [plot_data]


def get_stats_data(aln_obj, stats_idx, active_file_set, active_taxa_set,
                     additional_args, ns=None):
    """
    Given an aln_obj, this function will execute the according method to
    generate plot data

    :param aln_obj: AlignmentObject
    :param stats_idx: string, identifier that maps to an AlignmentObject method
    :return: data for plot production
    """

    # Creating deepcopy to perform changes without impacting the main attribute
    main_aln = deepcopy(aln_obj)
    main_aln.set_database_connections(aln_obj.cur, aln_obj.con)

    # Update alignment object according to active file and taxa sets
    main_aln.update_active_alignments(active_file_set)
    main_aln.update_taxa_names(active_taxa_set)

    # Check if active data sets are not empty. If so, raise an exception
    if main_aln.alignments == OrderedDict() or not main_aln.taxa_names:
        return [EmptyAlignment("Active alignment is empty")]

    # List of gene specific idx. These plots only have one gene for the footer
    gene_specific = ["Pairwise sequence similarity gn"]

    if stats_idx in gene_specific:
        footer = [1, len(main_aln.taxa_names)]
    else:
        footer = [len(main_aln.alignments), len(main_aln.taxa_names)]

    methods = {"Gene occupancy": main_aln.gene_occupancy,
               "Distribution of missing data sp":
                   main_aln.missing_data_per_species,
               "Distribution of missing data":
                   main_aln.missing_data_distribution,
               "Distribution of missing orthologs":
                   main_aln.missing_genes_per_species,
               "Distribution of missing orthologs avg":
                   main_aln.missing_genes_average,
               "Distribution of sequence size":
                   main_aln.average_seqsize_per_species,
               "Distribution of sequence size all": main_aln.average_seqsize,
               "Cumulative distribution of missing genes":
                   main_aln.cumulative_missing_genes,
               "Proportion of nucleotides or residues":
                   main_aln.characters_proportion,
               "Proportion of nucleotides or residues sp":
                   main_aln.characters_proportion_per_species,
               "Pairwise sequence similarity": main_aln.sequence_similarity,
               "Pairwise sequence similarity sp":
                   main_aln.sequence_similarity_per_species,
               "Pairwise sequence similarity gn":
                   main_aln.sequence_similarity_gene,
               "Segregating sites": main_aln.sequence_segregation,
               "Segregating sites sp": main_aln.sequence_segregation_per_species,
               "Segregating sites gn": main_aln.sequence_segregation_gene,
               "Segregating sites prop": main_aln.sequence_segregation,
               "Alignment length/Polymorphism correlation":
                   main_aln.length_polymorphism_correlation,
               "Distribution of taxa frequency":
                   main_aln.taxa_distribution,
               "Allele Frequency Spectrum":
                   main_aln.allele_frequency_spectrum,
               "Allele Frequency Spectrum prop":
                   main_aln.allele_frequency_spectrum,
               "Allele Frequency Spectrum gn":
                   main_aln.allele_frequency_spectrum_gene,
               "Missing data outliers": main_aln.outlier_missing_data,
               "Missing data outliers sp": main_aln.outlier_missing_data_sp,
               "Segregating sites outliers": main_aln.outlier_segregating,
               "Segregating sites outliers sp": main_aln.outlier_segregating_sp,
               "Sequence size outliers sp": main_aln.outlier_sequence_size_sp,
               "Sequence size outliers": main_aln.outlier_sequence_size}

    if additional_args:
        plot_data = methods[stats_idx](ns=ns, **additional_args)
    else:
        plot_data = methods[stats_idx](ns)

    return [plot_data, footer]


__author__ = "Diogo N. Silva"
