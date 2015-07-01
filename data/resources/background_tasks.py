__author__ = 'diogo'

from process import data
from process.sequence import AlignmentList
from ortho import orthomcl_pipeline as ortho_pipe
from ortho import OrthomclToolbox as OrthoTool

from os.path import join, basename
from copy import deepcopy
import logging


def load_proc(aln_list, file_list, nm):
    try:
        if aln_list:
            aln_list.add_alignment_files(file_list, shared_namespace=nm)
            aln_obj = aln_list
        else:
            aln_obj = AlignmentList(file_list, shared_namespace=nm)
        nm.alns = aln_obj
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


def process_dispatch(nm, temp_dir, proteome_files, protein_min_len, 
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
                                    bin_path=join(cur_dir, "ortho",
                                                  "orthomclBlastParser"))
        if nm.k:
            ortho_pipe.remove_duplicate_entries()
        if nm.k:
            nm.t = "Loading USEARCH output to database"
            nm.c = 6
            ortho_pipe.load_blast(temp_dir)
        if nm.k:
            nm.t = "Obtaining Pairs"
            nm.c = 7
            ortho_pipe.pairs(temp_dir)
        if nm.k:
            ortho_pipe.dump_pairs(temp_dir)
        if nm.k:
            nm.t = "Running MCL"
            nm.c = 8
            ortho_pipe.mcl(mcl_inflation)
        if nm.k:
            nm.t = "Dumping groups"
            nm.c = 9
            ortho_pipe.mcl_groups(mcl_inflation, ortholog_prefix, "1000",
                                  group_prefix, bin_path=join(cur_dir, "ortho",
                                                "orthomclMclToGroups"))
        if nm.k:
            nm.t = "Filtering group files"
            nm.c = 10
            stats, groups_obj = ortho_pipe.export_filtered_groups(
                                                    mcl_inflation,
                                                    group_prefix,
                                                    orto_max_gene,
                                                    orto_min_sp, sqldb,
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


def process_execution(**kwargs):
    """
    Process execution function
    :param ns: Namespace object
    """

    try:
        write_aln = {}

        #####
        # Perform operations
        #####
        # Setting the alignment to use.
        # Update active file set of the alignment object
        aln_object = update_active_fileset(kwargs["aln_list"],
                                           kwargs["file_set_name"],
                                           kwargs["file_list"],
                                           kwargs["file_groups"])

        # Update active taxa set of the alignment object
        aln_object = update_active_taxaset(aln_object, kwargs["taxa_set_name"],
                                           kwargs["active_taxa_list"],
                                           kwargs["taxa_groups"])

        kwargs["ns"].proc_files = len(aln_object.alignments)

        # Filtering - This should be the first operation to be performed
        if kwargs["secondary_operations"]["filter"]:
            kwargs["ns"].msg = "Filtering alignment(s)"
            if kwargs["secondary_options"]["filter_file"]:
                filtered_aln_obj = deepcopy(aln_object)
            # Check if a minimum taxa representation was specified
            if kwargs["secondary_options"]["gap_filter"]:
                if kwargs["missing_filter_settings"][2]:
                    try:
                        filtered_aln_obj.filter_min_taxa(
                        kwargs["missing_filter_settings"][2])
                    except NameError:
                        aln_object.filter_min_taxa(
                            kwargs["missing_filter_settings"][2])
            # Filter by taxa
            if kwargs["secondary_options"]["taxa_filter"]:
                # Get taxa list from taxa groups
                taxa_list = kwargs["taxa_groups"]\
                    [kwargs["taxa_filter_settings"][1]]
                try:
                    filtered_aln_obj.filter_by_taxa(
                        kwargs["taxa_filter_settings"][0], taxa_list)
                except NameError:
                    aln_object.filter_by_taxa(kwargs["taxa_filter_settings"][0],
                                              taxa_list)
            # Filter codon positions
            if kwargs["secondary_options"]["codon_filter"]:
                try:
                    filtered_aln_obj.filter_codon_positions(
                    kwargs["codon_filter_settings"])
                except NameError:
                    aln_object.filter_codon_positions(
                        kwargs["codon_filter_settings"])
            # Filter missing data
            if kwargs["secondary_options"]["gap_filter"]:
                try:
                    filtered_aln_obj.filter_missing_data(
                        kwargs["missing_filter_settings"][0],
                        kwargs["missing_filter_settings"][1])
                except NameError:
                    aln_object.filter_missing_data(
                    kwargs["missing_filter_settings"][0],
                    kwargs["missing_filter_settings"][1])
            try:
                write_aln[kwargs["output_file"] + "_filtered"] = \
                    filtered_aln_obj
            except NameError:
                pass

        # Concatenation
        if kwargs["main_operations"]["concatenation"]:
            kwargs["ns"].msg = "Concatenating"
            aln_object = aln_object.concatenate()
            # In case filtered alignments are going to be saved in a different
            # file, concatenate them as well
            if write_aln:
                for name, aln in write_aln.items():
                    write_aln[name] = aln.concatenate()
            # Concatenation of ZORRO files
            if kwargs["secondary_options"]["zorro"]:
                kwargs["ns"].msg = "Concatenating ZORRO files"
                zorro_data = data.Zorro(aln_object, kwargs["zorro_suffix"])
                zorro_data.write_to_file(kwargs["output_file"])

        # Reverse concatenation
        if kwargs["main_operations"]["reverse_concatenation"]:
            kwargs["ns"].msg = "Reverse concatenating"
            if not kwargs["use_app_partitions"]:
                partition_obj = data.Partitions()
                # In case the partitions file is badly formatted or invalid, the
                # exception will be returned by the read_from_file method.
                er = partition_obj.read_from_file(kwargs["partitions_file"])
                aln_object = aln_object.retrieve_alignment(kwargs["rev_infile"])
                aln_object.set_partitions(partition_obj)
            aln_object = aln_object.reverse_concatenate()

        # Collapsing
        if kwargs["secondary_operations"]["collapse"]:
            kwargs["ns"].msg = "Collapsing alignment(s)"
            if kwargs["secondary_options"]["collapse_file"]:
                collapsed_aln_obj = deepcopy(aln_object)
                collapsed_aln_obj.collapse(haplotype_name=kwargs["hap_prefix"],
                                           haplotypes_file="_collapsed")
                write_aln[kwargs["output_file"] + "_collapsed"] = \
                    collapsed_aln_obj
            else:
                aln_object.collapse(haplotype_name=kwargs["hap_prefix"])

        # Gcoder
        if kwargs["secondary_operations"]["gcoder"]:
            kwargs["ns"].msg = "Coding gaps"
            if kwargs["secondary_options"]["gcoder_file"]:
                gcoded_aln_obj = deepcopy(aln_object)
                gcoded_aln_obj.code_gaps()
                write_aln[kwargs["output_file"] + "_coded"] = gcoded_aln_obj
            else:
                aln_object.code_gaps()
                
        # The output file(s) will only be written after all the required
        # operations have been concluded. The reason why there are two "if"
        # statement for "concatenation" is that the input alignments must be
        # concatenated before any other additional operatiokwargs["ns"]. If the
        # first if statement did not exist, then all additional options would
        # have to be manually written for both "conversion" and "concatenation".
        #  As it is, when "concatenation", the aln_obj is firstly converted
        # into the concatenated alignment, and then all additional
        # operations are conducted in the same aln_obj
        write_aln[kwargs["output_file"]] = aln_object
        kwargs["ns"].msg = "Writhing output"
        
        if kwargs["main_operations"]["concatenation"]:
            for name, obj in write_aln.items():
                obj.write_to_file(kwargs["output_formats"], name,
                        interleave=kwargs["secondary_options"]["interleave"],
                        partition_file=kwargs["create_partfile"],
                        use_charset=kwargs["use_nexus_partitions"],
                        phy_truncate_names=kwargs["use_nexus_partitions"])
        else:
            for name, obj in write_aln.items():
                name = name.replace(kwargs["output_file"], "")
                obj.write_to_file(kwargs["output_formats"],
                        output_suffix=name,
                        interleave=kwargs["secondary_options"]["interleave"],
                        partition_file=kwargs["create_partfile"],
                        output_dir=kwargs["output_dir"],
                        use_charset=kwargs["use_nexus_partitions"],
                        phy_truncate_names=kwargs["use_nexus_partitions"])

    except:
        # Log traceback in case any unexpected error occurs. See
        # self.log_file for whereabouts of the traceback
        logging.exception("Unexpected exit in Process execution")
        kwargs["ns"].exception = True


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
