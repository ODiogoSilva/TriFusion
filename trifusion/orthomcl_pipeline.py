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


import warnings

# Suppress import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    import os
    import sys
    import time
    import codecs
    import subprocess
    import shutil
    import traceback
    import argparse
    from os.path import abspath, join, basename

    try:
        from process.base import print_col, GREEN, RED
        from ortho import OrthomclToolbox as OT
        import ortho.orthomclInstallSchema as install_sqlite
        import ortho.orthomclPairs as make_pairs_sqlite
        import ortho.orthomclDumpPairsFiles as dump_pairs_sqlite
        import ortho.orthomclFilterFasta as FilterFasta
        import ortho.orthomclBlastParser as BlastParser
        import ortho.orthomclMclToGroups as MclGroups
        from ortho.error_handling import *
        from process.error_handling import KillByUser
    except ImportError:
        from trifusion.process.base import print_col, GREEN, RED
        from trifusion.ortho import OrthomclToolbox as OT
        import trifusion.ortho.orthomclInstallSchema as install_sqlite
        import trifusion.ortho.orthomclPairs as make_pairs_sqlite
        import trifusion.ortho.orthomclDumpPairsFiles as dump_pairs_sqlite
        import trifusion.ortho.orthomclFilterFasta as FilterFasta
        import trifusion.ortho.orthomclBlastParser as BlastParser
        import trifusion.ortho.orthomclMclToGroups as MclGroups
        from trifusion.ortho.error_handling import *
        from trifusion.process.error_handling import KillByUser


def install_schema(db_dir):
    """
    Install the schema for the mySQL database

    :param db_dir: string, directory for the sqlite database
    """

    print_col("Creating sqlite database", GREEN, 1)
    install_sqlite.execute(db_dir)


def check_unique_field(proteome_file, verbose=False, nm=None):
    """
    Checks the original proteome file for a field in the fasta header
    that is unique to all sequences
    """

    # Some files may have utf8 encoding problems so I used codecs here
    file_handle = codecs.open(proteome_file, "r", "cp1252")
    header_list = []

    header = ""
    for line in file_handle:

        if nm:
            if nm.stop:
                raise KillByUser("")

        if line.startswith(">"):
            header = line[1:].strip()
            # Store header in list format
            header_list.append(header.split("|"))

    # Get the size of the header fields
    header_field_size = len(header.split("|"))

    for i in range(header_field_size):

        if nm:
            if nm.stop:
                raise KillByUser("")

        temp_list = []
        for header in header_list:
            temp_list.append(header[i])

        if len(temp_list) == len(set(temp_list)) and len(set(temp_list)) ==\
                len(header_list):

            # The orthoMCL program uses an index starting from 1, so the +1 is
            #  a necessary adjustment
            if verbose:
                print_col("\t Using unique header field {}".format(i), GREEN, 1)
            return i

    # Ideally, a unique field should be found before this code. If not, raise
    #  exception
    raise NoUniqueField("The proteome file {} has no unique field".format(
        os.path.basename(proteome_file)))


def prep_fasta(proteome_file, code, unique_id, verbose=False, nm=None):

    if verbose:
        print_col("\t Preparing file for USEARCH", GREEN, 1)

    # Storing header list to check for duplicates
    header_list = []

    # Storing dictionary with header and sequence for later use
    seq_storage = {}

    # Will prevent writing
    lock = True

    # File handles
    file_in = open(proteome_file)
    file_out = open(proteome_file.split(".")[0] + "_mod.fas", "w")

    for line in file_in:

        if nm:
            if nm.stop:
                raise KillByUser("")

        if line.startswith(">"):
            if line not in header_list:
                fields = line.split("|")
                unique_str = fields[unique_id].replace(" ", "_")
                seq_storage["%s|%s" % (code, unique_str)] = ""
                header_list.append(line)
                file_out.write(">%s|%s\n" % (code, unique_str))
                lock = True
            else:
                lock = False
        elif lock:
            seq_storage["%s|%s" % (code, unique_str)] += line.strip()
            file_out.write(line)

    # Close file handles:
    file_in.close()
    file_out.close()

    return seq_storage


def adjust_fasta(file_list, dest, nm=None):

    print_col("Adjusting proteome files", GREEN, 1)

    # Create compliant fasta directory
    cf_dir = join(dest, "backstage_files", "compliantFasta")
    if not os.path.exists(cf_dir):
        os.makedirs(cf_dir)
    else:
        for f in os.listdir(cf_dir):
            os.remove(join(cf_dir, f))

    # Setup progress information
    if nm:
        if nm.stop:
            KillByUser("")
        # Get total number of files for total progress
        nm.total = len(file_list)
        nm.counter = 0

    for proteome in file_list:
        # Get code for proteome
        code_name = proteome.split(os.path.sep)[-1].split(".")[0]

        if nm:
            if nm.stop:
                raise KillByUser("")
            nm.counter += 1
            nm.msg = "Adjusting file {}".format(basename(proteome))

        # Check the unique ID field
        unique_id = check_unique_field(proteome, True, nm)

        # Adjust fasta
        # stg = prep_fasta(proteome, code_name, unique_id)
        prep_fasta(proteome, code_name, unique_id, nm)

        protome_file_name = proteome.split(os.path.sep)[-1].split(".")[0] + \
                            ".fasta"

        shutil.move(proteome.split(".")[0] + "_mod.fas",
                    join(cf_dir, protome_file_name))


def filter_fasta(min_len, max_stop, db, dest, nm=None):

    print_col("Filtering proteome files", GREEN, 1)

    cp_dir = join(dest, "backstage_files", "compliantFasta")

    FilterFasta.orthomcl_filter_fasta(cp_dir, min_len, max_stop, db, dest, nm)


def allvsall_usearch(goodproteins, evalue, dest, cpus, usearch_outfile,
                     usearch_bin="usearch", nm=None):

    print_col("Perfoming USEARCH All-vs-All (may take a while...)", GREEN, 1)

    # FNULL = open(os.devnull, "w")
    usearch_cmd = [usearch_bin,
                   "-ublast",
                   join(dest, "backstage_files", goodproteins),
                   "-db",
                   join(dest, "backstage_files", goodproteins),
                   "-blast6out",
                   join(dest, "backstage_files", usearch_outfile),
                   "-evalue", str(evalue),
                   "--maxaccepts",
                   "0",
                   "-threads",
                   str(cpus)]

    if nm:
        # The subprocess.Popen handler cannot be passed directly in Windows
        # due to pickling issues. So I pass the pid of the process instead.
        subp = subprocess.Popen(usearch_cmd)
        nm.subp = subp.pid
        subp.wait()
        nm.subp = None
    else:
        _ = subprocess.Popen(usearch_cmd).wait()


def blast_parser(usearch_ouput, dest, db_dir, nm):

    print_col("Parsing BLAST output", GREEN, 1)

    BlastParser.orthomcl_blast_parser(
        join(dest, "backstage_files", usearch_ouput),
        join(dest, "backstage_files", "compliantFasta"),
        db_dir,
        nm)


def pairs(db_dir, nm=None):

    print_col("Finding pairs for orthoMCL", GREEN, 1)

    make_pairs_sqlite.execute(db_dir, nm=nm)


def dump_pairs(db_dir, dest, nm=None):

    print_col("Dump files from the database produced by the orthomclPairs "
              "program", GREEN, 1)

    dump_pairs_sqlite.execute(db_dir, dest, nm=nm)


def mcl(inflation_list, dest, mcl_file="mcl", nm=None):

    print_col("Running mcl algorithm", GREEN, 1)
    mcl_input = join(dest, "backstage_files", "mclInput")
    mcl_output = join(dest, "backstage_files", "mclOutput_")

    for val in inflation_list:

        mcl_cmd = [mcl_file,
                   mcl_input,
                   "--abc",
                   "-I",
                   val,
                   "-o",
                   mcl_output + val.replace(".", "")]

        if nm:
            # The subprocess.Popen handler cannot be passed directly in Windows
            # due to pickling issues. So I pass the pid of the process instead.
            subp = subprocess.Popen(mcl_cmd)
            nm.subp = subp.pid
            subp.wait()
            nm.subp = None
        else:
            _ = subprocess.Popen(mcl_cmd).wait()


def mcl_groups(inflation_list, mcl_prefix, start_id, group_file, dest,
                nm=None):

    print_col("Dumping groups", GREEN, 1)

    # Create a results directory
    results_dir = join(dest, "Orthology_results")
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    mcl_output = join(dest, "backstage_files", "mclOutput_")

    if nm:
        if nm.stop:
            raise KillByUser("")
        nm.total = len(inflation_list)
        nm.counter = 0

    for val in inflation_list:

        if nm:
            if nm.stop:
                raise KillByUser("")
            nm.counter += 1

        MclGroups.mcl_to_groups(
            mcl_prefix,
            start_id,
            mcl_output + val.replace(".", ""),
            os.path.join(results_dir, group_file + "_" + str(val) + ".txt"),
            nm=nm)


def export_filtered_groups(inflation_list, group_prefix, gene_t, sp_t, sqldb,
                           db, tmp_dir, dest, nm=None):

    print_col("Exporting filtered groups to protein sequence files", GREEN, 1)

    stats_storage = {}
    groups_obj = OT.MultiGroupsLight(tmp_dir)

    if nm:
        if nm.stop:
            raise KillByUser("")

    for val in inflation_list:
        # Create a directory that will store the results for the current
        # inflation value
        inflation_dir = join(dest, "Orthology_results", "Inflation%s" % val)
        if not os.path.exists(inflation_dir):
            os.makedirs(inflation_dir)

        group_file = join(dest, "Orthology_results",
                          group_prefix + "_%s.txt" % val)

        # Create Group object
        group_obj = OT.GroupLight(group_file, gene_t, sp_t)
        # Add group to the MultiGroups object
        groups_obj.add_group(group_obj)
        # Export filtered groups and return stats to present in the app
        stats = group_obj.basic_group_statistics()
        # Retrieve fasta sequences from the filtered groups
        group_obj.retrieve_sequences(sqldb, db,
                                     dest=join(inflation_dir, "Orthologs"),
                                     shared_namespace=nm)
        # os.remove(sqldb)
        stats_storage[val] = stats

    return stats_storage, groups_obj


def main():

    # The inclusion of the argument definition in main, makes it possible to
    # import this file as a module and not triggering argparse. The
    # alternative of using a if __name__ == "__main__" statement does not
    # work well with the entry_points parameter of setup.py, since they call
    # the main function but do nothing inside said statement.
    parser = argparse.ArgumentParser(description="Command line interface for "
        "TriFusion Orthology search module")

    parser.add_argument("-in", dest="infile", type=str,
                        required=True, help="Provide the path "
                        "to the directory containing the proteome files")

    # Execution modes
    exec_modes = parser.add_argument_group("Execution modes")
    exec_modes.add_argument("-n", action="store_const", const=True,
                            dest="normal",
                            help="Complete run of the pipeline")
    exec_modes.add_argument("-a", action="store_const", const=True,
                            dest="adjust",
                            help="Only adjust proteome fasta files")
    exec_modes.add_argument("-na", action="store_const", const=True,
                            dest="no_adjust",
                            help="Complete run of the pipeline without "
                                 "adjusting fasta files")

    # Input formatting
    input_format = parser.add_argument_group("Input formatting")
    input_format.add_argument("-d", action="store_const", const=True,
                              dest="code", help="Do not convert input proteome"
                              " file names because the file names are already "
                              "in code (e.g. Homo_sapiens.fas -> HoSap.fas")
    input_format.add_argument("-sep", dest="separator", help="Specify the "
                              "separator in the input files (e.g. '_' is the"
                              " separator in 'Homo_sapiens.fas'). This "
                              "parameter is ignored if the '-d' option is set")

    # Search options
    search_opts = parser.add_argument_group("Ortholog search options")
    search_opts.add_argument("--min-length", dest="min_length", type=int,
                             default=10, help="Set minimum length allowed "
                             "for protein sequences (default is '%(default)s')")
    search_opts.add_argument("--max-stop", dest="max_stop", type=int,
                             default=20, help="Set maximum percentage of "
                             "stop codons in protein sequences (default is "
                             "'%(default)s')")
    search_opts.add_argument("--db", dest="database",
                             default="goodProteins", help="Name of search "
                             "database (default is '%(default)s')")
    search_opts.add_argument("--search-out", dest="search_out",
                             default="AllVsAll.out", help="Name of the "
                             "search output file containing the All-vs-All "
                             "protein comparisons")
    search_opts.add_argument("-evalue", dest="evalue", default=1E-5,
                             help="Set the e-value cut off for search "
                             "operation (default is '%(default)s')")
    search_opts.add_argument("-inflation", dest="inflation", nargs="+",
                             default=["3"],
                             choices=[str(x) for x in xrange(1, 6)],
                             help="Set inflation values for ortholog group"
                             " clustering. Multiple values may be provided "
                             "but values are limited to the range [1, 5]")

    # Output options
    output_opts = parser.add_argument_group("Output options")
    output_opts.add_argument("-o", dest="output_dir", default=os.getcwd(),
                             help="Output directory")
    output_opts.add_argument("-prefix", dest="prefix", default="Ortholog",
                             help="Set the prefix name for each ortholog "
                             "cluster (default is '%(default)s')")
    output_opts.add_argument("-id", dest="id_num", type=int, default=1,
                             help="Set the starting number for the ortholog "
                             "clusters (default is '%(default)s')")
    output_opts.add_argument("--groups-file", dest="groups_file",
                             default="groups", help="Set the name of the "
                             "group files from the output of MCL (default is "
                             "'%(default)s')")
    output_opts.add_argument("--min-species", dest="min_sp", default=1,
                             type=float, help="Set the minimum number of "
                             "species required for an ortholog cluster to be "
                             "converted into protein sequence. This option "
                             "will only affect the protein sequence files, "
                             "not the group file output.")
    output_opts.add_argument("--max-gene-copy", dest="max_gn", default=100,
                             type=int, help="Set the maximum number of gene "
                             "copies from the same taxon for each ortholog "
                             "cluster. This option will only affect the "
                             "protein sequence files, not the group file "
                             "output.")

    # Miscellaneous options
    misc_options = parser.add_argument_group("Miscellaneous options")
    misc_options.add_argument("-np", dest="cpus", default=1, help="Number of "
                              "CPUs to be used during search operation ("
                              "default is '%(default)s')")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    arg = parser.parse_args()

    # Crete temp directory
    tmp_dir = join(os.getcwd(), ".tmp")
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    print_col("Executing OrthoMCL pipeline at %s %s" % (
        time.strftime("%d/%m/%Y"), time.strftime("%I:%M:%S")), GREEN, 1)

    try:
        start_time = time.time()

        # Arguments
        input_dir = arg.infile
        output_dir = arg.output_dir
        # name_separator = arg.separator
        min_length = arg.min_length
        max_percent_stop = arg.max_stop
        database_name = join(os.getcwd(), output_dir, "backstage_files",
                             arg.database)
        usearch_out_name = arg.search_out
        evalue_cutoff = arg.evalue
        cpus = arg.cpus
        inflation = arg.inflation
        prefix = arg.prefix
        start_id = arg.id_num
        groups_file = arg.groups_file
        min_sp = arg.min_sp
        max_gn = arg.max_gn

        # Get proteome files
        if not os.path.exists(input_dir):
            print_col("The input directory %s does not exist. Exiting." %
                      input_dir, RED, 1)

        proteome_files = [abspath(join(input_dir, x)) for x in os.listdir(
            input_dir)]

        # Create and change working directory
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        # os.chdir(output_dir)

        # Create directory that will store intermediate files during orthology
        # search
        int_dir = join(output_dir, "backstage_files")
        if not os.path.exists(int_dir):
            os.makedirs(int_dir)
        os.chdir(int_dir)

        if arg.normal:
            install_schema(tmp_dir)
            adjust_fasta(proteome_files, output_dir)
            filter_fasta(min_length, max_percent_stop, database_name,
                         output_dir)
            allvsall_usearch(database_name, evalue_cutoff, output_dir, cpus,
                             usearch_out_name)
            blast_parser(usearch_out_name, output_dir, tmp_dir, None)
            pairs(tmp_dir)
            dump_pairs(tmp_dir, output_dir)
            mcl(inflation, output_dir)
            mcl_groups(inflation, prefix, start_id, groups_file, output_dir)
            export_filtered_groups(inflation, groups_file, max_gn, min_sp,
                                   "tmp.sql3", database_name, tmp_dir,
                                   output_dir)

        elif arg.adjust:
            adjust_fasta(proteome_files, output_dir)

        elif arg.no_adjust:
            install_schema(tmp_dir)
            filter_fasta(min_length, max_percent_stop, database_name,
                         output_dir)
            allvsall_usearch(database_name, evalue_cutoff, output_dir, cpus,
                             usearch_out_name)
            blast_parser(usearch_out_name, output_dir, tmp_dir, None)
            pairs(tmp_dir)
            dump_pairs(tmp_dir, output_dir)
            mcl(inflation, output_dir)
            mcl_groups(inflation, prefix, start_id, groups_file, output_dir)
            export_filtered_groups(inflation, groups_file, max_gn, min_sp,
                                   "tmp.sql3", database_name, tmp_dir,
                                   output_dir)

        print_col("OrthoMCL pipeline execution successfully completed in %s "
                  "seconds" % (round(time.time() - start_time, 2)), GREEN, 1)

        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
    except Exception as e:
        print(e.message)
        traceback.print_exc()
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)
        print_col("Program exited with errors!", RED, 1)


if __name__ == "__main__":

    main()


__author__ = "Diogo N. Silva"
