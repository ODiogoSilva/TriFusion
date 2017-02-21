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
#

"""
This module deals with the conversion of protein sequences into their
corresponding nucleotide sequences. Since the conversion from protein to DNA
cannot be made without knowing the nucleotide sequence, this module contains
functions that compile and store DNA sequences, convert them into amino acid
sequences and then tries to match them to the original protein sequences.
"""

try:
    from process.error_handling import KillByUser
except ImportError:
    from trifusion.process.error_handling import KillByUser

from os.path import join
import subprocess
import os

dna_map = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
    'TGC': 'C', 'TGT': 'C', 'TGA': '', 'TGG': 'W'}


def translate(sequence):
    """
    Translates a DNA string into an amino acid sequence
    :param sequence: string. DNA sequence
    :return: String. Protein sequence
    """

    sequence = sequence.lower().replace("-", "")
    aa_sequence = ""

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]

        # Check if codon is multiple of 3. If not, place an X (missing data)
        if len(codon) == 3:
            # Check if there is missing data on the codon. If so, place an X
            if "n" in codon:
                aa = ""
            else:
                try:
                    aa = dna_map[codon.upper()]
                except KeyError:
                    aa = ""
        else:
            aa = ""

        aa_sequence += aa

    return aa_sequence


def create_db(f_list, dest="./", ns=None):
    """
    Creates a fasta database file containing the translated protein sequences
    from the cds files. The final transcripts.fas file will be use
    by USEARCH to get matches between the original protein sequences and their
    nucleotide counterparts. A dictionary database will also be created where
    the transcript headers will be associated with the original DNA sequence,
    so that they will be later retrieved
    :param f_list. List, containing the file names of the transcript files
    """

    output_handle = open(join(dest, "transcripts.fas"), "w")
    id_dic = {}

    if ns:
        if ns.stop:
            raise KillByUser("")

        ns.progress = 0
        ns.max_pb = len(f_list)

    for f in f_list:
        handle = open(f)
        seq = ""
        header = ""

        if ns:
            if ns.stop:
                raise KillByUser("")
            ns.progress += 1

        for line in handle:

            if ns:
                if ns.stop:
                    raise KillByUser("")

            if line.startswith(">"):
                if seq != "":
                    aa_seq = translate(seq)
                    output_handle.write(">%s\n%s\n" % (header, aa_seq))
                    id_dic[header] = seq

                header = line.strip()[1:].replace(" ", ";;")
                seq = ""
            else:
                seq += line.strip()

    output_handle.close()

    return id_dic


def create_query(input_list, dest="./"):
    """
    To speed things up, all sequences in the input protein files will be
    concatenated into a single file, which will be used as query in USEARCH.
    :param input_list: List, with file names of the protein files to convert
    """

    f_handle = open(join(dest, "query.fas"), "w")
    query_db = {}

    for f in input_list:
        handle = open(f)
        for line in handle:
            if line.startswith(">"):
                query_db[f].append(line.strip()[1:].replace(" ", ";;"))
                f_handle.write(line.replace(" ", ";;"))
            else:
                f_handle.write(line)

        handle.close()

    f_handle.close()

    return query_db


def create_query_from_dict(protein_dict):
    """
    Analogous to create_query, but begins from a dictionary provided by
    processed group files
    :param protein_dict: dictionary
    """

    query_db = {"group": []}
    f_handle = open("query", "w")

    for cl in protein_dict:
        for seq_id, seq in cl:
            query_db[cl].append(seq_id.replace(" ", ";;"))
            f_handle.write(">%s\n%s\n" % (seq_id.replace(" ", ";;"), seq))

    f_handle.close()

    return query_db


def pair_search(usearch_bin, dest="./"):
    """
    This will use USEARCH to search for translated transcript sequences
    identical to the original protein files
    """

    query_path = join(dest, "query.fas")
    db_path = join(dest, "transcripts.fas")
    out_path = join(dest, "pairs.out")

    subprocess.Popen([usearch_bin,
                      "-usearch_global",
                      query_path,
                      "-db",
                      db_path,
                      "-id",
                      "1",
                      "-maxaccepts",
                      ".9",
                      "-blast6out",
                      out_path]).wait()


def get_pairs(dest="./", ns=None):
    """
    Parses the output of USEARCH and creates a dictionary with the header
    pairs between original protein and transcripts
    """

    file_h = open(join(dest, "pairs.out"))
    pair_db = {}

    if ns:
        if ns.stop:
            raise KillByUser("")
        p = 0
        with open(join(dest, "pairs.out")) as f:
            for p, _ in enumerate(f):
                pass
        ns.max_pb = p + 1
        ns.progress = 0

    for l in file_h:

        if ns:
            if ns.stop:
                raise KillByUser("")
            ns.progress += 1

        fields = l.split("\t")
        pair_db[fields[0]] = fields[1]

    file_h.close()
    return pair_db


def convert_protein_file(pairs, group_obj, id_db, output_dir, shared_ns):
    """
    A given protein file will be converted into their corresponding nucleotide
    sequences using a previously set database using the create_db function
    :return:
    """

    # Create handle for file storing bad sequence headers.
    bad_file = open(join(output_dir, "missed_sequences.log"), "w")

    for line, cl in zip(group_obj.groups(),
                        group_obj.iter_species_frequency()):

        if shared_ns:
            if shared_ns.stop:
                raise KillByUser("")

        if group_obj._get_compliance(cl) == (1, 1):

            line = group_obj._remove_tx(line)

            fields = line.split(":")
            orto_name = fields[0]
            seq_headers = fields[-1].split()

            f_handle = open(join(output_dir, orto_name) + ".fas", "w")

            for h in seq_headers:
                if h in pairs:
                    seq = id_db[pairs[h]]
                    shared_ns.good += 1
                    f_handle.write(">%s\n%s\n" % (h.replace(";;", " "), seq))
                else:
                    shared_ns.missed += 1
                    bad_file.write("{}\t{}\n".format(orto_name, h))


def convert_group(sqldb, cds_file_list, protein_db, group_sequences,
                usearch_bin, output_dir, shared_namespace=None):
    """
    Convenience function that wraps all required operations to convert protein
    to nucleotide files from a Group object
    """

    if shared_namespace:
        shared_namespace.act = "Creating database"
        shared_namespace.missed = 0
        shared_namespace.good = 0
    # Create database
    id_db = create_db(cds_file_list, output_dir, shared_namespace)

    if shared_namespace:
        shared_namespace.act = "Creating query"

        # Kill switch
        if shared_namespace.stop:
            raise KillByUser("")

    # Create query for USEARCH
    group_sequences.retrieve_sequences(sqldb, protein_db, output_dir,
                                       outfile="query.fas",
                                       shared_namespace=shared_namespace)

    if shared_namespace:
        # Kill switch
        if shared_namespace.stop:
            raise KillByUser("")

    # Execute search
    if shared_namespace:
        shared_namespace.act = "Performing search"
    pair_search(usearch_bin, output_dir)

    if shared_namespace:
        # Kill switch
        if shared_namespace.stop:
            raise KillByUser("")

    pair_db = get_pairs(output_dir, ns=shared_namespace)
    # Convert files

    if shared_namespace:
        shared_namespace.act = "Converting to nucleotide"
    convert_protein_file(pair_db, group_sequences, id_db, output_dir,
                         shared_namespace)

    # Remove temporary files
    temp_files = [join(output_dir, "query.fas"),
                  join(output_dir, "transcripts.fas"),
                  join(output_dir, "pairs.out")]

    for f in temp_files:
        os.remove(f)

__author__ = "Diogo N. Silva"
