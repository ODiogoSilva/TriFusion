#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
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
#  Author: Diogo N. Silva
#  Version: 0.1
#  Last update: 11/02/14

"""
This module deals with the conversion of protein sequences into their
corresponding nucleotide sequences. Since the conversion from protein to DNA
cannot be made without knowing the nucleotide sequence, this module contains
functions that compile and store DNA sequences, convert them into amino acid
sequences and then tries to match them to the original protein sequences.
"""

from process.sequence import Alignment
from process.base import Base
import os
import pickle

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


def create_db(f_list):
    """
    This function is meant to be used with DNA fasta files. It creates a
    dictionary database with the protein sequence as a key and the original DNA
    sequence as value
    :param f_list: list, containing the file names of DNA/CDS/Transcript files
    :return: dictionary. Translated protein as key, original DNA as value
    """

    dna_db = {}
    id_db = {}

    for f in f_list:
        print("\rCreating database (Processing file %s)" % f, end="")
        aln = Alignment(f)
        for header, nucl_seq in aln.alignment.items():

            protein_seq = translate(nucl_seq)
            dna_db[protein_seq] = [nucl_seq, header]
            id_db[header] = nucl_seq

    #pickle.dump(dna_db, open("dna_db.key", "wb"))
    pickle.dump(id_db, open("dna_db.key", "wb"))

    return dna_db, id_db


def convert_protein_file(p_file, db, id_db, db_files, outfile_suffix="_dna"):
    """
    A given protein file will be converted into their corresponding nucleotide
    sequences using a previously set database using the create_db function
    :param p_file: string. File name of the protein file
    :param db: dictionary. Database of original DNA sequences and translated
    sequences
    :param outfile_suffix: string. Suffix to append at the end of the p_file
    name. This output file will contained the converted DNA sequences
    :return:
    """

    bad = 0

    for p in p_file:
        aln = Alignment(p)
        output_handle = open(p.split(".")[0] + outfile_suffix + ".fa", "w")

        for key, seq in aln.alignment.items():
            seq = seq.upper().replace("X", "")
            if seq in db:
                output_handle.write(">%s\n%s\n" % (key, db[seq][0]))
            else:
                unique_id = key.split("|")[-1].split()[0]
                res = os.popen("grep %s %s" % (unique_id, " ".join(db_files)))
                try:
                    s = res.read()
                    header = s.split(":")[1].split("\n")[0].strip()[1:]
                    bs = Base()
                    header = bs.rm_illegal(header)
                    output_handle.write(">%s\n%s\n" % (key, id_db[header]))
                except IndexError:
                    bad += 1
                except UnicodeError:
                    bad += 1
    else:
        print("%s sequences not retrieved" % bad)

__author__ = 'diogo'
