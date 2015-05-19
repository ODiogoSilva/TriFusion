#!/usr/bin/python3

"""
Made by Diogo Silva, Fernando Alves
"""

import sqlite3 as lite


def parse_seqs(file_name):

    file_handle = open(file_name)

    lines = []

    for line in file_handle:
        if line.strip() != "":
            f = line.split("\t")
            l = (f[0], f[1], f[2], f[3], float(f[4]), int(f[5]), float(f[6]),
                 float(f[7]))
            lines.append(l)

    lines = tuple(lines)

    return lines

##############################################################


def execute(similar_seqs_file):
    con = lite.connect("orthoDB.db")

    with con:

        cur = con.cursor()

        lines = parse_seqs(similar_seqs_file)

        cur.executemany("INSERT INTO SimilarSequences VALUES(?, ?, ?, ?, "
                        "?, ? ,?, ?)", lines)

