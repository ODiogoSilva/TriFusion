#!/usr/bin/python2

"""
Made by Diogo Silva, Fernando Alves
"""

import sqlite3 as lite
import os


def printInparalogsFile (cur, filename):

    cur.execute("select taxon_id, sequence_id_a, sequence_id_b, normalized_score\
        from InParalog\
        order by taxon_id, sequence_id_a, sequence_id_b asc")

    file = open(filename, "w")

    with file:
        while True:

            row = cur.fetchone()
            if row == None:
                break

            file.write(row[1]+"\t"+row[2]+"\t"+str((float(row[3]) * 1000 + .5)/1000)+"\n")

################################################################

def printOrthologsFile (cur, filename):

    cur.execute("select taxon_id_a, taxon_id_b, sequence_id_a, sequence_id_b, normalized_score\
        from Ortholog\
        order by taxon_id_a, taxon_id_b, sequence_id_a, sequence_id_b asc")

    file = open(filename, "w")

    with file:
        while True:

            row = cur.fetchone()
            if row == None:
                break

            file.write(row[2]+"\t"+row[3]+"\t"+str((float(row[4]) * 1000 + .5)/1000)+"\n")

################################################################

def printCoOrthologsFile (cur, filename):

    cur.execute("select taxon_id_a, taxon_id_b, sequence_id_a, sequence_id_b, normalized_score\
        from CoOrtholog\
        order by taxon_id_a, taxon_id_b, sequence_id_a, sequence_id_b asc")

    file = open(filename, "w")

    with file:
        while True:

            row = cur.fetchone()
            if row == None:
                break

            file.write(row[2]+"\t"+row[3]+"\t"+str((float(row[4]) * 1000 + .5)/1000)+"\n")

################################################################

def printMclAbcFile (cur, filename):

    cur.execute("select sequence_id_a, sequence_id_b, normalized_score\
        from InParalog\
        union\
        select sequence_id_a, sequence_id_b, normalized_score\
        from Ortholog\
        union\
        select sequence_id_a, sequence_id_b, normalized_score\
        from CoOrtholog")

    file = open(filename, "w")

    with file:
        while True:

            row = cur.fetchone()
            if row == None:
                break

            file.write(row[0]+"\t"+row[1]+"\t"+str((float(row[2]) * 1000 + .5)/100)+"\n")


def execute(db_dir):
    con = lite.connect(os.path.join(db_dir, "orthoDB.db"))

    with con:

        cur = con.cursor()

        printOrthologsFile(cur, "orthologs.txt")

        printInparalogsFile(cur, "inparalogs.txt")

        printOrthologsFile(cur, "coorthologs.txt")

        printMclAbcFile(cur, "mclInput")

