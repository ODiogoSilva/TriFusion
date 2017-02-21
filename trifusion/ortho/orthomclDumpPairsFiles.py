#!/usr/bin/python2
# -*- coding: utf-8 -*-

import sqlite3 as lite
import os

try:
    from process.error_handling import KillByUser
except ImportError:
    from trifusion.process.error_handling import KillByUser


def printInparalogsFile (cur, filename, nm=None):

    cur.execute("select taxon_id, sequence_id_a, sequence_id_b, normalized_score\
        from InParalog\
        order by taxon_id, sequence_id_a, sequence_id_b asc")

    file_fh = open(filename, "w")

    with file_fh:
        while True:

            if nm:
                if nm.stop:
                    raise KillByUser("")

            row = cur.fetchone()
            if row is None:
                break

            file_fh.write("{}\t{}\t{}\n".format(row[1],
                                                row[2],
                                                str((float(row[3]) * 1000 + .5) / 1000)))


################################################################


def printOrthologsFile (cur, filename, nm=None):

    cur.execute("select taxon_id_a, taxon_id_b, sequence_id_a, sequence_id_b, normalized_score\
        from Ortholog\
        order by taxon_id_a, taxon_id_b, sequence_id_a, sequence_id_b asc")

    file_fh = open(filename, "w")

    with file_fh:
        while True:

            if nm:
                if nm.stop:
                    raise KillByUser("")

            row = cur.fetchone()
            if row is None:
                break

            file_fh.write("{}\t{}\t{}\n".format(row[2],
                                                row[3],
                                                str((float(row[4]) * 1000 + .5) / 1000)))

################################################################


def printCoOrthologsFile (cur, filename, nm=None):

    cur.execute("select taxon_id_a, taxon_id_b, sequence_id_a, sequence_id_b, normalized_score\
        from CoOrtholog\
        order by taxon_id_a, taxon_id_b, sequence_id_a, sequence_id_b asc")

    file_fh = open(filename, "w")

    with file_fh:
        while True:

            if nm:
                if nm.stop:
                    raise KillByUser("")

            row = cur.fetchone()
            if row is None:
                break

            file_fh.write("{}\t{}\t{}\n".format(row[2],
                                                row[3],
                                                str((float(row[4]) * 1000 + .5) / 1000)))

################################################################


def printMclAbcFile (cur, filename, nm=None):

    cur.execute("select sequence_id_a, sequence_id_b, normalized_score\
        from InParalog\
        union\
        select sequence_id_a, sequence_id_b, normalized_score\
        from Ortholog\
        union\
        select sequence_id_a, sequence_id_b, normalized_score\
        from CoOrtholog")

    file_fh = open(filename, "w")

    with file_fh:
        while True:

            if nm:
                if nm.stop:
                    raise KillByUser("")

            row = cur.fetchone()
            if row is None:
                break

            file_fh.write("{}\t{}\t{}\n".format(row[0],
                                                row[1],
                                                str((float(row[2]) * 1000 + .5) / 1000)))


def execute(db_dir, dest, nm=None):
    con = lite.connect(os.path.join(db_dir, "orthoDB.db"))

    # Set up progression information
    if nm:
        if nm.stop:
            raise KillByUser("")
        nm.total = 4
        nm.counter = 0

    with con:

        cur = con.cursor()

        if nm:
            if nm.stop:
                raise KillByUser("")
            nm.counter = 1

        printOrthologsFile(cur, os.path.join(dest, "backstage_files",
                                             "orthologs.txt"), nm=nm)
        if nm:
            if nm.stop:
                raise KillByUser("")
            nm.counter = 2

        printInparalogsFile(cur, os.path.join(dest, "backstage_files",
                                              "inparalogs.txt"), nm=nm)
        if nm:
            if nm.stop:
                raise KillByUser("")
            nm.counter = 3
        printCoOrthologsFile(cur, os.path.join(dest, "backstage_files",
                                               "coorthologs.txt"), nm=nm)
        if nm:
            if nm.stop:
                raise KillByUser("")
            nm.counter = 4
        printMclAbcFile(cur, os.path.join(dest, "backstage_files",
                                          "mclInput"), nm=nm)

if __name__ == "__main__":
    execute(".", ".")

__author__ = "Fernando Alves and Diogo N. Silva"
