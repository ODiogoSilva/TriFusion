#!/usr/bin/python2
# -*- coding: utf-8 -*-

import sqlite3 as lite
import os


def execute(db_dir, similar_seqs_file):
    con = lite.connect(os.path.join(db_dir, "orthoDB.db"))

    with con:

        cur = con.cursor()

        file_handle = open(similar_seqs_file)

        for line in file_handle:
            if line.strip() != "" and line:
                f = line.split("\t")
                
                l = (f[0], f[1], f[2], f[3], float(f[4]), float(f[5]),
                    float(f[6]), float(f[7]))
                
                cur.execute("INSERT INTO SimilarSequences VALUES(?, ?, ?, ?, "
                            "?, ? ,?, ?)", l)


if __name__ == "__main__":
    execute(".", "sss_nodups.txt")

__author__ = "Fernando Alves and Diogo N. Silva"
