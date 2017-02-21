#!/usr/bin/python2
# -*- coding: utf-8 -*-

import sqlite3 as lite
import os

##############################################################


def createSimilarSequencesTable(cur):

    cur.execute("CREATE TABLE SimilarSequences (\
        QUERY_ID VARCHAR(60),\
        SUBJECT_ID VARCHAR(60),\
        QUERY_TAXON_ID VARCHAR(40),\
        SUBJECT_TAXON_ID VARCHAR(40),\
        EVALUE_MANT FLOAT,\
        EVALUE_EXP INT,\
        PERCENT_IDENTITY FLOAT,\
        PERCENT_MATCH FLOAT,\
        PRIMARY KEY(QUERY_ID, SUBJECT_ID)\
        )")

    cur.execute("CREATE INDEX ss_qtaxexp_ix\
        ON SimilarSequences(query_id, subject_taxon_id,\
        evalue_exp, evalue_mant,\
        query_taxon_id, subject_id)")

    cur.execute("CREATE INDEX ss_seqs_ix\
        ON SimilarSequences(query_id, subject_id,\
        evalue_exp, evalue_mant, percent_match)")

##############################################################


def createInParalogTable (cur):

    cur.execute("CREATE TABLE InParalog (\
        SEQUENCE_ID_A VARCHAR(60),\
        SEQUENCE_ID_B VARCHAR(60),\
        TAXON_ID VARCHAR(40),\
        UNNORMALIZED_SCORE FLOAT,\
        NORMALIZED_SCORE FLOAT)")

    cur.execute("CREATE INDEX inparalog_seqa_ix\
        ON InParalog(sequence_id_a)")

    cur.execute("CREATE INDEX inparalog_seqb_ix\
        ON InParalog(sequence_id_b)")

##############################################################


def createOrthologTable(cur):

    cur.execute("CREATE TABLE Ortholog (\
        SEQUENCE_ID_A VARCHAR(60),\
        SEQUENCE_ID_B VARCHAR(60),\
        TAXON_ID_A VARCHAR(40),\
        TAXON_ID_B VARCHAR(40),\
        UNNORMALIZED_SCORE FLOAT,\
        NORMALIZED_SCORE FLOAT)")

    cur.execute("CREATE INDEX ortholog_seq_a_ix\
        ON Ortholog(sequence_id_a)")

    cur.execute("CREATE INDEX ortholog_seq_b_ix\
        ON Ortholog (sequence_id_b)")

##############################################################


def createCoOrthologTable(cur):

    cur.execute("CREATE TABLE CoOrtholog (\
        SEQUENCE_ID_A VARCHAR(60),\
        SEQUENCE_ID_B VARCHAR(60),\
        TAXON_ID_A VARCHAR(40),\
        TAXON_ID_B VARCHAR(40),\
        UNNORMALIZED_SCORE FLOAT,\
        NORMALIZED_SCORE FLOAT)")

    cur.execute("CREATE INDEX coortholog_seq_a_ix\
        ON CoOrtholog(sequence_id_a)")

    cur.execute("CREATE INDEX coortholog_seq_b_ix\
        ON CoOrtholog (sequence_id_b)")

##############################################################


def createInterTaxonMatchView(cur):

    cur.execute("DROP VIEW IF EXISTS InterTaxonMatch")

    cur.execute("CREATE VIEW InterTaxonMatch\
        AS SELECT ss.query_id, ss.subject_id, ss.subject_taxon_id,\
        ss.evalue_mant, ss.evalue_exp\
        FROM SimilarSequences ss\
        WHERE ss.subject_taxon_id != ss.query_taxon_id")

##############################################################


def execute(out_dir):

    # Remove any previous DB
    if os.path.exists(os.path.join(out_dir, "orthoDB.db")):
        os.remove(os.path.join(out_dir, "orthoDB.db"))

    con = lite.connect(os.path.join(out_dir, "orthoDB.db"))

    with con:

        cur = con.cursor()

        createSimilarSequencesTable(cur)
        createInParalogTable(cur)
        createOrthologTable(cur)
        createCoOrthologTable(cur)
        createInterTaxonMatchView(cur)

if __name__ == "__main__":
    execute("./")


__author__ = "Fernando Alves and Diogo N. Silva"
