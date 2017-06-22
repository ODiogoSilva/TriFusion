#!/usr/bin/python2
# -*- coding: utf-8 -*-

import sqlite3 as lite
import os
import math

try:
    from process.error_handling import KillByUser
except ImportError:
    from trifusion.process.error_handling import KillByUser


"""my @steps = ( # Common
       ['updateMinimumEvalueExponent'],
       ['bestQueryTaxonScore'],
       ['qtscore_ix'],

       # Ortholog
       ['bestHit'],
       ['best_hit_ix'],
       ['ortholog', ["drop table BestHit"]],
       ['orthologTaxon'],
       ['orthologAvg'],
       ['orthologAvgIndex'],
       ['orthologsNormalization', ["drop table OrthologAvgScore", "drop table OrthologTaxon", "drop table OrthologTemp"]],

       # InParalog
       ['bestInterTaxonScore', ["drop table BestQueryTaxonScore"]],
       ['bis_uids_ix'],
       ['uniqueSimSeqsQueryId'],
       ['ust_qids_ix'],
       ['betterHit', ["drop table BestInterTaxonScore", "drop table UniqSimSeqsQueryId"]],
       ['better_hit_ix'],
       ['inParalog', ["drop table BetterHit"]],
       ['inParalogTaxonAvg'],
       ['orthologUniqueId'],
       ['orthologUniqueIdIndex'],
       ['inplgOrthTaxonAvg', ["drop table OrthologUniqueId"]],
       ['inParalogAvg',["drop table InParalogTaxonAvg", "drop table InplgOrthTaxonAvg"]],
       ['inParalogAvgIndex'],
       ['inParalogsNormalization', ["drop table InParalogAvgScore", "drop table InParalogTemp"]],

       # CoOrtholog
       ['inParalog2Way'],
       ['in2a_ix'],
       ['in2b_ix'],
       ['ortholog2Way'],
       ['ortholog2WayIndex'],
       ['inplgOrthoInplg'],
       ['inParalogOrtholog'],
       ['coOrthologCandidate', ["drop table Ortholog2Way", "drop table InParalog2Way", "drop table InplgOrthoInplg", "drop table InParalogOrtholog"]],
       ['coOrthologNotOrtholog', ["drop table CoOrthologCandidate"]],
       ['coOrthologNotOrthologIndex'],
       ['coOrtholog', ["drop table CoOrthNotOrtholog"]],
       ['coOrthologTaxon'],
       ['coOrthologAvg'],
       ['coOrthologAvgIndex'],
       ['coOrthologsNormalization', ["drop table CoOrthologAvgScore", "drop table CoOrthologTaxon", "drop table CoOrthologTemp"]],
       ['cleanall', ["truncate table InParalog", "truncate table Ortholog", "truncate table CoOrtholog"]],
      );
"""
################################################################################
############################### Auxiliar    #################################
################################################################################

def log(value):
    return math.log10(value)

def orthologTaxonSub (cur, co):

    #assuming in perl a var is true if not ""
    coCaps = "" if co == "" else "Co"
    t1 = coCaps + 'OrthologTaxon'
    t2 = coCaps + 'OrthologTemp'

    cur.execute("create table %s as\
        select case\
        when taxon_id_a < taxon_id_b\
        then taxon_id_a\
        else taxon_id_b\
        end as smaller_tax_id,\
        case\
        when taxon_id_a < taxon_id_b\
        then taxon_id_b\
        else taxon_id_a\
        end as bigger_tax_id,\
        unnormalized_score\
        from %s " % (t1, t2))

#from (?) ", (coCaps + 'OrthologTaxon', coCaps + 'OrthologTemp'))

def normalizeOrthologsSub (cur, co, table):

    #assuming in perl a var is true if not ""
    coCaps = "" if co == "" else "Co"
    co = "o" if co == "" else "coO"
    t1 = coCaps + 'OrthologAvgScore'
    t2 = coCaps + 'OrthologTaxon'

    cur.execute("create table %s as\
        select smaller_tax_id, bigger_tax_id, avg(unnormalized_score) avg_score\
        from %s\
        group by smaller_tax_id, bigger_tax_id" % (t1, t2))

################################################################
    t1 = co+'orthoAvg_ix'
    t2 = coCaps+'OrthologAvgScore'

    cur.execute("create unique index %s on %s(smaller_tax_id,bigger_tax_id,avg_score)" % (t1, t2))

################################################################

    t1 = coCaps + 'OrthologTemp'
    t2 = coCaps + 'OrthologAvgScore'

    cur.execute("insert into %s (sequence_id_a, sequence_id_b, taxon_id_a, taxon_id_b, unnormalized_score, normalized_score)\
        select ot.sequence_id_a, ot.sequence_id_b, ot.taxon_id_a, ot.taxon_id_b, ot.unnormalized_score, ot.unnormalized_score/a.avg_score\
        from %s ot, %s a\
        where min(ot.taxon_id_a, ot.taxon_id_b) = a.smaller_tax_id\
        and max(ot.taxon_id_a, ot.taxon_id_b) = a.bigger_tax_id" % (table, t1, t2))


################################################################################
############################### Common tables #################################
################################################################################
def commonTempTables (cur):


    cur.execute("select min(evalue_exp)\
        from SimilarSequences\
        where evalue_mant != 0")

    tup = cur.fetchone()
    minEvalueExp = tup[0] - 1

    cur.execute("update SimilarSequences\
        set evalue_exp = ?\
        where evalue_exp = 0 and evalue_mant = 0", (minEvalueExp,))

##########################################################################

    cur.execute("create table BestQueryTaxonScore as\
        select im.query_id as query_id, im.subject_taxon_id as subject_taxon_id, low_exp.evalue_exp as evalue_exp, min(im.evalue_mant) as evalue_mant\
        from InterTaxonMatch im,\
        (select query_id, subject_taxon_id, min(evalue_exp) as evalue_exp\
        from InterTaxonMatch\
        group by query_id, subject_taxon_id) low_exp\
        where im.query_id = low_exp.query_id\
        and im.subject_taxon_id = low_exp.subject_taxon_id\
        and im.evalue_exp = low_exp.evalue_exp\
        group by im.query_id, im.subject_taxon_id, low_exp.evalue_exp")

################################################################################

    cur.execute("create unique index qtscore_ix on BestQueryTaxonScore(query_id, subject_taxon_id, evalue_exp, evalue_mant)")

################################################################################
############################### Orthologs #####################################
################################################################################
def orthologs(cur):

    cur.execute("create table BestHit as\
        select s.query_id, s.subject_id,\
        s.query_taxon_id, s.subject_taxon_id,\
        s.evalue_exp, s.evalue_mant\
        from SimilarSequences s, BestQueryTaxonScore cutoff\
        where s.query_id = cutoff.query_id\
        and s.subject_taxon_id = cutoff.subject_taxon_id\
        and s.query_taxon_id != s.subject_taxon_id\
        and s.evalue_exp <= -5\
        and s.percent_match >= 50\
        and (s.evalue_mant < 0.01\
        or s.evalue_exp = cutoff.evalue_exp\
        and s.evalue_mant = cutoff.evalue_mant)")

    cur.execute("create unique index best_hit_ix on BestHit(query_id,subject_id)")

######################################################################

    cur.execute("create table OrthologTemp as\
        select bh1.query_id as sequence_id_a, bh1.subject_id as sequence_id_b,\
        bh1.query_taxon_id as taxon_id_a, bh1.subject_taxon_id as taxon_id_b,\
        case\
        when bh1.evalue_mant < 0.01 or bh2.evalue_mant < 0.01\
        then (bh1.evalue_exp + bh2.evalue_exp) / -2\
        else\
        (log(bh1.evalue_mant * bh2.evalue_mant)\
        + bh1.evalue_exp + bh2.evalue_exp) / -2\
        end as unnormalized_score\
        from BestHit bh1, BestHit bh2\
        where bh1.query_id < bh1.subject_id\
        and bh1.query_id = bh2.subject_id\
        and bh1.subject_id = bh2.query_id")

######################################################################

    orthologTaxonSub(cur, '')

######################################################################

    normalizeOrthologsSub(cur, '', "Ortholog")


################################################################################
############################### InParalogs ####################################
################################################################################
def inparalogs (cur):

    cur.execute("create table BestInterTaxonScore as\
        select im.query_id as query_id, low_exp.evalue_exp as evalue_exp, min(im.evalue_mant) as evalue_mant\
        from BestQueryTaxonScore im,\
        (select query_id, min(evalue_exp) as evalue_exp\
        from BestQueryTaxonScore\
        group by query_id) low_exp\
        where im.query_id = low_exp.query_id\
        and im.evalue_exp = low_exp.evalue_exp\
        group by im.query_id, low_exp.evalue_exp")

###########################################################################

    cur.execute("create unique index bis_uids_ix on BestInterTaxonScore(query_id)")

###########################################################################

    cur.execute("create table UniqSimSeqsQueryId as\
        select distinct s.query_id from SimilarSequences s")

###########################################################################

    cur.execute("create unique index ust_qids_ix on UniqSimSeqsQueryId (query_id)")

###########################################################################

    cur.execute("CREATE TABLE BetterHit as\
        select s.query_id, s.subject_id,\
        s.query_taxon_id as taxon_id,\
        s.evalue_exp, s.evalue_mant\
        from SimilarSequences s, BestInterTaxonScore bis\
        where s.query_id != s.subject_id \
        and s.query_taxon_id = s.subject_taxon_id\
        and s.query_id = bis.query_id\
        and s.evalue_exp <= -5\
        and s.percent_match >= 50\
        and (s.evalue_mant < 0.001\
        or s.evalue_exp < bis.evalue_exp\
        or (s.evalue_exp = bis.evalue_exp and s.evalue_mant <= bis.evalue_mant))\
        union\
        select s.query_id, s.subject_id, s.query_taxon_id as taxon_id, s.evalue_exp, s.evalue_mant\
        from SimilarSequences s\
        where s.query_taxon_id = s.subject_taxon_id \
        and s.evalue_exp <= -5\
        and s.percent_match >= 50\
        and s.query_id in \
        (SELECT distinct ust.query_id\
        from UniqSimSeqsQueryId ust\
        LEFT OUTER JOIN BestInterTaxonScore bis ON bis.query_id = ust.query_id\
        WHERE bis.query_id IS NULL)")

###########################################################################

    cur.execute("create index better_hit_ix on BetterHit (query_id,subject_id)")

###########################################################################

    cur.execute("create table InParalogTemp as\
        select bh1.query_id as sequence_id_a, bh1.subject_id as sequence_id_b,\
        bh1.taxon_id,\
        case\
        when bh1.evalue_mant < 0.01 or bh2.evalue_mant < 0.01\
        then (bh1.evalue_exp + bh2.evalue_exp) / -2\
        else\
        (log(bh1.evalue_mant * bh2.evalue_mant)\
        + bh1.evalue_exp + bh2.evalue_exp) / -2\
        end as unnormalized_score\
        from BetterHit bh1, BetterHit bh2\
        where bh1.query_id < bh1.subject_id\
        and bh1.query_id = bh2.subject_id\
        and bh1.subject_id = bh2.query_id")

################################################################

    cur.execute("create table InParalogTaxonAvg as\
        select avg(i.unnormalized_score) average, i.taxon_id as taxon_id\
        from InParalogTemp i\
        group by i.taxon_id")

################################################################

    cur.execute("create table OrthologUniqueId as\
        select distinct(sequence_id) from (\
        select sequence_id_a as sequence_id from Ortholog\
        union\
        select sequence_id_b as sequence_id from Ortholog) i")

################################################################

    cur.execute("create unique index ortho_uniq_id_ix on OrthologUniqueId (sequence_id)")

################################################################

    cur.execute(" create table InplgOrthTaxonAvg as\
        select avg(i.unnormalized_score) average, i.taxon_id as taxon_id\
        from InParalogTemp i\
        where i.sequence_id_a in\
        (select sequence_id from OrthologUniqueId)\
        or i.sequence_id_b in\
        (select sequence_id from OrthologUniqueId)\
        group by i.taxon_id")

################################################################

    cur.execute("create table InParalogAvgScore as\
        select case\
        when orth_i.average is NULL\
        then all_i.average\
        else orth_i.average\
        end as avg_score,\
        all_i.taxon_id\
        from InParalogTaxonAvg all_i LEFT OUTER JOIN InplgOrthTaxonAvg orth_i\
        ON all_i.taxon_id = orth_i.taxon_id")

################################################################

    cur.execute("create unique index inparalog_avg_ix on InParalogAvgScore(taxon_id,avg_score)")

################################################################

    cur.execute(" insert into InParalog(sequence_id_a, sequence_id_b, taxon_id, unnormalized_score, normalized_score)\
        select it.sequence_id_a, it.sequence_id_b, it.taxon_id, it.unnormalized_score, it.unnormalized_score/a.avg_score\
        from InParalogTemp it, InParalogAvgScore a\
        where it.taxon_id = a.taxon_id")

################################################################################
############################### CoOrthologs ###################################
################################################################################
def coorthologs (cur):

    cur.execute("create table InParalog2Way as\
        select sequence_id_a, sequence_id_b from InParalog\
        union\
        select sequence_id_b as sequence_id_a, sequence_id_a as sequence_id_b from InParalog")

######################################################################

    cur.execute("create unique index in2a_ix on InParalog2Way(sequence_id_a, sequence_id_b)")

######################################################################

    cur.execute("create unique index in2b_ix on InParalog2Way(sequence_id_b, sequence_id_a)")

######################################################################

    cur.execute("create table Ortholog2Way as\
        select sequence_id_a, sequence_id_b from Ortholog\
        union\
        select sequence_id_b as sequence_id_a, sequence_id_a as sequence_id_b from Ortholog")

######################################################################

    cur.execute("create unique index ortholog2way_ix on Ortholog2Way(sequence_id_a, sequence_id_b)")

######################################################################

    cur.execute("create table InplgOrthoInplg as\
        select ip1.sequence_id_a, ip2.sequence_id_b\
        from Ortholog2Way o, InParalog2Way ip2, InParalog2Way ip1\
        where ip1.sequence_id_b = o.sequence_id_a\
        and o.sequence_id_b = ip2.sequence_id_a")

##################################################################

    cur.execute("create table InParalogOrtholog as\
        select ip.sequence_id_a, o.sequence_id_b\
        from InParalog2Way ip, Ortholog2Way o\
        where ip.sequence_id_b = o.sequence_id_a")

##################################################################

    cur.execute("create table CoOrthologCandidate as\
        select distinct\
        min(sequence_id_a, sequence_id_b) as sequence_id_a,\
        max(sequence_id_a, sequence_id_b) as sequence_id_b\
        from (select sequence_id_a, sequence_id_b from InplgOrthoInplg\
        union\
        select sequence_id_a, sequence_id_b from InParalogOrtholog) t")
 
######################################################################

    cur.execute("create table CoOrthNotOrtholog as\
        SELECT cc.sequence_id_a, cc.sequence_id_b\
        FROM CoOrthologCandidate cc\
        LEFT OUTER JOIN Ortholog o\
        ON cc.sequence_id_a = o.sequence_id_a\
        AND cc.sequence_id_b = o.sequence_id_b\
        WHERE o.sequence_id_a IS NULL")

#####################################################################

    cur.execute("create index cno_ix on CoOrthNotOrtholog(sequence_id_a,sequence_id_b)")

######################################################################

    cur.execute("create table CoOrthologTemp as\
        select candidate.sequence_id_a, candidate.sequence_id_b,\
        ab.query_taxon_id as taxon_id_a, ab.subject_taxon_id as taxon_id_b,\
        case\
        when ab.evalue_mant < 0.00001 or ba.evalue_mant < 0.00001\
        then (ab.evalue_exp + ba.evalue_exp) / -2\
        else\
        (log(ab.evalue_mant * ba.evalue_mant)\
        + ab.evalue_exp + ba.evalue_exp) / -2\
        end as unnormalized_score\
        from SimilarSequences ab, SimilarSequences ba, CoOrthNotOrtholog candidate\
        where ab.query_id = candidate.sequence_id_a\
        and ab.subject_id = candidate.sequence_id_b\
        and ab.evalue_exp <= -5\
        and ab.percent_match >= 50\
        and ba.query_id = candidate.sequence_id_b\
        and ba.subject_id = candidate.sequence_id_a\
        and ba.evalue_exp <= -5\
        and ba.percent_match >= 50")
    
######################################################################

    orthologTaxonSub(cur, 'co')

######################################################################

    normalizeOrthologsSub(cur, "Co", "CoOrtholog")


def execute(db_dir, nm=None):
    con = lite.connect(os.path.join(db_dir, "orthoDB.db"))

    with con:

        if nm:
            if nm.stop:
                raise KillByUser("")
            nm.total = 4
            nm.counter = 0
            nm.msg = None

        con.create_function("log", 1, log)

        cur = con.cursor()

        for func in [commonTempTables, orthologs, inparalogs, coorthologs]:

            if nm:
                if nm.stop:
                    raise KillByUser("")
                nm.counter += 1

            func(cur)

if __name__ == '__main__':
    execute(".")


__author__ = "Fernando Alves and Diogo N. Silva"
