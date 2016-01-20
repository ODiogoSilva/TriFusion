#!/usr/bin/python2

import os
import re


"""
Read all fasta files from a folder, placing the genes present on those fasta into a single variable
"""
#################################################################################################################
def getGenesFromFasta(fastaFilesDir) :

    genes = {}

    for fasta in os.listdir(fastaFilesDir):
        #skip hidden files and not fasta files
        if fasta.startswith('.') or not fasta.endswith(".fasta"):
            continue

        splitted = fasta.split(".")
        taxon = splitted[len(splitted) - 1]
        fastaFile = open(os.path.join(fastaFilesDir, fasta), "r")
        gene = ''
        length = 0
        for line in fastaFile:
            if not line:
                continue
            line = re.sub('\>(\S+)', '', line)
            if not line:
                if gene:
                    genes["gene"]["length"] = length
                    genes["gene"]["taxon"] = taxon
                gene = line
                length = 0
            else:#TODO not sure if this is correct
                length += len(line)

        if gene:
            genes["gene"]["length"] = length
            genes["gene"]["taxon"] = taxon
        fastaFile.close()

    return genes

#################################################################################################################
def getTaxonAndLength (subject, genes):

    subject["queryTaxon"] = genes[subject["queryId"]]["taxon"];
    subject["subjectTaxon"] = genes[subject["subjectId"]]["taxon"]
    subject["queryLength"] = genes[subject["queryId"]]["length"]
    subject["subjectLength"] = genes[subject["subjectId"]]["length"]
    try:
        subject["subjectTaxon"]
    except KeyError:
        print "couldn't find taxon for gene " + subject["subjectId"]
    try:
        subject["queryTaxon"]
    except KeyError:
        print "couldn't find taxon for gene " + subject["queryId"]

    return subject["queryLength"] < subject["subjectLength"]

#################################################################################################################
def printPreviousSubject(subject):
    nonOverlapMatchLen = computeNonOverlappingMatchLength(subject)

    percentIdent = int(subject["totalIdentities"] / subject["totalLength"] * 10 + .5)/10;
    shorterLength = subject["queryLength"] if subject["queryShorter"] else subject["subjectLength"]
    percentMatch = int(nonOverlapMatchLen / shorterLength * 1000 + .5) / 10;
    print subject["queryId"] + "\t" + subject["subjectId"] + "\t" + subject["queryTaxon"] + "\t" + subject["subjectTaxon"] + "\t" + subject["evalueMant"] + "\t" + subject["evalueExp"] + "\t" + percentIdent + "\t" + percentMatch + "\n"

#################################################################################################################
# this (corrected) version of formatEvalue provided by Robson de Souza
def formatEvalue (evalue):
    if evalue.startswith('e'):
        evalue = '1'.evalue
    return [round(float(x), 2) for x in evalue.split("e")]

#################################################################################################################
def computeNonOverlappingMatchLength (subject):

    hsps = [str(x) for x in sorted([int(x) for x in subject["hspspans"]])]
    first = hsps.pop(0)

    if not first:
        return 0

    startEnd = getStartEnd(first)
    start = startEnd[0]
    end = startEnd[1]
    len = 0
    for h in hsps:
        hspStartEnd = getStartEnd(h)
        hspStart = hspStartEnd[0]
        hspEnd = hspStartEnd[0]

        if hspEnd <= end:
            continue
        if hspStart <= end:
            end = hspEnd
        else:
            len += end - start + 1
            start = hspStart
            end = hspEnd

    len += end - start + 1

    return len

#################################################################################################################
# flip orientation if nec.
def getStartEnd (h):
    start = h[0]
    end = h[1]
    if start > end:
        end = h[0]
        start = h[1]
 
    return(start,end)

#################################################################################################################

def orthomclBlastParser(blastFileName, fastaFilesDir):

    genes = getGenesFromFasta(fastaFilesDir);
    blastFile = open(blastFileName, "r")

    prevSubjectId = ''
    prevQueryId = ''
    subject = {} # hash to hold subject info
    # queryShorter not used
    queryShorter = ''

    for line in blastFile:
        splitted = line.split()

        queryId = splitted[0]
        subjectId = splitted[1]
        percentIdentity = splitted[2]
        length = int(splitted[3])
        # mismatches not used
        mismatches = splitted[4]
        # ngaps not used
        ngaps = splitted[5]
        queryStart = splitted[6]
        queryEnd = splitted[7]
        subjectStart = splitted[8]
        subjectEnd = splitted[9]
        evalue = splitted[10]
        # bits not used
        bits = splitted[11]

        if queryId != prevQueryId or subjectId != prevSubjectId:

            # print previous subject
            if subject:
                printPreviousSubject(subject)

            # initialize new one from first HSP
            prevSubjectId = subjectId
            prevQueryId = queryId

            subject = {}
            subject["queryId"] = queryId
            subject["subjectId"] = subjectId
            subject["queryShorter"] = getTaxonAndLength(subject, genes)

            tup = formatEvalue(evalue) # from first hsp
            subject["evalueMant"] = tup[0]
            subject["evalueExp"] = tup[1]

        # get additional info from subsequent HSPs
        hspspan = {subjectStart, subjectEnd}
        if subject and subject["queryShorter"] :
            hspspan = {queryStart, queryEnd}
        subject["hspspans"] = hspspan
        subject["totalIdentities"] += percentIdentity * length
        subject["totalLength"] += length

    printPreviousSubject(subject)

orthomclBlastParser("../test/AllVsAll.out", "../test/compliantFasta")