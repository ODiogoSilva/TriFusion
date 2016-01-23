#!/usr/bin/python2

import os
import re
import sys
import pprint
from operator import itemgetter

VAR_LENGTH=0
VAR_TAXON=1

"""
Read all fasta files from a folder, placing the genes present on those fasta into a single variable
"""
#################################################################################################################
def getGenesFromFasta(fastaFilesDir) :

    genes = {}

    #for all fast files in fasta directory
    for fasta in os.listdir(fastaFilesDir):
        #skip hidden files and not fasta files
        if fasta.startswith('.') or not fasta.endswith(".fasta"):
            continue

        #get taxon from file name
        splitted = fasta.split(".")
        taxon = splitted[0]

        #open file
        fastaFile = open(os.path.join(fastaFilesDir, fasta), "r")

        gene = ''
        length = 0
        newGene = False
        for line in fastaFile:
            #clean '\n'
            line = line.strip()
            #ignore empty lines and stop codons
            if not line :# or line.endswith("*"):
                continue

            newGene = True if line.startswith(">") else False

            if newGene:
                #save previous gene info
                if gene:
                    genes[gene][VAR_LENGTH] = length

                #save new gene info
                gene = line[1:]
                genes[gene] = [None, taxon]

                #reset vars
                length = 0
            else:
                length += len(line)

        genes[gene][VAR_LENGTH] = length

        fastaFile.close()

#        pprint.pprint(genes, width=-1)
#        exit()

    return genes

#################################################################################################################
def getTaxonAndLength (subject, genes):

    subject["queryTaxon"] = genes[subject["queryId"]][VAR_TAXON];
    subject["subjectTaxon"] = genes[subject["subjectId"]][VAR_TAXON]
    subject["queryLength"] = genes[subject["queryId"]][VAR_LENGTH]
    subject["subjectLength"] = genes[subject["subjectId"]][VAR_LENGTH]

#    print(subject)

    try:
        subject["subjectTaxon"]
    except KeyError:
        print ("couldn't find taxon for gene " + subject["subjectId"])
    try:
        subject["queryTaxon"]
    except KeyError:
        print ("couldn't find taxon for gene " + subject["queryId"])

    return int(subject["queryLength"]) < int(subject["subjectLength"])

#################################################################################################################
def printPreviousSubject(subject):
    nonOverlapMatchLen = computeNonOverlappingMatchLength(subject)

    percentIdent = int(subject["totalIdentities"] / subject["totalLength"] * 10 + .5)/10;
    shorterLength = subject["queryLength"] if subject["queryShorter"] else subject["subjectLength"]

    print(nonOverlapMatchLen)
    print(shorterLength)

    percentMatch = int(nonOverlapMatchLen / shorterLength * 1000 + .5) / 10;
    print (subject["queryId"] + "\t" + subject["subjectId"] + "\t" + subject["queryTaxon"] + "\t" + subject["subjectTaxon"] + "\t" + str(subject["evalueMant"]) + "\t" + str(subject["evalueExp"]) + "\t" + str(percentIdent) + "\t" + str(percentMatch))

#################################################################################################################
# this (corrected) version of formatEvalue provided by Robson de Souza
def formatEvalue (evalue):
    if evalue == '0':
        return (0,0)
    if evalue.startswith('e'):
        evalue = '1'.evalue
    return [round(float(x), 2) for x in evalue.split("e")]

#################################################################################################################
def computeNonOverlappingMatchLength (subject):
    #flatten lists
    hsps = subject["hspspans"]
    hsps.sort(key=lambda x: x[0])
    original = hsps.pop(0)

    original = (int(original[0]), int(original[1]))

    start = 0
    end = 1
    
    if not original:
        return 0

    original = getStartEnd(original)
    length = 0

    for h in hsps:
        h = getStartEnd(h)
        if h[end] < original[end]: #does not extend
            continue 
        if h[start] <= original[end]: #overlaps
            original = (original[start], h[end])  #extend end ... already dealt with if new end is less
        else:  #there is a gap in between
            length = original[end] - original[start] + 1
            original = (h[start], h[end])

    length += original[end] - original[start] + 1 # deal with the last one 
    return length

#################################################################################################################
# flip orientation if nec.
def getStartEnd (h):
    start = h[0]
    end = h[1]
    if start > end:
        end = h[0]
        start = h[1]
 
    return(int(start),int(end))

#################################################################################################################

def orthomclBlastParser(blastFileName, fastaFilesDir):

    genes = getGenesFromFasta(fastaFilesDir);
    blastFile = open(blastFileName, "r")

    prevSubjectId = ''
    prevQueryId = ''
    subject = {} # hash to hold subject info

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
            subject["totalIdentities"] = 0
            subject["totalLength"] = 0
            subject["hspspans"] = []

        # get additional info from subsequent HSPs
        hspspan = (subjectStart, subjectEnd)
        if subject and subject["queryShorter"] :
            hspspan = (queryStart, queryEnd)
        subject["hspspans"].append(hspspan)
        subject["totalIdentities"] += float(percentIdentity) * length
        subject["totalLength"] += length

    printPreviousSubject(subject)

orthomclBlastParser("../test/AllVsAll.out", "../test/compliantFasta")
