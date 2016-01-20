#!/usr/bin/python

import os
import re
from EmptyFolder import *

minLength = 10
maxStopPercent = 20
good = open("goodProteins.fasta", "w")
bad = open("poorProteins.fasta", "w")

def handleSeq(seq, length, stopCnt):
	isBad = 0;
	stopPercent = ((length - stopCnt)/length)* 100;

	if length < minLength or stopPercent > maxStopPercent:
		bad.write(seq + "\n")
		isBad = 1
	else:
		good.write(seq + "\n")
	
	return isBad

def orthomclFilterFasta(inputDir, minLength, maxStopPercent):

	minLength = minLength
	maxStopPercent = maxStopPercent

	filenames = [os.path.join(inputDir, x) for x in os.listdir(inputDir)]
	#TODO
	if not filenames:
		raise EmptyFolder("The provided folder is empty.")

	rejectRates = []


	for fileName in filenames:
		if fileName.startswith('.'):
			continue

		inputFile = open(fileName, 'r')
		seqCount = 0
		rejectSeqCount = 0
		currentSeq = ""
		currentLen = 0
		currentStopCnt = 0

		# process lines of one file
		for line in inputFile:
			if line.startswith('>'):
				if (currentSeq):
					seqCount += 1
					rejectSeqCount += handleSeq(currentSeq, currentLen, currentStopCnt)
					currentSeq = "";
					currentLen = 0;
					currentStopCnt = 0;
			else:
				lineLen = len(line)
				currentLen += lineLen
				line = re.sub('[^A-Za-z]', '', line)
				currentStopCnt += lineLen - len(line) # this removes the stop codon from line
			
			currentSeq += line
	
		rejectSeqCount += handleSeq(currentSeq, currentLen, currentStopCnt)
		seqCount += 1;

		# add file stats to reject count if it qualifies
		if (rejectSeqCount):
			pct = rejectSeqCount/seqCount * 100;
			if (pct > 10):
				rejectRates.append([inputFile, pct])

		inputFile.close()


orthomclFilterFasta("/home/fernando/Dropbox/Diogo_Fernando/perl_scripts/cena" ,minLength, maxStopPercent)
