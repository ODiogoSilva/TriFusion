#!/usr/bin/python

def orthomclMCLToGroups(prefix, startId, inFileName, outFileName):

	if not isinstance(startId, int)
		raise TypeError("StartId is not a number")

	inputFile = open(inFileName, "r")
	outFile = open(outFileName, "w")

	for line in inputFile:
		out.write(prefix+startId+": "+line)
		startId++

