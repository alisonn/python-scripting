#!/usr/bin/env python

# Find mean coverage:
# input: genomecov file (gene - start - end - coverage) where each gene and positions are grouped together in list format, desired output name
# output: tsv file with mean coverage 
# fencepost - if end of file, take current name --> average the associated coverage

import sys

filename = sys.argv[1]
fh = open(filename)
line = fh.readline()

outname = sys.argv[2]
of = open(outname, 'w')

# housekeeping
TECov = {}
prevTE = ""
currTE = ""
count = 0

while line:
	parseLine = line.split()
	currTE = parseLine[0]
	currCov = float(parseLine[3])
	# add first TE manually to the dict
	if (bool(TECov) == False):
		print "we are adding the first gene"
		TECov[currTE] = currCov

	# option 1: TE is in the dict --> add to the value
	# option 2: TE is not in the dict --> a) mean cov for prev, b) new key pair
	if currTE in TECov: 
		TECov[currTE] = TECov[currTE]+currCov 

	# we are starting new transposon, so calc the prev TE's mean coverage, write to file, and reset counter for next TE
	else:
		TECov[prevTE] = TECov[prevTE] / count
		of.write( prevTE + '\t' + str(TECov[prevTE]) + '\n' )
		count = 0  
		TECov[currTE] = currCov

	# keep track of the TE we last saw
	prevTE = currTE	
	count += 1
	line = fh.readline()

# find mean of last guy manually and write
TECov[prevTE] = TECov[prevTE] / count
of.write( prevTE + '\t' + str(TECov[prevTE]) + '\n' )

fh.close()
of.close()

#EOF
