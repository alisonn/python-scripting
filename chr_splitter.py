#!/usr/bin/env python

# Chromosome Splitter
# alisonn - adapted from ebrown
# Input: a bed file (resulting from coverageBed)  
# Output: temp coverage files for each chromosome 
# Method: Create a dictionary [chr --> chr_file.txt] and write to each file given [chr]

import sys
import os

infh = sys.argv[1]
infile = open(infh, 'r')

# list of chromosomes
chromo = ["2L", "2R", "3L", "3R", "X", "Y", "4"]
# dict [ chr --> chr.tmp file for writing] 
chrFileDict = {
    c: open(c + '.tmp', 'w') for c in chromo 
}

for line in infile:
    if line == '': break
    line = line.strip()
    split = line.split('\t')
    chrom = split[0]
    if chrom not in chromo:
        continue
    if chrom in chromo:
        chrFileDict[chrom].write(str(line) + '\n')

infile.close()
# close all temp files
for f in chrFileDict.values():
    f.close()

# Done :-)
