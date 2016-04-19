#!/usr/bin/env python

# Finds genes under same GO accession in a sleuth output file.
# Then creates a file to make a 2x2 contingency table based on b-value.
# Load the sleuth results into a dictionary
# (key = gene/transcript, value = b-value)

import sys

# prompt for genes to look in dictionary
# @param string of genes separated by spaces
# @return nothing (everything is printed to screen)
def getQuery(argv[2]):
    query = raw_input("Input comma-separated gene queries: ")
    return query

# load the results
# heuristic: I know I don't have duplicate gene names in sleuth files
def fillRepos (fileName)
    sleuth_repos = {}
    sf = open(fileName)
    line = sf.readline()

    while line:
        # have to parse each line at a time
        parsed = line.split()
        name = parsed[0]
        pval = float(parsed[1])
        qval = float(parsed[2])
        sign = float(parsed[3])
        values = (pval, qval, sign) 

        sleuth_repos[name] = values

        line = sf.readline()
    done

# prompt for sleuth file to use
# @param sleuth file to open
# @return nothing, dictionary is updated with sleuth file

def getSleuth():
    sleuth_file = raw_input("python-compatible sleuth results file: ")
    fillRepos (sleuth_file)
