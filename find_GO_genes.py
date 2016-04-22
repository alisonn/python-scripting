#!/usr/bin/env python

# Finds genes under same GO accession in a sleuth output file.
# Then creates a file to make a 2x2 contingency table based on b-value.
# Load the sleuth results into a dictionary
# (key = gene/transcript, value = pval, qval, bval)

import sys

# prompt for genes to look in dictionary
# @param string of genes separated by spaces
# @return nothing (everything is printed to screen)
def getQuery():
    query = raw_input("Input comma-separated gene queries: ")
    return query

# load the results
# heuristic: no  duplicate gene names in sleuth files
def updateRepos(fileName):
    sleuth_repos = {}
    sf = open(fileName)
    line = sf.readline()

    while line:
        # parse each line at a time
        parsed = line.split()
        name = parsed[0]
        pval = float(parsed[1])
        qval = float(parsed[2])
        sign = float(parsed[3])
        values = (pval, qval, sign) 
        # add to the dictionary and grab the next line
        sleuth_repos[name] = values
        line = sf.readline()
    print "Dictionary fully loaded with genes and associated values! :-)" 

# prompt for sleuth file to use
# @param sleuth file to open
# @return nothing, dictionary is updated with sleuth file
def getSleuth():
    sleuth_file = raw_input("Input python-compatible sleuth results file: ")
    updateRepos(sleuth_file)

def retrieveQuery(sleuth_repos, query)
    geneList =  query.split(',')
    # check for loop syntax ...
    for i in 0 to length(geneList-1):
        currVal = sleuth_repos[geneList[i]]
        if currVal[0] > 0.05:
            if currVal[2] > 0: 
                # write it is downreg in winter
            else:
                #write that it is upreg in winter
        #done
    print "Done with processing your query."
