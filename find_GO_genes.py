#!/usr/bin/env python

# Finds genes under same GO accession in a sleuth output file.
# Then creates a file to make a 2x2 contingency table based on b-value.
# Load the sleuth results into a dictionary
# (key = gene/transcript, value = pval, qval, bval)

import sys
import re

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
        pval = parsed[1]
        if pval != "NA":
            pval = float(pval)
            qval = float(parsed[2])
            sign = float(parsed[3])
            values = (pval, qval, sign) 
            # add to the dictionary and grab the next line
            sleuth_repos[name] = values
            line = sf.readline()
        # If we hit NA values, then we are at end of significant transcripts
        else: 
            print "Dictionary fully loaded with genes and associated significant values! :-)" 
            break 
    return sleuth_repos

# prompt for sleuth file to use
# @param sleuth file to open
# @return nothing, dictionary is updated with sleuth file
def getSleuth():
    sleuth_file = raw_input("Input python-compatible sleuth results file: ")
    return updateRepos(sleuth_file)

def getIsoforms(sleuth_repos, query):
    # search in the dictionary keys all possible transcripts 
    allGenes = sleuth_repos.keys()
    query = "^" + query
    for i in xrange(len(allGenes)):
        re.search(query, allGenes[i])

def retrieveQuery(sleuth_repos, query):
    geneList =  query.split(',')
    p_threshold = 0.05

    # check that the gene/transcript exists
    # add support for regular expression
    for i in xrange(len(geneList)):
        currGene = geneList[i]
        if currGene in sleuth_repos:
            currVal = sleuth_repos[geneList[i]]
            if currVal[0] <= p_threshold:
                if currVal[2] < 0: 
                    print currGene + " is down-regulated in starvation condition."
                else:
                    print currGene + " is up-regulated in starvation condition."
        else: 
            print currGene + " was not found or was not significant."

    print "Done with processing your query!"

# interactive prompting for user input
def main():
    repository = {}
    while True:
        # What do you want to do?
        print "What do you want to do?"
        choice = raw_input("S: load dictionary Q: query search or X: exit... ")
        if choice == "S":
            repository = getSleuth()
        elif choice == "Q" and bool(repository):
            query = str(getQuery())
            retrieveQuery(repository, query)
        elif choice == "Q" and not bool(repository):
            print "Sorry, no dictionary has been loaded."
        elif choice == "X":     
            # deal with this later
            print "Okay, closing gene finder."
            break
        else:
            print "Sorry, not an option."
# close
main()
