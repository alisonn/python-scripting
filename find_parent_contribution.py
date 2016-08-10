#!/usr/bin/env python

# Find Parent Contribution: This script figures out how much each parent contributes to allelic imbalance in NGS data
# Input: VCF File with Parents, GATK-output File with Ref and Alt allele columns
# User Input: Regions of interest to asses parental contribution
# Output: If you prompt - a .csv file with regions of interest and the main contributing parent per region/SNP

import sys

# load an initial VCF file
def loadVCF():
    vcf_file = raw_input("Input the path to your VCF File: ")
    return updateVCF(vcf_file)

# update VCF dictionary: yup, this works! tuple --> tuple
# @param: path to vcf file
# @return: a dictionary containing chr, pos, and par genotypes
def updateVCF(fileName):

    #initiate a dictionary of dictionaries (tuple --> tuple)
    vcf_repos = {} 
    vcf = open(fileName)
    line = vcf.readline()

    # populate the dictionary with parent-SNP information
    while line:

        if line[0] == "c":
            # the line contains information about SNPs at particular sites - so process
            # I just want to know if parent is ref or alt easily
            parsed = line.split()
            currChr = parsed[0]
            currSite = parsed[1]
            currPos = (currChr, currSite)

            #start processing parental genotypes here
            parent1 = ""
            parent2 = ""
            if parsed[9][0] == "1":
                parent1 = "alt"
                parent2 = "ref"
            else:
                parent1 = "ref"
                parent2 = "alt"

            vcf_repos[currPos] = (parent1, parent2) 

        line = vcf.readline() 

    vcf.close() 
    return vcf_repos

# load a GATK output file
def loadGATK():
    GATK_file_path = raw_input("Input the path to your GATK output file (.csv): ")
    return updateGATK(GATK_file_path)

## update GATK dictionary - yup, this works fine as well
# GATK_repos is in format of (chr, snp) --> (altFreq, rawDepth)
def updateGATK(fileName):
    gatk = open(fileName)
    GATK_repos = {}
    line = gatk.readline()

    while line:
        if line[1] == "h":
            parsed = line.split()
            currChr = parsed[0]
            currSNP = parsed[1]
            currPos = (currChr, currSNP)
            currAltCt = float(parsed[6])
            currTotalCt = float(parsed[7])
            currAltFreq = float(0)
            if currTotalCt != float(0):
                currAltFreq = currAltCt/currTotalCt
            currRawDepth = parsed[10]
            GATK_repos[currPos] = (currAltFreq, currRawDepth)
        line = gatk.readline()

    gatk.close()
    return GATK_repos

## get the user's query
# this should be in format of chr2L:11111111 (split on commas)
# @return a string of queries separated by commas
def getQuery():
    query = raw_input("Input your regions of interest (enter ALL if you want to find ALL possible sites): ")
    return query

# function that takes in a string of queries in a particular format
# @param a string of queries in format of chr:position,chr:position,...,chr:position
# @return a set of tuples to be used as keys for both vcf and gatk dictionaries
def queryToTuple(query):
    parsedQuery = query.split(",")
    # find out how many times to loop and create enough tuples as keys 
    tupleSet = set()
    # we are indexed at 0, so edge case
    for i in range (0, len(parsedQuery)-1):
        preTuple = parsedQuery[i].split(":") 
        currChr = preTuple[0]
        currSNP = preTuple[1]
        currPos = (currChr, currSNP) 
        tupleSet.add(currPos)
    return tupleSet

# special case of retrieving queries when you want all the sites
# @param GATK repository that has been filled already
# @return a list of keys as tuples from the GATK dictionary
def getAllSites(gatk_repos):
    # the max valid sites you can use are those from the GATK-out file    
    tupleSet = set() 
    for key in gatk_repos:
        tupleSet.add(key)
    return tupleSet

# Search for the query and save to a dictionary?
# @param set of tuples, vcf dictionary, and gatk dictionary (I think these are passed by reference)
# @return a results dictionary of tuples (chr, snp) --> tuple (altFreq, parent, rawDepth)
def getQueryResults(queryTupleSet, vcf_repos, GATK_repos):

    results_repos = {}

    #queryTupleSet contains tuples that are keys compatible with both vcf and GATK dictionaries
    # iterate over query Tuple Set - for each element do the following
    for posit in queryTupleSet:

        if posit in vcf_repos.keys() and posit in gatk_repos.keys():
            
            currParAlleles = vcf_repos[posit]
            currAltFreq = gatk_repos[posit][0]
            currRawDepth = gatk_repos[posit][1] 

            # 1 = parent 1 , 2 = parent 2, 1.5 = no imbalance, 0.0 = no information
            if currAltFreq == float(0.5):
                results_repos[posit] = (float(1.5), float(currAltFreq), float(currRawDepth))

            elif (currAltFreq > float(0.5) and gatk_repos[0] == "alt") or (currAltFreq < 0.5 and gatk_repos[0] == "ref"):
                results_repos[posit] = (float(1.0), float(currAltFreq), float(currRawDepth))

            elif (currAltFreq > float(0.5) and gatk_repos[0] == "ref") or (currAltFreq < 0.5 and gatk_repos[0] == "alt"):
                results_repos[posit] = (float(2.0), float(currAltFreq), float(currRawDepth))

            else:
                results_repos[posit] = (float(0), float(currAltFreq), float(currRawDepth))

        else:
            print "This SNP does not exist in either vcf or gatk files:  " + posit[0] + ": " + posit[1]

    return results_repos

# okay go test this alison wtf
# write the query results to a file
def queryToFile(results_repos):
    file_name = raw_input("Input the name of the output file to write: ")

    print contig + "\t" + pos + "\t" + parent + "\t" + altFreq + "\t" + rawDepth + "\n"

    for key, value in results_repos:
        print key[0] + "\t" + key[1] + "\t" + results_repos[key][0] + "\t" + results_repos[key][1] + "\t" + results_repos[key][2] + "\n"
    # okay find the function for writing to file.....

# main - yup this is it!!!
def main():
    # initialize the variables
    vcf_repos = {}
    gatk_repos = {}
    results_repos = {}
    query = ""
    querySet = set()
    
    print "\n"
    print "#####################################"
    print "WELCOME TO FIND PARENTAL CONTRIBUTION"
    print "#####################################" + "\n"
    print "This program calculates parental contributions to allelic imbalance in NGS data."
    print "You will need VCF and GATK files for each sample." + "\n"

    while True:

        print "What do you want to do? Enter one of the options"
        print "Upload VCF: V ...  Upload GATK: G ... Calculate Parent Contribution: P ... Write Results: W ... Exit: X ..."
        choice = raw_input("Your choice: ")
        if choice == "V":
            vcf_repos = loadVcf()
        elif choice == "G":
            gatk_repos = loadGATK() 
        elif choice == "P":
            query = getQuery()
            if query == "ALL":
                querySet = getAllSites()
            else:
                querySet == getQueryTuples
        elif choice == "X":
            print "Okay, closing FIND PARENTAL CONTRIBUTION"
            break
        else:
            print "Sorry, not an option. Try again."
# close

main()
