#!/usr/bin/env python

# Find Parent Contribution: This script figures out how much each parent contributes to allelic imbalance in NGS data
# Input: VCF File with Parents, GATK-output File with Ref and Alt allele columns
# User Input: Regions of interest to asses parental contribution
# Output: If you prompt - a .csv file with regions of interest and the main contributing parent per region/SNP


# Need binomial from scipy and argparse for parameter import
import sys
import argparse
import numpy as np
from scipy import stats

# load an initial VCF file
def loadVCF(VCF_file_path):
    return updateVCF(VCF_file_path)

# update VCF dictionary: yup, this works! tuple --> tuple
# @param: path to vcf file
# @return: a dictionary containing chr, pos, and par genotypes
def updateVCF(fileName):
    vcf_repos = {}

    try:
        #initiate a dictionary of dictionaries (tuple --> tuple)
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

    except IOError:
        print "IOError: cannot open " + fileName + "."

    return vcf_repos

# load a GATK output file
def loadGATK(GATK_file_path):
    return updateGATK(GATK_file_path)

## update GATK dictionary - yup, this works fine as well
# GATK_repos is in format of (chr, snp) --> (altFreq, rawDepth)
def updateGATK(fileName):
    GATK_repos = {}

    try:
        gatk = open(fileName)
        line = gatk.readline()


# what you want is chr, pos --> currAltFreq (in case), altAlleleCt, readDepth
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
                GATK_repos[currPos] = (currAltFreq, currAltCt, currRawDepth)
            line = gatk.readline()

        gatk.close()

    except IOError as e:
        print "I/O error: cannot open " + fileName + "."

    return GATK_repos 

## get the user's query
# this should be in format of chr2L:11111111 (split on commas)
# @return a string of queries separated by commas
def getQuery():
    query = raw_input("Input your regions of interest (enter ALL if you want to find ALL possible sites): ")
    return query

# function that takes in a string of queries in a particular format - yup, this works! :-)
# @param a string of queries in format of chr:position,chr:position,...,chr:position
# @return a set of tuples to be used as keys for both vcf and gatk dictionaries
def queryToTuple(query):
    parsedQuery = query.split(",")
    # find out how many times to loop and create enough tuples as keys 
    tupleSet = set()
    # we are indexed at 0, so edge case
    for i in range (0, len(parsedQuery)):
        preTuple = parsedQuery[i].split(":") 
        currChr = preTuple[0]
        currSNP = preTuple[1]
        currPos = (currChr, currSNP) 
        tupleSet.add(currPos)
    return tupleSet

# special case of retrieving queries when you want all the sites - yup, this works! :-)
# @param GATK repository that has been filled already
# @return a list of keys as tuples from the GATK dictionary
def getAllSites(gatk_repos):
    # the max valid sites you can use are those from the GATK-out file    
    tupleSet = set() 
    for key in gatk_repos:
        tupleSet.add(key)
    print "There are " + str(len(tupleSet)) + " elements in your query search!"
    return tupleSet

# @param set of tuples, vcf dictionary, and gatk dictionary (I think these are passed by reference)
# @return a results dictionary of tuples (chr, snp) --> tuple (altFreq, parent, rawDepth)
def getQueryResults(queryTupleSet, vcf_repos, gatk_repos):
    results_repos = {}
    #queryTupleSet contains tuples that are keys compatible with both vcf and GATK dictionaries
    # iterate over query Tuple Set - for each element do the following
    for posit in queryTupleSet:
        if (posit in vcf_repos.keys()) and (posit in gatk_repos.keys()):
            
            currParAlleles = vcf_repos[posit]
            currAltFreq = gatk_repos[posit][0]
            currAltCt = gatk_repos[posit][1]
            currRawDepth = gatk_repos[posit][2] 

            #TODO: Test the scipy package with binom_test()
            # do binomial test per snp site
            if currRawDepth == float(0):
                currPval = float(1)
            else:
                currPval =  stats.binom_test(currAltCt, currRawDepth, p = 0.5)

            # only report who is the parent with alternate allele
            if vcf_repos[posit][0] == "alt":
                currAltParent = "1"
            else:
                currAltParent = "2"
            results_repos[posit] = (currAltFreq, currPval, currAltCt, currRawDepth, currAltParent)

        else:
            print "This SNP does not exist in either vcf or gatk files:  " + posit[0] + ": " + posit[1]

    print "DONE!"
    return results_repos

# write the query results to a file
def queryToFile(results_repos, file_name):
#    file_name = raw_input("Input the name of the output file to write: ")
    f = open(file_name, 'w')
    f.write("# parent: (1.0=leftmost parent from VCF, 2.0=rightmost parent from VCF, 1.5=both parents, 0.0=no information)\n")
    f.write("contig\tpos\taltFreq\tp-value\taltCount\treadDepth\taltParent\n")
    for key in results_repos:
        line = str(key[0]) + "\t" + str(key[1]) + "\t" + str(results_repos[key][0]) + "\t" + str(results_repos[key][1]) + "\t" + str(results_repos[key][2]) + "\t" + str(results_repos[key][3]) + "\t" + str(results_repos[key][4]) +  "\n"
        f.write(line)
    f.close()

# prints out the header
def printHeader():
    #header
    print "\n\n"
    print "#####################################"
    print "     FIND PARENTAL CONTRIBUTION      "
    print "#####################################\n\n"
    print "This program calculates parental contributions to allelic imbalance in NGS data."
    print "You will need VCF and GATK files for each sample\n"

# specifies required inputs for processing
# @return: object whose attributes are the different inputs for main 
def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_file")
    parser.add_argument("gatk_file")
    parser.add_argument("--query", "-q", type=str, default="ALL", help="Default is all sites from GATK")
    parser.add_argument("out_file")
    args = parser.parse_args()
    return args 

# need 4 arguments: VCF_path, GATK_path, query, file_out
# easier for command line handling
def mainCmdLine():
    #initialize the variables again...
    vcf_repos = {}
    gatk_repos = {}
    results_repos = {}
    query = ""
    querySet = set()

    printHeader()
    args = getArgs()

    vcf_repos = loadVCF(args.vcf_file)
    gatk_repos = loadGATK(args.gatk_file)
    query = args.query
    if query == "ALL":
        querySet = getAllSites(gatk_repos)
    else:
        querySet = queryToTuple(query)
    results_repos = getQueryResults(querySet, vcf_repos, gatk_repos)
    queryToFile(results_repos, args.out_file)

# close
mainCmdLine()
