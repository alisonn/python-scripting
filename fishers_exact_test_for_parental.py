#!/usr/bin/env python

# Fishers Exact Test for Parental Alleles: Creates 2x2 contingency tables with parental_contribution_v2.py output and
# calculates p-value and oddsratio when data are sorted into Quartile1: Summer+Parent1, Q2: Summer+Parent2; Q3: Winter+Parent1, Q4: Winter+Parent2
# positional arguments -1 original vcf file -2 non-starved rep1 -3 non-starved rep2 -4 starved rep1 -5 starved rep2 -6 output in format "chr, pos, pvalue, oddsratio, Q1-Q4"

# TODO - figure out how to normalize readDepth across all samples to prevent overestimates/bias
# calculate readDepth of each sample?

# we need argparse for positional arguments
# need numpy and scipy for stats.fisher_exact([table])
import sys
import argparse
import numpy
from scipy import stats

# There are 2 dictionaries/mappings: (1) Treatment 1 Parental AI (2) Treatment 2 Parental AI 

# so that I don't have to rewrite parsed line shit over again
def getInfoFromLine (line):
    parsed = line.split()
    chrom = parsed[0]
    snp = parsed[1]
    altCt = parsed[4]
    refCt = parsed[5] - altCt
    altPar = parsed[6]
    numReps = [1] # default
    pos = (chrom, snp)
    posInfo = (altCt, refCt, altPar, numReps)
    result = (pos, posInfo)

    return result

# read in the parental alleles into dictionary and return dictionary to user
def loadParentCont (aiFile):
    print "Initiating a repository with parental allelic imbalance information\n"
    ai_repos = {}
    try:

        ai = open(aiFile)
        line = ai.readline()
        while line:
            if line[0] == "c":
                # tuple to tuple of tuples i.e. [chr, pos] -->  [altCount, readDepth, AltParent]
                currInfo = getInfoFromLine(line)
                ai_repos[currInfo[0]] = currInfo[1] 
            line = ai.readline()

        ai.close()

    except IOError as e:
        print "I/O error: cannot open " + aiFile + "."

    return ai_repos

# update the repository with replicate information
def updateRepos (ai_repos, repFile):
    print "Updating a repository with replicate information \n"
    try:
        rep = open(repFile)
        line = ai.readline()
        while line:

            repInfo = getInfoFromLine(line)
            currPos = repInfo[0]
            currInfo = repInfo[1]
            # add to the altCt with the replicate info
            # if the snp is identified in both reps
            if currPos in repos:
                ai_repos[currPos][0] += currInfo[0]
                ai_repos[currPos][1] += currInfo[1]
                ai_repos[currPos][2] += 1 # this snp site contains info from 2 replicates
            # otherwise tag on another key --> info
            else:
                ai_repos[currPos] = currInfo
        rep.close()
    except IOError as e:
        print "I/O error: cannot open " + repFile + "."
    return repos

# conducts fisher's exact test and returns a line to print
def getResults(currSNP, nonStarv_repos, starv_repos):
    result = ""
    # assuming everything works perfectly (i.e. snp exists in both gatk-parent files)
    if currSNP in nonStarv_repos and currSNP in starv_repos:
        fishers_row1 = [nonStarv_repos[currSNP][0], nonStarv_repos[currSNP][1]]
        fishers_row2 = [starv_repos[currSNP][0], starv_repos[currSNP][1]]
        oddsRatio, pval = stats.fishers_exact([fishers_row1, fishers_row2]) 
        result = currSNP[0] + "\t" + currSNP[1] + "\t" + oddsRatio + "\t" + pval + "\t" + fishers_row1[0] + "\t" + fishers_row1[1] + "\t" + fishers_row2[0] + "\t" + fishers_row2[1] + "\t" + nonStarv_repos[currSNP][3] + "\t" + starv_repos[currSNP][3] + "\n"
    # if the snp doesn't exist in one of the dictionaries
    elif currSNP not in nonStarv_repos:
        result = currSNP[0] + "\t" + currSNP[1] + "\tNA\tNA\t" + nonStarv_repos[currSNP][0] + "\t" + nonStarv_repos[currSNP][1] + "\tNA\tNA\t" + nonStarv_repos[currSNP][3] + "\tNA\n"
    elif currSNP not in starv_repos:
        result = currSNP[0] + "\t" + currSNP[1] + "\tNA\tNA\tNA\tNA\t" + starv_repos[currSNP][0] + "\t" + starv_repos[currSNP][1] + "\tNA\t" + starv_repos[currSNP][3] + "\n" 

    # if the snp does not exist in either dictionary
    else:
        result = currSNP[0] + "\t" + currSNP[1] + "\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\n" 

    return result

def printResults(vcf, nonStarv_repos, starv_repos, output):

    print "Computing fisher's exact test and printing results to file . . . \n"
    vcf = open(vcf_file)
    # read the vcf line by line
    out = open(output, 'w')
    out.write("# parent 1.0 = summer; parent 2.0 = winter\n")
    out.write("# contig\tpos\toddsRatio\tpval\tNS_Par1\tNS_Par2\tS_Par1\tS_Par2\tnumRepsNS\tnumRepsS\n")
    line = vcf.readline() 
    
    while line:
        vcfPosition = getInfoFromLine(line)
        currSNP = vcfPosition[0]
        resultToPrint = getResults(currSNP, nonStarv_repos, starv_repos) 
        out.write(resultToPrint)
        line = vcf.readline()

    vcf.close()
    out.close()
    #done - nothing to return

def getArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_file")
    parser.add_argument("envA_rep1")
    parser.add_argument("envA_rep2")
    parser.add_argument("envB_rep1")
    parser.add_argument("envB_rep2")
    parser.add_argument("output_file")
    args = parser.parse_args()
    return args

def printHeader():
    print "\n\n#########################################\n" 
    print "FISHER'S EXACT TEST FOR ALLELIC IMBALANCE\n"
    print "#########################################\n\n"

# main
def main():
    
    printHeader()
    args = getArgs()

    #TODO test if this works...
    starv_repos = loadParentCont(args.envA_rep1)
    starv_repos = updateRepos(starv_repos,args.envA_rep2)
    nonStarv_repos = loadParent(args.envB_rep1)
    nonStarv_repos = updateRepos(nonStarv_repos, args.envB_rep2)
    
    printResults(args.vcf_file, nonStarv_repos, starv_repos, args.output_file)

main()
# close
