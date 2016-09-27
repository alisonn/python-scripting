#!/usr/bin/env python

# Fishers Exact Test for Parental Alleles: Creates 2x2 contingency tables with parental_contribution_v2.py output and
# calculates p-value and oddsratio when data are sorted into Quartile1: Summer+Parent1, Q2: Summer+Parent2; Q3: Winter+Parent1, Q4: Winter+Parent2
# positional arguments -1 original vcf file -2 starved rep1 -3 starved rep2 -4 non-starved rep1 -5 non-starved rep2 -6 output in format "chr, pos, pvalue, oddsratio, Q1-Q4"

# we need argparse for positional arguments
# need numpy and scipy for stats.fisher_exact([table])
import sys
import argparse
import numpy
from scipy import stats

# There are 2 dictionaries/mappings: (1) Treatment 1 Parental AI (2) Treatment 2 Parental AI 
# readDepth in result output is NOT SCALED... ergo, you can figure out the scaling factor and how much information is lost when normalizing across replicates
# TODO: create a set of objects instead of dictionary of tuples --> tuples
def getInfoFromLine (line):
    parsed = line.split()
    chrom = parsed[0]
    snp = parsed[1]
    altCt = float(parsed[4])
    readDepth = float(parsed[5])
    refCt = readDepth - altCt
    altPar = float(parsed[6])
    numReps = float(1) # default
    pos = (chrom, snp)

    posInfo = (altCt, refCt, altPar, numReps, readDepth)
#    posInfo = (altCt, refCt, altPar, numReps)
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
            if line[1] == "h":
                # tuple to tuple of tuples i.e. [chr, pos] -->  [altCount, readDepth, AltParent, numReps, readDepth(unscaled)]
                currInfo = getInfoFromLine(line)
                ai_repos[currInfo[0]] = currInfo[1] 
            line = ai.readline()

        ai.close()

    except IOError as e:
        print "I/O error: cannot open " + aiFile + "."

    return ai_repos

def scaleReplicates (factor, repInfo):
    # scaling factor = we gotta scale each of the counts by the scaling factor and remake a new tuple
    scaleAltCt = int(factor*repInfo[0])
    scaleRefCt = int(factor*repInfo[1])
    scaleRepInfo = (scaleAltCt, scaleRefCt, repInfo[2], repInfo[3], repInfo[5]) 

    return scaleRepInfo 

# update the repository with replicate information
def updateRepos (ai_repos, repFile):
    print "Updating a repository with replicate information \n"
    try:
        rep = open(repFile)
        line = rep.readline()
        while line:
            if line[1] == "h":
                repInfo = getInfoFromLine(line)
                currPos = repInfo[0]
                currInfo = repInfo[1]
                # add to the altCt with the replicate info
                # if the snp is identified in both reps
                # tuples are immutable so create a new tuple with entries as sum of prev and new info
                if currPos in ai_repos:
                    # TODO: normalize replicates to rep with lowest coverage 
                    if ai_repos[currPos][4] < currInfo[4]:
                        scaleFactor = float(ai_repos[currPos][4]/currInfo[4])
                        scaledInfo = scaleReplicates(scaleFactor, currInfo)
                        updatedTuple = (ai_repos[currPos][0]+currInfo[0], ai_repos[currPos][1]+currInfo[1], ai_repos[currPos][2], ai_repos[currPos][3]+float(1))
                        ai_repos[currPos] = updatedTuple

                    elif ai_repos[currPos][4] > currInfo[4]:
                        scaleFactor = float(currInfo[4]/ai_repos[currPos][4])
                        scaledInfo = scaleReplicates(scaleFactor, ai_repos[currPos][4])
                        updatedTuple = (scaledInfo[0]+currInfo[0], scaledInfo[1]+currInfo[1], currInfo[2], scaledInfo[3]+float(1), scaledInfo[4]+currInfo[4])
                        ai_repos[currPos] = updatedTuple

                    else: 
                        updatedTuple = (ai_repos[currPos][0]+currInfo[0], ai_repos[currPos][1]+currInfo[1], ai_repos[currPos][2]+currInfo[2], ai_repos[currPos][3]+float(1), ai_repos[currPos][4]+currInfo[4])
                        ai_repos[currPos = updatedTuple

                # otherwise tag on another key --> info
                else:
                    ai_repos[currPos] = currInfo
            line = rep.readline()
        rep.close()
    except IOError as e:
        print "I/O error: cannot open " + repFile + "."
    return ai_repos

# conducts fisher's exact test and returns a line to print
def getResults(currSNP, nonStarv_repos, starv_repos):
    result = ""

    # TODO: normalize across the rows so that there is equi amts from both starv and non-starv environments - gotta normalize by the smallest readDepth/Environment
    # TODO: may be helpful to normalize the readDepth from "updatewRepInfo"
    # assuming everything works perfectly (i.e. snp exists in both gatk-parent files)
    if currSNP in nonStarv_repos and currSNP in starv_repos:
        fishers_row1 = [nonStarv_repos[currSNP][0], nonStarv_repos[currSNP][1]]
        fishers_row2 = [starv_repos[currSNP][0], starv_repos[currSNP][1]]
        oddsRatio, pval = stats.fisher_exact([fishers_row1, fishers_row2]) 
        result = currSNP[0] + "\t" + currSNP[1] + "\t" + str(oddsRatio) + "\t" + str(pval) + "\t" + str(fishers_row1[0]) + "\t" + str(fishers_row1[1]) + "\t" + str(fishers_row2[0]) + "\t" + str(fishers_row2[1]) + "\t" + str(nonStarv_repos[currSNP][3]) + "\t" + str(starv_repos[currSNP][3]) + "\n"
    # if the snp doesn't exist in one of the dictionaries

    elif currSNP not in starv_repos and currSNP not in nonStarv_repos:
        result = currSNP[0] + "\t" + currSNP[1] + "\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\t0\n"

    elif currSNP not in starv_repos:
        result = str(currSNP[0]) + "\t" + str(currSNP[1]) + "\tNA\tNA\t" + str(nonStarv_repos[currSNP][0]) + "\t" + str(nonStarv_repos[currSNP][1]) + "\tNA\tNA\t" + str(nonStarv_repos[currSNP][3]) + "\tNA\n"

    elif currSNP not in nonStarv_repos:
        result = str(currSNP[0]) + "\t" + str(currSNP[1]) + "\tNA\tNA\tNA\tNA\t" + str(starv_repos[currSNP][0]) + "\t" + str(starv_repos[currSNP][1]) + "\tNA\t" + str(starv_repos[currSNP][3]) + "\n" 

    return result

def printResults(vcf_file, nonStarv_repos, starv_repos, output):

    print "Computing fisher's exact test and printing results to file . . . \n"
    vcf = open(vcf_file)
    # read the vcf line by line
    out = open(output, 'w')
    out.write("# parent 1.0 = summer; parent 2.0 = winter\n")
    out.write("# readDepth is not scaled so you can estimate the amount of data lost due to normalization across replicates") 
    out.write("# contig\tpos\toddsRatio\tpval\tNS_Par1\tNS_Par2\tS_Par1\tS_Par2\trDepthNS\trDepthS\tnumRepsNS\tnumRepsS\n")
    print "Output header written\n"
    line = vcf.readline() 
    
    while line:
        if line[0] != "#" or line[1] == "h":
            currSNP = (line.split()[0], line.split()[1])
            resultToPrint = getResults(currSNP, nonStarv_repos, starv_repos) 
            out.write(resultToPrint)
        line = vcf.readline()

    print "Finished obtaining and writing the results!\n"
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
    print "\n\n\n\n\n#########################################\n" 
    print "FISHER'S EXACT TEST FOR ALLELIC IMBALANCE\n"
    print "#########################################\n\n"

# main
def main():
    
    printHeader()
    args = getArgs()

    starv_repos = loadParentCont(args.envA_rep1)
    starv_repos = updateRepos(starv_repos,args.envA_rep2)
    nonStarv_repos = loadParentCont(args.envB_rep1)
    nonStarv_repos = updateRepos(nonStarv_repos, args.envB_rep2)
    # write the run log
#    output = args.output_file(open, "w")
#    output.write("# run arguments: " + args.vcf_file + " " + args.envA_rep1 + " " + args.envA_rep2 + " " + args.envB_rep1 + " " + args.envB_rep2 + " " + args.output_file + ". \n")
#    output.close()
    printResults(args.vcf_file, nonStarv_repos, starv_repos, args.output_file)

main()
# close
