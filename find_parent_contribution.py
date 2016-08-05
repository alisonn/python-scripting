#!/usr/bin/env python

# Find Parent Contribution: This script figures out how much each parent contributes to allelic imbalance in NGS data
# Input: VCF File with Parents, GATK-output File with Ref and Alt allele columns
# User Input: Regions of interest to asses parental contribution
# Output: A .csv file with regions of interest and the main contributing parent per region/SNP

import sys
from collections import namedtuple

# load an initial VCF file
def loadVCF():
    vcf_file = raw_input("Input the path to your VCF File: ")
    return updateVCF(vcf_file)

# update VCF dictionary 
# @param: path to vcf file
# @return: a dictionary containing chr, pos, and par genotypes
def updateVCF(fileName):

    #initiate a dictionary of dictionaries (tuple --> tuple)
    vcf_repos = namedtuple("") 
    vcf = open(fileName)
    line = vcf.readline()

    # populate the dictionary with parent-SNP information
    while line:
        if line[0] != "c":
            line = vcf.readline()
        else:
            # the line contains information about SNPs at particular sites - so process
            # I just want to know if parent is ref or alt easily
            parsed = line.split()
            currChr = parsed[0]
            currSite = parsed[1]
            currPos = (currChr, currSite)
            vcf_repos[currPos]
            # process line here
            if parsed[9][0] == "1":
                parent1 = "alt"
                parent2 = "ref"
            else:
                parent1 = "ref"
                parent2 = "alt"
            print "At " + currChr + " " + currSite + " Parent 1 is " + parent1 + " and Parent 2 is " + parent2
            line = vcf.readline() 
    vcf.close() # once we're done populating dictionary, close file
    return vcf_repos

loadVCF()

## load a GATK output file
#def loadGATK():
#    GATK_file = raw_input("Input the path to your GATK output file (.csv): ")
#    return updateGATK(GATK_file)
#
## update GATK dictionary 
#def updateGATK():
#    GATK_repos = {}
#    # I think this is dictionary of chr --> dictionary [int --> tuple of ints? let's have default be the first parent 
#    # written in the VCF file...]
#    # PROCESS THE GATK FILE INTO DICTIONARY FORMAT
#    return GATK_repos
#
## get the user's query
#def getQuery():
#    query = raw_input("Input your regions of interest: ")
#    return query
#
## write the currently stored query to a file
#def queryToFile(query):
#    file_name = raw_input("Input the name of the output file to write: ")
#    # WRITE THE RESULTS OF QUERY...
#    print "Wrote " + file_name + " with your results!"
#
### main
##def main():
##    input=""
##    # while loop this shit until input = "X"
##    while True:
##        print "#####################################"
##        print "WELCOME TO FIND PARENTAL CONTRIBUTION"
##        print "#####################################"
##        print "What do you want to do? Enter one of the options"
##        choice = raw_input("Upload VCF: V  Upload GATK: G  Calculate Parent Contribution: P  Exit: X")
##        if choice == "V":
##            # Do stuff with the VCF
##        elif choice == "G": # need some conditional statements... 
##            # Do stuff with GATK
##        elif choice == "X":
##            print "Okay, closing FIND PARENTAL CONTRIBUTION"
##            break
##        else:
##            print "Sorry, not an option. Try again."
### close
### call main becuase nothing with happen otherwise
##    #main()
