#!/usr/bin/env python

## take in tsv file from CAT --> tabulate the counts for genes found in YL, YM, YS etc
import sys
import os

inFile = open(sys.argv[1], 'r')
ys_cov = open(sys.argv[2], 'r')
ym_cov = open(sys.argv[3], 'r')
yl_cov = open(sys.argv[4], 'r')

line = inFile.readline()

## @ param: file name for sizes
## @ return: map of contig --> contig occurrences
def populateCovMap(covFile): 
	covData = open(covFile, 'r')
	curr_line = covData.readline()
	contig = curr_line.split()[0]
	copies = curr_line.split()[7]

	cov_map = {}

	while curr_line:
		cov_map[contig] = copies
		covData.readline()	

	covData.close()
	return cov_map

# initiate a dictionary for mod gff (location,gene) --> (alternative locations)
# @param: geneFile
# @return: a map (geneName --> (locations of the alt copies, as a list))
def populateGeneMap(geneFile):
	genes = open(geneFile, 'r')
	curr_line = genes.readline()

	geneMap = {}
	while curr_line:
		contig = curr_line.split()[0]
		src_gene = curr_line.split()[1]
		alt_loc = curr_line.split()[2].split(",") # retrive these locations as a list
		cat_name = curr_line.split()[3]

		keyStr = contig + "," + cat_name + ","  + src_gene
		geneMap[keyStr] = alt_loc

		curr_line = genes.readline()
	
	genes.close()
	return mainMap


def tabulateCopies(geneMap, copiesMap):
	
	finalMap = {}
	ys_count = 0
	ym_count = 0
	yl_count = 0
	auto_count = 0
	x_count = 0
	autosome_set = {"Muller_C", "Muller_E", "Muller_F", "Muller_B"}
	for key in geneMap:
	## KEY: contig,catName,srcGene 
	## VALUE: [location1,location2,location3....etc]
		## get the list of the YL/YM/YS locations, and tabulate occurrences	
		curr_locations = geneMap[key]
		for i in curr_locations:
			new_key = curr_locations[i]## key to check	
			chromosome = new_key.split(":")[0]

			if chromosome == "Muller_A-AD":
				x_count += 1

			elif autosome_set.contains(chromosome):
				auto_count += 1

			else:  # we have one of the Y chromosome contigs
				y_strain = chromosome.split("_")[0]
				y_contig = chromosome.split("_")[1] + "_" + chromosome.split("_")[2]
				relevantMap = copiesMap[y_strain] ## extract the map that is relevant
				copies = relevantMap[y_contig]
				if y_strain == "YL":
					yl_count += 1
				elif y_strain == "YM":
					ym_count += 1
				else:
					ys_count += 1	

		## final check - check if the key is YS/YM/YL, then check also

	return finalMap

def checkCopies(strain, copiesMap):

	return copyNumber

def main():
	mapYL = populateCovMap("LarY_cffc_asm_v2_maskedAutoX.copyNumber.tsv")
	mapYM = populateCovMap("MedY_cffc_asm_v2_maskedAutoX.copyNumber.tsv")
	mapYS = populateCovMap("SmaY_cffc_asm_v2_maskedAutoX.copyNumber.tsv")
	copiesMap = {"YL":mapYL, "YM":mapYM, "YS":mapYS}
	## ADD FILE NAME HERE 
	geneMap = populateGeneMap()
	tabulateCopies(geneMap, copiesMap)
