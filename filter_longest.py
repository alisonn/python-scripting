#!/usr/bin/env python 
#  
# Alison Hanh Nguyen: alisonn@berkeley.edu 
# Bachtrog Lab at UCBerkeley  
#
# FILTER_LONGEST: Takes in Trinity.fasta headers (must create a separate file) 
# and returns genes/contigs with their longest isoform in a .txt file. 
#
# Useful for filtering through many Trinity isoforms
#

inputFile = "/home/alisonn/01_assemblies/temp.txt"

# initiated nested repositionary: gene --> length --> associated isoform  
reposit = {}

# first open the file 
f = open(inputFile)

# read the first line of the fasta header file 
line = f.readline() 

# while file is not empty, keep reading line by line until file (i.e. line) is empty
while line: 		
		
	# parse the fasta header to get only a unique gene name 
	header = line[1:] # cut off > in headers
	splitHeader = header.split() #"TRfdfdfd|fdf_fdfd_fdf /split/ len=INT" 
	currLength = int(splitHeader[1][4:]) # takes "len=INT" and retrieves casted INT
	
	tempGene = splitHeader[0]
	splitGene = tempGene.split('_')
	currGene = splitGene[0]+"_"+splitGene[1] # retrieves TR*{0-9}|c0_g*{0-9}
	currIsoform = splitGene[2] # retrieves "i*{0-9}"	
	
	# if gene does not exist, add it and add length and assoc isoform
	if currGene not in reposit:
		reposit[currGene] = {}
		reposit[currGene][currLength] = currIsoform					
		
	# if gene exists, compare sequence lengths
	else: 		
		# if new isoform length is longer than one in repositionary, replace
		for storedLength in reposit[currGene]:		
			if (currLength > storedLength):
				reposit[currGene].popitem() # pop the only item 
				reposit[currGene][currLength] = currIsoform # add new value :-) # 
				
	line = f.readline() # get the next line, ahh you infinite loop.

# wahoo! we've reached the end of the file 	
f.close()


# write to a file called "outputList.txt" in current directory
f = open('outputList.txt', 'w')

for gene, value in reposit.iteritems():
	for length, isoform in reposit[gene]: 
		contig = gene+"_"+reposit[gene][length]
		f.write(contig+'\n')
		
f.close()