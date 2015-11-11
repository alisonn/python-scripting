#! /usr/bin/env python 

## FILTER_LONGEST:  @param: temp.txt file containing only fasta headers from Trinity 
## @return: unique genes/contigs with longest isoform only in Trinity format
inputFile="/home/alisonn/01_assemblies/temp.txt"

## initiated nested dictionary: gene --> length --> associated isoform  
dict = {}; 

## first open the file 
f = open(inputFile)

## read the first line of the fasta header file 
line = f.readline() 

## while file is not empty, keep reading line by line until file (i.e. line) is empty
while line: 		
		
	## parse the fasta header to get only a unique gene name 
	tempGene = line.split()[0][1:] ## cut off ">" in fasta headers
	currGene = tempGene.split('_')[0]+"_"+tempGene.split('_')[1] ## retrieves TR*{0-9}|c0_g*{0-9}
	currIsoform = tempGene.split('_')[2] ## retrieves "i*{0-9}"
	currLength = line.split()[1][4:] ## takes "len=INT" and retrieves INT

	
	## if gene does not exist, add it and add length and assoc isoform
	if not (currGene in dict):
		dict[currGene]={}
		dict[currGene][currLength]=currIsoform					
		
	## if gene already exists, compare sequence lengths
	else : 
		## if new isoform length is longer than one in dictionary, replace
		for storedLength in dict[currGene]:		
			if (currLength > storedLength):
				dict[currGene].popitem() ## pop the only item 
				dict[currGene][currLength] = currIsoform ## add new value :-) ## 
				
	line = f.readline() ## get the next line, ahh you infinite loop.

## wahoo! we've reached the end of the file 	
f.close()


## write to a file called "outputList.txt" in current directory

contig="" 
f=open('outputList.txt', 'w')
for gene, value in dict.iteritems():
	for length in dict[gene]: 
		contig=gene+"_"+dict[gene][length]
		print contig
		f.write(contig+'\n')
f.close()