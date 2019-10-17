# ComFunc
A module with the common functions that we use in CS 116: BIOINFORMATICS
Current Version = 1.0 

# Installation
1. Download Comfunc.py and put it in the same directory as your Jupyter Notebook files for BioInformatics
2. import Comfunc at the top of your program 

# Usage
## parse_fasta(FileName)
Takes the name of a fasta file and returns the DNA sequence in it as an uppercase String 

	import Comfunc
	DNA_String = Comfunc.parse_fasta('CFTR.fasta')
	print(DNA_String)
	
	output = 'AATTGGAAGCAAATGACATCACAGCAGGTCAGAGAAAAAGGGTTGAG...'

## parse_gb(FileName)
Takes the name of a Genbank file and returns the DNA sequence in it as an uppercase String
	
	import Comfunc
	DNA_String = Comfunc.parse_gb('hras.gb')
	print(DNA_String)

	output = 'CCACCCCGAGCCTAGAGAAGGCTCCTAGCTTGGCTGGGCCTCATGGGGCCCTTGTCTGCTT...'

## new_fasta(newFileName, body, heading)
Takes the name of the file you want to create, the body of the file(the DNA sequence), and the heading of the file creates a new fasta file
	
	import Comfunc
	DNA_String = 'CCACCCCGAGCCTAGAGAAGGCTCCTAGCTTGGCTGGGCCTCATGGGGCCCTTGTCTGCTT...'
	Comfunc.new_fasta('CFTR.fasta', DNA_String, '>CFTR coding')
	
## get_CDS(FileName)
Takes a genbank file and returns the CDS list as a python formatted list
	
	import Comfunc
	python_CDS_list = Comfunc.get_CDS('hras.gb')

## find_exons(FileName, return_list = True)
return_list:bool, optional - True will return a list of exons, false will return the full exon string
Takes a genbank file and returns either a list of the exon sequences if return_list == True
	
