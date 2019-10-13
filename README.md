# ComFunc
A module with the common functions that we use in CS 116: BIOINFORMATICS

# Installation
Download Comfunc.py and put it in the same directory as your Notebook files for BioInformatics

#Usage
parse_fasta(FileName)
	Takes the name of a fasta file and returns the DNA sequence in it as an uppercase String 

	'DNA_String = parse_fasta('CFTR.fasta')
	print(DNA_String)
	
	output = 'AATTGGAAGCAAATGACATCACAGCAGGTCAGAGAAAAAGGGTTGAG...''

parse_gb(FileName)
	Takes the name of a Genbank file and returns the DNA sequence in it as an uppercase String
	EX.
	DNA_String = parse_gb('hras.gb')
	print(DNA_String)

	output = 'CCACCCCGAGCCTAGAGAAGGCTCCTAGCTTGGCTGGGCCTCATGGGGCCCTTGTCTGCTT...'

new_fasta(newFileName, body, heading
