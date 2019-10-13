'''
TODO:
A function that translates DNA to Amino Acidse
Find coordinates of open reading frames

DONE:
parse fasta
parse gb
File into fasta
find CDS list
Finds exons
'''
import re

def parse_fasta(FileName):
	#open the file
    file = open(FileName, 'r')
    file.close()
    #Split the content of the file into a list based newline characters and remove the descriptor line
    DNA_sequence = file.readlines()
    DNA_sequence = DNA_sequence[1:]
    #Join the list into a string
    DNA_String = ''.join(DNA_sequence)
    #remove the newline chars
    DNA_String = DNA_String.replace("\n", "")
    #Makes it uppercase
    DNA_String = DNA_String.upper()
    
    return DNA_String

def parse_gb(FileName):
	#open file
	file = open(FileName, 'r')
	gb_file = file.read()
	file.close()

	#extract the DNA sequence
	origin = gb_file.find("ORIGIN")
	end = gb_file.find('//')
	file_Chunk = gb_file[origin:end]
	Chunk_lines = file_Chunk.split('\n')
	
	DNA_seq = ''
	#for each line in the relevant seciton, split it into a list based on whitespace and then append the DNA portion of the line to DNA_seq
	for i in Chunk_lines:
		lst = i.split()
		DNA_seq += ''.join(lst[1:])
	return DNA_seq.upper()

def new_fasta(newFileName, body, heading):
	newfile = open(newFileName, 'w+')
	heading += '\n'
	newfile.write(heading)
	newfile.write(body)
	newfile.close()

def convert_CDS(FileName):
	#open file
	f = open(FileName, 'r')
	File = f.read()
	f.close()

	#Extract the CDS list
	CDS_begin = File.find("CDS")
	CDS_end = File.find('/gene', CDS_begin)
	CDS_chunk = File[CDS_begin:CDS_end]

	#Convert the CDS list into a Python formatted list 
	CDS_chunk = re.sub("\s", "", CDS_chunk)
	CDS_chunk = re.sub("CDSjoin", "", CDS_chunk)
	CDS_chunk = re.sub("\.\.", ",", CDS_chunk)
	CDS_chunk = re.sub("[\(\)]", "", CDS_chunk) 
	CDS_list = CDS_chunk.split(",")

	#Convert the Strings in the list to ints
	for i in range(0,len(CDS_list)):
		CDS_list[i] = int(CDS_list[i])

	return CDS_list

def find_exons(FileName, list = True):
	CDS_list = convert_CDS(FileName)

	exon_sequence = ""
	exon_sequence_list = []
	DNA_seq = parse_gb(FileName)
	# now pull it out
	for index in range(0,len(CDS_list),2):
	    begin = int(CDS_list[index]) - 1
	    end = int(CDS_list[index+1])
	    exon_sequence += DNA_seq[begin:end]
	    exon_sequence_list.append(DNA_seq[begin:end])
	if(list):
		return exon_sequence_list
	else:
		return exon_sequence

def find_introns(FileName, list = True):

	#Gets the CDS list and the DNA sequence
	CDS_list = convert_CDS(FileName)
	DNA_seq = parse_gb(FileName)

	#variable declaration
	intron_coor = []
	intron_seq = ""
	intron_list = []

	#
	for i in range(0,len(CDS_list)-1, 2):
		intron_coor.append(int(CDS_list[i+1]))
		intron_coor.append(int(CDS_list[i+1])-1)

	for index in range(0,len(intron_coor),2):
	    begin = int(intron_coor[index]) - 1
	    end = int(intron_coor[index+1])
	    intron_seq += DNA_seq[begin:end]
	    intron_list.append(DNA_seq[begin:end])
	if(list):
		return intron_list
	else:
		return intron_seq

def translate(DNA):
	genecode = {
	    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
	    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

	prot_seq = ''
	for i in range(0,len(DNA), 3):
		codon = DNA[i:i+3]
		AA = genecode[codon]
		prot_seq += AA

	return prot_seq


