###### Jenny James, last edited October 3rd 2019

###### this set of functions were written to return translated strings of amino acids, from an inputed nucleotide string.
###### the functions have self explanatory names, but codon_codes is required as setup. If a different translation table is needed, it can be specified.

import re
from collections import defaultdict
from collections import OrderedDict



### required setup to create the translation dictionary

universal = """
TTT F Phe      TCT S Ser      TAT Y Tyr      TGT C Cys  
TTC F Phe      TCC S Ser      TAC Y Tyr      TGC C Cys  
TTA L Leu      TCA S Ser      TAA * Ter      TGA * Ter  
TTG L Leu i    TCG S Ser      TAG * Ter      TGG W Trp  

CTT L Leu      CCT P Pro      CAT H His      CGT R Arg  
CTC L Leu      CCC P Pro      CAC H His      CGC R Arg  
CTA L Leu      CCA P Pro      CAA Q Gln      CGA R Arg  
CTG L Leu i    CCG P Pro      CAG Q Gln      CGG R Arg  

ATT I Ile      ACT T Thr      AAT N Asn      AGT S Ser  
ATC I Ile      ACC T Thr      AAC N Asn      AGC S Ser  
ATA I Ile      ACA T Thr      AAA K Lys      AGA R Arg  
ATG M Met i    ACG T Thr      AAG K Lys      AGG R Arg  

GTT V Val      GCT A Ala      GAT D Asp      GGT G Gly  
GTC V Val      GCC A Ala      GAC D Asp      GGC G Gly  
GTA V Val      GCA A Ala      GAA E Glu      GGA G Gly  
GTG V Val      GCG A Ala      GAG E Glu      GGG G Gly  
"""


def codon_translations(table):
	translation_dict = defaultdict(list)
	for v, k in re.findall('([A-Z][A-Z][A-Z]) (.) [A-Z][a-z][a-z]', table):
		translation_dict[k].append(v)
	return translation_dict

codon_codes = codon_translations(universal)



### translation functions

def translate_string_nogaps(string):
	string = re.sub('-', '', string)
	protein = []
# 	print(string)
	if '&' not in string:
		frame = int(string.split('\t')[5].split('\n')[0])
	else:
		frame = int(string.split('\t')[5].split('&')[0])
#	print(int(frame))
	codon_string = ([string.split('\n')[1][frame:][i:i+3] for i in range(0, len(string.split('\n')[1][frame:]), 3)])
#	print(codon_string[0])
	for codon in codon_string:
		if 'P' in codon: ### an identified polymorphic site!
			protein.append('?')
		elif 'N' in codon: ### an unknown nt in the sequencing data, replace with gaps.
			protein.append('-')		
		else:
			for k, v in codon_codes.items():
				if codon in v:
					protein.append(k)
	return(protein)


def translate_string(string):
	protein = []
	if '&' not in string:
		frame = int(string.split('\t')[5].split('\n')[0])
	else:
		frame = int(string.split('\t')[5].split('&')[0])
	codon_string = ([string.split('\n')[1][frame:][i:i+3] for i in range(0, len(string.split('\n')[1][frame:]), 3)])
	for codon in codon_string:
		if '-' in codon:
			protein.append('-')	
		if 'P' in codon: ### an identified polymorphic site!
			protein.append('?')
		elif 'N' in codon: ### an unknown nt in the sequencing data, replace with gaps.
			protein.append('-')		
		else:
			for k, v in codon_codes.items():
				if codon in v:
					protein.append(k)
	return(protein)


def translate_string_frame_checker(string):
	string = re.sub('-', '', string)
	protein = []
	count = 0
# 	print(string)
	codon_string = ([string.split('\n')[1][i:i+3] for i in range(0, len(string.split('\n')[1]), 3)])
	for codon in codon_string:
		if 'P' in codon: ### include other unknown nucleotides as exceptions here if required
			protein.append('?')
		else:
			for k, v in codon_codes.items():
				if codon in v:
					protein.append(k)											
	if '*' in protein[:-1]: ### the protein is not in frame.
		count = count + 1
		protein.clear()	
		codon_string1 = ([string.split('\n')[1][i:i+3] for i in range(1, len(string[1:].split('\n')[1]), 3)])
		for codon1 in codon_string1:
			if 'P' in codon1: ### include other unknown nucleotides as exceptions here if required
				protein.append('?')
			else:
				for k, v in codon_codes.items():
					if codon1 in v:
						protein.append(k)	
	if '*' in protein[:-1]: ### the protein is not in frame.
		count = count + 1
		protein.clear()	
		codon_string2 = ([string.split('\n')[1][i:i+3] for i in range(2, len(string[2:].split('\n')[1]), 3)])
		for codon2 in codon_string2:
			if 'P' in codon2: ### include other unknown nucleotides as exceptions here if required
				protein.append('?')
			else:
				for k, v in codon_codes.items():
					if codon2 in v:
						protein.append(k)
	if '*' in protein[:-1]:
		print('still not in frame. Problems!')
		Problems = "Stops"
			
	else:
		Problems = "No stops"
	return(Problems, protein, count) ### where count is frame	
