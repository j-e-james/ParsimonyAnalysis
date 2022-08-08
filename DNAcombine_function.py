import re

from translation_table_function import universal ### translation table to use
from translation_table_function import codon_translations ### write a defaultdict list of codons Ter is stop
from translation_table_function import translate_string 

codon_codes = codon_translations(universal)


### for exons, in list in correct order, return combined string for each species.
def combine_exon(DNA_list, species):### this function assumes the sequence data has already been reverse translated if required 
	species_string = []				### expecting fasta format, requires some directional (reading frame) information within the header
	for DNA_region in DNA_list:		
		if species in DNA_region:		
			species_string.append(DNA_region)
# 	print(species_string)
	title = '&'.join([x.strip('<').split('\n')[0] for x in species_string])
	seq = ''.join([x.split('\n')[1] for x in species_string]) 
	full_region = title + '\n' + seq
	return full_region


# 
# ### implicit assumption- there will only ever be 2 DNA regions? This might not be true...
# 
# def combine_DNA(DNA_list, species):### this function assumes the sequence data has already been reverse translated if required 
# 	species_string = []				### expecting fasta format, requires some directional (reading frame) information within the header
# 	for DNA_region in DNA_list:		
# 		if species in DNA_region:		
# 			species_string.append(DNA_region)
# 	title = '&'.join([x.strip('<').split('\n')[0] for x in species_string])
# 	seq = ''.join([x.split('\n')[1] for x in species_string])
# 	full_region = title + '\n' + seq 
# 	revtitle = '&'.join([x.strip('<').split('\n')[0] for x in reversed(species_string)])
# 	revseq = ''.join([x.split('\n')[1] for x in reversed(species_string)]) 
# 	full_rev_region = revtitle + '\n' + revseq
# 	return(full_region, full_rev_region)


def combine_DNA(DNA_list, species, order):### this function assumes the sequence data has already been reverse translated if required 

# 	print(order)
	
	species_string = []				### expecting fasta format, requires some directional (reading frame) information within the header
	for DNA_region in DNA_list:		
		if species in DNA_region:		
			species_string.append(DNA_region)
	title = '&'.join([x.strip('<').split('\n')[0] for x in species_string])
#	seq = re.sub('-','',''.join([x.split('\n')[1] for x in species_string]))
	seq = ''.join([x.split('\n')[1] for x in species_string])
	full_region = title + '\n' + seq 
	revtitle = '&'.join([x.strip('<').split('\n')[0] for x in reversed(species_string)])
#	revseq = re.sub('-','',''.join([x.split('\n')[1] for x in reversed(species_string)]))
	revseq = ''.join([x.split('\n')[1] for x in reversed(species_string)])
	full_rev_region = revtitle + '\n' + revseq
	if order == 'none' or order == 'forward':
		return(full_region)
	if order == 'backward':
		return(full_rev_region)



