
### for a given list of polymorphic sites, and a list of fasta format sequences, mask polymorphic sites with a P. 
### both lists organised in exon order.

def mask_pol_sites(pol_sites, protein_seq, species_):
	protein_pol_masked = []
	for x, b in (zip(pol_sites[species_],[x for x in protein_seq if species_ in x])):
		variant_sites = []

		for a in x: ### a = variants per exon
			if len(a) == 0:
				b_replace = b.split('\n')[1]
			elif len(a) != 0:
				a = [int(x) for x in a]
				if '&' not in b:
					species_coding_range = (int(b.split('_')[3]), int(b.split('_')[4].split('[')[0]))
					species_coding_range = sorted(species_coding_range)
				else:
					coding_range = (int(b.split('&')[0].split('_')[3]), int(b.split('&')[1].split('_')[3]), int(b.split('&')[0].split('_')[4].split('[')[0]), int(b.split('&')[1].split('_')[4].split('[')[0]))				
					species_coding_range = (sorted(coding_range)[0], sorted(coding_range)[-1])	
# 				print(species_coding_range)

			for c in a:
				c_place = (range(species_coding_range[0], species_coding_range[1]+1).index(c))
# 				print(c_place)
				variant_sites.append(c_place)	
										
		direction = b.split('\t')[4]
		b_replace = b.split('\n')[1]
		count = 0 ### cannot perform easy string splicing due to gap characters in alignments (-). Count positional information, ignoring these.

		if direction == '+':
			for num, char in enumerate(b_replace):
				if char != '-':
					count = count + 1
					if any([x == count for x in variant_sites]):
						### have checked manually. Due to inclusive counting of nt positions
						b_replace = b_replace[:int(num+1)] + 'P' + b_replace[int(num+2):] 
									
		elif direction == '-':
			for num, char in enumerate(reversed(b_replace)):
				if char != '-':
					count = count + 1
					if any([x == count for x in variant_sites]):
						### have checked manually: not intuitive. Reversed string indexing and inclusive nt counting 
						b_replace = b_replace[:int(num+2)*(-1)] + 'P' + b_replace[(int(num+1))*(-1):]

# 		print(b_replace)
		protein_pol_masked.append(b.split('\n')[0]+ '\n' + b_replace)

	return protein_pol_masked
	
	