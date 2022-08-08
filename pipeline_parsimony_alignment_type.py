### Created October 3rd. This new version will also identify, and write to file, polymorphism variant data

### Note on calling the script
### The script takes 2 arguments: chromosome number, as an integer, and regions in bp, as integers seperated by '-'.
### Scientific numbers are expected. An example: 1 10e6-50e6. User arguments are also used to label output files.
###
### The reference genome used is mouse, which is the following:
### chr_regions = {1:'1-196e6', 2:'1-182.5e6', 3:'1-160.5e6', 4:'1-157e6', 5:'1-152e6', chrX':'1-171.5e6'
### 6:'1-150e6', 7:'1-145.5e6', 8:'1-129.5e6', 9:'1-125e6', 10:'1-131e6', 11:'1-122.5e6',
### 12:'1-120.5e6', 13:'1-120.5e6', 14:'1-125e6',
### 15:'1-104.5e6', 16:'1-98.5e6', 17:'1-95e6', 18:'1-91e6', 19:'1-61.5e6', 'chrY':'1-92e6'}
###
### this code is specific to the quadruplet mouse-rat-human-chimp. To generalise would take a bit of work, alternatively these species need to be substituted 
### with others throughout the code- including latin names.
### written using python 3 - (not my native-running through pipenv shell)

### packages
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from collections import defaultdict
import re
import requests, sys
import time, datetime

### own functions - code requires the following general functions:
from translation_table_function import universal ### translation table to use
from translation_table_function import codon_translations
from translation_table_function import translate_string
from DNAcombine_function import combine_DNA
from DNAcombine_function import combine_exon
from mask_polymorphisms_function import mask_pol_sites
from gettriplets_function import get_triplets

### a quick note on reading this code: throughout, and in variable names, I use exons and CDS regions somewhat interchangeably.
### I always mean CDS- CDS regions are straightforward to translate into protein, and tend to include trailing stop codons.

# verbose setup- will print some additional output if true, helpful for tracking progress of programme.
# if set to false, the programme takes a while and prints very little.
verbose = False
verboseprint = print if verbose else lambda *a, **k: None

# calling function to get a dictionary of amino acids as keys, with their codons as values.
codon_codes = codon_translations(universal)

dir = "/Users/jennyjames/Dropbox/MASEL/Parsimony_analysis/"
server = "http://rest.ensembl.org" ### the ensembl database used for calls!

chr = str(sys.argv[1])
vals = str(sys.argv[2])


####################################### CHECKPOINT ########################################
################################# setup over- get started #################################
###########################################################################################


for i in range(int(float(vals.split('-')[0]))-1+5000000, int(float(vals.split('-')[1]))-1+5000000, 5000000):

	with requests.Session() as s:
		region = str(i-4999999)+'-'+str(i)

		genes = defaultdict(list) ### dictionary of genes in region, key as protein name and exons as list
								  ### exons tab delimited start, end, score (sometimes), strand, phase

		startTime = datetime.datetime.now()
		print("It is now : %s"%startTime)
		print(str(chr) +' '+ region+' starting...')

		### get exon positions in the chromosome (defined by mouse) and region defined above
		### ensembl sometimes needs time to process a request- if the problem is server-side, wait and try again.
		ext = "/overlap/region/mus_musculus/"+str(chr)+":"+region+"?feature=cds"
		while True:
			try:
				verboseprint('trying CDS call')
				verboseprint(server+ext)
				r = s.get(server+ext, headers={"Content-Type" : "text/x-gff3"})
				r.raise_for_status()
			except requests.exceptions.HTTPError as err:
				if r.status_code == 503:
					verboseprint('waiting 10 seconds and trying again...')
					time.sleep(11)
					continue
				else:
					r.raise_for_status()
					sys.exit()
			else:
				break

		exons = [x for x in r.text.split('\n')[2:] if len(x) > 1] ###
		for x in exons:
			genes[x.split('protein_id=')[1]].append(x.split('\tCDS\t')[1].split('\tID=CDS')[0])

		### removing gene duplicates: sequences with identical ORFs that are listed under multiple names
		### retaining splice variants that result in different CDS (such as genes that appear to be missing exons)
		
		genes_to_remove = []
		unique_genes = ['something']
		for gene, exons in genes.items():
			if exons == unique_genes[-1]:
				genes_to_remove.append(gene)
			else:
				unique_genes.append(exons)
		for gene in genes_to_remove:
			del genes[gene]

		verboseprint(genes.items())	
		with open('mouse_reference_gene_list.txt', 'a') as notefile:
			for gene, exons in genes.items():
				notefile.write(gene + ', ')
				notefile.write(', '.join(exons) + '\n')
		

		gene_sequences = defaultdict(list) ### dictionary with key of protein name as in genes, with aligned sequences:
										   ### aligned sequences organised in their exons nested list of lists.
										   ### some exons are extracted as alignments in multiple sections.
										   ### Alignment sequences are reverse complimented if needed.
										   ### sequence names end with '_exon_info' mus musculus positions start, end, score, strand, phase for reference
										   ### some genes have very similar CDSs, e.g. if there are splice variants, but will not be identical.



		gene_sequences_variants_masked = defaultdict(list)	### as above, however, polymorphic sites in any species have been masked with a P.

		for gene, exons in genes.items():

			for exon in exons:
				Get_exon = "/alignment/region/mus_musculus/"+str(chr)+":"+exon.split('\t')[0]+'-'+exon.split('\t')[1]+":1?species_set_group=amniotes;display_species_set=mus_musculus;display_species_set=rattus_norvegicus;display_species_set=homo_sapiens;display_species_set=pan_troglodytes;method=PECAN"

				while True:
					verboseprint('trying alignment call')
					exon_r = s.get(server+Get_exon, headers={"Content-Type" : "text/x-phyloxml"})
					if exon_r.status_code == 200:
						verboseprint('good request ' + server + Get_exon)
						break
					elif exon_r.status_code == 503: ### server error
						verboseprint(exon_r.status_code)
						verboseprint('waiting 10 seconds and trying again...')
						time.sleep(11)
					else:
						verboseprint(exon_r.status_code) ### likely client error
						verboseprint('bad request ' + server + Get_exon)
						break

				if not exon_r.ok: ### loop broke on response, this CDS region does not exist in the alignment. Break to next exon
					break

				exon_seq = exon_r.text
				exon_alignments_fasta = []

				### identify aligned regions from our focal species. These must have an aligned region for all 4 focal species, hence the following filters.
				pattern = re.compile('<name>.*?</mol_seq>', re.DOTALL)
				align_regions = re.findall(pattern, exon_seq)
				align_regions = [x for x in align_regions if 'Ancestor' not in x]

				if len(align_regions) % 4 == 0 and len(align_regions) >= 4: ### all 4 species present for every aligned region present.
					species = ('mus_musculus', 'rattus_norvegicus', 'homo_sapiens', 'pan_troglodytes')
					evaluator = True
					for sp in species:
						if any(sp in x for x in align_regions):
							pass
						else:
							evaluator = False 
					if evaluator == True:
				
				
						for x in align_regions:
							verboseprint(x)
							x = x.split('\n')
							for y in x:
								if '<name>' in y:
									name = '>'+y.split('<name>')[1].split('</name>')[0]
								if  '<mol_seq' in y:
									seq = y.split('>')[1].split('<')[0] ### extract just DNA sequence
									seq = Seq(seq, generic_dna)	### biopython Seq object
									if exon.split('\t')[3] == '-': ### orientation, from mouse reference.
										seq_format = str(seq.reverse_complement()) ### use biopython to reverse complement strands
									else:
										seq_format = str(seq)

									exon_fasta = name + '_mouse_ref\t' + exon + '\n' + seq_format + '\n'

									exon_alignments_fasta.append(exon_fasta)


					gene_sequences[gene].append([x for x in exon_alignments_fasta])


####################################### CHECKPOINT ########################################
###### have identified all the sequences present in all species in the alignment file #####
###########################################################################################


		for gene, protein in gene_sequences.items():
		
			SNP_identifiers = defaultdict(list) ### all sites polymorphic in the exon in at least one species, organised as species, identifiers

			SNP_freq_passed = defaultdict(list)	### all polymorphic sites that pass a MAF threshold of 5%

			SNP_details_passed = defaultdict(list) ### details of polymorphisms, as a nested dictionary: gene, organised with species as keys, and variants detail rows to write as values.

			new_protein = [] ### combining sequences that make up one exon, for each species for each gene

			for exon_alignment in protein:
				### call function to combine partial exon sequences to make up complete exons.

				exon_positions = 'none'
				if len(exon_alignment) > 4:
					direction = [seq.split('\t')[4] for seq in exon_alignment]
					mouse_refs = [(int(seq.split('_')[3]),int(seq.split('_')[4].split('[')[0]))  for seq in exon_alignment if 'mus_musculus' in seq]
					
					reverse = (mouse_refs[0][0]-mouse_refs[1][1])
					forward = (mouse_refs[0][1]-mouse_refs[1][0])
					
					if direction[0] == '+':					 
						if reverse == 1 or reverse == -1:
							exon_positions = 'backward'
						if forward == 1 or forward == -1:
							exon_positions = 'forward'
					elif direction[0] == '-':
						if reverse == 1 or reverse == -1:
							exon_positions = 'forward'
						if forward == 1 or forward == -1:
							exon_positions = 'backward'
												
				### sometimes returning nothing?
				mouse_string = combine_DNA(exon_alignment, 'mus_musculus', exon_positions)
				rat_string = combine_DNA(exon_alignment, 'rattus_norvegicus', exon_positions)
				chimp_string = combine_DNA(exon_alignment, 'pan_troglodytes', exon_positions)
				human_string = combine_DNA(exon_alignment, 'homo_sapiens', exon_positions)
				
				
				
				if len(mouse_string) > 5:
					new_protein.append(mouse_string)
					new_protein.append(rat_string)
					new_protein.append(chimp_string)
					new_protein.append(human_string)
			
			for species in new_protein:
							
				species_name = species.split('_')[0].split('>')[1]+'_'+species.split('_')[1]
				
				species_chr = species.split('_')[2]
				if '&' not in species:
					species_coding_range = species.split('_')[3]+'-'+species.split('_')[4].split('[')[0]
				else:
					coding_range = (int(species.split('&')[0].split('_')[3]), int(species.split('&')[1].split('_')[3]), int(species.split('&')[0].split('_')[4].split('[')[0]), int(species.split('&')[1].split('_')[4].split('[')[0]))
					species_coding_range = str((sorted(coding_range))[0])+'-'+str((sorted(coding_range))[-1])

				### find all variation that overlaps CDS, per species
				ext = "/overlap/region/"+species_name+"/"+species_chr+":"+species_coding_range+"?feature=variation"

				while True:
					try:
						verboseprint('trying variant call')
						variant_r = s.get(server + ext, headers={"Content-Type" : "text/x-gff3"}) ### return as tab delimited line as in VCF
						variant_r.raise_for_status()
					except requests.exceptions.HTTPError as err:
						if variant_r.status_code == 503 or variant_r.status_code == 504 or variant_r.status_code == 500:
							verboseprint('variants server trouble, waiting 10 seconds and trying again...') # server error- we can wait it out like before
							time.sleep(11)
							continue
						else:
							variant_r.raise_for_status() # something else has gone wrong- this is our problem...
							sys.exit()
					else:
						break

				SNP_name = re.findall('ID=sequence_variant:(.+);alleles.+;consequence_type=missense_variant',variant_r.text)

				SNP_identifiers[species_name].append(SNP_name)
				

####################################### CHECKPOINT ########################################
###############  identify all variant sites in the exon in any of the species  ############
###########################################################################################

			for species, variants_per_exon in SNP_identifiers.items():
				
				ext = "/variation/"+species+"?pops=1;population_genotypes=0;genotypes=0;phenotypes=0"
				headers={"Content-Type" : "text/x-xml"}
				for z in variants_per_exon:
					
					# holder lists for variants that occur in an exon. SNP_passed= variant position, SNP_details= variant details
					SNP_passed = []
					SNP_details = []
					
					if len(z)== 0:
						SNP_passed.append([])
												
					elif len(z) < 200:
						str_variants = '"' + '", "'.join([x for x in z]) + '"'
						str_variants = '{ "ids" : [' + str_variants + ' ] }'

						### get details for all variants that overlap the exon. There is a possible error to throw:
						### there is a maximum size limit for this request of 200 variants
						### store results in the SNP_freq_passed dictionary

						while True:
							try:
								verboseprint('trying variant details call')
								variant_details_r = s.post(server+ext, headers=headers, data=str_variants)
								variant_details_r.raise_for_status()
							except requests.exceptions.HTTPError as err:
								if variant_details_r.status_code == 503: # wait it out
									verboseprint('waiting 10 seconds and trying again...')
									time.sleep(11)
									continue
								else:
									variant_details_r.raise_for_status() # client error
									sys.exit()
							else:
								break

						variant_details_r = variant_details_r.text
						count = 0
				
						
						for var in variant_details_r.split('var_class: SNP\n'):
							frequency = re.findall('frequency: &#39;(.+)&#39;', var)
							
							## if the variant reaches 5% threshold in ANY population
							if any([0.5 > float(x) >= 0.05 for x in frequency]):
							
								count = count + 1
								
								var_freq = re.findall('frequency: &#39;(.+)&#39;', var)
								var_nt = re.findall('allele_string: (./.)', var)
								var_position = re.findall('      end: (.+)\n', var)
								var_name = re.findall('(.+): \n  MAF', var)
											
								polymorphism_string = ';'.join(var_freq) + ',' + ';'.join(var_nt) + ',' + ';'.join(var_position) + ',' + ';'.join(var_name)		
								
								SNP_details.append(polymorphism_string)
								
								for x in var_position:
									SNP_passed.append([x for x in var_position])
									
								count = count + 1
								
						if count == 0:
							SNP_passed.append([])
		#					SNP_details.append([])


					else: ### count of variants is over 200, break up into groups.
						variant_details_list = []
						z_groups = float(len(z))/float(190)
						z_groups = round(z_groups + 0.4)
						for x in range(190, len(z)+1, 190):
							xgroup = z[x-190:x]
							str_variants = '"' + '", "'.join([x for x in xgroup]) + '"'
							str_variants = '{ "ids" : [' + str_variants + ' ] }'

							while True:
								try:
									verboseprint('trying variant details call')
									variant_details_r = s.post(server+ext, headers=headers, data=str_variants)
									variant_details_r.raise_for_status()
								except requests.exceptions.HTTPError as err:
									if variant_details_r.status_code == 503:
										verboseprint('waiting 10 seconds and trying again...')
										time.sleep(11)
										continue
									else:
										variant_details_r.raise_for_status()
										sys.exit()
								else:
									break

							variant_details_list.append(variant_details_r.text)

						variant_details_r = ''.join([x for x in variant_details_list])						
						count = 0
					
						
						for var in variant_details_r.split('var_class: SNP\n'):
							frequency = re.findall('frequency: &#39;(.+)&#39;', var)
							## if the variant reaches 5% threshold in ANY population
							if any([0.5 > float(x) >= 0.05 for x in frequency]):
							
								count = count + 1
								
								var_freq = re.findall('frequency: &#39;(.+)&#39;', var)
								var_nt = re.findall('allele_string: (./.)', var)
								var_position = re.findall('      end: (.+)\n', var)
								var_name = re.findall('(.+): \n  MAF', var)
											
								polymorphism_string = ';'.join(var_freq) + ',' + ';'.join(var_nt) + ',' + ';'.join(var_position) + ',' + ';'.join(var_name)	
								
								SNP_details.append(polymorphism_string)
								
								for x in var_position:
									SNP_passed.append([x for x in var_position])
								count = count + 1
								
						if count == 0:
							SNP_passed.append([])
		
					SNP_freq_passed[species].append(SNP_passed)
					SNP_details_passed[species].append(SNP_details)
					
									
			mouse_pols = mask_pol_sites(SNP_freq_passed, new_protein, 'mus_musculus')
			rat_pols = mask_pol_sites(SNP_freq_passed, new_protein, 'rattus_norvegicus')
			chimp_pols = mask_pol_sites(SNP_freq_passed, new_protein, 'pan_troglodytes')
			human_pols = mask_pol_sites(SNP_freq_passed, new_protein, 'homo_sapiens')

			masked_chunk_list = []
			for i in range(len(mouse_pols)):
				masked_chunk_list.append([mouse_pols[i], rat_pols[i], chimp_pols[i], human_pols[i]])
			gene_sequences_variants_masked[gene] = masked_chunk_list
									
			chunked_list = []
			for i in range(0, len(new_protein), 4):
				chunked_list.append(new_protein[i:i + 4])
			gene_sequences[gene] = chunked_list
			
####################################### CHECKPOINT ########################################
############# write details of polymorphic sites to file, but no translation  #############
###########################################################################################
									
			with open(dir + str(chr)+':'+vals+'_'+'rat_polymorphisms.csv', 'a') as rat_polfile, open(dir + str(chr)+':'+vals+'_'+'mouse_polymorphisms.csv', 'a') as mouse_polfile, open(dir + str(chr)+':'+vals+'_'+'chimp_polymorphisms.csv', 'a') as chimp_polfile, open(dir + str(chr)+':'+vals+'_'+'human_polymorphisms.csv', 'a') as human_polfile:
					
				### zip together SNP variant details and their corresponding exons			
				mouse_info = zip(SNP_details_passed['mus_musculus'], mouse_pols)
				rat_info = zip(SNP_details_passed['rattus_norvegicus'], rat_pols)
				chimp_info = zip(SNP_details_passed['pan_troglodytes'], chimp_pols)
				human_info = zip(SNP_details_passed['homo_sapiens'], human_pols)
				
				
				### write variant details to file, along with their context: one row per variant
				for line in mouse_info:
					polymorphisms, sequence = line[0], line[1]
					if len(polymorphisms) >= 1:
						for y in polymorphisms:
							mouse_polfile.write(gene + ',' + y + ',' + re.sub('\n', '\t', sequence) + '\n')

				for line in rat_info:
					polymorphisms, sequence = line[0], line[1]
					if len(polymorphisms) >= 1:
						for y in polymorphisms:
							rat_polfile.write(gene + ',' + y + ',' + re.sub('\n', '\t', sequence) + '\n')

				for line in chimp_info:
					polymorphisms, sequence = line[0], line[1]
					if len(polymorphisms) >= 1:
						for y in polymorphisms:
							chimp_polfile.write(gene + ',' + y + ',' + re.sub('\n', '\t', sequence) + '\n')

				for line in human_info:
					polymorphisms, sequence = line[0], line[1]
					if len(polymorphisms) >= 1:
						for y in polymorphisms:
							human_polfile.write(gene + ',' + y + ',' + re.sub('\n', '\t', sequence) + '\n')


####################################### CHECKPOINT ########################################
############## identify useful substitution sites, and write results to file  #############
###########################################################################################

		with open(dir + 'amniotePECAN' + 'mouse_reference_alignment_gene_list.txt', 'a') as ref_file, open(dir + 'amniotePECAN'+ str(chr)+':'+vals+'_'+'rat_results.csv', 'a') as rat_file, open(dir+ 'amniotePECAN' + str(chr)+':'+vals+'_'+'mouse_results.csv', 'a') as mouse_file, open(dir+ 'amniotePECAN' + str(chr)+':'+vals+'_'+'chimp_results.csv', 'a') as chimp_file, open(dir+ 'amniotePECAN' + str(chr)+':'+vals+'_'+'human_results.csv', 'a') as human_file:

			species = ('mus_musculus', 'rattus_norvegicus', 'pan_troglodytes', 'homo_sapiens')
	
			for gene, CDS in gene_sequences_variants_masked.items():
			
				ref_file.write(gene + ', ')
				for exon in CDS:
					ref_file.write(exon[0].split('\n')[0] + ', ')
				ref_file.write('\n')
			
				for exons, exon in zip(gene_sequences[gene], CDS):
				
					mouse_exons = []
					rat_exons = []
					chimp_exons = []
					human_exons = []
						
					mouse_exon = []
					rat_exon = []
					chimp_exon = []
					human_exon = []	
								
					sequence = (x.split('\n')[1] for x in exons)
					for x in sequence:
						print(x)
					sequences = zip(*sequence)
					
					mouse_exons.append(exons[0].split('\n')[0] + '\n')
					rat_exons.append(exons[1].split('\n')[0] + '\n')
					chimp_exons.append(exons[2].split('\n')[0] + '\n')
					human_exons.append(exons[3].split('\n')[0] + '\n')
										
# 					sequence = zip(x for x in sequence)
									
					for position in sequences:
						if position != ('-', '-', '-', '-'): 
							mouse_exons.append(position[0])
							rat_exons.append(position[1])
							chimp_exons.append(position[2])
							human_exons.append(position[3])
					
					mouse_exons = ''.join(mouse_exons)				
					rat_exons = ''.join(rat_exons)
					chimp_exons = ''.join(chimp_exons)
					human_exons = ''.join(human_exons)
					
					mouse_protein = ''.join(translate_string(mouse_exons))
					rat_protein = ''.join(translate_string(rat_exons))
					chimp_protein = ''.join(translate_string(chimp_exons))
					human_protein = ''.join(translate_string(human_exons))

					sequence = (x.split('\n')[1] for x in exon)
					sequences = zip(*sequence)
					
					mouse_exon.append(exon[0].split('\n')[0] + '\n')
					rat_exon.append(exon[1].split('\n')[0] + '\n')
					chimp_exon.append(exon[2].split('\n')[0] + '\n')
					human_exon.append(exon[3].split('\n')[0] + '\n')
										
# 					sequence = zip(x for x in sequence)
									
					for position in sequences:
						if position != ('-', '-', '-', '-'): 
							mouse_exon.append(position[0])
							rat_exon.append(position[1])
							chimp_exon.append(position[2])
							human_exon.append(position[3])
					
					mouse_exon = ''.join(mouse_exon)				
					rat_exon = ''.join(rat_exon)
					chimp_exon = ''.join(chimp_exon)
					human_exon = ''.join(human_exon)
								
					mouse_masked_protein = ''.join(translate_string(mouse_exon))
					rat_masked_protein = ''.join(translate_string(rat_exon))
					chimp_masked_protein = ''.join(translate_string(chimp_exon))
					human_masked_protein = ''.join(translate_string(human_exon))
			
					print(mouse_masked_protein)
					print(rat_masked_protein)
					print(chimp_masked_protein)
					print(human_masked_protein)

					zip_seqs = zip(mouse_masked_protein, rat_masked_protein, chimp_masked_protein, human_masked_protein)
					
					### retain a count- so relative position of substitution is known.
					count = 0
					for aas in zip_seqs:
						count = count + 1
						if len(set(aas)) != 1: 
						
							if '?' not in aas and '-' not in aas:
						
						
								if aas[0] != aas[1] and aas[0] != aas[2] and aas[0] != aas[3]:
									mouse_new = (aas[0])
									mouse_anc = aas[3]
									mouse_pos = count
									mouse_file.write(gene + ',')
									mouse_file.write(mouse_new + ',')
									mouse_file.write(mouse_anc + ',')
									mouse_file.write(mouse_protein + ',')
									mouse_file.write(str(mouse_pos) + ',')
									mouse_file.write(mouse_masked_protein + ',')
									mouse_file.write(exon[0].split('\n')[1] + ',')
									mouse_file.write(exon[0].split('\n')[0] + '\n')


								if aas[1] != aas[0] and aas[1] != aas[2] and aas[1] != aas[3]:
									rat_new = (aas[1])
									rat_anc = aas[3]
									rat_pos = count
									rat_file.write(gene + ',')
									rat_file.write(rat_new + ',')
									rat_file.write(rat_anc + ',')
									rat_file.write(rat_protein + ',')
									rat_file.write(str(rat_pos) + ',')
									rat_file.write(rat_masked_protein + ',')
									rat_file.write(exon[1].split('\n')[1] + ',')
									rat_file.write(exon[1].split('\n')[0] + '\n')


								if aas[2] !=  aas[0] and aas[2] != aas[1] and aas[2] != aas[3]:
									chimp_new = aas[2]
									chimp_anc = aas[0]
									chimp_pos = count
									chimp_file.write(gene + ',')
									chimp_file.write(chimp_new + ',')
									chimp_file.write(chimp_anc + ',')
									chimp_file.write(chimp_protein + ',')
									chimp_file.write(str(chimp_pos) + ',')
									chimp_file.write(chimp_masked_protein + ',')
									chimp_file.write(exon[2].split('\n')[1] + ',')
									chimp_file.write(exon[2].split('\n')[0] + '\n')

								if aas[3] !=  aas[0] and aas[3] != aas[1] and aas[3] != aas[2]:
									human_new = aas[3]
									human_anc = aas[0]
									human_pos = count
									human_file.write(gene + ',')
									human_file.write(human_new + ',')
									human_file.write(human_anc + ',')
									human_file.write(human_protein + ',')
									human_file.write(str(human_pos) + ',')
									human_file.write(human_masked_protein + ',')
									human_file.write(exon[3].split('\n')[1] + ',')
									human_file.write(exon[3].split('\n')[0] + '\n')
						
								
		print('We got there! 5mb took: %s'%(datetime.datetime.now()-startTime))								
								
							
