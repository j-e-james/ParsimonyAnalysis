# -*- coding: utf-8 -*-
import os, re


dir = "/Users/jennyjames/Dropbox/MASEL/Parsimony_analysis_Frameshift/"

def combine_files(speciesIdentifier, outfile):
	for file in os.listdir(dir):
		if speciesIdentifier in file:
			fileref = file.split('_')[0] 
			with open(file, 'r') as resfile:
				resfile= resfile.readlines()
				resfile = [re.sub('\n','',x) for x in resfile]
				resfile = [x+','+fileref+'\n' for x in resfile]
				for x in resfile:
					outfile.write(x)


with open('Mouse_substitutions.csv', 'a') as mousefile, open('Rat_substitutions.csv', 'a') as ratfile, open('Human_substitutions.csv', 'a') as humanfile, open('Chimp_substitutions.csv', 'a') as chimpfile:

	combine_files('mouse_results', mousefile)
	combine_files('rat_results', ratfile)
	combine_files('human_results', humanfile)
	combine_files('chimp_results', chimpfile)


	
with open('Mouse_polymorphisms.csv', 'a') as mousefile, open('Rat_polymorphisms.csv', 'a') as ratfile, open('Human_polymorphisms.csv', 'a') as humanfile, open('Chimp_polymorphisms.csv', 'a') as chimpfile:

	combine_files('mouse_polymorphisms', mousefile)
	combine_files('rat_polymorphisms', ratfile)
	combine_files('human_polymorphisms', humanfile)
	combine_files('chimp_polymorphisms', chimpfile)
				
				
