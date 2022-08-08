# ParsimonyAnalysis
All Python3 scripts used to identify substitutions unique to a single species, in the species quadruplet of human-chimp-mouse-rat. 
To be included in the results, the substitution had to be unique in one of the species, and identical in all other species. Other sites were ignored in this analysis.

Scripts also identify polymorphic sites in the species quadruplet and masks them from substitution results.

Main script: pipeline_parsimony_alignment_type.py

Requires the provided additional python functions, in addition to packages Bio.Seq, Bio.Alphabet, collections, re, requests, sys, time and datetime. 
This script pulls alignments from ensembl by chromosome number and position, which are arguments provided by the user. 
This script was run in chunks of 5mb (50e6 bp), for speed and to  prevent issues with the ensembl server, using the complete mouse genome as a reference for chromosome and position arguments. 

The script writes relatively large output files, with 1 produced per species in the quadruplet.
Each line of the results file gives the gene name, the species-specific substitution, the ancestral state as inferred from parsimony, the protein, the position in the protein of the subsitution, the protein with any ancestral sites masked, the reference of the site in the alignment, and exon information.


Post-processing of files can be done using combine_files.py.
