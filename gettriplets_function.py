
## prepare the alignment fasta format strings to ID substitutions
def get_triplets(mouse, rat, human, chimp):
	zip_mask_seqs = zip(mouse.split('\n')[1], rat.split('\n')[1], chimp.split('\n')[1], human.split('\n')[1])
	zip_mask_seqs = (x for x in zip_mask_seqs if x != ('-', '-', '-', '-')) ### remove instances where '-' is present at a point in every species.
	zip_mask_seqs = ([x for x in zip_mask_seqs]) ### replacing the generator in order to index
	aligned_mask_triplets = ([zip_mask_seqs[i:i+3] for i in range(0, sum(1 for i in zip_mask_seqs), 3)])
	return aligned_mask_triplets