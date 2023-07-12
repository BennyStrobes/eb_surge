import numpy as np 
import os
import sys
import pdb




def extract_hvg(hvg_file):
	f = open(hvg_file)
	hvg = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		hvg[line] = 1
	f.close()
	return hvg 


def extract_array_of_gene_and_top_variants_according_to_standard_eqtl_results(standard_eqtl_results_file, hvg):
	f = open(standard_eqtl_results_file)
	head_count = 0
	gene_to_variant_and_pvalue = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# extract relevent fields
		gene_id = data[0]
		variant_id = data[4]
		pvalue = float(data[5])

		# Ignore genes not in HVG
		if gene_id not in hvg:
			continue

		# Add to dictionary
		if gene_id not in gene_to_variant_and_pvalue:
			gene_to_variant_and_pvalue[gene_id] = (variant_id, pvalue)
		else:
			old_pvalue = gene_to_variant_and_pvalue[gene_id][1]
			if pvalue < old_pvalue:
				gene_to_variant_and_pvalue[gene_id] = (variant_id, pvalue)
	f.close()

	tuple_arr = []
	for gene_id in [*gene_to_variant_and_pvalue]:
		tuple_arr.append((gene_id, gene_to_variant_and_pvalue[gene_id][0], gene_to_variant_and_pvalue[gene_id][1]))

	tuple_arr_sorted = sorted(tuple_arr, key=lambda tup: tup[2])

	kk = 1
	sig_vec = []
	num_genes = len(tuple_arr_sorted)
	fdr_thresh = .05
	for gene_tuple in tuple_arr_sorted:
		bf_pvalue = gene_tuple[2]
		fdr = num_genes*bf_pvalue/kk 
		kk = kk + 1
		if fdr > fdr_thresh:
			sig = False
		else:
			sig = True
		sig_vec.append(sig)

	pdb.set_trace()









######################
# Command line args
######################
test_names_file = sys.argv[1]  # output file
standard_eqtl_results_file = sys.argv[2]  # File containing standard eqtl results for all variant-gene pairs
hvg_file = sys.argv[3]  # Highly variable genes file


# First extract highly variable genes
hvg = extract_hvg(hvg_file)

# Extract array of gene and the gene's top variant and the corresponding pvalue
all_variant_gene_pairs = extract_array_of_gene_and_top_variants_according_to_standard_eqtl_results(standard_eqtl_results_file, hvg)