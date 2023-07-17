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


def extract_array_of_gene_and_top_variants_according_to_standard_eqtl_results(standard_eqtl_results_file, hvg, gene_set, valid_snps):
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
		gene_id = data[1]
		variant_id = data[0]
		pvalue = float(data[4])

		if variant_id not in valid_snps:
			# Quick error checking
			if pvalue != 1.0:
				print('assumption eororor')
				pdb.set_trace()
			continue

		# Ignore genes not in HVG
		if gene_id not in hvg and gene_set == 'hv_genes':
			continue

		# Add to dictionary
		if gene_id not in gene_to_variant_and_pvalue:
			gene_to_variant_and_pvalue[gene_id] = (variant_id, pvalue, 1)
		else:
			new_count = gene_to_variant_and_pvalue[gene_id][2] + 1
			old_pvalue = gene_to_variant_and_pvalue[gene_id][1]
			old_variant_id = gene_to_variant_and_pvalue[gene_id][0]
			if pvalue < old_pvalue:
				gene_to_variant_and_pvalue[gene_id] = (variant_id, pvalue, new_count)
			else:
				gene_to_variant_and_pvalue[gene_id] = (old_variant_id, old_pvalue, new_count)
	f.close()

	tuple_arr = []
	for gene_id in [*gene_to_variant_and_pvalue]:
		tuple_arr.append((gene_id, gene_to_variant_and_pvalue[gene_id][0], gene_to_variant_and_pvalue[gene_id][1]*gene_to_variant_and_pvalue[gene_id][2]))

	tuple_arr_sorted = sorted(tuple_arr, key=lambda tup: tup[2])

	return tuple_arr_sorted





def extract_valid_snps(genotype_input_file):
	f = open(genotype_input_file)
	valid_snps = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		valid_snps[data[1]] = 1
	f.close()
	return valid_snps



######################
# Command line args
######################
test_names_file = sys.argv[1]  # output file
standard_eqtl_results_file = sys.argv[2]  # File containing standard eqtl results for all variant-gene pairs
hvg_file = sys.argv[3]  # Highly variable genes file
gene_set = sys.argv[4]  # Whether or not to include HVG
num_genes = int(sys.argv[5])  # Threshold for significance of variant-gene pairs
genotype_input_file =sys.argv[6]

# Extract valid snps (that we have genotype for)
valid_snps = extract_valid_snps(genotype_input_file)

# First extract highly variable genes
hvg = extract_hvg(hvg_file)

# Extract array of gene and the gene's top variant and the corresponding pvalue
all_variant_gene_pairs = extract_array_of_gene_and_top_variants_according_to_standard_eqtl_results(standard_eqtl_results_file, hvg, gene_set, valid_snps)

# Open output file handle
t = open(test_names_file,'w')
t.write('variant_id\tgene_name\n')

# Don't allow variant to be used twice
used_variants = {}

# Loop through genes
for tuple_id in all_variant_gene_pairs[:num_genes]:
	variant_name = tuple_id[1]
	gene_name = tuple_id[0]

	# Don't allow for repeat variants
	if variant_name in used_variants:
		continue
	used_variants[variant_name] = 1

	# print to output
	t.write(variant_name + '\t' + gene_name + '\n')

# Close file handle
t.close()

