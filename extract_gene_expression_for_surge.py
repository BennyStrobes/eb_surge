import numpy as np 
import os
import sys
import pdb


def extract_genes_for_surge(test_names_file):
	f = open(test_names_file)
	dicti = {}
	arr = []

	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[1]
		arr.append(gene_name)
		dicti[gene_name] = 1

	return np.asarray(arr), dicti




def create_mapping_from_gene_name_to_expression_array(expression_input_file, genes_dicti):
	# Initialize output mapping
	mapping = {}

	head_count = 0
	f = open(expression_input_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_name = data[0]
		if gene_name in genes_dicti:
			if gene_name in mapping:
				print('assumption erorro')
				pdb.set_trace()
			mapping[gene_name] = np.asarray(data[1:]).astype(float)
	f.close()
	return mapping








######################
# Command line args
######################
test_names_file = sys.argv[1]
expression_input_file = sys.argv[2]
surge_expression_file = sys.argv[3]

# Extract genes
ordered_genes, genes_dicti = extract_genes_for_surge(test_names_file)

# Quick error check
if len(ordered_genes) != len(genes_dicti):
	print('asssumtpion oeroro')
	pdb.set_trace()


# Create mapping from gene name to expression array
gene_name_to_expression = create_mapping_from_gene_name_to_expression_array(expression_input_file, genes_dicti)

# Create matrix of gene expression
expr_mat = []
for gene_name in ordered_genes:
	expr_mat.append(gene_name_to_expression[gene_name])
expr_mat = np.asarray(expr_mat)


# Save to output
np.save(surge_expression_file, expr_mat)


