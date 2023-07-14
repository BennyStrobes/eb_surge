import numpy as np 
import os
import sys
import pdb
import h5py
import scanpy as sc
from anndata import AnnData
import pandas as pd




def standardize_gene_expression_and_cap_outliers(un_normalized_expression, capper=10.0):
	temp_expr = (un_normalized_expression - np.mean(un_normalized_expression))/np.std(un_normalized_expression)
	temp_expr[temp_expr > capper] = capper
	temp_expr[temp_expr < (-1.0*capper)] = (-1.0*capper)
	temp_expr = temp_expr - np.mean(temp_expr)
	return temp_expr

def extract_dictionary_list_of_tested_genes(qtl_test_names_file):
	f = open(qtl_test_names_file)
	dicti = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		dicti[data[0]] = 1
	f.close()
	return dicti

##############################
# Command line args
##############################
pseudocell_pseudobulk_adata_file = sys.argv[1]
standard_eqtl_input_data_dir = sys.argv[2]
qtl_test_names_file = sys.argv[3]

#  extract_dictionary_list_of_tested_genes
tested_genes = extract_dictionary_list_of_tested_genes(qtl_test_names_file)


# Load in adata file
adata = sc.read_h5ad(pseudocell_pseudobulk_adata_file)

# Extract ordered sample names
pseudocell_sample_names = np.asarray(adata.obs.index)

# Extract sample covariate mat
pseudocell_sample_covs = np.asarray(adata.obs)
pseudocell_sample_cov_names = np.asarray(adata.obs.keys())

if pseudocell_sample_covs.shape[0] != len(pseudocell_sample_names):
	print('assumption eroror')
	pdb.set_trace()

# Extract number of pseudocell samples
n_samples = len(pseudocell_sample_names)

# Print pseudocell sample covariates to output file
t = open(standard_eqtl_input_data_dir+ 'pseudobulk_sample_covariates.txt','w')
# print header
t.write('pseudocell_name\t' + '\t'.join(pseudocell_sample_cov_names) + '\n')
for sample_iter in range(n_samples):
	t.write(pseudocell_sample_names[sample_iter] + '\t' + '\t'.join(pseudocell_sample_covs[sample_iter, :].astype(str)) + '\n')
t.close()


# Extract ordered gene names
ordered_gene_symbol_names =  np.asarray(adata.var.index)
ordered_ensamble_ids = np.asarray(adata.var['gene_ids'])

# Quick error checking
if len(ordered_ensamble_ids) != len(ordered_gene_symbol_names):
	print('assumption eroror')
	pdb.set_trace()

# Extract number of genes
n_genes = adata.X.shape[1]
# Quick error checking
if n_genes != len(ordered_gene_symbol_names):
	print('assumption error')
	pdb.set_trace()



# Print expression matrix to tab delimited file
tab_delimited_expression_matrix_file = standard_eqtl_input_data_dir + 'pseudobulk_expression_standardized.txt'
# Print header
t = open(tab_delimited_expression_matrix_file,'w')
t.write('gene_id\t' + '\t'.join(pseudocell_sample_names) + '\n')
# Loop through genes
for gene_iter in range(n_genes):
	# Extract relevent fields
	gene_name = ordered_gene_symbol_names[gene_iter]
	un_normalized_expression = adata.X[:, gene_iter]

	gene_name = ordered_gene_symbol_names[gene_iter]
	if gene_name not in tested_genes:
		continue 

	# Standardize gene expression
	normalized_expression = standardize_gene_expression_and_cap_outliers(un_normalized_expression)

	# Print to output file
	t.write(gene_name + '\t' + '\t'.join(normalized_expression.astype(str)) + '\n')
t.close()


# Print Gene id file
gene_id_file = standard_eqtl_input_data_dir + 'pseudobulk_expression_gene_ids.txt'
t = open(gene_id_file,'w')
t.write('gene_name\tensamble_id\n')
for gene_iter in range(n_genes):
	gene_name = ordered_gene_symbol_names[gene_iter]
	if gene_name not in tested_genes:
		continue 
	t.write(ordered_gene_symbol_names[gene_iter] + '\t' + ordered_ensamble_ids[gene_iter] + '\n')
t.close()






