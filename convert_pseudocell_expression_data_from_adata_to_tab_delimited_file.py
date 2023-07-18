import numpy as np 
import os
import sys
import pdb
import h5py
import scanpy as sc
from anndata import AnnData
import pandas as pd
from sklearn.decomposition import PCA
import rnaseqnorm



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

# Generate expression PC loadings and variance explained of those expression PCs
def generate_pca_scores_and_variance_explained(X, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file):
	# Load in data
	#X = np.loadtxt(filtered_standardized_sc_expression_file)

	# Run PCA (via SVD)
	#uuu, sss, vh = np.linalg.svd(np.transpose(X), full_matrices=False)
	#svd_loadings = np.transpose(vh)[:,:num_pcs]
	#ve = (np.square(sss)/np.sum(np.square(sss)))[:num_pcs]

	# Faster in sklearn
	_pca = PCA(n_components=num_pcs, svd_solver='arpack')
	svd_loadings = _pca.fit_transform(X)
	ve = _pca.explained_variance_ratio_

	# Save to output file
	np.savetxt(filtered_cells_pca_file, svd_loadings, fmt="%s", delimiter='\t')

	# Compute variance explained
	np.savetxt(filtered_cells_pca_ve_file, ve, fmt="%s", delimiter='\n')

def normalize_expression_and_generate_expression_pcs(raw_pseudobulk_expression, gene_names, sample_names, sample_level_normalization, gene_level_normalization, num_pcs, pseudobulk_expression_file, expression_pc_output_stem):
	# Initialize output normalized expression matrix
	normalized_expression = np.zeros(raw_pseudobulk_expression.shape)

	##################################
	# Perform sample level normalization
	##################################
	if sample_level_normalization == 'qn':
		print('quantile normalizing')
		df = pd.DataFrame(np.transpose(raw_pseudobulk_expression))
		temp_out = rnaseqnorm.normalize_quantiles(df)
		raw_pseudobulk_expression = np.transpose(np.asarray(temp_out))

	##################################
	# Perform gene level normalization
	##################################
	if gene_level_normalization == 'zscore':
		for gene_num in range(normalized_expression.shape[1]):
			temp_expr = (raw_pseudobulk_expression[:, gene_num] - np.mean(raw_pseudobulk_expression[:, gene_num]))/np.std(raw_pseudobulk_expression[:, gene_num])
			temp_expr[temp_expr > 10.0] = 10.0
			temp_expr[temp_expr < -10.0] = -10.0
			temp_expr = temp_expr - np.mean(temp_expr)
			normalized_expression[:, gene_num] = temp_expr
	elif gene_level_normalization == 'ign':
		# Code from GTEx v8
		# Project each gene onto a gaussian
		df = pd.DataFrame(np.transpose(raw_pseudobulk_expression))
		norm_df = rnaseqnorm.inverse_normal_transform(df)
		normalized_expression = np.transpose(np.asarray(norm_df))
	else:
		print(gene_level_normalization + ' gene level normalization method currently not implemented')
		pdb.set_trace()

	# Save normalized pseudobulk gene expression to output file
	t = open(tab_delimited_expression_matrix_file,'w')
	t.write('gene_id\t' + '\t'.join(sample_names) + '\n')
	# Loop through genes
	for gene_iter in range(len(gene_names)):
		t.write(gene_names[gene_iter] + '\t' + '\t'.join(normalized_expression[:, gene_iter].astype(str)) + '\n')
	t.close()	

	# Run PCA on pseudobulk data
	pca_file = expression_pc_output_stem + 'pca_scores.txt'
	pca_ve_file = expression_pc_output_stem + 'pca_pve.txt'
	generate_pca_scores_and_variance_explained(normalized_expression, num_pcs, pca_file, pca_ve_file)

	return


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


log_transformed_expression_mat = []
gene_names = []
# Loop through genes
for gene_iter in range(n_genes):
	# Extract relevent fields
	gene_name = ordered_gene_symbol_names[gene_iter]
	un_normalized_expression = adata.X[:, gene_iter]

	gene_name = ordered_gene_symbol_names[gene_iter]
	if gene_name not in tested_genes:
		continue 
	log_transformed_expression_mat.append(un_normalized_expression)
	gene_names.append(gene_name)
log_transformed_expression_mat = np.asarray(log_transformed_expression_mat)
gene_names = np.asarray(gene_names)


#####################
# Normalize expression and generate expression pcs
#####################	
# Options for sample level normalization are currently 'none'
sample_level_normalization = 'qn'
# Options for gene level normalization are 'zscore' and 'ign'
gene_level_normalization = 'zscore'
# number of pcs
num_pcs = 200
# output root
tab_delimited_expression_matrix_file = standard_eqtl_input_data_dir + 'pseudobulk_expression_standardized.txt'
expression_pc_output_stem = standard_eqtl_input_data_dir + 'pseudobulk_expression_pcs_'

normalize_expression_and_generate_expression_pcs(np.transpose(log_transformed_expression_mat), gene_names, pseudocell_sample_names, sample_level_normalization, gene_level_normalization, num_pcs, tab_delimited_expression_matrix_file, expression_pc_output_stem)




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






















############
#OLD
############
'''
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



'''




