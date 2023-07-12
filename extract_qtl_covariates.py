import numpy as np 
import os
import sys
import pdb




def create_numeric_covariate_matrix_from_categorical_cov(input_vec):
	unique_categories = np.unique(input_vec)
	n_categories = len(unique_categories) - 1
	n_samples = len(input_vec)

	return_mat = np.zeros((n_samples, n_categories))

	for category_iter in range(n_categories):
		return_mat[:, category_iter] = 1.0*(input_vec == unique_categories[category_iter])

	return return_mat







#######################
# command line args
#######################
psuedocell_covariate_file = sys.argv[1]
pseudocell_expression_pc_file = sys.argv[2]
pseudocell_technical_cov_file = sys.argv[3]
n_pcs = int(sys.argv[4])
qtl_covariate_file = sys.argv[5]  # Output file

# Extract pseudocell sample names according to gene expression (need to make sure covariate data matches this)
tmp_data = np.loadtxt(psuedocell_covariate_file,dtype=str,delimiter='\t')
pseudocell_sample_names = tmp_data[1:,0]


# Load in expression pc data
expr_pc_tmp = np.loadtxt(pseudocell_expression_pc_file, dtype=str,delimiter='\t')
# Error checking to make sure sample names agree
if np.array_equal(expr_pc_tmp[1:,0], pseudocell_sample_names) == False:
	print('asssumption error')
	pdb.set_trace()

expr_pc_mat = expr_pc_tmp[1:, 1:(n_pcs+1)]
expr_pc_col_names = expr_pc_tmp[0, 1:(n_pcs+1)]


# Create technical cov mat for first covarate
tech_cov_1_mat = create_numeric_covariate_matrix_from_categorical_cov(tmp_data[1:,6])  # Sex 
tech_cov_2_mat = create_numeric_covariate_matrix_from_categorical_cov(tmp_data[1:,7])  # lib.prep.batch

# Concatenate 
tech_cov_mat = np.hstack((tech_cov_1_mat, tech_cov_2_mat))
tech_cov_mat_column_names = []
tech_cov_mat_column_names.append('sex')
for itera in range(tech_cov_2_mat.shape[1]):
	tech_cov_mat_column_names.append('lib_prep_batch_' + str(itera))
tech_cov_mat_column_names = np.asarray(tech_cov_mat_column_names)

# Concatenate everything
qtl_cov_mat = np.hstack((expr_pc_mat, tech_cov_mat))
qtl_cov_names = np.hstack((expr_pc_col_names, tech_cov_mat_column_names))


# Print to output file
t = open(qtl_covariate_file,'w')
# Header
t.write('pseudocell_name\t' + '\t'.join(qtl_cov_names) + '\n')
# loop through sample names
n_samples = qtl_cov_mat.shape[0]
for sample_iter in range(n_samples):
	t.write(pseudocell_sample_names[sample_iter] + '\t' + '\t'.join(qtl_cov_mat[sample_iter, :].astype(str)) + '\n')
t.close()



