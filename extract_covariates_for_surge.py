import numpy as np 
import os
import sys
import pdb







test_names_file = sys.argv[1]
covariate_input_file = sys.argv[2]
surge_covariate_file = sys.argv[3]


# Load in data
cov = np.loadtxt(covariate_input_file,dtype=str,delimiter='\t')

# Remove column labels and row labels
cov_values = cov[1:,1:].astype(float)

# Standardize cov
for cov_iter in range(cov_values.shape[1]):
	# Standardize covariates
	cov_values[:, cov_iter] = (cov_values[:, cov_iter] - np.mean(cov_values[:, cov_iter]))/np.std(cov_values[:, cov_iter])


# Add intercept 
cov_values_plus_intercept = np.hstack((np.ones((cov_values.shape[0],1)), cov_values))

# Save to output
np.savetxt(surge_covariate_file, cov_values_plus_intercept.astype(str), fmt="%s", delimiter='\t')