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
cov_values = cov[1:,1:]

# Save to output
np.savetxt(surge_covariate_file, cov_values, fmt="%s", delimiter='\t')