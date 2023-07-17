import numpy as np 
import os
import sys
import pdb



def extract_variants_for_surge(test_names_file):
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
		variant_name = data[0]
		arr.append(variant_name)
		dicti[variant_name] = 1

	return np.asarray(arr), dicti



def create_mapping_from_variant_name_to_individual_genotype_array(genotype_input_file, variant_dicti):
	mapping = {}
	f = open(genotype_input_file)
	head_count = 0 
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[1]
		if variant_id in variant_dicti:
			mapping[variant_id] = np.asarray(data[6:]).astype(float)
	f.close()
	return mapping






########################
# Command line args
########################
test_names_file = sys.argv[1]
genotype_input_file = sys.argv[2]
genotype_individual_to_sample_mapping_file = sys.argv[3]
surge_genotype_file = sys.argv[4]


# Extract variants
ordered_variants, variant_dicti = extract_variants_for_surge(test_names_file)

# Genotype individual to sample mapping
genotype_individual_to_sample_mapping = np.loadtxt(genotype_individual_to_sample_mapping_file) - 1 #  Minus 1 cause in python
genotype_individual_to_sample_mapping = genotype_individual_to_sample_mapping.astype(int)

# Create mapping from variant name to individual-level genotype array
variant_name_to_individual_genotype = create_mapping_from_variant_name_to_individual_genotype_array(genotype_input_file, variant_dicti)


# Create matrix of genotype
variant_mat = []
for variant_name in ordered_variants:
	individual_level_genotype = variant_name_to_individual_genotype[variant_name]

	# Convert to sample level
	sample_level_genotype = individual_level_genotype[genotype_individual_to_sample_mapping]
	variant_mat.append(sample_level_genotype)
variant_mat = np.asarray(variant_mat)



# Save to output
np.save(surge_genotype_file, variant_mat)



