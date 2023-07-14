import numpy as np 
import os
import sys
import pdb






#######################
# Command line args
#######################
genotype_data_file = sys.argv[1]  # Just used to get ordered genotype names
qtl_covariate_file = sys.argv[2]  # Used to get ordered rna sample names
genotype_individual_to_sample_mapping_file = sys.argv[3]  # Output file

# Extract ordered array rna samples individuals
rna_samples_indis = []
f = open(qtl_covariate_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count = head_count + 1
		continue
	sample_name = data[0]
	individual_name = sample_name.split(':')[0]
	rna_samples_indis.append(individual_name)
f.close()
rna_samples_indis = np.asarray(rna_samples_indis)

# Extract genotype individual anmes
geno_indis = []
f = open(genotype_data_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	raw = data[6:]
	for raw_name in raw:
		geno_indis.append(raw_name.split('_')[0])
	break
f.close()
geno_indis = np.asarray(geno_indis)

# Create mapping from indi name to genotype position
indi_name_to_genotyped_position = {}
for ii, val in enumerate(geno_indis):
	if val in indi_name_to_genotyped_position:
		print('assumption error')
		pdb.set_trace()
	indi_name_to_genotyped_position[val] = ii + 1


# Print to mapping file
tmp = []
t = open(genotype_individual_to_sample_mapping_file,'w')
for rna_sample_indi in rna_samples_indis:
	t.write(str(indi_name_to_genotyped_position[rna_sample_indi]) + '\n')
	tmp.append(indi_name_to_genotyped_position[rna_sample_indi])
t.close()

# Quick error checking
tmp = np.asarray(tmp) - 1 # -1 is for python. this is built for R
if np.array_equal(geno_indis[tmp], rna_samples_indis) == False:
	print('assumption erorro')
	pdb.set_trace()

