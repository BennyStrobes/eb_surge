import numpy as np 
import os
import sys
import pdb









######################
# Command line args
#######################
old_standard_eqtl_results_file = sys.argv[1]  # INput file
output_file = sys.argv[2] # output file



f = open(old_standard_eqtl_results_file)
t = open(output_file,'w')
head_count = 0
used_tests = {}
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	gene = data[0]
	snp_id = data[4]
	chrom_num = data[2]
	test_name = gene + ':' + snp_id
	if test_name in used_tests:
		print('assumption oeroror')
		pdb.set_trace()

	t.write(gene + '\t' + snp_id + '\t' + chrom_num + '\n')
	used_tests[test_name] = 1

f.close()
t.close()