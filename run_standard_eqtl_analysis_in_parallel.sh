#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB




qtl_test_names_file="$1"
qtl_expression_file="$2"
qtl_covariate_file="$3"
qtl_genotype_file="$4"
qtl_sample_overlap_file="$5"
qtl_genotype_samples_to_rna_samples_mapping_file="$6"
qtl_output_root="$7"
job_number="$8"
num_jobs="$9"


source ~/.bash_profile
module load r/3.6.3


date
echo $qtl_output_root

Rscript run_standard_eqtl_analysis_in_parallel.R $qtl_test_names_file $qtl_expression_file $qtl_covariate_file $qtl_genotype_file $qtl_sample_overlap_file $qtl_genotype_samples_to_rna_samples_mapping_file $qtl_output_root $job_number $num_jobs
date