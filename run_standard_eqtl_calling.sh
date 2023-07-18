#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB





standard_eqtl_input_data_dir="$1"
standard_eqtl_results_dir="$2"



###################
# Parameters
num_jobs="100"

###################
# Input data
qtl_expression_file=${standard_eqtl_input_data_dir}"pseudobulk_expression_standardized.txt"
qtl_covariate_file=${standard_eqtl_input_data_dir}"pseudocell_sample_eqtl_covariates.txt"
qtl_genotype_file=${standard_eqtl_input_data_dir}"eqtl_genotype.traw"
qtl_sample_overlap_file=${standard_eqtl_input_data_dir}"eqtl_sample_repeat_vector.txt"
qtl_genotype_samples_to_rna_samples_mapping_file=${standard_eqtl_input_data_dir}"eqtl_genotype_individual_to_sample_mapping.txt"
qtl_test_names_file=${standard_eqtl_input_data_dir}"standard_eqtl_test_names.txt"


job_number="0"
qtl_output_root=$standard_eqtl_results_dir"standard_eqtl_results_"$job_number"_"$num_jobs"_"
if false; then
sbatch run_standard_eqtl_analysis_in_parallel.sh $qtl_test_names_file $qtl_expression_file $qtl_covariate_file $qtl_genotype_file $qtl_sample_overlap_file $qtl_genotype_samples_to_rna_samples_mapping_file $qtl_output_root $job_number $num_jobs
fi


if false; then
for job_number in $(seq 1 $(($num_jobs-1))); do 
	qtl_output_root=$standard_eqtl_results_dir"standard_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_standard_eqtl_analysis_in_parallel.sh $qtl_test_names_file $qtl_expression_file $qtl_covariate_file $qtl_genotype_file $qtl_sample_overlap_file $qtl_genotype_samples_to_rna_samples_mapping_file $qtl_output_root $job_number $num_jobs
done
fi

if false; then
python merge_parallelized_standard_eqtl_calls.py $standard_eqtl_results_dir"standard_eqtl_results_" $num_jobs
fi
