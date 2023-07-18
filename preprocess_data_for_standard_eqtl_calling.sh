#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB


standard_eqtl_input_data_dir="$1"
pseudocell_pseudobulk_adata_file="$2"
plink_eb_genotype_file_stem="$3"
pseudocell_technical_cov_file="$4"
pseudocell_to_donor_id_mapping_file="$5"
n_pcs="$6"
old_standard_eqtl_results_file="$7"



source ~/.bash_profile



##########################
# Extract genotype data
##########################
eqtl_genotype_data_stem=${standard_eqtl_input_data_dir}"eqtl_genotype"
plink --bfile $plink_eb_genotype_file_stem --recode A-transpose --out $eqtl_genotype_data_stem


##########################
# Extract eqtl variant-gene pairs (tests) to test
# Limit to variants we have genotype data for
##########################
qtl_test_names_file=${standard_eqtl_input_data_dir}"standard_eqtl_test_names.txt"
python extract_variant_gene_pairs_to_test_for_standard_eqtl_analysis.py $old_standard_eqtl_results_file $qtl_test_names_file $eqtl_genotype_data_stem".traw"


##########################
# Extract expression data
###########################
python convert_pseudocell_expression_data_from_adata_to_tab_delimited_file.py $pseudocell_pseudobulk_adata_file $standard_eqtl_input_data_dir $qtl_test_names_file


##########################
# Extract covariate data
##########################
psuedocell_covariate_file=${standard_eqtl_input_data_dir}"pseudobulk_sample_covariates.txt"
qtl_covariate_file=${standard_eqtl_input_data_dir}"pseudocell_sample_eqtl_covariates.txt"
qtl_sample_repeat_file=${standard_eqtl_input_data_dir}"eqtl_sample_repeat_vector.txt"
pseudocell_expression_pc_file=${standard_eqtl_input_data_dir}"pseudobulk_expression_pcs_pca_scores.txt"
python extract_qtl_covariates.py $psuedocell_covariate_file $pseudocell_expression_pc_file $pseudocell_technical_cov_file $n_pcs $qtl_covariate_file $qtl_sample_repeat_file

##########################
# Create vector mapping from individual level genotype to sample level genotype
##########################
genotype_individual_to_sample_mapping_file=${standard_eqtl_input_data_dir}"eqtl_genotype_individual_to_sample_mapping.txt"
python create_mapping_from_genotype_individual_order_to_rna_sample_individual_order.py ${eqtl_genotype_data_stem}".traw" $qtl_covariate_file $genotype_individual_to_sample_mapping_file




