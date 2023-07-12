#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB


standard_eqtl_input_data_dir="$1"
pseudocell_pseudobulk_adata_file="$2"
plink_eb_genotype_file_stem="$3"
pseudocell_expression_pc_file="$4"
pseudocell_technical_cov_file="$5"
pseudocell_to_donor_id_mapping_file="$6"
n_pcs="$7"



source ~/.bash_profile


if false; then
python convert_pseudocell_expression_data_from_adata_to_tab_delimited_file.py $pseudocell_pseudobulk_adata_file $standard_eqtl_input_data_dir
fi

psuedocell_covariate_file=${standard_eqtl_input_data_dir}"pseudobulk_sample_covariates.txt"
qtl_covariate_file=${standard_eqtl_input_data_dir}"pseudocell_sample_eqtl_covariates.txt"
python extract_qtl_covariates.py $psuedocell_covariate_file $pseudocell_expression_pc_file $pseudocell_technical_cov_file $n_pcs $qtl_covariate_file
