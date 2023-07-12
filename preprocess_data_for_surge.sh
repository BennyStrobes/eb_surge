


surge_input_data_dir="$1"
pseudocell_technical_cov_file="$2"
pseudocell_expression_pc_file="$3"
pseudocell_pseudobulk_adata_file="$4"
pseudocell_to_donor_id_mapping_file="$5"
standard_eqtl_results_file="$6"
plink_eb_genotype_file_stem="$7"
eb_hvg_file="$8"

source ~/.bash_profile

# Create file stem where all output files are written to
file_stem=${surge_input_data_dir}"eb_surge_standard_eqtl_input_"


# Extract names of tests (variant-gene pairs) to run eqtl factorization on
test_names_file=${file_stem}"test_names.txt"
python extract_tests_for_surge.py $test_names_file $standard_eqtl_results_file $eb_hvg_file