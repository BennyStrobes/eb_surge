


surge_input_data_dir="$1"
pseudocell_technical_cov_file="$2"
pseudocell_pseudobulk_adata_file="$3"
pseudocell_to_donor_id_mapping_file="$4"
standard_eqtl_results_file="$5"
plink_eb_genotype_file_stem="$6"
eb_hvg_file="$7"
standard_eqtl_input_data_dir="$8"

source ~/.bash_profile

# Create file stem where all output files are written to
global_file_stem=${surge_input_data_dir}"eb_surge_input_standard_eqtl_"

################################
# Iterate over parameters deciding which tests to include in analysis
gene_set_arr=( "all_genes" "hv_genes")
num_genes_arr=( "500" "1000" "1500" "2000")
################################

# Loop through covariate methods
for gene_set in "${gene_set_arr[@]}"; do
for num_genes in "${num_genes_arr[@]}"; do

	echo $gene_set"_"$num_genes

	# Create new file stem for these parameters
	file_stem=${global_file_stem}${num_genes}"_"${gene_set}"_"


	# Extract names of tests (variant-gene pairs) to run eqtl factorization on
	test_names_file=${file_stem}"test_names.txt"
	genotype_input_file=${standard_eqtl_input_data_dir}"eqtl_genotype.traw"
	python extract_tests_for_surge.py $test_names_file $standard_eqtl_results_file $eb_hvg_file $gene_set $num_genes $genotype_input_file

	# Extract gene expression
	expression_input_file=${standard_eqtl_input_data_dir}"pseudobulk_expression_standardized.txt"
	surge_expression_file=${file_stem}"expression.npy"
	python extract_gene_expression_for_surge.py $test_names_file $expression_input_file $surge_expression_file

	# Extract genotype for SURGE
	genotype_input_file=${standard_eqtl_input_data_dir}"eqtl_genotype.traw"
	genotype_individual_to_sample_mapping_file=${standard_eqtl_input_data_dir}"eqtl_genotype_individual_to_sample_mapping.txt"
	surge_genotype_file=${file_stem}"genotype.npy"
	python extract_genotype_for_surge.py $test_names_file $genotype_input_file $genotype_individual_to_sample_mapping_file $surge_genotype_file

	# Extract covariates for SURGE
	covariate_input_file=${standard_eqtl_input_data_dir}"pseudocell_sample_eqtl_covariates.txt"
	surge_covariate_file=${file_stem}"covariates.txt"
	python extract_covariates_for_surge.py $test_names_file $covariate_input_file $surge_covariate_file

	# Extract sample repeat for SURGE
	cp ${standard_eqtl_input_data_dir}"eqtl_sample_repeat_vector.txt" ${file_stem}"sample_repeat.txt"


	# Extract sample covariates for surge analysis
	cp ${standard_eqtl_input_data_dir}"pseudobulk_sample_covariates.txt" ${file_stem}"sample_information.txt"

done
done
