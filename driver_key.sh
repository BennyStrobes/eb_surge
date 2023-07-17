##################################
# Input data
##################################

# File containing pseudocell_technical covariates
# Shared by Merlin Li
pseudocell_technical_cov_file="/data/abattle4/merlin/eb_project/surge/data_ben/cov.csv"

# File containing pseudocell expression pcs
# Contains 100 expression pcs
# Shared by Merlin Li
pseudocell_expression_pc_file="/data/abattle4/merlin/eb_project/surge/data_ben/expr_pc100.csv"

# File containing pseudocell pseudobulk adata file
# Shared by Merlin Li
pseudocell_pseudobulk_adata_file="/data/abattle4/merlin/eb_project/surge/data_ben/eb_prlog1pf_pseudobulk.h5ad"

# File containing mapping from pseudocell name to donor id
# Shared by Merlin Li
pseudocell_to_donor_id_mapping_file="/data/abattle4/merlin/eb_project/surge/data_ben/pseudocell_donor_map.csv"

# File contianing standard eQTL results
# Shared by Merlin Li
old_standard_eqtl_results_file="/data/abattle4/merlin/eb_project/surge/data_ben/association_summary_bonf.tsv"

# Plink file stem containing eb genotype data
# Shared by Merlin Li
plink_eb_genotype_file_stem="/data/abattle4/merlin/eb_project/surge/data_ben/genotype/yri_maf0.1_all.hg38"

# File containing list of 5000 highly variable genes identified in the EB data
# Shared by Josh Popp
eb_hvg_file="/scratch16/abattle4/bstrober/eb_surge/input_data/eb_hvgs.txt"



##################################
# Output data
##################################
# Output root
output_root="/scratch16/abattle4/bstrober/eb_surge/"

# Directory containing input data for standard eqtl calling
standard_eqtl_input_data_dir=${output_root}"standard_eqtl_input/"

# Directory containing standard eqtl calling results
standard_eqtl_results_dir=${output_root}"standard_eqtl_results/"

# Directory containing Processed input data for surge
surge_input_data_dir=${output_root}"surge_input/"

# Directory containing Surge results
surge_results_dir=${output_root}"surge_results/"

# Directory containing visualizations of Surge results
visualization_dir=${output_root}"visualization/"



######################################
# Analysis parameters
######################################
n_pcs="20"


############################
# Preprocess data for standard eQTL calling
if false; then
sh preprocess_data_for_standard_eqtl_calling.sh $standard_eqtl_input_data_dir $pseudocell_pseudobulk_adata_file $plink_eb_genotype_file_stem $pseudocell_expression_pc_file $pseudocell_technical_cov_file $pseudocell_to_donor_id_mapping_file $n_pcs $old_standard_eqtl_results_file
fi

############################
# Run standard eqtl calling
if false; then
sh run_standard_eqtl_calling.sh $standard_eqtl_input_data_dir $standard_eqtl_results_dir
fi


############################
# Preprocess data for SURGE
standard_eqtl_results_file=$standard_eqtl_results_dir"standard_eqtl_results_merged.txt"
if false; then
sh preprocess_data_for_surge.sh $surge_input_data_dir $pseudocell_technical_cov_file $pseudocell_expression_pc_file $pseudocell_pseudobulk_adata_file $pseudocell_to_donor_id_mapping_file $standard_eqtl_results_file $plink_eb_genotype_file_stem $eb_hvg_file $standard_eqtl_input_data_dir
fi



############################
# Run SURGE
# SURGE parameters
permutation_type="False"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="20"
data_filter="none"
delta_elbo_thresh="1e-2"
lambda_v="1"
seed="1"
variance_param="1e-3"
ard_variance_param="1e-3"
# SURGE input data
gene_set_arr=( "all_genes" "hv_genes")
num_genes_arr=( "500" "1000" "1500" "2000")
if false; then
for gene_set in "${gene_set_arr[@]}"; do
for num_genes in "${num_genes_arr[@]}"; do
	input_data_stem="eb_surge_input_standard_eqtl_"${num_genes}"_"${gene_set}
	test_names_file=${surge_input_data_dir}${input_data_stem}"_test_names.txt"
	expression_training_file=${surge_input_data_dir}${input_data_stem}"_expression.npy"
	genotype_training_file=${surge_input_data_dir}${input_data_stem}"_genotype.npy"
	covariate_file=${surge_input_data_dir}${input_data_stem}"_covariates.txt"
	sample_overlap_file=${surge_input_data_dir}${input_data_stem}"_sample_repeat.txt"


	output_stem=$surge_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_"
	sbatch run_surge.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh
done
done
fi




