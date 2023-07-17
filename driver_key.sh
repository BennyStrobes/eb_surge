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



