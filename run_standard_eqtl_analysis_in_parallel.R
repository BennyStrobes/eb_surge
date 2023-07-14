args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)
library(reshape2)
library(hash)




create_gene_name_to_expression_vector_hash <- function(qtl_expression_file, unique_genes) {
	# Initialize mapping from gene id to expression vec
	gene_id_to_expression_vector = hash()

	# Stream files
	stop = FALSE
	count = 0

	f = file(qtl_expression_file, "r")


	line = readLines(f, n=1)

	while(!stop) {
		count = count + 1

		line_info <- strsplit(line,'\t')[[1]]
		gene_id <- line_info[1]
		if (gene_id %in% unique_genes) {
			expression_vec = as.numeric(line_info[2:length(line_info)])
			gene_id_to_expression_vector[[gene_id]] = expression_vec
		}

		line = readLines(f, n=1)
		if (length(line) == 0) {
			stop=TRUE
			close(f)
		}
	}

	return(gene_id_to_expression_vector)
}

create_snp_to_genotype_vector_hash <- function(qtl_genotype_file, unique_snps) {
	# Initialize mapping from snp id to genotype vec
	snp_id_to_genotype_vector = hash()

	# Stream files
	stop = FALSE
	count = 0

	f = file(qtl_genotype_file, "r")


	line = readLines(f, n=1)

	while(!stop) {
		count = count + 1

		line_info <- strsplit(line,'\t')[[1]]
		snp_id <- line_info[2]
		if (snp_id %in% unique_snps) {
			genotype_vec = as.numeric(line_info[7:length(line_info)])
			snp_id_to_genotype_vector[[snp_id]] = genotype_vec
		}

		line = readLines(f, n=1)
		if (length(line) == 0) {
			stop=TRUE
			close(f)
		}
	}

	return(snp_id_to_genotype_vector)
}




run_eqtl_lmm <- function(expr, geno, covariates, groups) {
	fit_full <- lmer(expr ~ geno + covariates + (1|groups), REML=FALSE)
	fit_null <- lmer(expr ~ covariates + (1|groups), REML=FALSE)
	lrt <- anova(fit_null,fit_full)
	coefs <- data.frame(coef(summary(fit_full)))
	beta <- coefs[2,1]
	std_err <- coefs[2,2]
	pvalue <- lrt[[8]][2]
	return(list(eqtl_pvalue=pvalue, eqtl_beta=beta, eqtl_std_err=std_err))
}






#####################
# Command line args
#####################
qtl_test_names_file = args[1]
qtl_expression_file = args[2]
qtl_covariate_file = args[3]
qtl_genotype_file = args[4]
qtl_sample_overlap_file = args[5]
qtl_genotype_samples_to_rna_samples_mapping_file = args[6]
qtl_output_root = args[7]
job_number = as.numeric(args[8])
num_jobs = as.numeric(args[9])



#######################
# Load in test names
test_names = read.table(qtl_test_names_file, header=FALSE)
# Get number of tests
total_lines = length(test_names$V1)
print(paste0(total_lines, " total tests"))

#######################
# Determine number of lines each parrallelized job will complete
lines_per_job = ceiling(total_lines/num_jobs)
print(paste0(lines_per_job, " tests on this parallel job"))
start_num = job_number*lines_per_job
end_num = (job_number + 1)*lines_per_job
if (end_num > total_lines) {
	end_num = total_lines
}
# Subset test names (for this parallel run)
test_names_subset = test_names[(start_num+1):end_num,]

#######################
# Load in data
# Covariates
covariates <- read.table(qtl_covariate_file, header=TRUE)
covariates = as.matrix(covariates[,2:dim(covariates)[2]])
# Sample repeat vector
groups <- read.table(qtl_sample_overlap_file, header=FALSE)$V1
# Genotype individual to rna samples mapping
genotype_indi_to_rna_samples <- read.table(qtl_genotype_samples_to_rna_samples_mapping_file, header=FALSE)$V1

#######################
# Extract genotype data for all snps in test_names_subset
print("Loading in genotype data")
ordered_snps = as.character(test_names_subset$V2)
unique_snps = unique(ordered_snps)
snp_id_to_genotype_vector = create_snp_to_genotype_vector_hash(qtl_genotype_file, unique_snps)

#######################
# Extract expression data for all gene in test_names_subset
print("Loading in Expression data")
ordered_genes = as.character(test_names_subset$V1)
unique_genes = unique(ordered_genes)
gene_name_to_expression_vector = create_gene_name_to_expression_vector_hash(qtl_expression_file, unique_genes)


#######################
# Run eQTL analysis
print('Loading complete.')
print('Running eQTL analysis')

# Open output file handle
output_file <- paste0(qtl_output_root, "results.txt")
print(output_file)
sink(output_file)

# Get number of tests
n_tests = dim(test_names_subset)[1]

# Loop through tests
for (test_iter in 1:n_tests) {
	gene_id <- as.character(test_names_subset$V1[test_iter])
	rs_id <- as.character(test_names_subset$V2[test_iter])

	# Load in expression and genotype data
	expr = gene_name_to_expression_vector[[gene_id]]
	geno_ind = snp_id_to_genotype_vector[[rs_id]]
	geno_samp = geno_ind[genotype_indi_to_rna_samples]
	norm_geno = (geno_samp - mean(geno_samp))/(sd(geno_samp))

	new_line <- paste0(rs_id, "\t", gene_id)
	tryCatch(
		{
			lmm_results = run_eqtl_lmm(expr, norm_geno, covariates, groups)
			new_line <- paste0(new_line, "\t", lmm_results$eqtl_beta, "\t", lmm_results$eqtl_std_err, "\t", lmm_results$eqtl_pvalue)
			cat(paste0(new_line, "\n"))
        },
        error = function(e) {
        	new_line <- paste0(new_line, "\t", 0.0 ,"\t", 1.0, "\t", 1.0)
			cat(paste0(new_line, "\n"))
        }
    )


}
sink()



