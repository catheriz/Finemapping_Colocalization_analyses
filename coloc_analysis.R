# Load required libraries
library(coloc)
library(stringr)
library(data.table)
set.seed(8)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
variant <- args[1]
summ_file_1 <- args[2]  # User-specified GWAS summary statistics file
summ_file_2 <- args[3]  # User-specified eQTL or other summary statistics file
ld_file_1 <- args[4]     # User-specified GWAS LD matrix file
ld_file_2 <- args[5]     # User-specified eQTL LD matrix file
lower <- as.numeric(args[6])  # User-specified lower bound
upper <- as.numeric(args[7])  # User-specified upper bound
chromosome <- as.numeric(args[8]) # User-specified chromosome
coloc_prob <- as.numeric(args[9]) # User-specified colocalization probability

# Read and process GWAS summary statistics
gwas_data <- fread(summ_file_1)
gwas_data$position <- as.numeric(str_split_fixed(gwas_data$SNP, ':', 4)[,2])
gwas_data$chr <- as.numeric(gsub('chr', '', str_split_fixed(gwas_data$SNP, ':', 4)[,1]))

gwas_data <- gwas_data[gwas_data$position >= lower & gwas_data$position <= upper & gwas_data$chr == chromosome, ]
gwas_data$varbeta <- (gwas_data$se) ** 2
gwas_data$z <- gwas_data$b / gwas_data$se
gwas_data <- gwas_data[, .(snp = SNP, chromosome = chr, position, allele1 = A1, allele2 = A2, beta = b, se, z, N, varbeta)]
setorder(gwas_data, chromosome, position)
gwas_data <- as.list(gwas_data)

# Load GWAS LD matrix
ldfile_1 <- as.matrix(fread(ld_file_1))
colnames(ldfile_1) <- gwas_data$snp
rownames(ldfile_1) <- gwas_data$snp
gwas_data$LD <- ldfile_1
gwas_data$type <- 'cc'

# Run SuSiE analysis
gwas_S <- runsusie(gwas_data, n = unique(gwas_data$N), L = 5, max_iter = 2500, verbose = TRUE, coverage = 0.95, min_abs_corr = sqrt(0.4), refine = TRUE)
summary(gwas_S)

# Extract credible sets
L <- length(gwas_S$sets$cs)
for (i in 1:L) {
  name_list <- names(gwas_S$sets$cs)[i]
  name_list_num <- as.numeric(gsub('L', '', name_list))
  q <- gwas_data$snp[gwas_data$snp %in% gsub("^L\\d+\\.", "", names(unlist(gwas_S$sets$cs[name_list])))]
  g <- as.data.frame(t(gwas_S$alpha[name_list_num, colnames(gwas_S$alpha) %in% q]))
  colnames(g) <- paste0(name_list, '_pip')
  g$name <- colnames(gwas_S$alpha)[colnames(gwas_S$alpha) %in% q]
  g$coverage <- gwas_S$sets$coverage[i]
  g <- g[, c(3,1,2)]
  fwrite(g, paste0(SNP, '_', name_list, '_credible_set_gwas.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Read and process eQTL summary statistics
eqtl_data <- fread(summ_file_2)
eqtl_data$position <- as.numeric(str_split_fixed(eqtl_data$SNP, ':', 4)[,2])
eqtl_data$chr <- as.numeric(gsub('chr', '', str_split_fixed(eqtl_data$SNP, ':', 4)[,1]))

eqtl_data <- eqtl_data[eqtl_data$position >= lower & eqtl_data$position <= upper & eqtl_data$chr == chromosome, ]
eqtl_data$varbeta <- (eqtl_data$se) ** 2
eqtl_data$z <- eqtl_data$b / eqtl_data$se
eqtl_data <- eqtl_data[, .(snp = SNP, chromosome = chr, position, allele1 = A1, allele2 = A2, beta = b, se, z, N, varbeta)]
setorder(eqtl_data, chromosome, position)
eqtl_data <- as.list(eqtl_data)

# Load eQTL LD matrix
ldfile_2 <- as.matrix(fread(ld_file_2))
colnames(ldfile_2) <- eqtl_data$snp
rownames(ldfile_2) <- eqtl_data$snp
eqtl_data$LD <- ldfile_2
eqtl_data$type <- 'cc'

# Run SuSiE analysis
eqtl_S <- runsusie(eqtl_data, n = unique(eqtl_data$N), L = 5, max_iter = 2500, verbose = TRUE, coverage = 0.95, min_abs_corr = sqrt(0.4), refine = TRUE)
summary(eqtl_S)

# Extract credible sets
L <- length(eqtl_S$sets$cs)
for (i in 1:L) {
  name_list <- names(eqtl_S$sets$cs)[i]
  name_list_num <- as.numeric(gsub('L', '', name_list))
  q <- gwas_data$snp[gwas_data$snp %in% gsub("^L\\d+\\.", "", names(unlist(gwas_S$sets$cs[name_list])))]
  g <- as.data.frame(t(eqtl_S$alpha[name_list_num, colnames(eqtl_S$alpha) %in% q]))
  colnames(g) <- paste0(name_list, '_pip')
  g$name <- colnames(eqtl_S$alpha)[colnames(eqtl_S$alpha) %in% q]
  g$coverage <- eqtl_S$sets$coverage[i]
  g <- g[, c(3,1,2)]
  fwrite(g, paste0(SNP, '_', name_list, '_credible_set_eqtl.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE)
}

# Perform colocalization analysis
susie.res <- coloc.susie(gwas_S, eqtl_S)
susie.res_summary <- as.data.frame(susie.res$summary)
write.table(susie.res_summary, paste0(variant, '_coloc_credible_set_summary.txt'), quote = FALSE, row.names = TRUE, col.names = TRUE)

# Extract credible sets with high posterior probability of colocalization and output list of SNPs within the credible sets
rows_with_high_PP_H4 <- which(susie.res_summary$PP.H4.abf > coloc_prob)
susie.res_result <- as.data.frame(susie.res$results)

if (length(rows_with_high_PP_H4) > 0) {
  susie.res_result <- lapply(rows_with_high_PP_H4, function(row) {
    list(
      snp = susie.res_result$snp, # Extract the first column (SNP ID)
      SNP.PP.H4.abf = susie.res_result[, row + 1] # Extract the (row + 1)th column
    )
  })
}
susie.res_result <- as.data.frame(susie.res_result)
susie.res_result <- susie.res_result[order(susie.res_result$SNP.PP.H4.abf, decreasing = TRUE),]
write.table(susie.res_result, paste0(variant, '_coloc_all_SNPs.txt'), quote = FALSE, row.names = TRUE, col.names = TRUE)

