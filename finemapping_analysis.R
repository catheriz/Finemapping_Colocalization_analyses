# Load required libraries
library(susieR)
library(stringr)
library(data.table)
set.seed(8)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
variant <- args[1]
summ_file <- args[2]  # User-specified summary statistics file
ld_file <- args[3]     # User-specified LD matrix file
lower <- as.numeric(args[4])  # User-specified lower bound
upper <- as.numeric(args[5])  # User-specified upper bound
chromosome <- as.numeric(args[6]) # User-specified chromosome

# Read input summary statistics
gwas_data <- fread(summ_file)
gwas_data$position <- as.numeric(str_split_fixed(gwas_data$SNP, ':', 4)[,2])


# Filter GWAS data
gwas_data <- gwas_data[gwas_data$position >= lower & gwas_data$position <= upper & gwas_data$chr == chromosome, ]
gwas_data$varbeta <- (gwas_data$se) ** 2
gwas_data$z <- gwas_data$b / gwas_data$se
gwas_data <- gwas_data[, .(snp = SNP, chromosome = chr, position, allele1 = A1, allele2 = A2, beta = b, se, z, N, varbeta)]
setorder(gwas_data, chromosome, position)

# Load LD matrix
ldfile <- as.matrix(fread(ld_file))
colnames(ldfile) <- gwas_data$snp
rownames(ldfile) <- gwas_data$snp
gwas_data$LD <- ldfile
gwas_data$type <- 'cc'

# Run SuSiE analysis
gwas_S <- runsusie(gwas_data, n = unique(gwas_data$N), L = 5, max_iter = 2500, verbose = TRUE, coverage = 0.95, min_abs_corr = sqrt(0.4), refine = TRUE)
summary(gwas_S)

# Extract credible sets
L <- length(gwas_S$sets$cs)
for (i in 1:L) {
  name_list <- names(gwas_S$sets$cs)[i]
  name_list_num <- as.numeric(gsub('L', '', name_list))
  q <- gwas_data$snp[as.numeric(rownames(gwas_data)) %in% unlist(gwas_S$sets$cs[name_list])]
  g <- as.data.frame(t(gwas_S$alpha[name_list_num, colnames(gwas_S$alpha) %in% q]))
  colnames(g) <- paste0(name_list, '_pip')
  g$name <- colnames(gwas_S$alpha)[colnames(gwas_S$alpha) %in% q]
  g$coverage <- gwas_S$sets$coverage[i]
  g <- g[, c(3,1,2)]
  fwrite(g, paste0(SNP, '_', name_list, '_credible_set_susieR_results.txt'), quote = FALSE, row.names = FALSE, col.names = TRUE)
}

