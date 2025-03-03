library(data.table)
library(dplyr)
library(locuszoomr)
library(AnnotationHub)
library(stringr)
library(ggplot2)
library(ggrepel)
library(cowplot)
ah <- AnnotationHub()
ensDb_v113 <- ah[["AH119325"]]
ensDb_v111 <- ah[["AH116291"]]

args <- commandArgs(trailingOnly = TRUE)
disease = args[1]
variants <- unlist(strsplit(args[2], ","))  # Accept multiple variants, comma-separated
gene <- ifelse(args[3] == "NULL", NULL, unlist(strsplit(args[3], ";")))
index_snp_1 = args[4]
index_snp_2 = args[5]
summ_file_1 <- args[6]  # User-specified GWAS summary statistics file
summ_file_2 <- args[7]  # User-specified eQTL or other summary statistics file
ld_file_1 <- args[8]     # User-specified GWAS LD matrix file
ld_file_2 <- args[9]     # User-specified eQTL LD matrix file

# Function to read and preprocess summary statistics
read_summ_data <- function(file) {
  summ <- fread(file,head=T)
  summ[, SNP := gsub("chr", "", SNP)]
  summ[, `:=` (Z = b / se,
               position = as.numeric(str_split_fixed(SNP, ":", 4)[,2]),
               chromosome = as.numeric(gsub("chr", "", str_split_fixed(SNP, ":", 4)[,1])))]
  summ <- summ[, .(rsid = SNP, chromosome, position, Z, p)]
  return(summ)
}

# Function to read and preprocess LD matrix
read_ld_data <- function(file) {
  ld <- fread(file, select = c("SNP_B", "R2"))
  ld[, SNP_B := gsub("chr", "", SNP_B)]
  setnames(ld, c("SNP_B", "R2"), c("rsid", "r2"))
  return(ld)
}


# Read summary statistics and LD data
summ_1 <- read_summ_data(summ_file_1)
ld_1 <- read_ld_data(ld_file_1)
indat_1 <- merge(summ_1, ld_1, by = "rsid", all.y = TRUE)

summ_2 <- read_summ_data(summ_file_2)
ld_2 <- read_ld_data(ld_file_2)
indat_2 <- merge(summ_2, ld_2, by = "rsid", all.y = TRUE)

# Construct labels_graph with all variants
variants_clean <- gsub("chr", "", variants)  # Clean variant names
labels_graph <- c(index_snp_1, index_snp_2, variants_clean)

# Create a variant string
variant_str <- gsub("[:,]", "_", paste(variants, collapse = "_"))

# Output PDF file
pdf(paste0(disease, "_", variant_str, "_coloc_SNP.pdf"), width = 9, height = 6)

multi_layout(
  nrow = 1,
  ncol = 2,
  heights = c(9, 3),

  plots = {
    # First locus plot
    {
      locus_plot(
        loc_1, 
        use_layout = FALSE, 
        legend_pos = NULL, 
        labels = unique(labels_graph), 
        LD_scheme = adjustcolor(c("grey", "royalblue", "cyan2", "green3", "orange", "red", "purple"), alpha.f = 0.65),
        maxrows = 2, 
        cex.text = 1.2, 
        cex = 1.2, 
        cex.axis = 1.2, 
        cex.lab = 1.2,
        label_x = c(-4, 3),
        label_y = c(3, 3),
        point.args = list( col = loc_1$data$col )
      )
    }
    # Second locus plot
    {
      locus_plot(
        loc_2, 
        use_layout = FALSE, 
        legend_pos = NULL, 
        labels = unique(labels_graph), 
        LD_scheme = adjustcolor(c("grey", "royalblue", "cyan2", "green3", "orange", "red", "purple"), alpha.f = 0.65),

        maxrows = 2, 
        cex.text = 1.2, 
        cex = 1.2, 
        cex.axis = 1.2, 
        cex.lab = 1.2,
        label_x = c(-4, 3),
        label_y = c(3, 3),
        point.args = list( col = loc_2$data$col )
      )
    }
  }
)
dev.off()
