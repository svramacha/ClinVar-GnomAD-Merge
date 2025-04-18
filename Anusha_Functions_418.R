library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(GenomicRanges)
library(Rsamtools)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)


# ------------------------
# Extract Gnomad -- need to concatenate all files 
#-------------------------

extract_gnomad_gene <- function(vcf_path, gene_name, chrom, start, end) {
  region <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end, names = gene_name))
  tabix <- TabixFile(vcf_path)
  vcf <- readVcf(tabix, "hg38", param = region)
  
  df <- as.data.frame(rowRanges(vcf))
  df$FAFmax_faf95_max_joint <- unlist(info(vcf)$fafmax_faf95_max_joint)
  df$alt_values <- sapply(df$ALT, as.character)
  df$merge <- paste(df$start, df$REF, df$alt_values, sep = " ")
  df <- df %>% filter(FILTER == "PASS")
  
  return(df)
}

# ------------------------
# Extract Clinvar
# ------------------------

extract_clinvar_gene <- function(clinvar_path, gene_name) {
  vcf <- readVcf(clinvar_path, "hg38")
  geneinfo <- info(vcf)$GENEINFO
  gene_variants <- grepl(gene_name, geneinfo, ignore.case = TRUE)
  vcf_gene <- vcf[gene_variants, ]
  
  info_df <- as.data.frame(info(vcf_gene))
  gr_df <- as.data.frame(rowRanges(vcf_gene))
  df <- cbind(info_df, gr_df)
  df$alt_values <- sapply(df$ALT, as.character)
  df$merge <- paste(df$start, df$REF, df$alt_values, sep = " ")
  
  return(df)
}


# ------------------------
# Extract CADD
# ------------------------

extract_cadd_gene <- function(cadd_path, chrom, start, end) {
  region <- GRanges(seqnames = chrom, ranges = IRanges(start = start, end = end))
  tabix <- TabixFile(cadd_path)
  raw_data <- scanTabix(tabix, param = region)
  
  df <- read.table(text = unlist(raw_data), header = FALSE, sep = "\t")
  colnames(df) <- c("CHROM", "start", "REF", "ALT", "RawScore", "PHRED")
  df$merge <- paste(df$start, df$REF, df$ALT, sep = " ")
  
  return(df)
}
