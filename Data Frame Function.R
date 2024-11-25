library(BiocManager)
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)
library(stringr)
library(tidyr)
library(Biostrings)


# Original code for CHEK2
rng_chek2 <- GRanges(seqnames = "chr22", ranges = IRanges(
  start = 28687820,
  end = 28742014,
  names = "CHEK2"
))

gnomad_vcf_chr22 <- "/Volumes/Seagate/Capstone/gnomad.joint.v4.1.sites.chr22.vcf.bgz"

tab_gnomad <- TabixFile(gnomad_vcf_chr22)
vcf_rng_gnomad_chek2 <- readVcf(tab_gnomad, "hg38", param=rng_chek2)

gnomad_fixed <- as.data.frame(rowRanges(vcf_rng_gnomad_chek2))

info <- info(vcf_rng_gnomad_chek2)
faf <- unlist(info$fafmax_faf95_max_joint)

gnomad_fixed$FAFmax_faf95_max_joint <- faf

gnomad_fixed$alt_values <- sapply(gnomad_fixed$ALT, function(x) as.character(x))

unique_alt_values <- unique(gnomad_fixed$alt_values)

gnomad_fixed$merge <- paste(gnomad_fixed$start, gnomad_fixed$REF, gnomad_fixed$alt_values, sep = " ")

clinvar_vcf <- "/Volumes/Seagate/Capstone/clinvar.vcf.gz"
vcf_clinvar <- readVcf(clinvar_vcf, "hg38") 

geneinfo_data <- info(vcf_clinvar)$GENEINFO

contains_chek2 <- grep("CHEK2", geneinfo_data, ignore.case = TRUE, value = TRUE)

chek2_variants <- grepl("CHEK2", geneinfo_data, ignore.case = TRUE)

vcf_chek2 <- vcf_clinvar[chek2_variants, ]

info_chek2 <- info(vcf_chek2)
chek2_info <- as.data.frame(info_chek2)

gr_chek2 <- rowRanges(vcf_chek2)
chek2_gr <- as.data.frame(gr_chek2)

chek2_df <- cbind(chek2_info, chek2_gr)

chek2_df_filtered <- chek2_df %>%
  select(start, REF, ALT, AF_ESP, ALLELEID, CLNHGVS, CLNSIG, MC, GENEINFO, RS)

chek2_df_filtered$alt_values <- sapply(chek2_df_filtered$ALT, function(x) as.character(x))

chek2_df_filtered$merge <- paste(chek2_df_filtered$start, chek2_df_filtered$REF, chek2_df_filtered$alt_values, sep = " ")

combined_inner_join <- inner_join(chek2_df_filtered, gnomad_fixed, by = "merge")



# Converting original code into a function 

library(GenomicRanges)
library(VariantAnnotation)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(org.Hs.eg.db)

process_gene_variants <- function(gene_name, gnomad_vcf_path, clinvar_vcf_path) {

# 1. Fetch genomic coordinates for the gene
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  gene_id <- mapIds(org.Hs.eg.db, keys = gene_name, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  
  if (is.na(gene_id)) {
    stop(paste("Gene", gene_name, "not found in annotation database."))
  }
  
  gene_coords <- genes(txdb, filter = list(gene_id = gene_id))
  
  if (length(gene_coords) == 0) {
    stop(paste("No genomic coordinates found for gene", gene_name))
  }
  
  gene_range <- gene_coords[1]  # Use the first result if there are multiple
  
# 2. Load and filter gnomAD data
  tab_gnomad <- TabixFile(gnomad_vcf_path)
  vcf_gnomad <- readVcf(tab_gnomad, "hg38", param = gene_range)
  
  gnomad_fixed <- as.data.frame(rowRanges(vcf_gnomad))
  gnomad_info <- info(vcf_gnomad)
  faf_values <- unlist(gnomad_info$fafmax_faf95_max_joint)
  
  gnomad_fixed$FAFmax_faf95_max_joint <- faf_values
  gnomad_fixed$alt_values <- sapply(gnomad_fixed$ALT, as.character)
  gnomad_fixed$merge <- paste(gnomad_fixed$start, gnomad_fixed$REF, gnomad_fixed$alt_values, sep = " ")
  
# 3. Load ClinVar data
  vcf_clinvar <- readVcf(clinvar_vcf_path, "hg38")
  geneinfo_data <- info(vcf_clinvar)$GENEINFO
  
  gene_variants <- grepl(gene_name, geneinfo_data, ignore.case = TRUE)
  vcf_gene <- vcf_clinvar[gene_variants, ]
  
  clinvar_info <- info(vcf_gene)
  clinvar_df <- as.data.frame(clinvar_info)
  clinvar_gr <- as.data.frame(rowRanges(vcf_gene))
  
  clinvar_combined <- cbind(clinvar_df, clinvar_gr)
  
  clinvar_filtered <- clinvar_combined %>% dplyr::select(start, REF, ALT, AF_ESP, ALLELEID, CLNHGVS, CLNSIG, MC, GENEINFO, RS) %>%
    mutate(
      alt_values = sapply(ALT, as.character),
      merge = paste(start, REF, alt_values, sep = " ")
    )
  
# 4. Combine gnomAD and ClinVar data using inner join
  combined_inner_join <- inner_join(clinvar_filtered, gnomad_fixed, by = "merge")
  
  return(combined_inner_join)
}

## Example usage 
# If you want multiple genes -> add them and separate them by commas (e.g. genes <- c("CHEK2", "BRCA1", "TP53"))
genes <- c("CHEK2")

results <- lapply(genes, function(gene) {
  process_gene_variants(
    gene_name = gene,
    gnomad_vcf_path = "/Volumes/Seagate/Capstone/gnomad.joint.v4.1.sites.chr22.vcf.bgz",
    clinvar_vcf_path = "/Volumes/Seagate/Capstone/clinvar.vcf.gz"
  )
})

# Access results for individual genes
chek2_data <- results[[1]]
# brca1_data <- results[[2]]
































