# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
BiocManager::install(c("VariantAnnotation", "GenomicRanges", "Rsamtools"))
library(VariantAnnotation)
library(GenomicRanges)
library(dplyr)

# Define regions of interest
define_regions <- function(seqname, start, end, names) {
  GRanges(seqnames = seqname, ranges = IRanges(start = start, end = end, names = names))
}

# Define region for CHEK2
region_chek2 <- define_regions("chr22", 28687820, 28742014, "CHEK2")
region_chek2_clinvar <- define_regions("22", 28687820, 28742014, "CHEK2")  # No 'chr' for ClinVar

# Extract gnomAD data
gnomad_vcf_file <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr22.vcf.bgz"
gnomad_data <- extract_vcf_data(gnomad_vcf_file, region_chek2) %>%
  adjust_gnomad_cadd_chr() %>%
  mutate(merge = paste(POS, REF, ALT, sep = " "))

# Extract ClinVar data
clinvar_vcf_file <- "F:/Capstone/Resources/ClinVar/clinvar.vcf.gz"
clinvar_data <- extract_vcf_data(clinvar_vcf_file, region_chek2_clinvar) %>%
  adjust_clinvar_chr() %>%
  mutate(merge = paste(POS, REF, ALT, sep = " "))

# Merge gnomAD and ClinVar data
merged_gnomad_clinvar <- inner_join(clinvar_data, gnomad_data, by = "merge")

# Extract CADD data
cadd_file <- "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"
cadd_data <- extract_cadd_data(cadd_file, region_chek2) %>%
  adjust_gnomad_cadd_chr() %>%
  mutate(merge = paste(POS, REF, ALT, sep = " "))

# Merge all data
final_combined_data <- merged_gnomad_clinvar %>%
  inner_join(cadd_data, by = "merge")

# Check final combined data
head(final_combined_data)

# Save the final combined dataset
write.table(
  final_combined_data,
  file = "F:/Capstone/Results/final_combined_data.txt",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)
