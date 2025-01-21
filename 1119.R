## SET UP 
exonfile <- 'C:/Users/sramachandran/Documents/GitHub/Capstone_2024/Resources/CodingLocs_hg38.rds'
readRDS(exonfile)
exonpositions <- readRDS(exonfile)

# Install BiocManager to handle Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install VariantAnnotation for reading VCF files
BiocManager::install("VariantAnnotation")

# Install GenomicRanges for handling genomic regions
BiocManager::install("GenomicRanges")

# Install dplyr for data manipulation
install.packages("dplyr")

# Load the VariantAnnotation package for VCF file handling
library(VariantAnnotation)

# Load GenomicRanges for defining and handling regions of interest
library(GenomicRanges)

# Load dplyr for data manipulation
library(dplyr)

# Define region of interest
# Specific for CHR 22
rng <- GRanges(seqnames="chr22", ranges=IRanges(
  start=c(50301422, 50989541), 
  end=c(50312106, 51001328),
  names=c("gene_79087", "gene_644186")))

# Define region of CHEK2: 
rng_chek2 <- GRanges(seqnames = "chr22", ranges = IRanges(
  start = 28687820,
  end = 28742014,
  names = "CHEK2"
))

# Path to gnomAD VCF for chromosome 22
gnomad_vcf_chr22 <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr22.vcf.bgz"

# Open the VCF with TabixFile for subsetting
tab_gnomad <- TabixFile(gnomad_vcf_chr22)
vcf_rng_gnomad <- readVcf(tab_gnomad, "hg38", param=rng)

# Try for just CHEK2 
vcf_rng_gnomad_chek2 <- readVcf(tab_gnomad, "hg38", param=rng_chek2)

# Extract fixed fields for CHROM, POS, REF, ALT
gnomad_fixed <- as.data.frame(rowRanges(vcf_rng_gnomad_chek2))

# Rename columns for clarity
gnomad_fixed <- gnomad_fixed %>%
  rename(
    CHROM = seqnames,
    POS = start
  )

# Check contents
head(gnomad_fixed)

### Having trouble addin this fafmax

##
info <- info(vcf_rng_gnomad_chek2)
faf <- unlist(info$fafmax_faf95_max_joint)

# Extract the FAFmax_faf95_max_joint allele frequency from INFO fields
gnomad_fixed$FAFmax_faf95_max_joint <- faf

# Check resulting gnomAD data frame
head(gnomad_fixed)

print(gnomad_fixed$ALT)

library(stringr)

install.packages("tidyr")
library(tidyr)

# Load the Biostrings package
library(Biostrings)

# Convert the ALT column to character
gnomad_fixed$alt_values <- sapply(gnomad_fixed$ALT, function(x) as.character(x))

# Print the first few rows to verify
head(gnomad_fixed)

# Get unique values from the new column
unique_alt_values <- unique(gnomad_fixed$alt_values)

# Print the unique values
print(unique_alt_values)

gnomad_fixed$merge <- paste(gnomad_fixed$start, gnomad_fixed$REF, gnomad_fixed$alt_values, sep = " ")


## Choose one of the large pop 
# Should NA be considered rare or is it number of samples
# how many sampels from each pop and what is the count 
# maybe get rid of small indels 

### MERGE CLINVAR AND GNOMAD 
# Convert the ALT column to character for ClinVar
chek2_df_filtered$alt_values <- sapply(chek2_df_filtered$ALT, function(x) as.character(x))

# Create new column 
chek2_df_filtered$merge <- paste(chek2_df_filtered$start, chek2_df_filtered$REF, chek2_df_filtered$alt_values, sep = " ")

# Merge Clivar and Gnomad 
combined_inner_join <- inner_join(chek2_df_filtered, gnomad_fixed, by = "merge")

## DOuble checks 

# Count the number of duplicate rows
num_duplicates <- sum(duplicated(combined_df))

# Print the number of duplicate rows
print(num_duplicates)

# Count the number of unique rows
num_unique_rows <- nrow(unique(combined_df))

# Print the number of unique rows
print(num_unique_rows)

# Find the lowest value in the 'start' column
lowest_value <- min(chek2_df_filtered$start, na.rm = TRUE)

# Print the lowest value
print(lowest_value)

# Find the highest value in the 'start' column
highest_value <- max(chek2_df_filtered$start, na.rm = TRUE)

# Print the highest value
print(highest_value)

# Progress since 10/29

#######

# ClinVar Extraction 
clinvar_vcf <- "F:/Capstone/Resources/ClinVar/clinvar.vcf.gz"
vcf_clinvar <- readVcf(clinvar_vcf, "hg38")

# Open the VCF with TabixFile for subsetting
tab_clinvar <- TabixFile(clinvar_vcf)
vcf_rng_clinvar <- readVcf(tab_clinvar, "hg38", param=rng)

# Extract fixed fields for CHROM, POS, REF, ALT
clinvar_fixed <- as.data.frame(rowRanges(vcf_rng_clinvar))[, c("seqnames", "start", "REF", "ALT")]

## Looking at distributions 

# Find percentage of clinically significant variants for CHEK2 
total_variants <- nrow(combined_inner_join)
clinically_significant_variants <- sum(combined_inner_join$CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"))
percentage_clinically_significant <- (clinically_significant_variants / total_variants) * 100
percentage_clinically_significant

# show quantiles for fafmax
hist(combined_inner_join$FAFmax_faf95_max_joint, breaks = 30, main = "Distribution of FAFmax_faf95_max_joint",
     xlab = "FAFmax_faf95_max_joint", ylab = "Frequency")
quantiles_fafmax <- quantile(combined_inner_join$FAFmax_faf95_max_joint, probs = seq(0, 1, 0.25), na.rm = TRUE)
quantiles_fafmax

# look at allele freq distribution for clinically sig variants 

# Filter the data for clinically significant variants
clinically_significant_data <- subset(combined_inner_join, CLNSIG %in% c("Pathogenic", "Likely_pathogenic","Pathogenic/Likely_pathogenic" ))

# Plot histogram of FAFmax_faf95_max_joint for clinically significant variants
hist(clinically_significant_data$FAFmax_faf95_max_joint, breaks = 30, main = "FAFmax_faf95_max_joint for Clinically Significant Variants",
     xlab = "FAFmax_faf95_max_joint", ylab = "Frequency")

# Calculate quantiles for FAFmax_faf95_max_joint for clinically significant variants
quantiles_fafmax_significant <- quantile(clinically_significant_data$FAFmax_faf95_max_joint, probs = seq(0, 1, 0.25), na.rm = TRUE)
quantiles_fafmax_significant

## Check to see if one is skews data 
# For "Pathogenic" variants
pathogenic_data <- subset(combined_inner_join, CLNSIG == "Pathogenic")
median_pathogenic <- quantile(pathogenic_data$FAFmax_faf95_max_joint, probs = 0.5, na.rm = TRUE)
median_pathogenic

# For "Likely pathogenic" variants
likely_pathogenic_data <- subset(combined_inner_join, CLNSIG == "Likely_pathogenic")
median_likely_pathogenic <- quantile(likely_pathogenic_data$FAFmax_faf95_max_joint, probs = 0.5, na.rm = TRUE)
median_likely_pathogenic

# For "Pathogenic/Likely pathogenic" variants
pathogenic_likely_pathogenic_data <- subset(combined_inner_join, CLNSIG == "Pathogenic/Likely_pathogenic")
median_pathogenic_likely_pathogenic <- quantile(pathogenic_likely_pathogenic_data$FAFmax_faf95_max_joint, probs = 0.5, na.rm = TRUE)
median_pathogenic_likely_pathogenic

A <- sapply(chek2_df$CLNSIG, length)
table(A)

table(chek2_df$CLNSIG)

A <- sapply(chek2_df$CLNSIGCONF, length)
table(A)

I <- which(A==2)
chek2_df$CLNSIGCONF[I[1]]

I <- which(A==1)
chek2_df$CLNSIGCONF[I[1]]

## Downloading CADD 

if (!requireNamespace("Rsamtools", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rsamtools")

library(Rsamtools)

cadd_file <- "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"
# Open the CADD file with TabixFile
tabix_file <- TabixFile(cadd_file)

# Define manual regions of interest (update as necessary)
rng_manual <- GRanges(seqnames = "22", ranges = IRanges(
  start = c(28687820, 28700000), 
  end = c(28690000, 28705000),
  names = c("Region_1", "Region_2")
))

# Query the CADD file for the manually defined regions
cadd_data <- scanTabix(tabix_file, param = rng_manual)

# Convert the extracted data into a dataframe
cadd_df <- read.table(text = unlist(cadd_data), header = FALSE, sep = "\t")

# Checking the header 
# Specify the path to your CADD gzipped file
file_path <- "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"

# Open the gzipped file in text mode
con <- gzfile(file_path, "rt")

# Read the first few lines to check the header
header_lines <- readLines(con, n = 10)  # Adjust 'n' to read more lines if needed

# Close the file connection
close(con)

# Print the header lines
print(header_lines)

## Change the column names so that merge is possible 
colnames(cadd_df) <- c("CHROM", "start", "REF", "ALT", "RawScore", "PHRED")
head(cadd_df)

# Create merge column for CADD 
cadd_df$merge <- paste(cadd_df$start, cadd_df$REF, cadd_df$ALT, sep = " ")

# Combined with inner join and CADD 
final_combined_data <- merge(
  combined_inner_join, 
  cadd_df, 
  by = "merge", 
  all.x = TRUE  # Keep all rows from combined_inner_join
)

