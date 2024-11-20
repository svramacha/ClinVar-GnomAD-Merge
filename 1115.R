## SET UP 

exonfile <- 'C:/Users/sramachandran/Documents/GitHub/Capstone_2024/Resources/CodingLocs_hg38.rds'
exonfile <- '/Users/anushaakhtar/Documents/GitHub/Capstone_2024/Resources/CodingLocs_hg38.rds'

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
library(BiocManager)

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
gnomad_vcf_chr22 <- "/Volumes/Seagate/Capstone/gnomad.joint.v4.1.sites.chr22.vcf.bgz"

# Open the VCF with TabixFile for subsetting
tab_gnomad <- TabixFile(gnomad_vcf_chr22)
vcf_rng_gnomad <- readVcf(tab_gnomad, "hg38", param=rng)

# Try for just CHEK2 
vcf_rng_gnomad_chek2 <- readVcf(tab_gnomad, "hg38", param=rng_chek2)

# Extract fixed fields for CHROM, POS, REF, ALT
gnomad_fixed <- as.data.frame(rowRanges(vcf_rng_gnomad_chek2))


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

# Load in clinvar data 
# Define Clinvar path to VCF 
clinvar_vcf <- "/Volumes/Seagate/Capstone/clinvar.vcf.gz"
vcf_clinvar <- readVcf(clinvar_vcf, "hg38") 

# Check and filter for gene chek2
# Extract the GENEINFO field from the INFO column
geneinfo_data <- info(vcf_clinvar)$GENEINFO

# Check if any entries in GENEINFO contain "CHEK2"
contains_chek2 <- grep("CHEK2", geneinfo_data, ignore.case = TRUE, value = TRUE)

# Filter the variants where GENEINFO contains "CHEK2"
chek2_variants <- grepl("CHEK2", geneinfo_data, ignore.case = TRUE)

# Subset the VCF to keep only rows where GENEINFO contains "CHEK2"
vcf_chek2 <- vcf_clinvar[chek2_variants, ]

# Check that the file is filtered only for chek2
geneinfo_values <- info(vcf_chek2)$GENEINFO
unique_geneinfo_values <- unique(geneinfo_values)

# Converting vcf_chek2 to a data frame
info_chek2 <- info(vcf_chek2)
chek2_info <- as.data.frame(info_chek2)

# Extract the genomic ranges from chek2 vcf and convert to data frame
gr_chek2 <- rowRanges(vcf_chek2)
chek2_gr <- as.data.frame(gr_chek2)

# Combine both data frames 
chek2_df <- cbind(chek2_info, chek2_gr)

# Filter the chek2 df
chek2_df_filtered <- chek2_df %>%
  select(start, REF, ALT, AF_ESP, ALLELEID, CLNHGVS, CLNSIG, MC, GENEINFO, RS)

# Convert the ALT column to character for ClinVar
chek2_df_filtered$alt_values <- sapply(chek2_df_filtered$ALT, function(x) as.character(x))

# Create new column 
chek2_df_filtered$merge <- paste(chek2_df_filtered$start, chek2_df_filtered$REF, chek2_df_filtered$alt_values, sep = " ")

### MERGE CLINVAR AND GNOMAD 
# Merge Clivar and Gnomad 
combined_inner_join <- inner_join(chek2_df_filtered, gnomad_fixed, by = "merge")


# Classifying pathogenicity for values labeled as "conflicting"
unique(combined_inner_join$CLNSIG)

# Cleaning the CLNSIG column 
combined_inner_join$CLNSIG <- trimws(combined_inner_join$CLNSIG)
combined_inner_join$CLNSIG <- as.factor(combined_inner_join$CLNSIG)
str(combined_inner_join$CLNSIG)
str(combined_inner_join$FAFmax_faf95_max_joint)

# Scatterplot to visualize the distribution of allele freq. and clingsig
library(ggplot2)

ggplot(combined_inner_join, aes(x = CLNSIG, y = FAFmax_faf95_max_joint, color = CLNSIG)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter points
  labs(
    title = "Scatterplot of Allele Frequency by Clinical Significance",
    x = "Clinical Significance",
    y = "Allele Frequency"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### Using k nearest neighbor to determine pathogenicity for conflicting values
library(FNN)
library(VIM)

# # Replace "conflicting" with NA in the 'CLNSIG' column
combined_inner_join$pathgen <- ifelse(grepl("conflicting", combined_inner_join$CLNSIG, ignore.case = TRUE), 
                                     NA, 
                                     combined_inner_join$CLNSIG)

# Check specific rows where CLNSIG is NA
combined_inner_join[is.na(combined_inner_join$pathgen), ]

# Perform k-NN imputation
imputed_data <- kNN(combined_inner_join, variable = "pathgen", dist_var = "FAFmax_faf95_max_joint", k = 5)

# Extract the imputed dataset
imputed_data_clean <- imputed_data[, 1:ncol(combined_inner_join)]

print("Imputed Dataset:")
print(head(imputed_data_clean))

# Verify imputed values
summary(imputed_data_clean$pathgen)

# Replace the original column with the imputed one
combined_inner_join$pathgen <- imputed_data$pathgen

# Compare original and imputed values
comparison <- data.frame(
  Original = combined_inner_join$CLNSIG,
  Imputed = imputed_data$pathgen
)

# Matching up imputed data with categorical variables
comparison_non_missing <- comparison[!is.na(comparison$Original), ]
print(comparison_non_missing)

category_mapping <- c("10" = "Uncertain_significance", 
                      "5" = "Likely_benign", 
                      "2" = "Benign/Likely_benign", 
                      "9" = "Pathogenic/Likely_pathogenic",
                      "6" = "Likely_pathogenic",
                      "1" = "Benign",
                      "7" = "not provided",
                      "8" = "Pathogenic")

# Covert imputed category "pathgen" into categories 
combined_inner_join$pathgen <- as.character(category_mapping[as.character(imputed_data$pathgen)])

# Scatterplot for visualization again 
ggplot(combined_inner_join, aes(x = pathgen, y = FAFmax_faf95_max_joint, color = pathgen)) +
  geom_point(size = 3, alpha = 0.7) +  # Scatter points
  labs(
    title = "Scatterplot of Allele Frequency by Clinical Significance",
    x = "Clinical Significance",
    y = "Allele Frequency"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Make x-axis labels vertical
    plot.title = element_text(hjust = 0.5)  # Center the title
  )

# Plotting what the conflicting values were assigned as 


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

