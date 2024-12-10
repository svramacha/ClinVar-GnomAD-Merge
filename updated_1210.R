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

library(stringr)

install.packages("tidyr")
library(tidyr)

# Load the Biostrings package
library(Biostrings)


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

# RUN THIS ONE!!!!!
# USE THIS FOR THE PURPOSE OF CHEK2
# Try for just CHEK2 
tab_gnomad <- TabixFile(gnomad_vcf_chr22)
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

## Add in allele frequency (the faf95_max_joint) 
info <- info(vcf_rng_gnomad_chek2)
faf <- unlist(info$fafmax_faf95_max_joint)
# Extract the FAFmax_faf95_max_joint allele frequency from INFO fields
gnomad_fixed$FAFmax_faf95_max_joint <- faf

### Now we have created a DF that has gnomad data with the allele freq 

# Convert the ALT column to character
gnomad_fixed$alt_values <- sapply(gnomad_fixed$ALT, function(x) as.character(x))

# Print the first few rows to verify
head(gnomad_fixed)

# Get unique values from the new column
unique_alt_values <- unique(gnomad_fixed$alt_values)

# Print the unique values
print(unique_alt_values)

# This was to change 
gnomad_fixed <- gnomad_fixed %>%
  rename(start = POS)

gnomad_fixed$merge <- paste(gnomad_fixed$start, gnomad_fixed$REF, gnomad_fixed$alt_values, sep = " ")

### MERGE CLINVAR AND GNOMAD 
# Convert the ALT column to character for ClinVar
chek2_df_filtered$alt_values <- sapply(chek2_df_filtered$ALT, function(x) as.character(x))

# Create new column 
chek2_df_filtered$merge <- paste(chek2_df_filtered$start, chek2_df_filtered$REF, chek2_df_filtered$alt_values, sep = " ")

# Merge Clivar and Gnomad 
combined_inner_join <- inner_join(chek2_df_filtered, gnomad_fixed, by = "merge")

print(combined_inner_join)

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
final_combined_data_2 <- merge(
  combined_inner_join, 
  cadd_df, 
  by = "merge", 
  all.x = TRUE  # Keep all rows from combined_inner_join
)

# Download 
write.csv(final_combined_data_2, "combined_df_2.csv", row.names = FALSE)

# Now adding in the logistic regression 

# Load required packages
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages("pROC")
library(pROC)

# Prepare data for modeling
# Recode pathogenicity (1 for Pathogenic/Likely pathogenic, 0 for others)
final_combined_data_2 <- final_combined_data_2 %>%
  mutate(
    pathogenic = ifelse(
      CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"),
      1,
      0
    )
  )

# Check data summary
table(final_combined_data_2$pathogenic)

# Build a logistic regression model using CADD PHRED score
logistic_model <- glm(pathogenic ~ PHRED, data = final_combined_data_2, family = binomial())

# Model summary
summary(logistic_model)

# Filter the dataset
filtered_data <- final_combined_data_2 %>%
  filter(!is.na(PHRED))

# Predict probabilities
filtered_data <- filtered_data %>%
  mutate(predicted_prob = predict(logistic_model, newdata = filtered_data, type = "response"))

# Derive Thresholds 

# First calculate the quantiles of PHRED scores 
quantiles_cadd <- quantile(filtered_data$PHRED, probs = seq(0, 1, 0.25), na.rm = TRUE)
print(quantiles_cadd)

# Generate ROC curve
roc_curve <- roc(filtered_data$pathogenic, filtered_data$PHRED)

# Plot the ROC curve
plot(roc_curve, main = "ROC Curve for CADD PHRED Scores")
abline(a = 0, b = 1, lty = 2, col = "gray")  # Diagonal line for random chance

# Optimal threshold using Youden's index
optimal_threshold <- coords(roc_curve, "best", ret = "threshold", best.method = "youden")
print(optimal_threshold)

# Visualize predictions 

# Ensure PHRED is numeric
filtered_data$PHRED <- as.numeric(filtered_data$PHRED)

# Convert pathogenic to factor - so that it is categorical data 
filtered_data$pathogenic <- as.factor(filtered_data$pathogenic)

# Extract the numeric value from the data frame
optimal_threshold <- optimal_threshold$threshold

# Visualize predicted probabilities vs PHRED scores
ggplot(filtered_data, aes(x = PHRED, y = predicted_prob, color = pathogenic)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = optimal_threshold, linetype = "dashed", color = "red") +
  labs(title = "Predicted Pathogenicity vs. CADD PHRED Score",
       x = "CADD PHRED Score", y = "Predicted Probability") +
  scale_color_manual(name = "Pathogenic", values = c("0" = "blue", "1" = "red"))

# ## FROM ANDREW'S CODE 

# Create dis_ben column - based on pathogenicity and is numeric
filtered_data <- filtered_data %>%
  mutate(dis_ben = ifelse(pathogenic == 1, "Pathogenic", "Not Pathogenic"))

# Create dis_ben_numeric column
filtered_data <- filtered_data %>%
  mutate(dis_ben_numeric = ifelse(dis_ben == "Pathogenic", 1, 0))

# Set probability of pathogenicity and number of imputations
prop_path <- 0.04  # 4% probability of pathogenicity
num_imputations <- 25  # Number of imputations

# Initialize storage for model estimates
model_estimates_list <- list()

# Run imputation
for (i in 1:num_imputations) {
  # Weight adjustment based on prop_path
  filtered_data <- filtered_data %>%
    mutate(weight = ifelse(dis_ben == "Pathogenic", prop_path / mean(dis_ben_numeric == 1), 1))
  
  # Logistic regression using weighted dataset
  logistic_model <- glm(dis_ben_numeric ~ PHRED, family = binomial(), data = filtered_data, weights = weight)
  
  # Store model estimates
  model_estimates_list[[i]] <- summary(logistic_model)$coefficients
}

# Summarize results
estimates_summary <- do.call(rbind, model_estimates_list)
estimates_summary <- data.frame(estimates_summary)
colnames(estimates_summary) <- c("Estimate", "StdError", "zValue", "pValue")

# Print summary
print(estimates_summary)


   


