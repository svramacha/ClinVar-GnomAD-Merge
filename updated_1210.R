## SET UP 
exonfile <- '~/GitHub/Capstone_2024/Resources/CodingLocs_hg38.rds'

readRDS(exonfile)
exonpositions <- readRDS(exonfile)

# Install BiocManager to handle Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

# Install VariantAnnotation for reading VCF files
BiocManager::install("VariantAnnotation")

# Install GenomicRanges for handling genomic regions
BiocManager::install("GenomicRanges")


library(BiocManager)

# Install dplyr for data manipulation
install.packages("dplyr")

library(dplyr)

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
# gnomad_fixed <- gnomad_fixed %>%
#   rename(
#     CHROM = seqnames,
#     POS = start
#   )

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
# gnomad_fixed <- gnomad_fixed %>%
#   rename(start = POS)

gnomad_fixed$merge <- paste(gnomad_fixed$start, gnomad_fixed$REF, gnomad_fixed$alt_values, sep = " ")

# Add in Clinvar data 
clinvar_vcf <- 'F:/Capstone/Resources/ClinVar/clinvar.vcf.gz'
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



### MERGE CLINVAR AND GNOMAD 
# Convert the ALT column to character for ClinVar
chek2_df$alt_values <- sapply(chek2_df$ALT, function(x) as.character(x))

# Create new column 
chek2_df$merge <- paste(chek2_df$start, chek2_df$REF, chek2_df$alt_values, sep = " ")



# Outer join to merge Clivar and Gnomad
combined_outer_join <- merge(chek2_df, gnomad_fixed, by = "merge", all = TRUE)

print(combined_outer_join)


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
  combined_outer_join, 
  cadd_df, 
  by = "merge", 
  all.x = TRUE  # Keep all rows from combined_inner_join
)

# Download 
write.csv(final_combined_data_2, "combined_df_2.csv", row.names = FALSE)



## Adding in the other annotations
# Add in BayesDel
library(R.utils)

# Decompress .log.gz file
gunzip("F:/Capstone/Resources/BayesDel/hg19/BayesDel_170824_addAF.log.gz", destname = "F:/Capstone/Resources/BayesDel/hg19/BayesDel_170824_addAF.log", remove = FALSE)
gunzip("F:/Capstone/Resources/BayesDel/hg19/BayesDel_170824_noAF.log.gz", destname = "F:/Capstone/Resources/BayesDel/hg19/BayesDel_170824_noAF.log", remove = FALSE)

# Extract .tgz file
untar("F:/Capstone/Resources/BayesDel/hg19/BayesDel_170824_addAF.tgz", exdir = "F:/Capstone/Resources/BayesDel/hg19")
untar("F:/Capstone/Resources/BayesDel/hg19/BayesDel_170824_noAF.tgz", exdir = "F:/Capstone/Resources/BayesDel/hg19")



# Now adding in the logistic regression 

# Load required packages
library(dplyr)
install.packages("ggplot2")
library(ggplot2)
install.packages("pROC")
library(pROC)
library(nnet) # For multinomial logistic regression

# Prepare data for modeling
unique(final_combined_data_2$CLNSIG)

# Recode pathogenicity (1 for Pathogenic/Likely pathogenic, -1 for Benign/Likely benign, and 0 for VUS/uncertain)
# Assigns "pathogenic" variants first, then "benign", and if neither then as "VUS"


#change 
faf_cutoff <- 10^-4
prop_path <- 0.4

# Recode pathogenicity (1 for Pathogenic/Likely pathogenic, 0 for others)

final_combined_data_2 <- final_combined_data_2 %>%
  mutate(pathogenic = 0,  
    pathogenic = ifelse(
      CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"), 
      1, 
      ifelse(
        CLNSIG %in% c("Benign", "Likely_benign", "Benign/Likely_benign"), 
        -1, 
        0
      )
    )
  )

# Check data summary
table(final_combined_data_2$pathogenic)
str(final_combined_data_2$pathogenic)
      CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic")
      1
      pathogenic 
    pathogenic = ifelse(
      CLNSIG %in% c("Benign", "Likely_benign", "Benign/Likely_benign"),
      -1,
      pathogenic) 
    pathogenic = ifelse(
      is.na(CLNSIG) & FAFmax_faf95_max_joint > faf_cutoff, 
      -1, 
      pathogenic 
    )
    
# Clinsig is NA - then it is benign 
# cutoff for gnomad freq will be an argument 
  

# Check data summary
table(final_combined_data_2$pathogenic)


# Weighted

# First regression preformed for begnin/pathogenic - have to create a seperate df for this regression 
# create lm model data 
# take out anything = 0 
lm_data <- final_combined_data_2 %>%
  filter(pathogenic != 0)

# create the weights - ratio of the frequency of 
# line 577/585 analysis_functions
# total number of observations 
# cases = pathogenic (1), control = bengin (-1)
# set prop_path = 0.4 - if i took a random position and changed the allele - prob that it would become pathogenic 
# also will have to do something similar to line 597 

# line 625 has the glm with the weight in it 

# Build a logistic regression model using CADD PHRED score
logistic_model <- glm(pathogenic ~ PHRED, data = final_combined_data_2, family = binomial())

# Model summary
summary(logistic_model)


# Filter the dataset - filter this before hand in the final combined data part (move this up)
filtered_data <- final_combined_data_2 %>%
  filter(!is.na(PHRED))

# Create new dataset (path/benign) filtering out 0's
path_benign <- filtered_data %>%
  filter(pathogenic != 0) %>% mutate(pathogenic=factor(pathogenic, levels = c("-1","1"))) 

# Use kernel density estimation to smooth the distributions of benign and pathogenic scores
# Kernel density estimation for PHRED scores
benign_density <- density(path_benign$PHRED[path_benign$pathogenic == -1])
pathogenic_density <- density(path_benign$PHRED[path_benign$pathogenic == 1])

# For reproducibility 
set.seed(123)

# Generate synthetic PHRED values
synthetic_benign <- sample(benign_density$x, size = 50, prob = benign_density$y / sum(benign_density$y), replace = TRUE)
synthetic_pathogenic <- sample(pathogenic_density$x, size = 50, prob = pathogenic_density$y / sum(pathogenic_density$y), replace = TRUE)

# Create synthetic dataframes
synthetic_data <- data.frame(
  PHRED = c(synthetic_benign, synthetic_pathogenic),
  pathogenic = factor(c(rep(-1, length(synthetic_benign)), rep(1, length(synthetic_pathogenic))),
                      levels = c("-1", "1"))
)

# Add missing columns from path_benign to synthetic_data
missing_cols <- setdiff(colnames(path_benign), colnames(synthetic_data))
for (col in missing_cols) {
  synthetic_data[[col]] <- NA # Assign NA for simplicity
}

# Reorder synthetic_data columns to match path_benign
synthetic_data <- synthetic_data[, colnames(path_benign)]

# Combine datasets
path_benign_wf <- rbind(path_benign, synthetic_data)

# Density plot to compare real vs synthetic data
ggplot() +
  geom_density(data = path_benign, aes(x = PHRED, fill = as.factor(pathogenic)), alpha = 0.4) +
  geom_density(data = synthetic_data, aes(x = PHRED, color = as.factor(pathogenic)), linetype = "dashed") +
  labs(title = "PHRED Density: Real vs Synthetic", fill = "Real Data", color = "Synthetic Data")


# Stopping point




# Sampling data
# data_fake = path_benign[sample(1:nrow(path_benign), size=15), ]
# data_fake$PHRED = data_fake$PHRED[sample(1:15)]
# path_benign_wf = rbind(path_benign, data_fake)

# Randomly shuffled PHRED
path_benign$fake_y = sample(path_benign$PHRED)

# Build a binomial logistic regression model using CADD PHRED score using new dataset (path/benign)
logistic_model <- glm(pathogenic ~ PHRED, data = path_benign_wf, family = binomial())

# Balancing distribution between regular and synthetic and data (DIDN'T WORK)
# library(ROSE) - use for binary imbalance
# install.packages("ROSE")
# path_benign_combined <- rbind(path_benign, synthetic_data)
# ggplot(path_benign_combined, aes(x = PHRED, fill = as.factor(pathogenic))) +
#   geom_density(alpha = 0.4) +
#   labs(title = "PHRED Score Distribution: Combined Data", x = "PHRED", fill = "Pathogenicity")
# balanced_data <- ROSE(pathogenic ~ PHRED, data = path_benign_combined, seed = 123)$data
# table(balanced_data$pathogenic)
# balanced_data$PHRED_scaled <- scale(balanced_data$PHRED)
# logistic_model <- glm(pathogenic ~ PHRED_scaled, data = balanced_data, family = binomial())
summary(logistic_model)



# Check on data distribution 
boxplot(PHRED ~ pathogenic, data = path_benign_wf, main = "PHRED by Pathogenicity")

library(glmnet)
X <- model.matrix(pathogenic ~ PHRED, data = path_benign_wf)[, -1]
y <- path_benign_wf$pathogenic

y <- as.numeric(as.factor(path_benign_wf$pathogenic)) - 1  # Convert to binary (0, 1)
length(y) == nrow(X)  # Must return TRUE

any(is.na(path_benign_wf$PHRED))  # Ensure there are no missing values

lasso_model <- cv.glmnet(X, y, family = "binomial", alpha = 1) # alpha = 1 for Lasso

#final_combined_data_2$pathogenic <- factor(final_combined_data_2$pathogenic, levels = c(-1, 0, 1))

#multinomial_model <- multinom(pathogenic ~ PHRED, data = final_combined_data_2)

# Model summary
summary(logistic_model)

# Using SMOTE to balance data
install.packages("smotefamily")
install.packages("caret")
install.packages("nnet")
library(smotefamily)
library(caret)
library(nnet)

# Convert data frame to numeric to insert into SMOTE
# path_benign_wf$pathogenic_num <- as.numeric(path_benign_wf$pathogenic)
#your_column <- factor(c(1, 0))
# Correct conversion
#numeric_column <- as.numeric(levels(your_column))[your_column]
#print(numeric_column)
path_benign_wf_num <- path_benign_wf[, sapply(path_benign_wf, is.numeric)] # convert all columns
#path_benign_wf_num$pathogenic_num <- as.numeric(levels(path_benign_wf$pathogenic))[path_benign_wf$pathogenic] # only convert pathogenic columns so levels stay correct

# if (is.factor(path_benign_wf$pathogenic)) {
#   # Convert the factor column `pathogenic` to numeric without altering its original numeric values
#   path_benign_wf$pathogenic_num <- as.numeric(levels(path_benign_wf$pathogenic))[path_benign_wf$pathogenic]
# } else {
#   # If `pathogenic` is already numeric, copy it directly
#   path_benign_wf$pathogenic_num <- path_benign_wf$pathogenic
# }

# # Step 1: Convert `pathogenic` factor to numeric without altering original numeric values
# path_benign_wf$pathogenic_num <- as.numeric(levels(path_benign_wf$pathogenic))[path_benign_wf$pathogenic]
# 
# # Step 2: Replace -1 values with 0
# path_benign_wf$pathogenic_num[path_benign_wf$pathogenic_num == -1] <- 0
# 
# # Verify the result
# table(path_benign_wf$pathogenic_num)  # Check the distribution

# Check for NAs after conversion
sum(is.na(path_benign_wf_num$pathogenic_num))
sum(is.na(path_benign_wf_num$PHRED))

# Keep only relevant columns
path_benign_wf_num = subset(path_benign_wf_num, select = c(pathogenic_num, PHRED))

smote_results <- SMOTE(X = path_benign_wf_num, target = path_benign_wf_num$pathogenic_num, K = 5, dup_size = 0)

# Combine results into one data frame to input into logistic regression 
oversampled_data <- data.frame(smote_results$data)
oversampled_data$pathogenic_num[oversampled_data$pathogenic_num == -1] <- 0

# Build a binomial logistic regression model using CADD PHRED score using new SMOTE dataset
logistic_model_SMOTE <- glm(pathogenic_num ~ PHRED, data = oversampled_data, family = binomial())

summary(logistic_model_SMOTE)
table(oversampled_data$pathogenic_num)



# Predict probabilities
filtered_data <- filtered_data %>%
  mutate(predicted_prob = predict(multinomial_model, newdata = filtered_data, type = "response"))

## Derive Thresholds 

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


   


