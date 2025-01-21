# Building logistic regression model with chromosome 10

library(dplyr)
library(ggplot2)
library(pROC)

# Define region of interest
# Specific for CHR 10
rng <- GRanges(seqnames="chr10", ranges=IRanges(
  start=c(89765432), 
  end=c(89774321),
  names=c("gene_123456")))

# Define region of PTEN: 
rng_pten <- GRanges(seqnames = "chr10", ranges = IRanges(
  start = 87863114,
  end = 87971941,
  names = "PTEN"
))

# Path to gnomAD VCF for chromosome 10
gnomad_vcf_chr10 <- "F:/Capstone/Resources/gnomAD/gnomad.joint.v4.1.sites.chr10.vcf.bgz"

# Open the VCF with TabixFile for subsetting
tab_gnomad <- TabixFile(gnomad_vcf_chr10)
vcf_rng_gnomad_pten <- readVcf(tab_gnomad, "hg38", param=rng_pten)

# Extract fixed fields for CHROM, POS, REF, ALT
gnomad_fixed <- as.data.frame(rowRanges(vcf_rng_gnomad_pten))

## Add in allele frequency (the faf95_max_joint) 
info <- info(vcf_rng_gnomad_pten)
faf <- unlist(info$fafmax_faf95_max_joint)

# Extract the FAFmax_faf95_max_joint allele frequency from INFO fields
gnomad_fixed$FAFmax_faf95_max_joint <- faf

# Convert the ALT column to character
gnomad_fixed$alt_values <- sapply(gnomad_fixed$ALT, function(x) as.character(x))

# Get unique values from the new column
unique_alt_values <- unique(gnomad_fixed$alt_values)

# Print the unique values
print(unique_alt_values)

gnomad_fixed$merge <- paste(gnomad_fixed$start, gnomad_fixed$REF, gnomad_fixed$alt_values, sep = " ")

# Add in Clinvar data 
clinvar_vcf <- 'F:/Capstone/Resources/ClinVar/clinvar.vcf.gz'
vcf_clinvar <- readVcf(clinvar_vcf, "hg38")

geneinfo_data <- info(vcf_clinvar)$GENEINFO

contains_pten <- grep("PTEN", geneinfo_data, ignore.case = TRUE, value = TRUE)

pten_variants <- grepl("PTEN", geneinfo_data, ignore.case = TRUE)

vcf_pten <- vcf_clinvar[pten_variants, ]

info_pten <- info(vcf_pten)
pten_info <- as.data.frame(info_pten)

gr_pten <- rowRanges(vcf_pten)
pten_gr <- as.data.frame(gr_pten)

pten_df <- cbind(pten_info, pten_gr)

### MERGE CLINVAR AND GNOMAD 
# Convert the ALT column to character for ClinVar
pten_df$alt_values <- sapply(pten_df$ALT, function(x) as.character(x))

# Create new column 
pten_df$merge <- paste(pten_df$start, pten_df$REF, pten_df$alt_values, sep = " ")

# Outer join to merge Clivar and Gnomad
combined_outer_join_pten <- merge(pten_df, gnomad_fixed, by = "merge", all = TRUE)

# Adding CADD
cadd_file <- "F:/Capstone/Resources/CADD/v1.7/whole_genome_SNVs.tsv.gz"

# Open the CADD file with TabixFile
tabix_file <- TabixFile(cadd_file)

# Query the CADD file for ch10
cadd_data <- scanTabix(tabix_file, param = rng)

# Convert the extracted data into a dataframe
cadd_df <- read.table(text = unlist(cadd_data), header = FALSE, sep = "\t")

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

# Adding in logistic regression 

# Recode pathogenicity (1 for Pathogenic/Likely pathogenic, -1 for Benign/Likely benign, and 0 for VUS/uncertain)
# Assigns "pathogenic" variants first, then "benign", and if neither then as "VUS"

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

# Same assingment as above with faf_cutoff (didn't run)
faf_cutoff <- 10^-4
prop_path <- 0.4

final_combined_data_2 <- final_combined_data_2 %>%
  mutate(
    pathogenic = ifelse(
      CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"), 
      1,  # Assign 1 for pathogenic classifications
      ifelse(
        CLNSIG %in% c("Benign", "Likely_benign", "Benign/Likely_benign"), 
        -1,  # Assign -1 for benign classifications
        ifelse(
          is.na(CLNSIG) & FAFmax_faf95_max_joint > faf_cutoff, 
          -1,  # Assign -1 if FAFmax_faf95_max_joint exceeds faf_cutoff and CLNSIG is NA
          0  # Default to 0 if none of the above conditions are met
        )
      )
    )
  )

# First regression preformed for begnin/pathogenic - have to create a seperate df for this regression 
# create lm model data 
# take out anything = 0 
lm_data <- final_combined_data_2 %>%
  filter(pathogenic != 0)

# Build a logistic regression model using CADD PHRED score
# Filter the dataset - filter this before hand in the final combined data part (move this up)
filtered_data <- final_combined_data_2 %>%
  filter(!is.na(PHRED))

# Create new dataset (path/benign) filtering out 0's
path_benign <- filtered_data %>%
  filter(pathogenic != 0) %>% mutate(pathogenic=factor(pathogenic, levels = c("-1","1"))) 

logistic_model <- glm(pathogenic ~ PHRED, data = path_benign, family = binomial())

# Model summary
summary(logistic_model)






