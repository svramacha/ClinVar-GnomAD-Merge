# Define Uncertain Variants from ClinVar
clinvar_uncertain_variants <- final_combined_data_2 %>%
  filter(CLNSIG == "Uncertain_significance", !is.na(PHRED))

# Define Uncertain Variants based on gnomAD frequency threshold (e.g., FAFmax_faf95_max_joint > 0.01)
gnomad_uncertain_variants <- final_combined_data_2 %>%
  filter(is.na(CLNSIG) & FAFmax_faf95_max_joint > 0.01, !is.na(PHRED))

# Predict probabilities using the logistic model for both groups
clinvar_uncertain_variants$prob_pathogenic <- predict(logistic_model, newdata = clinvar_uncertain_variants, type = "response")
gnomad_uncertain_variants$prob_pathogenic <- predict(logistic_model, newdata = gnomad_uncertain_variants, type = "response")

# Classify both sets of uncertain variants
clinvar_uncertain_variants <- clinvar_uncertain_variants %>%
  mutate(
    Source = "ClinVar",
    classification = case_when(
      prob_pathogenic < 0.2 ~ "Very Strong Benign",
      prob_pathogenic >= 0.2 & prob_pathogenic < 0.4 ~ "Strong Benign",
      prob_pathogenic >= 0.4 & prob_pathogenic < 0.55 ~ "Moderate Benign",
      prob_pathogenic >= 0.55 & prob_pathogenic < 0.65 ~ "Indeterminate",
      prob_pathogenic >= 0.65 & prob_pathogenic < 0.8 ~ "Supporting Pathogenic",
      prob_pathogenic >= 0.8 & prob_pathogenic < 0.9 ~ "Moderate Pathogenic",
      prob_pathogenic >= 0.9 ~ "Strong Pathogenic",
      TRUE ~ "Indeterminate"
    )
  )

gnomad_uncertain_variants <- gnomad_uncertain_variants %>%
  mutate(
    Source = "gnomAD",
    classification = case_when(
      prob_pathogenic < 0.2 ~ "Very Strong Benign",
      prob_pathogenic >= 0.2 & prob_pathogenic < 0.4 ~ "Strong Benign",
      prob_pathogenic >= 0.4 & prob_pathogenic < 0.55 ~ "Moderate Benign",
      prob_pathogenic >= 0.55 & prob_pathogenic < 0.65 ~ "Indeterminate",
      prob_pathogenic >= 0.65 & prob_pathogenic < 0.8 ~ "Supporting Pathogenic",
      prob_pathogenic >= 0.8 & prob_pathogenic < 0.9 ~ "Moderate Pathogenic",
      prob_pathogenic >= 0.9 ~ "Strong Pathogenic",
      TRUE ~ "Indeterminate"
    )
  )

# Combine both into a single table
all_uncertain_variants <- bind_rows(clinvar_uncertain_variants, gnomad_uncertain_variants)

# Count the number of variants in each classification category, split by source (ClinVar or gnomAD)
phred_uncertain_table <- table(all_uncertain_variants$Source, all_uncertain_variants$classification)

# Convert to DataFrame for better readability
phred_uncertain_df_2 <- as.data.frame(phred_uncertain_table)
colnames(phred_uncertain_df_2) <- c("Source", "Category", "Count")

# Print the new table
print(phred_uncertain_df_2)