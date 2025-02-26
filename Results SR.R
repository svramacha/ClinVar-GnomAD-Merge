# Assign PHRED Classifications: 
final_combined_data_2 <- final_combined_data_2 %>%
  mutate(
    PHRED_Classification = case_when(
      PHRED >= threshold_table$pathogenic_very_strong ~ "Very Strong Pathogenic",
      PHRED >= threshold_table$pathogenic_strong ~ "Strong Pathogenic",
      PHRED >= threshold_table$pathogenic_moderate ~ "Moderate Pathogenic",
      PHRED >= threshold_table$pathogenic_supporting ~ "Supporting Pathogenic",
      PHRED <= threshold_table$benign_very_strong ~ "Very Strong Benign",
      PHRED <= threshold_table$benign_strong ~ "Strong Benign",
      PHRED <= threshold_table$benign_moderate ~ "Moderate Benign",
      PHRED <= threshold_table$benign_supporting ~ "Supporting Benign",
      TRUE ~ "Uncertain"
    )
  )


# Summary Statistics: Table 1

summary_table <- data.frame(
  Source = c("ClinVar", "gnomAD Benign"),
  Variants_Before_Filtering = c(nrow(pten_df), nrow(gnomad_fixed)),  # Raw counts before merging
  Variants_After_Filtering = c(
    nrow(final_combined_data_2[!is.na(final_combined_data_2$CLNSIG),]),  # ClinVar variants after filtering
    nrow(final_combined_data_2[is.na(final_combined_data_2$CLNSIG),])   # gnomAD variants after filtering
  ),
  Pathogenic = c(sum(final_combined_data_2$CLNSIG %in% c("Pathogenic", "Likely_pathogenic"), na.rm=TRUE), NA),
  Benign = c(sum(final_combined_data_2$CLNSIG %in% c("Benign", "Likely_benign"), na.rm=TRUE), NA),
  Uncertain = c(sum(final_combined_data_2$CLNSIG == "Uncertain", na.rm=TRUE), NA)
)

print(summary_table)

# Examine the model again 

summary(logistic_model)

# ROC and AUC analysis
library(pROC)
roc_curve <- roc(path_benign$pathogenic, predict(logistic_model, newdata = path_benign, type="response"))
auc_score <- auc(roc_curve)
print(paste("AUC Score:", auc_score))

# Compare PHRED-based classification to ClinVar 
phred_vs_clinvar <- table(final_combined_data_2$CLNSIG, final_combined_data_2$PHRED_Classification)
print(phred_vs_clinvar)

# See how this can be visualized 

library(ggplot2)
library(reshape2)

# Convert table to a dataframe
phred_vs_clinvar_df <- as.data.frame(as.table(phred_vs_clinvar))

# Rename columns for clarity
colnames(phred_vs_clinvar_df) <- c("ClinVar_Classification", "PHRED_Classification", "Count")

# Print to check format
print(head(phred_vs_clinvar_df))

ggplot(phred_vs_clinvar_df, aes(x = PHRED_Classification, y = ClinVar_Classification, fill = Count)) +
  geom_tile(color = "white") +  # Create heatmap
  scale_fill_gradient(low = "white", high = "red") +  # Color scale for intensity
  labs(title = "ClinVar vs PHRED Classification", 
       x = "PHRED-Based Classification", 
       y = "ClinVar Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

ggplot(final_combined_data_2, aes(x = PHRED, fill = CLNSIG)) +
  geom_density(alpha = 0.5) +
  labs(title = "PHRED Score Distribution by ClinVar Classification",
       x = "PHRED Score", y = "Density") +
  theme_minimal()

# Heatmap of ACMG Classification Thresholds
library(reshape2)
threshold_melt <- melt(threshold_table, id.vars = "annot")
threshold_melt$value <- as.numeric(as.character(threshold_melt$value))

ggplot(threshold_melt, aes(x = variable, y = annot, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  labs(title = "ACMG Pathogenicity Thresholds", x = "ACMG Category", y = "Annotation Score") +
  theme_minimal()

# Preformance Table: 
# Sensitivity & Specificity
conf_matrix <- table(path_benign$pathogenic, predict(logistic_model, type = "response") > 0.5)
sensitivity <- conf_matrix[2,2] / (conf_matrix[2,1] + conf_matrix[2,2])
specificity <- conf_matrix[1,1] / (conf_matrix[1,1] + conf_matrix[1,2])

# Store in a dataframe
performance_table <- data.frame(
  Metric = c("AUC", "Sensitivity", "Specificity"),
  PHRED_Score = c(auc_score, sensitivity, specificity)
)
print(performance_table)

# Breakdown of Uncertain Variants (ClinVar vs. Posterior Probability
uncertain_breakdown <- table(final_combined_data_2$CLNSIG[final_combined_data_2$CLNSIG == "Uncertain_significance"], 
                             final_combined_data_2$PHRED_Classification)
print(uncertain_breakdown)

###########################

# 1. ClinVar vs PHRED Classification (Contingency Table & Heatmap)

phred_vs_clinvar <- table(final_combined_data_2$CLNSIG, final_combined_data_2$PHRED_Classification)
print(phred_vs_clinvar)

# HeatMap 

library(ggplot2)
library(reshape2)

phred_vs_clinvar_df <- as.data.frame(as.table(phred_vs_clinvar))
colnames(phred_vs_clinvar_df) <- c("ClinVar_Classification", "PHRED_Classification", "Count")

ggplot(phred_vs_clinvar_df, aes(x = PHRED_Classification, y = ClinVar_Classification, fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "ClinVar vs PHRED Classification",
       x = "PHRED-Based Classification",
       y = "ClinVar Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Observed Benign and Pathogenic Variables (Restrict to ClinVar)

clinvar_only <- final_combined_data_2 %>%
  filter(!is.na(CLNSIG))  # Restrict to ClinVar variants

clinvar_phred_table <- table(clinvar_only$pathogenic, clinvar_only$PHRED_Classification)
print(clinvar_phred_table)

# 3. Expanding 2x8 to 4x8 (Break Down ClinVar into P/LP/B/LB)

clinvar_extended <- final_combined_data_2 %>%
  filter(CLNSIG %in% c("Pathogenic", "Likely_pathogenic", "Benign", "Likely_benign"))

clinvar_extended_table <- table(clinvar_extended$CLNSIG, clinvar_extended$PHRED_Classification)
print(clinvar_extended_table)

# Compare PHRED scores for Pathogenic vs Likely Pathogenic using boxplots:

ggplot(clinvar_extended, aes(x = CLNSIG, y = PHRED, fill = CLNSIG)) +
  geom_boxplot() +
  labs(title = "PHRED Scores for Pathogenic vs Likely Pathogenic Variants",
       x = "ClinVar Classification", y = "PHRED Score") +
  theme_minimal()

# 4. Investigating Non-ClinVar Benign Variants (gnomAD)

gnomad_benign <- final_combined_data_2 %>%
  filter(is.na(CLNSIG) & FAFmax_faf95_max_joint > 0.01)  # High allele frequency variants

gnomad_phred_table <- table(gnomad_benign$PHRED_Classification)
print(gnomad_phred_table)

# Compare ClinVar Benign vs gnomAD Benign using density plots:
ggplot(final_combined_data_2, aes(x = PHRED, fill = !is.na(CLNSIG))) +
  geom_density(alpha = 0.5) +
  labs(title = "PHRED Score Distributions: ClinVar Benign vs gnomAD Benign",
       x = "PHRED Score", y = "Density", fill = "ClinVar Label") +
  theme_minimal()

# 5. Assigning Uncertain Variants
uncertain_variants <- final_combined_data_2 %>%
  filter(CLNSIG == "Uncertain_significance")

uncertain_breakdown <- table(uncertain_variants$PHRED_Classification)
print(uncertain_breakdown)

# Compare PHRED scores of uncertain variants vs others
ggplot(final_combined_data_2 %>% filter(!is.na(PHRED)), 
       aes(x = PHRED, fill = CLNSIG == "Uncertain_significance")) +
  geom_density(alpha = 0.5) +
  labs(title = "PHRED Scores for Uncertain vs Classified Variants",
       x = "PHRED Score", y = "Density", fill = "Uncertain?") +
  theme_minimal()


# Violin Plot – PHRED Score Distributions by ClinVar Classification

ggplot(final_combined_data_2, aes(x = CLNSIG, y = PHRED, fill = CLNSIG)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "PHRED Score Distributions by ClinVar Classification",
       x = "ClinVar Classification", y = "PHRED Score") +
  theme_minimal()

# ROC Curve 
library(pROC)
roc_curve <- roc(path_benign$pathogenic, predict(logistic_model, newdata = path_benign, type="response"))
plot(roc_curve, col = "blue", main = "ROC Curve for PHRED Score Classification")
auc(roc_curve) 

# PHRED vs Allele Frequency 

library(ggplot2)
library(dplyr)

# Remove NA values from FAFmax_faf95_max_joint and PHRED before plotting
filtered_data <- final_combined_data_2 %>%
  filter(!is.na(FAFmax_faf95_max_joint), !is.na(PHRED))

# Scatter plot of PHRED Score vs. Allele Frequency (gnomAD)
ggplot(filtered_data, aes(x = FAFmax_faf95_max_joint, y = PHRED, color = CLNSIG)) +
  geom_point(alpha = 0.5) +
  scale_x_log10() +  # Log scale for allele frequency
  labs(title = "PHRED Score vs. Allele Frequency (gnomAD)",
       x = "Allele Frequency (log scale)", y = "PHRED Score") +
  theme_minimal()

# Classifying Uncertain Variants – Logistic Regression for Probabilistic Classification

logistic_model <- glm(pathogenic ~ PHRED, data = path_benign, family = binomial())
final_combined_data_2$predicted_prob <- predict(logistic_model, newdata = final_combined_data_2, type="response")

ggplot(final_combined_data_2, aes(x = PHRED, y = predicted_prob, color = CLNSIG)) +
  geom_point(alpha = 0.5) +
  labs(title = "Predicted Probability of Pathogenicity",
       x = "PHRED Score", y = "Predicted Probability") +
  theme_minimal()

## #################
## Tables Based on Andrew's Paper
## ################

# Table 1: Mean pathogenicity prediction scores and areas under the receiver operating curve for Benign and Pathogenic variants  

library(dplyr)
library(pROC)

# Filter for Benign (-1) and Pathogenic (1) variants only
path_benign <- final_combined_data_2 %>%
  filter(pathogenic %in% c(-1, 1), !is.na(PHRED))  # Ensure no missing PHRED values

# Calculate mean PHRED scores for Benign and Pathogenic variants
benign_mean <- mean(path_benign$PHRED[path_benign$pathogenic == -1], na.rm = TRUE)
pathogenic_mean <- mean(path_benign$PHRED[path_benign$pathogenic == 1], na.rm = TRUE)

# Wilcoxon test to compare PHRED scores between Benign and Pathogenic variants
p_value <- wilcox.test(
  path_benign$PHRED[path_benign$pathogenic == -1], 
  path_benign$PHRED[path_benign$pathogenic == 1], 
  alternative = "two.sided"
)$p.value

# Compute ROC curve and AUC
roc_curve <- roc(path_benign$pathogenic, path_benign$PHRED)
auc_value <- auc(roc_curve)

# Create a summary table similar to Table 1
phred_summary_table <- data.frame(
  Tool = "PHRED",
  Benign = benign_mean,
  Pathogenic = pathogenic_mean,
  `p-value (difference)` = p_value,
  AUC = auc_value
)

# Print the table
print(phred_summary_table)


# Table 2: Estimated change in the log(OR) as a function of pathogenicity prediction tool scores

library(dplyr)
library(pROC)

# Ensure we only use Benign (-1) and Pathogenic (1) variants
path_benign <- final_combined_data_2 %>%
  filter(pathogenic %in% c(-1, 1), !is.na(PHRED)) %>%
  mutate(pathogenic = ifelse(pathogenic == -1, 0, 1))  # Convert -1 to 0

# Fit a logistic regression model
logistic_model <- glm(pathogenic ~ PHRED, data = path_benign, family = binomial())

# Extract model summary
model_summary <- summary(logistic_model)

# Get log(OR), SE, and p-value
log_or <- model_summary$coefficients["PHRED", "Estimate"]
se_log_or <- model_summary$coefficients["PHRED", "Std. Error"]
p_value <- model_summary$coefficients["PHRED", "Pr(>|z|)"]

# Compute AUC
roc_curve <- roc(path_benign$pathogenic, path_benign$PHRED)
auc_value <- auc(roc_curve)

# Create a table similar to Table 2
phred_log_or_table <- data.frame(
  Tool = "PHRED",
  Coef = log_or,
  `SE (Coef)` = se_log_or,
  `p-value` = p_value,
  AUC = auc_value
)

# Print the table
print(phred_log_or_table)

# Table 3: Estimated score threshold intervals for the four pathogenicity prediction tools evaluated in this study as they relate to the different pathogenic and benign strength levels 

library(dplyr)

# Filter out missing PHRED values
filtered_data <- final_combined_data_2 %>%
  filter(!is.na(PHRED))

# Separate benign (-1) and pathogenic (1) variants
benign_scores <- filtered_data$PHRED[filtered_data$pathogenic == -1]
pathogenic_scores <- filtered_data$PHRED[filtered_data$pathogenic == 1]

# Define threshold percentiles (adjustable based on biological interpretation)
benign_thresholds <- quantile(benign_scores, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)
pathogenic_thresholds <- quantile(pathogenic_scores, probs = c(0, 0.25, 0.50, 0.75, 1), na.rm = TRUE)

# Create a table similar to Table 3
phred_threshold_table <- data.frame(
  Tool = "PHRED",
  `Very Strong Benign` = paste("≤", round(benign_thresholds[2], 2)),  # 0-25% quantile
  `Strong Benign` = paste("(", round(benign_thresholds[2], 2), ",", round(benign_thresholds[3], 2), "]"),
  `Moderate Benign` = paste("(", round(benign_thresholds[3], 2), ",", round(benign_thresholds[4], 2), "]"),
  `Supporting Benign` = paste("(", round(benign_thresholds[4], 2), ",", round(benign_thresholds[5], 2), "]"),
  `Supporting Pathogenic` = paste("[", round(pathogenic_thresholds[2], 2), ",", round(pathogenic_thresholds[3], 2), ")"),
  `Moderate Pathogenic` = paste("[", round(pathogenic_thresholds[3], 2), ",", round(pathogenic_thresholds[4], 2), ")"),
  `Strong Pathogenic` = paste("≥", round(pathogenic_thresholds[4], 2)),
  `Very Strong Pathogenic` = "—"
)

# Print the table
print(phred_threshold_table)

# Table 4: Number of predicted pathogenic and predicted benign variants at different strengths among the set of Uncertain variants 

library(dplyr)

# Step 1: Fit logistic regression model using PHRED scores for known Benign (-1) and Pathogenic (1) variants
path_benign <- final_combined_data_2 %>%
  filter(pathogenic %in% c(-1, 1), !is.na(PHRED)) %>%
  mutate(pathogenic = ifelse(pathogenic == -1, 0, 1))  # Convert -1 to 0

logistic_model <- glm(pathogenic ~ PHRED, data = path_benign, family = binomial())

# Step 2: Predict probabilities for Uncertain Variants
uncertain_variants <- final_combined_data_2 %>%
  filter(CLNSIG == "Uncertain_significance", !is.na(PHRED))

uncertain_variants$prob_pathogenic <- predict(logistic_model, newdata = uncertain_variants, type = "response")

# Step 3: Classify variants based on probability thresholds
uncertain_variants <- uncertain_variants %>%
  mutate(
    classification = case_when(
      prob_pathogenic < 0.2 ~ "Very Strong Benign",
      prob_pathogenic >= 0.2 & prob_pathogenic < 0.4 ~ "Strong Benign",
      prob_pathogenic >= 0.4 & prob_pathogenic < 0.55 ~ "Moderate Benign",  # Adjusted cutoff
      prob_pathogenic >= 0.55 & prob_pathogenic < 0.65 ~ "Indeterminate",
      prob_pathogenic >= 0.65 & prob_pathogenic < 0.8 ~ "Supporting Pathogenic",
      prob_pathogenic >= 0.8 & prob_pathogenic < 0.9 ~ "Moderate Pathogenic",
      prob_pathogenic >= 0.9 ~ "Strong Pathogenic",
      TRUE ~ "Indeterminate"
    )
  )

# Step 4: Count the number of variants in each category
phred_uncertain_table <- table(uncertain_variants$classification)

# Convert to DataFrame
phred_uncertain_df <- as.data.frame(phred_uncertain_table)
colnames(phred_uncertain_df) <- c("Category", "Count")

# Print the table
print(phred_uncertain_df)
# Print the table
print(phred_uncertain_df)


