# Load necessary libraries
library(ggplot2)
library(dplyr)

# Ensure the pathogenic column is numeric
final_combined_data_2 <- final_combined_data_2 %>%
  mutate(pathogenic = as.numeric(pathogenic))

# Create a dataframe for variant counts BEFORE modeling
before_modeling <- data.frame(
  Source = rep(c("ClinVar", "gnomAD"), each = 3),
  Classification = rep(c("Pathogenic", "Benign", "VUS"), 2),
  Count = c(
    sum(pten_df$CLNSIG %in% c("Pathogenic", "Likely_pathogenic"), na.rm = TRUE),  # ClinVar Pathogenic
    sum(pten_df$CLNSIG %in% c("Benign", "Likely_benign"), na.rm = TRUE),          # ClinVar Benign
    sum(pten_df$CLNSIG %in% c("Uncertain significance", "VUS"), na.rm = TRUE),    # ClinVar VUS
    sum(is.na(gnomad_fixed$CLNSIG)),                                              # gnomAD Unknown
    sum(gnomad_fixed$FAFmax_faf95_max_joint > 0.0001, na.rm = TRUE),              # gnomAD Benign-like
    sum(gnomad_fixed$FAFmax_faf95_max_joint <= 0.0001, na.rm = TRUE)              # gnomAD VUS-like
  )
)

# Generate the bar plot BEFORE logistic regression
ggplot(before_modeling, aes(x = Classification, y = Count, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of Variants Before Logistic Regression", 
       x = "Variant Classification", 
       y = "Count") +
  theme_minimal()

# Ensure d_all is formatted correctly for AFTER modeling with ACMG bins
after_modeling <- d_all %>%
  mutate(ACMG_Classification = case_when(
    prob_disease >= 0.99 ~ "Very Strong Pathogenic",
    prob_disease >= 0.95 & prob_disease < 0.99 ~ "Strong Pathogenic",
    prob_disease >= 0.90 & prob_disease < 0.95 ~ "Moderate Pathogenic",
    prob_disease >= 0.50 & prob_disease < 0.90 ~ "Supporting Pathogenic",
    prob_disease > 0.10 & prob_disease < 0.50 ~ "Uncertain (VUS)",
    prob_disease <= 0.10 & prob_disease > 0.01 ~ "Supporting Benign",
    prob_disease <= 0.01 & prob_disease > 0.001 ~ "Moderate Benign",
    prob_disease <= 0.001 ~ "Strong Benign",
    TRUE ~ "Uncertain"
  )) %>%
  group_by(ACMG_Classification) %>%
  summarise(Count = n(), .groups = "drop")

# Generate the bar plot AFTER logistic regression with ACMG bins
ggplot(after_modeling, aes(x = ACMG_Classification, y = Count, fill = ACMG_Classification)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Reclassified Variants After Logistic Regression (ACMG Bins)", 
       x = "ACMG Classification", 
       y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for readahttp://127.0.0.1:47089/graphics/d114d61f-dc4d-4c78-bf2c-3ee6af1c238d.pngbility
