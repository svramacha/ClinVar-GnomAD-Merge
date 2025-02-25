### Figure 1: Density plot and historgram 

# Load required libraries
library(ggplot2)
library(dplyr)

# Ensure PHRED scores are numeric and classify source
final_combined_data_2 <- final_combined_data_2 %>%
  filter(!is.na(PHRED)) %>%
  mutate(Source = ifelse(is.na(CLNSIG), "gnomAD", "ClinVar"))

# Define soft color scheme
fill_colors <- c("ClinVar" = "#FFB6C1", "gnomAD" = "#ADD8E6")  # Light Pink & Light Blue
line_colors <- c("ClinVar" = "#E75480", "gnomAD" = "#4682B4")  # Slightly darker pink & blue

# Create the enhanced histogram + density plot
ggplot(final_combined_data_2, aes(x = PHRED, fill = Source, color = Source)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 40, 
                 alpha = 0.5, 
                 position = "identity", 
                 color = "black",  # Outline bars in black for visibility
                 linewidth = 0.3) +  # Thicker bar outlines
  geom_density(alpha = 0.7, linewidth = 1.5) +  # Slightly thicker density lines
  scale_fill_manual(values = fill_colors) +  # Apply pastel fill colors
  scale_color_manual(values = line_colors) +  # Apply slightly darker density lines
  labs(title = "Distribution of CADD (PHRED) Scores: gnomAD vs ClinVar",
       subtitle = "Comparing PHRED score distributions for gnomAD and ClinVar variants",
       x = "CADD Score",
       y = "Relative Frequency",
       fill = "Variant Source",
       color = "Density") +
  theme_minimal(base_size = 14) +  # Apply clean theme with larger font
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())  # Remove minor grid lines for a cleaner look

# 



### Table 1: Mean pathogenicity prediction scores and areas under the receiver operating curve for Benign and Pathogenic variants  

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


#### Table 3: Estimated score threshold intervals for the four pathogenicity prediction tools evaluated in this study as they relate to the different pathogenic and benign strength levels 

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

#### Table 4: Number of predicted pathogenic and predicted benign variants at different strengths among the set of Uncertain variants 

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





