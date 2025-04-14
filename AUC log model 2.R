# Load required libraries
library(dplyr)
library(pROC)

# Define estimated pathogenic proportion for weighting
prop_path_update <- 0.04  # Adjust as needed

# Run weighted logistic regression with imputation
tmp <- prob_imputation(filtered_data, prop_path = prop_path_update)

# Extract probabilistically imputed pathogenicity values
imputed_data <- tmp$plot_data  # This dataset contains weighted assignments

# Ensure data is formatted correctly for AUC calculation
imputed_data <- imputed_data %>%
  mutate(pathogenic = as.numeric(pathogenic),  # Convert factor to numeric
         pathogenic = ifelse(pathogenic == 0, 1, pathogenic))  # Adjust for ROC analysis

# Compute ROC curve on weighted data
roc_curve_weighted <- roc(imputed_data$pathogenic, imputed_data$PHRED)

# Calculate AUC
auc_weighted <- auc(roc_curve_weighted)

# Print the AUC value for the second model
cat("AUC for the weighted logistic regression model:", auc_weighted, "\n")

# Optionally, plot the ROC curve
plot(roc_curve_weighted, main = "ROC Curve - Weighted Logistic Model", col = "blue", lwd = 2)


######

# Ensure pathogenic column is binary (0 or 1)
imputed_data <- tmp$plot_data %>%
  mutate(pathogenic = as.numeric(pathogenic)) %>%  # Convert to numeric
  filter(!is.na(pathogenic) & pathogenic %in% c(0, 1))  # Remove NAs and ensure only 0 or 1

# Check unique values to confirm it's binary
unique_values <- unique(imputed_data$pathogenic)
cat("Unique values in pathogenic:", unique_values, "\n")

# Compute ROC curve only if pathogenic is binary
if (length(unique_values) == 2) {
  roc_curve_weighted <- roc(imputed_data$pathogenic, imputed_data$PHRED)
  auc_weighted <- auc(roc_curve_weighted)
  
  # Print AUC value
  cat("AUC for the weighted logistic regression model:", auc_weighted, "\n")
  
  # Plot ROC curve
  plot(roc_curve_weighted, main = "ROC Curve - Weighted Logistic Model", col = "blue", lwd = 2)
} else {
  cat("Error: Pathogenic column is not binary (contains values:", unique_values, ")\n")
}
