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
  labs(title = "Distribution of CADD Scores: gnomAD vs ClinVar",
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

### STEP 4: Compute ACMG Thresholds with Corrected Ranges
calculate_thresholds <- function(est_perm, filtered_data, prop_path, file="acmg_table_impute.txt") {
  
  rslts <- c()
  cutoffs <- c(2.406, 5.790, 33.53, 1124)  # ACMG cutoffs
  cutlabel <- c("supporting", "moderate", "strong", "very_strong") 
  pD <- prop_path
  odds_path <- pD / (1 - pD)
  
  min_func <- function(x, beta, v, cutoff, odds_path, var_type) {
    if (!var_type %in% c("path", "benign")) {
      stop("var_type must be 'path' or 'benign'")
    }
    
    var_adj <- ifelse(var_type == "path", 1, -1)
    se_log_odds <- sqrt(t(c(1, x)) %*% v %*% c(1, x)) 
    z <- qnorm(.95) 
    logodds_bound <- (beta[1] + beta[2] * x) - var_adj * z * se_log_odds
    odds_bound <- exp(var_adj * logodds_bound)
    odds_path <- odds_path^var_adj
    
    return(cutoff - odds_bound / odds_path)
  }
  
  for (annot in c("PHRED")) {
    
    betas <- est_perm[[annot]]$params
    v <- est_perm[[annot]]$cov
    int <- betas[1]   
    logor <- betas[2]  
    p <- est_perm[[annot]]$p[2] 
    se_int <- sqrt(v[1, 1])
    se_logor <- sqrt(v[2, 2])
    
    rng <- range(filtered_data[[annot]], na.rm = TRUE)  
    lower <- min(rng) 
    upper <- max(rng)  
    
    path_thresholds <- c()
    benign_thresholds <- c()
    
    for (j in seq_along(cutoffs)) {
      cutoff <- cutoffs[j]
      
      # Compute pathogenic thresholds
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta = betas, v = v,  
                cutoff = cutoff, odds_path = odds_path, var_type = "path")
      }, error = function(e) list(root = NA))
      
      path_thresholds[j] <- result$root
      
      # Compute benign thresholds
      result <- tryCatch({
        uniroot(min_func, lower = lower, upper = upper, beta = betas, v = v,
                cutoff = cutoff, odds_path = odds_path, var_type = "benign")
      }, error = function(e) list(root = NA))
      
      benign_thresholds[j] <- result$root
    }
    
    # Ensure benign thresholds are sorted in increasing order
    benign_thresholds <- sort(benign_thresholds)
    
    # Ensure pathogenic thresholds are sorted in increasing order
    path_thresholds <- sort(path_thresholds)
    
    # Format into ranges
    phred_threshold_table <- data.frame(
      Tool = "PHRED",
      `Very Strong Benign` = paste("≤", round(benign_thresholds[1], 2)),
      `Strong Benign` = paste("(", round(benign_thresholds[1], 2), ",", round(benign_thresholds[2], 2), "]"),
      `Moderate Benign` = paste("(", round(benign_thresholds[2], 2), ",", round(benign_thresholds[3], 2), "]"),
      `Supporting Benign` = paste("(", round(benign_thresholds[3], 2), ",", round(benign_thresholds[4], 2), "]"),
      `Supporting Pathogenic` = paste("[", round(path_thresholds[1], 2), ",", round(path_thresholds[2], 2), ")"),
      `Moderate Pathogenic` = paste("[", round(path_thresholds[2], 2), ",", round(path_thresholds[3], 2), ")"),
      `Strong Pathogenic` = paste("≥", round(path_thresholds[3], 2)),
      `Very Strong Pathogenic` = "—"
    )
    
    # Print and save table
    print(phred_threshold_table)
    write.table(file = file, phred_threshold_table, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }
  
  cat("Wrote ACMG threshold table to", file, "\n")
  
  return(phred_threshold_table)
}

### STEP 5: Run ACMG Threshold Calculation
threshold_table <- calculate_thresholds(est_perm, filtered_data, prop_path_update)

print(threshold_table)


### Threshold visual: 

# Load required libraries
library(ggplot2)
library(dplyr)

# Ensure thresholds are numeric
benign_thresholds <- c(6.44, 12.5, 15.52, 17.01)
path_thresholds <- c(20.69, 22.06, 24.93)

# Define custom color palette
custom_colors <- c("Benign" = "#99CCFF",  # Light Blue
                   "Pathogenic" = "#FFB6C1",  # Light Pink,
                   "Thresholds" = "#B2D8B2")  # Light Green for thresholds

# Create a structured dataset
plot_data <- data.frame(
  Classification = factor(c(
    "Very Strong Benign", "Strong Benign", "Moderate Benign", "Supporting Benign",
    "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic"
  ), levels = rev(c(  # Reverse order so that stronger pathogenicity appears on top
    "Very Strong Benign", "Strong Benign", "Moderate Benign", "Supporting Benign",
    "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic"
  ))),
  Min = c(
    0, benign_thresholds[1], benign_thresholds[2], benign_thresholds[3], 
    path_thresholds[1], path_thresholds[2], path_thresholds[3]
  ),
  Max = c(
    benign_thresholds[1], benign_thresholds[2], benign_thresholds[3], benign_thresholds[4], 
    path_thresholds[2], path_thresholds[3], 30  # Ensuring the last value isn't Inf
  ),
  Category = c(
    "Benign", "Benign", "Benign", "Benign", 
    "Pathogenic", "Pathogenic", "Pathogenic"
  )
)

# Format range labels correctly: ( Min , Max ]
plot_data$Label <- paste0("( ", round(plot_data$Min, 2), " , ", round(plot_data$Max, 2), " ]")

# Create the styled horizontal range plot
ggplot(plot_data, aes(y = Classification)) +
  # Add range bars (same color as dots)
  geom_segment(aes(x = Min, xend = Max, y = Classification, yend = Classification, color = Category),
               size = 3) +  # Increased size for better visibility
  # Add diverging dots at Min and Max
  geom_point(aes(x = Min, y = Classification, color = Category), size = 5) +
  geom_point(aes(x = Max, y = Classification, color = Category), size = 5) +
  # Add range labels next to the bars (slightly to the right of max value)
  geom_text(aes(x = Max + 1, y = Classification, label = Label), 
            size = 4, hjust = 0, color = "black") +  
  # Customize colors
  scale_color_manual(values = custom_colors) +
  # Extend x-axis limit to 35
  scale_x_continuous(limits = c(0, 35)) +  
  # Labels and theme
  labs(title = "ACMG Threshold Classification Ranges",
       x = "CADD Score",  
       y = "Classification") +
  theme_minimal() +
  theme(
    legend.position = "top", 
    legend.title = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),  
    axis.title.y = element_text(size = 14, face = "bold"),  
    axis.text.y = element_text(size = 12)
  )

####

### dont use this one !!!! 

# Load required libraries
library(ggplot2)
library(dplyr)

# Ensure thresholds are numeric
benign_thresholds <- c(6.44, 12.5, 15.52, 17.01)
path_thresholds <- c(20.69, 22.06, 24.93)

# Define custom color palette
custom_colors <- c("Benign" = "#99CCFF",  # Light Blue
                   "Pathogenic" = "#FFB6C1",  # Light Pink
                   "Thresholds" = "#B2D8B2")  # Light Green for thresholds

# Ensure thresholds are numeric
benign_thresholds <- c(6.44, 12.5, 15.52, 17.01)
path_thresholds <- c(20.69, 22.06, 24.93)

# Create a structured dataset
plot_data <- data.frame(
  Classification = factor(c(
    "Very Strong Benign", "Strong Benign", "Moderate Benign", "Supporting Benign",
    "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic"
  ), levels = rev(c(  # Reverse order so that stronger pathogenicity appears on top
    "Very Strong Benign", "Strong Benign", "Moderate Benign", "Supporting Benign",
    "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic"
  ))),
  Min = c(
    0, benign_thresholds[1], benign_thresholds[2], benign_thresholds[3], 
    path_thresholds[1], path_thresholds[2], path_thresholds[3]
  ),
  Max = c(
    benign_thresholds[1], benign_thresholds[2], benign_thresholds[3], benign_thresholds[4], 
    path_thresholds[2], path_thresholds[3], Inf
  ),
  Category = c(
    "Benign", "Benign", "Benign", "Benign", 
    "Pathogenic", "Pathogenic", "Pathogenic"
  )
)

# Format range labels correctly: ( Min , Max ]
plot_data$Label <- paste0("( ", round(plot_data$Min, 2), " , ", round(plot_data$Max, 2), " ]")

# Create the styled horizontal range plot
ggplot(plot_data, aes(y = Classification)) +
  # Add range bars (same color as dots)
  geom_segment(aes(x = Min, xend = Max, y = Classification, yend = Classification, color = Category),
               size = 2) +  
  # Add diverging dots at Min and Max
  geom_point(aes(x = Min, y = Classification, color = Category), size = 4) +
  geom_point(aes(x = Max, y = Classification, color = Category), size = 4) +
  # Add range labels next to the bars (slightly to the right of max value)
  geom_text(aes(x = Max + 2, y = Classification, label = Label), 
            size = 3, hjust = 0, color = "black") +  
  # Customize colors
  scale_color_manual(values = custom_colors) +
  # Labels and theme
  labs(title = "ACMG Threshold Classification Ranges",
       x = "CADD Score",  # Changed from PHRED to CADD
       y = "Classification") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())


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

### Stacked Bar Plot for table 4: 

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For number formatting

# Create a data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 4)),
  Subcategory = c("Very Strong", "Strong", "Moderate", "Supporting",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong", "Very Strong"),
  Count = c(494, 102, 62, 63, 
            82, 
            98, 98, 362, 0)  # Use 0 instead of "—"
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# **Explicitly fix the stacking order**: Supporting → Moderate → Strong → Very Strong
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Supporting", "Moderate", "Strong", "Very Strong", "Indeterminate"),
                                            ordered = TRUE)  

# **Combine Category and Subcategory to use in colors and labels**
classification_counts <- classification_counts %>%
  mutate(Combined = paste(Category, Subcategory, sep = " - "))

# **Define Correct Color Mapping (Benign = Blue, Pathogenic = Red, Indeterminate = Gray)**
custom_palette <- c(
  # **Benign (BP4)**
  "Benign (BP4) - Supporting"  = "#ADD8E6",  # Light Blue
  "Benign (BP4) - Moderate"    = "#66B2FF",  # Medium Blue
  "Benign (BP4) - Strong"      = "#3399FF",  # Darker Blue
  "Benign (BP4) - Very Strong" = "#003366",  # Darkest Blue
  
  # **Indeterminate**
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  
  # **Pathogenic (PP3)**
  "Pathogenic (PP3) - Supporting" = "#FFCCCC",  # Lightest Red
  "Pathogenic (PP3) - Moderate"   = "#FF9999",  # Medium Red
  "Pathogenic (PP3) - Strong"     = "#FF6666",  # Darker Red
  "Pathogenic (PP3) - Very Strong"= "#FF3333"   # Darkest Red
)

# **Sort the dataset to ensure correct stacking order**
classification_counts <- classification_counts %>%
  arrange(Category, Subcategory)  # This ensures ggplot2 stacks bars correctly

# **Create the vertical stacked bar plot with only counts inside the bars**
ggplot(classification_counts, aes(x = Category, y = Count, fill = factor(Combined, levels = rev(names(custom_palette))))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  geom_text(aes(label = ifelse(Count > 0, Count, "")),  # Label inside bars (only count)
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_y_continuous(labels = comma) +  # Use normal numeric format for y-axis
  scale_fill_manual(values = custom_palette, name = "Subcategory",
                    breaks = rev(names(custom_palette)), guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Count of Variants") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# **Save the corrected plot**
ggsave("Variant_Classification_StackedBar_CountsOnly.png", width = 8, height = 6, dpi = 300)





#### Different color scheme 

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For number formatting

# Create a data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 4)),
  Subcategory = c("Very Strong", "Strong", "Moderate", "Supporting",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong", "Very Strong"),
  Count = c(494, 102, 62, 63, 
            82, 
            98, 98, 362, 0)  # Use 0 instead of "—"
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# **Explicitly fix the stacking order**: Supporting → Moderate → Strong → Very Strong
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Supporting", "Moderate", "Strong", "Very Strong", "Indeterminate"),
                                            ordered = TRUE)  

# **Combine Category and Subcategory to use in colors and labels**
classification_counts <- classification_counts %>%
  mutate(Combined = paste(Category, Subcategory, sep = " - "))

custom_palette <- c(
  # **Benign (BP4) - Shades of Blue**
  "Benign (BP4) - Supporting"  = "#ADD8E6",  # Light Blue
  "Benign (BP4) - Moderate"    = "#66B2FF",  # Medium Blue
  "Benign (BP4) - Strong"      = "#3399FF",  # Darker Blue
  "Benign (BP4) - Very Strong" = "#003366",  # Darkest Blue
  
  # **Indeterminate - Neutral Gray**
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  
  # **Pathogenic (PP3) - Now Using Shades of Pink**
  "Pathogenic (PP3) - Supporting" = "#FFCCE5",  # Lightest Pink
  "Pathogenic (PP3) - Moderate"   = "#FF99CC",  # Medium Pink
  "Pathogenic (PP3) - Strong"     = "#FF6699",  # Darker Pink
  "Pathogenic (PP3) - Very Strong"= "#FF3385"   # Darkest Pink
)

# **Sort the dataset to ensure correct stacking order**
classification_counts <- classification_counts %>%
  arrange(Category, Subcategory)  # This ensures ggplot2 stacks bars correctly

# **Create the vertical stacked bar plot with only counts inside the bars**
ggplot(classification_counts, aes(x = Category, y = Count, fill = factor(Combined, levels = rev(names(custom_palette))))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  geom_text(aes(label = ifelse(Count > 0, Count, "")),  # Label inside bars (only count)
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  scale_y_continuous(labels = comma) +  # Use normal numeric format for y-axis
  scale_fill_manual(values = custom_palette, name = "Subcategory",
                    breaks = rev(names(custom_palette)), guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Count of Variants") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")


########
# Dont use for now 
######## 

# the other threshold table
# dont look at this

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

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For percent_format

# Create a data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 4)),
  Subcategory = c("Very Strong", "Strong", "Moderate", "Supporting",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong", "Very Strong"),
  Count = c(494, 102, 62, 63, 
            82, 
            98, 98, 362, 0)  # Use 0 in place of "—"
)

# Convert Category and Subcategory to factors with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# For the purposes of coloring, combine Category and Subcategory so that identical subcategory names in different main categories are distinct.
classification_counts <- classification_counts %>%
  mutate(Combined = ifelse(Category == "Indeterminate",
                           paste(Category, Subcategory, sep = " - "),
                           paste(Category, Subcategory, sep = " - ")))

# Calculate percentages within each Category
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# Define custom color palette for each Combined level
# (Benign: light blue shades, Indeterminate: gray, Pathogenic: light to dark red shades)
custom_palette <- c(
  "Benign (BP4) - Very Strong" = "#99CCFF",
  "Benign (BP4) - Strong"      = "#66B2FF",
  "Benign (BP4) - Moderate"    = "#3399FF",
  "Benign (BP4) - Supporting"  = "#ADD8E6",
  "Indeterminate - Indeterminate" = "#B2B2B2",
  "Pathogenic (PP3) - Supporting" = "#FFCCCC",
  "Pathogenic (PP3) - Moderate"   = "#FF9999",
  "Pathogenic (PP3) - Strong"     = "#FF6666",
  "Pathogenic (PP3) - Very Strong"= "#FF3333"
)

# Create the vertical stacked bar plot with counts & percentages
ggplot(classification_counts, aes(x = Category, y = Count, fill = Combined)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bar with counts
  geom_text(aes(label = paste0(Count, " (", round(Percent, 1), "%)")),  # Show both counts & percentages
            position = position_stack(vjust = 0.5), size = 4, color = "black") +  
  scale_y_continuous(labels = comma) +  # Use normal numeric format for y-axis
  scale_fill_manual(values = custom_palette) +
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Count of Variants",
       fill = "Subcategory") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Optionally, save the plot
ggsave("Variant_Classification_StackedBar.png", width = 8, height = 6, dpi = 300)

