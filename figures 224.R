library(pROC)
library(ggplot2)

# Compute ROC curve
roc_obj <- roc(path_benign$pathogenic, predict(logistic_model, type="response"))

# Plot ROC curve
ggplot(data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities), aes(x=FPR, y=TPR)) +
  geom_line(color = "blue", size = 1.2) +
  geom_abline(linetype="dashed", color="grey") +
  labs(title = "ROC Curve for CADD PHRED Score",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  annotate("text", x=0.6, y=0.2, label=paste("AUC =", round(auc(roc_obj),3)), size=5, color="red") +
  theme_minimal()

# CADD PHRED Score Distribution by Pathogenicity

ggplot(path_benign, aes(x=factor(pathogenic, levels=c("-1", "1")), y=PHRED, fill=pathogenic)) +
  geom_boxplot() +
  scale_fill_manual(values = c("blue", "red"), labels = c("Benign", "Pathogenic")) +
  labs(title="Distribution of CADD PHRED Scores by Pathogenicity",
       x="Variant Classification",
       y="CADD PHRED Score") +
  theme_minimal()

# Log model curve: 

# Create dataset for plotting logistic regression curve
logit_curve_data <- data.frame(
  PHRED = seq(min(path_benign$PHRED), max(path_benign$PHRED), length.out = 100)
)

# Predict probability of pathogenicity
logit_curve_data$prob_pathogenic <- predict(logistic_model, newdata=logit_curve_data, type="response")

# Plot logistic regression curve
ggplot(path_benign, aes(x=PHRED, y=as.numeric(pathogenic))) +
  geom_point(alpha=0.2) +
  geom_line(data=logit_curve_data, aes(x=PHRED, y=prob_pathogenic), color="red", size=1.2) +
  labs(title="Logistic Regression Model of Pathogenicity Probability",
       x="CADD PHRED Score",
       y="Predicted Probability of Pathogenicity") +
  theme_minimal()

# Classification of ClinVar Variants Using Pathogenicity Strength Estimates

# Count the number of variants in each classification level
clinvar_class_counts <- final_combined_data_2 %>%
  filter(!is.na(CLNSIG)) %>%
  group_by(CLNSIG) %>%
  summarise(Count = n())

# Bar plot
ggplot(clinvar_class_counts, aes(x=reorder(CLNSIG, Count), y=Count, fill=CLNSIG)) +
  geom_bar(stat="identity") +
  coord_flip() +
  labs(title="Classification of ClinVar Variants by Pathogenicity Level",
       x="Pathogenicity Classification",
       y="Number of Variants") +
  theme_minimal() +
  scale_fill_viridis_d()

# Violin plot: 

ggplot(final_combined_data_2, aes(x=factor(CLNSIG), y=PHRED, fill=CLNSIG)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  labs(title="CADD Score Distribution by ACMG Classification",
       x="ACMG Classification",
       y="CADD PHRED Score") +
  theme_minimal() +
  scale_fill_viridis_d()

# density plot:

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
