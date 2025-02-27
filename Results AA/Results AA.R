## Annotation vs Position Plot

# Plotting according to Chatgpt -- this is the one that works 
library(dplyr)
library(ggplot2)
library(purrr)  # Load purrr for map_dbl() and map_chr()


# Load gene info
gene_info <- read.table(file = "~/GitHub/ClinVar-GnomAD-Merge/pten_exon_postions", 
                        as.is = TRUE, header = FALSE, sep = "\t")
names(gene_info)[1:3] <- c("chr", "start", "end")
gene_info$exon <- 1:nrow(gene_info)

# Compute cumulative exon lengths
gene_info <- gene_info %>%
  mutate(exon_length = end - start + 1,
         cumulative_length = cumsum(exon_length) - exon_length)

# Function to compute c. position (using start and end)
find_cdot <- function(start_pos, end_pos) {
  pos <- ifelse(is.na(end_pos) | start_pos == end_pos, start_pos, floor((start_pos + end_pos) / 2))
  exon_found <- gene_info %>%
    filter(start <= pos & end >= pos) %>%
    mutate(pos_within_exon = pos - start + 1)
  
  if (nrow(exon_found) == 0) return(list(coding_pos = NA, exon = NA))
  coding_pos <- exon_found$cumulative_length + exon_found$pos_within_exon
  return(list(coding_pos = coding_pos, exon = exon_found$exon))
}

# Plot function for variants, exons, and pathogenicity
variant_plot <- function(data_all, prefix = "variant_plot") {
  if (!"start.x" %in% colnames(data_all)) stop("🚨 Missing 'start.x' column.")
  if (!"end.x" %in% colnames(data_all)) stop("🚨 Missing 'end.x' column.")
  if (!"PHRED" %in% colnames(data_all)) stop("🚨 Missing 'PHRED' column.")
  if (!"pathogenic" %in% colnames(data_all)) stop("🚨 Missing 'pathogenic' column.")  
  
  # Compute c. positions and map pathogenic values
  data_all <- data_all %>%
    mutate(start = as.numeric(start.x),
           end = as.numeric(end.x),
           pathogenic = factor(case_when(
             pathogenic == -1 ~ "Benign",
             pathogenic == 0 ~ "Uncertain",
             pathogenic == 1 ~ "Pathogenic",
             TRUE ~ "Unknown"  # 🛡️ Safety fallback
           ), levels = c("Benign", "Uncertain", "Pathogenic", "Unknown"))) %>%
    rowwise() %>%
    mutate(cdot_info = list(find_cdot(start, end))) %>%
    ungroup() %>%  
    mutate(c_position = map_dbl(cdot_info, ~ .x$coding_pos),
           exon = factor(map_chr(cdot_info, ~ as.character(.x$exon)))) %>%
    dplyr::select(-cdot_info)  # ✅ Now works without error
  
  # Plot with exon boundaries
  p <- ggplot(data_all, aes(x = c_position, y = PHRED, color = pathogenic)) +
    geom_point(size = 2, alpha = 0.5) +  # Decreased point size and increased transparency
    scale_color_manual(values = c("Benign" = "green", "Uncertain" = "blue", "Pathogenic" = "red", "Unknown" = "gray")) +
    geom_vline(xintercept = gene_info$cumulative_length + 1, color = "black", linetype = "dotted") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Variant Pathogenicity Across Coding Sequence",
      x = "c. Position (Coding DNA Position)",
      y = "PHRED Score",
      color = "Pathogenicity"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Title centered
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, face = "bold")
    )
  
  print(p)
  # Save plot
  ggsave(paste0(prefix, "_pathogenicity_plot.pdf"), p, width = 10, height = 6)
  cat("✅ Plot saved as:", paste0(prefix, "_pathogenicity_plot.pdf"), "\n")
}

# Usage Example
variant_plot(d_all, prefix = "final_variant_plot")






#### Visualization for table 1

# ROC Curve
library(pROC)
plot(roc_curve, col = "#1f78b4", main = paste("ROC Curve for PHRED Scores (AUC =", round(auc_value, 3), ")"))
abline(a = 0, b = 1, lty = 2, col = "gray")

# Density Plot
d <- ggplot(path_benign, aes(x = PHRED, fill = factor(pathogenic, labels = c("Benign", "Pathogenic")))) +
  geom_density(alpha = 0.7) +
  labs(title = "Density Plot of PHRED Scores by Variant Classification",
       x = "PHRED Score", y = "Density", fill = "Classification") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("#FFB6C1", "#ADD8E6")) +  # Apply pastel fill colors
  scale_color_manual(values = c("#E75480", "#4682B4")) +
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 13),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank()) 

# Save as a PDF
ggsave("density_plot.pdf", plot = d, device = "pdf", width = 10, height = 6)




# Boxplot/Violin Plot for Distribution
ggplot(path_benign, aes(x = factor(pathogenic, labels = c("Benign", "Pathogenic")), y = PHRED, fill = factor(pathogenic))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white") +
  labs(title = "PHRED Score Distribution by Variant Classification",
       x = "Variant Classification", y = "PHRED Score") +
  theme_minimal() +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +
  annotate("text", x = 1.5, y = max(path_benign$PHRED), label = paste("p-value:", signif(p_value, 3)))

# Bar plot for mean PHRED score
mean_data <- data.frame(
  Classification = c("Benign", "Pathogenic"),
  Mean_PHRED = c(benign_mean, pathogenic_mean)
)

ggplot(mean_data, aes(x = Classification, y = Mean_PHRED, fill = Classification)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = round(Mean_PHRED, 2)), vjust = -0.5) +
  labs(title = "Mean PHRED Scores by Variant Classification",
       y = "Mean PHRED Score") +
  theme_minimal() +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +
  annotate("text", x = 1.5, y = max(mean_data$Mean_PHRED) + 0.1, label = paste("p-value:", signif(p_value, 3)))






#### ACMG Threshold Visualization Plot (PHRED Score vs Probability)
library(ggplot2)

# Assuming threshold_table is already generated and contains PHRED thresholds
ggplot(d_all, aes(x = PHRED, y = prob_disease)) +
  geom_line(color = "blue", size = 1) +  # Probability curve
  geom_vline(data = threshold_table, aes(xintercept = as.numeric(pathogenic_supporting)), 
             linetype = "dashed", color = "orange", size = 1.2) +
  geom_vline(data = threshold_table, aes(xintercept = as.numeric(pathogenic_moderate)), 
             linetype = "dashed", color = "yellow", size = 1.2) +
  geom_vline(data = threshold_table, aes(xintercept = as.numeric(pathogenic_strong)), 
             linetype = "dashed", color = "red", size = 1.2) +
  geom_vline(data = threshold_table, aes(xintercept = as.numeric(pathogenic_very_strong)), 
             linetype = "dashed", color = "darkred", size = 1.2) +
  annotate("text", x = as.numeric(threshold_table$pathogenic_supporting), y = 0.1, label = "Supporting", angle = 90, vjust = -0.5) +
  annotate("text", x = as.numeric(threshold_table$pathogenic_moderate), y = 0.3, label = "Moderate", angle = 90, vjust = -0.5) +
  annotate("text", x = as.numeric(threshold_table$pathogenic_strong), y = 0.5, label = "Strong", angle = 90, vjust = -0.5) +
  annotate("text", x = as.numeric(threshold_table$pathogenic_very_strong), y = 0.7, label = "Very Strong", angle = 90, vjust = -0.5) +
  labs(title = "ACMG Threshold Visualization for PHRED Score",
       x = "PHRED Score",
       y = "Predicted Probability of Pathogenicity") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))






#### Distribution of Variants Before and Aftare Log Regression 
library(ggplot2)
library(dplyr)

library(ggplot2)
library(dplyr)
library(tidyr)

# Assuming 'filtered_data' and 'tmp$predicted_probabilities' are defined
filtered_data$predicted_prob <- d_all$prob_disease  # Replace with actual prediction column

# Remove non-finite values before plotting
variants_df <- filtered_data %>%
  dplyr::select(PHRED, predicted_prob) %>%
  filter(is.finite(PHRED), is.finite(predicted_prob)) %>%  # ✅ Filter out non-finite values
  pivot_longer(cols = c(PHRED, predicted_prob), names_to = "stage", values_to = "value")

# Plot without non-finite warnings
ggplot(variants_df, aes(x = value, fill = stage)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribution of Variants Before and After Logistic Regression",
    x = "Score / Predicted Probability",
    y = "Density"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Before Regression (PHRED)", "After Regression (Predicted Prob)")) +
  theme(legend.title = element_blank())






#### Distribution of weights
# Assuming 'tmp$weights' exists
weights_df <- data.frame(weights = d_all$weight)

ggplot(weights_df, aes(x = weights)) +
  geom_histogram(bins = 30, fill = "#619CFF", alpha = 0.7, color = "black") +
  labs(
    title = "Distribution of Weights Used in Logistic Regression",
    x = "Weight Value",
    y = "Frequency"
  ) +
  theme_minimal()

