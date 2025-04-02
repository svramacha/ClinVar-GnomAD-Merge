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
    scale_color_manual(values = c("Benign" = "#99CCFF", "Uncertain" = "#B2D8B2", "Pathogenic" = "#FFB6C1", "Unknown" = "gray")) +
    geom_vline(xintercept = gene_info$cumulative_length + 1, color = "black", linetype = "dotted") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Variant Pathogenicity Across Coding Sequence",
      x = "c. Position (Coding DNA Position)",
      y = "CADD Score",
      color = "Pathogenicity"
    ) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # Title centered
      strip.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 14, face = "bold")
    )
  
  print(p)
  # Save plot
  ggsave(paste0(prefix, "_pathogenicity_plot.png"), p, width = 10, height = 6, dpi = 300, device = "png")
  cat("✅ Plot saved as:", paste0(prefix, "_pathogenicity_plot.pdf"), "\n")
}

# Usage Example
variant_plot(d_all, prefix = "final_variant_plot")






#### Visualization for table 1

# ROC Curve
library(pROC)

# Open a PDF file to save the plot
pdf("roc_plot.pdf", width = 7, height = 5)

# Generate the ROC plot
plot(roc_curve, col = "#1f78b4", main = paste("ROC Curve for CADD Scores (AUC =", round(auc_value, 3), ")"))
abline(a = 0, b = 1, lty = 2, col = "gray")

# Close the PDF file
dev.off()


# Open a PNG file to save the plot
png("roc_plot.png", width = 7, height = 5, units = "in", res = 300)

# Generate the ROC plot
plot(roc_curve, col = "#1f78b4", main = paste("ROC Curve for CADD Scores (AUC =", round(auc_value, 3), ")"))
abline(a = 0, b = 1, lty = 2, col = "gray")

# Close the PNG file
dev.off()





# Density Plot
d <- ggplot(filtered_data, aes(x = PHRED, fill = factor(pathogenic, labels = c("Benign", "Pathogenic", "Uncertain")))) +
  geom_density(alpha = 0.7) +
  labs(title = "Density Plot of CADD Scores by Variant Classification",
       x = "CADD Score", y = "Relative Frequency", fill = "Classification") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Benign" = "#FFB6C1", "Pathogenic" = "#ADD8E6", "Uncertain" = "#B2D8B2")) +  # Apply specific colors for each classification
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 13),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())

# Save as a PDF
ggsave("density_plot.pdf", plot = d, device = "pdf", width = 10, height = 6)
ggsave("density_plot.png", plot = d, device = "png", width = 15, height = 8, dpi = 300)

# d_observed
d_ob <- ggplot(d_observed, aes(x = PHRED, fill = factor(pathogenic, labels = c("Benign", "Pathogenic")))) +
  geom_density(alpha = 0.7) +
  labs(title = "Density Plot of CADD Scores by Variant Classification",
       x = "CADD Score", y = "Relative Frequency", fill = "Classification") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c("Benign" = "#FFB6C1", "Pathogenic" = "#ADD8E6")) +  # Apply specific colors for each classification
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        plot.subtitle = element_text(size = 13),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())

# Boxplot/Violin Plot for Distribution

# Convert 'pathogenic' to a factor with meaningful labels
filtered_data$pathogenic <- factor(filtered_data$pathogenic, 
                                   levels = c(-1, 0, 1), 
                                   labels = c("Benign", "Uncertain", "Pathogenic"))

violin <- ggplot(filtered_data, aes(x = pathogenic, y = PHRED, fill = pathogenic)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot with transparency
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 1) +  # Boxplot overlay
  labs(title = "CADD Score Distribution by Variant Classification",
       x = "Variant Classification", y = "CADD Score") +
  theme_minimal() +
  scale_fill_manual(values = c("#FFB6C1", "#B2D8B2", "#ADD8E6")) +  # Colors for each category
  annotate("text", x = 2, y = max(filtered_data$PHRED, na.rm = TRUE) + 1, 
           label = paste("p-value:", signif(p_value, 3))) +  # Add p-value annotation
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5))

print(violin)

# Stack violin and denisty plots
library(gridExtra)
library(grid)

# Adjust margins to give more space (but not too much)
violin_adjusted <- violin + theme(plot.margin = margin(t = 50, r = 10, b = 10, l = 10))
d_adjusted      <- d + theme(plot.margin = margin(t = 50, r = 10, b = 10, l = 10))

# Create the labeled plots with enough space for the titles
violin_labeled <- arrangeGrob(violin_adjusted, top = textGrob("A", x = 0.05, y = 1.05, vjust = 1, gp = gpar(fontsize = 14, fontface = "bold")))
d_labeled      <- arrangeGrob(d_adjusted, top = textGrob("B", x = 0.05, y = 1.05, vjust = 1, gp = gpar(fontsize = 14, fontface = "bold")))

# Increase the size of the combined plot and adjust height/width for better visibility
combined_plot <- grid.arrange(
  violin_labeled, d_labeled,
  ncol = 1, 
  heights = unit(c(1.5, 1.5), "null")  # Adjust heights for bigger plots
)

# Show the combined plot
combined_plot




# Show the combined plot
print(combined_plot)
ggsave("combined_plot.png", combined_plot, width = 12, height = 10, dpi = 300)


# Save the combined plot as a PDF
ggsave("CADD_Score_Distribution.pdf", plot = combined_plot, width = 8, height = 10)
ggsave("CADD_Score_Distribution.png", plot = combined_plot, device = "png", width = 15, height = 10, dpi = 300)


# # Bar plot for mean PHRED score
# mean_data <- data.frame(
#   Classification = c("Benign", "Pathogenic"),
#   Mean_PHRED = c(benign_mean, pathogenic_mean)
# )
# 
# mean <- ggplot(mean_data, aes(x = Classification, y = Mean_PHRED, fill = Classification)) +
#   geom_col(width = 0.5) +
#   geom_text(aes(label = round(Mean_PHRED, 2)), vjust = -0.5) +
#   labs(title = "Mean CADD Scores by Variant Classification",
#        y = "Mean CADD Score") +
#   theme_minimal() +
#   scale_fill_manual(values = c("#FFB6C1", "#ADD8E6")) +
#   annotate("text", x = 1.5, y = max(mean_data$Mean_PHRED) + 0.1, label = paste("p-value:", signif(p_value, 3))) +
#   theme(legend.position = "none",
#         plot.title = element_text(hjust = 0.5))

#print(mean)

# Box plot to display means
# Assuming 'filtered_data' contains 'Classification' and 'PHRED' columns
# Convert 'pathogenic' to a factor with meaningful labels
filtered_data$pathogenic <- factor(filtered_data$pathogenic, 
                                   levels = c(-1, 0, 1), 
                                   labels = c("Benign", "Uncertain", "Pathogenic"))

# Box plot
mean <- ggplot(filtered_data, aes(x = pathogenic, y = PHRED, fill = pathogenic)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 2, alpha = 0.6) +  # Box plot with outliers
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +  # Add mean as a point
  labs(title = "CADD Score Distribution by Variant Classification",
       x = "Classification",
       y = "CADD Score") +
  theme_minimal() +
  scale_fill_manual(values = c("#FFB6C1", "#B2D8B2", "#ADD8E6")) +  # Colors for Benign, Uncertain, Pathogenic
  annotate("text", x = 2, y = max(filtered_data$PHRED, na.rm = TRUE) + 1, 
           label = paste("p-value:", signif(p_value, 3))) +  # Add p-value annotation
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

print(mean)



# Show violin, bar, and density plot side by side
library(ggplot2)
library(gridExtra)
library(grid)

# Adjust margins to prevent cutting off titles
violin_adjusted <- violin + theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10))
mean_adjusted   <- mean + theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10))
d_adjusted      <- d + theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))

# Create labels with enough space
violin_labeled <- arrangeGrob(violin_adjusted, top = textGrob("A", x = 0.05, y = 1.1, vjust = 1.5, gp = gpar(fontsize = 14, fontface = "bold")))
mean_labeled   <- arrangeGrob(mean_adjusted, top = textGrob("B", x = 0.05, y = 1.1, vjust = 1.5, gp = gpar(fontsize = 14, fontface = "bold")))
d_labeled      <- arrangeGrob(d_adjusted, top = textGrob("C", x = 0.05, y = 1.1, vjust = 1.5, gp = gpar(fontsize = 14, fontface = "bold")))

# Arrange the plots (2 on top, 1 on bottom)
combined_plot <- arrangeGrob(
  violin_labeled, mean_labeled, d_labeled,
  layout_matrix = rbind(c(1, 2), c(3, 3)), 
  heights = c(1, 1.3)  # Give more space to bottom plot
)

# Save as PDF with increased height to prevent squishing
ggsave("combined_plots.pdf", plot = combined_plot, width = 10, height = 10, dpi = 300)
ggsave("combined_plots.png", plot = combined_plot, device = "png", width = 15, height = 8, dpi = 300)




#### ACMG Threshold Visualization Plot (PHRED Score vs Probability)
library(ggplot2)

# Assuming threshold_table is already generated and contains PHRED thresholds
c <- ggplot(d_all, aes(x = PHRED, y = prob_disease)) +
  geom_line(color = "black", size = 1) +  # Probability curve
  geom_vline(data = threshold_table, aes(xintercept = as.numeric(pathogenic_supporting)), 
             linetype = "dashed", color = "#99CCFF", size = 1.2) +
  geom_vline(data = threshold_table, aes(xintercept = as.numeric(pathogenic_moderate)), 
             linetype = "dashed", color = "#FFB6C1", size = 1.2) +
  geom_vline(data = threshold_table, aes(xintercept = as.numeric(pathogenic_strong)), 
             linetype = "dashed", color = "#B2D8B2", size = 1.2) +
  geom_vline(data = threshold_table, aes(xintercept = as.numeric(pathogenic_very_strong)), 
             linetype = "dashed", color = "#D7BDE2" , size = 1.2) +
  annotate("text", x = as.numeric(threshold_table$pathogenic_supporting), y = 0.1, label = "Supporting", angle = 90, vjust = -0.5) +
  annotate("text", x = as.numeric(threshold_table$pathogenic_moderate), y = 0.3, label = "Moderate", angle = 90, vjust = -0.5) +
  annotate("text", x = as.numeric(threshold_table$pathogenic_strong), y = 0.5, label = "Strong", angle = 90, vjust = -0.5) +
  annotate("text", x = as.numeric(threshold_table$pathogenic_very_strong), y = 0.7, label = "Very Strong", angle = 90, vjust = -0.5) +
  labs(title = "Classification Threshold Visualization for CADD Score",
       x = "CADD Score",
       y = "Pathogenicity Probability") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

print(c)

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
a <- ggplot(plot_data, aes(y = Classification)) +
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
  labs(title = "Classification Threshold Ranges",
       x = "CADD Score",  
       y = "Classification") +
  theme_minimal() +
  theme(
    legend.position = "top", 
    legend.title = element_blank(),
    axis.title.x = element_text(size = 14),  
    axis.title.y = element_text(size = 14),  
    axis.text.y = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  


# Show ACMG plots side by side
library(gridExtra)
library(grid)


# Adjust margins to prevent cutting off titles
a_adjusted <- a + theme(plot.margin = margin(t = 30, r = 10, b = 10, l = 10))
c_adjusted <- c + theme(plot.margin = margin(t = 20, r = 10, b = 10, l = 10))

# Create labels with enough space
a_labeled <- arrangeGrob(a_adjusted, top = textGrob("A", x = 0.05, y = 1.1, vjust = 1.5, gp = gpar(fontsize = 14, fontface = "bold")))
c_labeled <- arrangeGrob(c_adjusted, top = textGrob("B", x = 0.05, y = 1.1, vjust = 1.5, gp = gpar(fontsize = 14, fontface = "bold")))


ac <- grid.arrange(a_labeled, c_labeled, ncol = 2)

ggsave("acmg_plot.pdf", plot = ac, device = "pdf", width = 15, height = 8)
ggsave("acmg_plot.png", plot = ac, device = "png", width = 15, height = 8, dpi = 300)








#### Distribution of Variants Before and After Log Regression 
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
  geom_density(aes(y = ..scaled..), alpha = 0.5) +  # Scale density to sum to 1
  labs(
    title = "Proportional Distribution of Variants Before and After Logistic Regression",
    x = "Score / Predicted Probability",
    y = "Proportion"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("#00BFC4", "#F8766D"), labels = c("Before Regression (PHRED)", "After Regression (Predicted Prob)")) +
  theme(legend.title = element_blank())




#### Distribution of weights
# Filter out weights equal to 0 or 1 and remove NAs
df_weights <- data.frame(weights = na.omit(d_all$weight[!(d_all$weight %in% c(0, 1))]))

# Create histogram of weights
weight_plot <- ggplot(df_weights, aes(x = weights)) +
  geom_histogram(binwidth = 0.05, fill = "#99CCFF", alpha = 0.7, color = "black") +
  labs(
    title = "Distribution of Weights Before Logistic Regression",
    x = "Weight Value",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

# Create histogram of weights for d_imp (d_imp is d_all so its the same)

library(dplyr)

# Group by pathogenicity and summarize weight statistics
d_all %>%
  group_by(pathogenic) %>%
  summarize(
    min_weight = min(weight, na.rm = TRUE),
    max_weight = max(weight, na.rm = TRUE),
    avg_weight = mean(weight, na.rm = TRUE),
    count = n() # this is showing that some uncertain variants have a weight of 0 
  )

d_all %>%
  filter(pathogenic == 0, weight == 1) %>%
  nrow() # showing that there are not uncertain variants with weight = 0

d_all %>%
  filter(pathogenic == 0, abs(weight - 1) < 1e-6) %>%
  nrow()



weight_imp <- ggplot(d_imp, aes(x = weight)) +
  geom_histogram(binwidth = 0.05, fill = "#99CCFF", alpha = 0.7, color = "black") +
  labs(
    title = "Distribution of Weights Before Logistic Regression",
    x = "Weight Value",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))


weight_bar_plot <- ggplot(d_observed, aes(x = as.factor(pathogenic), y = w, fill = as.factor(pathogenic))) +
  geom_bar(stat = "identity", position = "dodge", show.legend = FALSE) +
  scale_fill_manual(values = c("red", "blue")) + # Custom colors for the categories
  labs(
    title = "Bar Plot of Weights by Pathogenicity",
    x = "Pathogenic Category",
    y = "Weight Value"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

ggplot(d_all, aes(x = prob_disease)) +
  geom_histogram(binwidth = 0.05, fill = "#99CCFF", alpha = 0.7, color = "black") +
  labs(
    title = "Distribution of Predicted Probabilities for VUS",
    x = "Predicted Probability",
    y = "Count"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))


summary_table <- d_observed %>%
  group_by(pathogenic) %>%
  summarise(
    count = n(),                       # Frequency count
    min_weight = min(w, na.rm = TRUE),  # Minimum weight
    max_weight = max(w, na.rm = TRUE),  # Maximum weight
    avg_weight = mean(w, na.rm = TRUE)  # Average weight
  )

print(summary_table)

# Save the plot as PDF and PNG
ggsave("weight_distribution_plot.pdf", plot = weight_plot, device = "pdf", width = 15, height = 8)
ggsave("weight_distribution_plot.png", plot = weight_plot, device = "png", width = 15, height = 8, dpi = 300)
