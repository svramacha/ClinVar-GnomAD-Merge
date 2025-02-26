# Load required libraries
library(ggplot2)
library(dplyr)

# Load required libraries
library(ggplot2)
library(dplyr)

# Define a nicer color palette
custom_colors <- c("Benign" = "#99CCFF",  # Light Blue
                   "Pathogenic" = "#FFB6C1",  # Light Pink
                   "Thresholds" = "#B2D8B2")  # Light Green for thresholds

# Ensure thresholds are numeric
benign_thresholds <- c(6.44, 12.5, 15.52, 17.01)
path_thresholds <- c(20.69, 22.06, 24.93)

# Create structured dataset for plotting
plot_data <- data.frame(
  Classification = factor(c(
    "Very Strong Benign", "Strong Benign", "Moderate Benign", "Supporting Benign",
    "Supporting Pathogenic", "Moderate Pathogenic", "Strong Pathogenic"
  ), levels = rev(c(  # Reverse order so that Strong Pathogenic appears on top
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

# Now define Label separately after Min and Max exist
plot_data$Label <- paste0("[", round(plot_data$Min, 2), ", ", round(plot_data$Max, 2), "]")

# Create the horizontal stacked range plot
ggplot(plot_data, aes(y = Classification, xmin = Min, xmax = Max, fill = Category)) +
  geom_rect(aes(xmin = Min, xmax = Max, ymin = as.numeric(Classification) - 0.4, ymax = as.numeric(Classification) + 0.4), 
            color = "black", alpha = 0.8) +
  geom_text(aes(x = (Min + Max) / 2, y = Classification, label = Label), 
            size = 3, vjust = 0, color = "black") +  # Small font size for range labels
  scale_fill_manual(values = custom_colors) +
  labs(title = "ACMG Threshold Classification Ranges", x = "PHRED Score", y = "Classification") +
  theme_minimal()

# Range plot 

# Load required libraries
library(ggplot2)
library(dplyr)

# Define a nice color palette
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

# Create the styled horizontal range plot
ggplot(plot_data, aes(y = Classification)) +
  # Add range bars (horizontal)
  geom_segment(aes(x = Min, xend = Max, y = Classification, yend = Classification),
               color = "#DDCDBA", size = 2) +  # Light beige bars
  # Add diverging dots at Min and Max
  geom_point(aes(x = Min, y = Classification, color = Category), size = 4) +
  geom_point(aes(x = Max, y = Classification, color = Category), size = 4) +
  # Customize colors
  scale_color_manual(values = custom_colors) +
  # Labels and theme
  labs(title = "ACMG Threshold Classification Ranges",
       x = "PHRED Score",
       y = "Classification") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())

# Load required libraries
library(ggplot2)
library(dplyr)

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

# Add label for ranges
plot_data$Label <- paste0("[", round(plot_data$Min, 2), ", ", round(plot_data$Max, 2), "]")

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

######

# Load required libraries
library(ggplot2)
library(dplyr)

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


#### 

# Load required libraries
library(ggplot2)
library(dplyr)

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
p <- ggplot(plot_data, aes(y = Classification)) +
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
  # Extend x-axis to prevent cutoff
  scale_x_continuous(limits = c(0, max(plot_data$Max, na.rm = TRUE) + 5)) +  
  # Labels and theme
  labs(title = "ACMG Threshold Classification Ranges",
       x = "CADD Score",
       y = "Classification") +
  theme_minimal() +
  theme(legend.position = "top", legend.title = element_blank())

# Print the plot
print(p)

# Save the plot with adjusted dimensions to avoid cutoff issues
ggsave("CADD_threshold_plot.png", plot = p, width = 10, height = 5, dpi = 300)

#####

# Load required libraries
library(ggplot2)
library(dplyr)

# Define custom color palette
custom_colors <- c("Benign" = "#99CCFF",  # Light Blue
                   "Pathogenic" = "#FFB6C1",  # Light Pink,
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
    path_thresholds[2], path_thresholds[3], 30  # Ensure last value isn't Inf
  ),
  Category = c(
    "Benign", "Benign", "Benign", "Benign", 
    "Pathogenic", "Pathogenic", "Pathogenic"
  )
)

# Format range labels correctly: ( Min , Max ]
plot_data$Label <- paste0("( ", round(plot_data$Min, 2), " , ", round(plot_data$Max, 2), " ]")

# Create the improved horizontal range plot
p <- ggplot(plot_data, aes(y = Classification)) +
  # Add range bars (fix missing bars)
  geom_segment(aes(x = Min, xend = Max, y = Classification, yend = Classification, color = Category),
               size = 4, lineend = "round") +  # Increase thickness for visibility
  # Add diverging dots at Min and Max
  geom_point(aes(x = Min, y = Classification, color = Category), size = 5) +
  geom_point(aes(x = Max, y = Classification, color = Category), size = 5) +
  # Add range labels properly aligned (shifted to avoid overlap)
  geom_text(aes(x = Max + 2, y = Classification, label = Label), 
            size = 5, hjust = 0, vjust = 0.5, color = "black", fontface = "bold") +  
  # Customize colors
  scale_color_manual(values = custom_colors) +
  # Extend x-axis to prevent cutoff
  scale_x_continuous(limits = c(0, max(plot_data$Max, na.rm = TRUE) + 5)) +  
  # Labels and theme adjustments
  labs(title = "ACMG Threshold Classification Ranges",
       x = "CADD Score",
       y = "Classification") +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),  # Use a clean font
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Center title, bold
    axis.title.x = element_text(size = 14, face = "bold"),  # Bold x-axis title
    axis.title.y = element_text(size = 14, face = "bold"),  # Bold y-axis title
    axis.text.y = element_text(size = 12),  # Increase classification text size
    legend.position = "top", 
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  )

# Print the plot
print(p)

# Save the improved plot with adjusted dimensions
ggsave("CADD_threshold_plot_fixed.png", plot = p, width = 12, height = 6, dpi = 300)

####

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

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)  # For statistical annotation

# Filter for Benign (-1) and Pathogenic (1) variants only
path_benign <- final_combined_data_2 %>%
  filter(pathogenic %in% c(-1, 1), !is.na(PHRED))  # Ensure no missing PHRED values

# Define color scheme
box_colors <- c("Benign" = "#ADD8E6", "Pathogenic" = "#FFB6C1")  # Light Blue & Light Pink

# Perform Wilcoxon test for significance
wilcox_test <- wilcox.test(
  path_benign$PHRED[path_benign$pathogenic == -1], 
  path_benign$PHRED[path_benign$pathogenic == 1], 
  alternative = "two.sided"
)
p_value <- wilcox_test$p.value  # Extract p-value

# Create a boxplot with jitter to show PHRED score distributions
ggplot(path_benign, aes(x = factor(pathogenic, labels = c("Benign", "Pathogenic")), 
                        y = PHRED, fill = factor(pathogenic, labels = c("Benign", "Pathogenic")))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "black") +  # Add jittered points
  scale_fill_manual(values = box_colors) +  # Apply color scheme
  labs(title = "Distribution of CADD (PHRED) Scores: Benign vs. Pathogenic",
       subtitle = paste("Wilcoxon test p-value:", signif(p_value, 3)),
       x = "Variant Type",
       y = "CADD (PHRED) Score",
       fill = "Variant Classification") +
  theme_minimal(base_size = 14) +  # Apply clean theme
  theme(legend.position = "none",  # Remove redundant legend
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())  # Remove minor grid lines

####

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)  # For statistical annotation

# Filter for Benign (-1) and Pathogenic (1) variants only
path_benign <- final_combined_data_2 %>%
  filter(pathogenic %in% c(-1, 1), !is.na(PHRED))  # Ensure no missing PHRED values

# Define color scheme
box_colors <- c("Benign" = "#ADD8E6", "Pathogenic" = "#FFB6C1")  # Light Blue & Light Pink

# Perform Wilcoxon test for significance
wilcox_test <- wilcox.test(
  path_benign$PHRED[path_benign$pathogenic == -1], 
  path_benign$PHRED[path_benign$pathogenic == 1], 
  alternative = "two.sided"
)
p_value <- wilcox_test$p.value  # Extract p-value

# Create a boxplot with jitter to show PHRED score distributions
ggplot(path_benign, aes(x = factor(pathogenic, labels = c("Benign", "Pathogenic")), 
                        y = PHRED, fill = factor(pathogenic, labels = c("Benign", "Pathogenic")))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Transparent boxplot without outliers
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "black") +  # Add jittered points
  scale_fill_manual(values = box_colors) +  # Apply color scheme
  labs(title = "Distribution of CADD (PHRED) Scores: Benign vs. Pathogenic",
       subtitle = paste("Wilcoxon test p-value:", signif(p_value, 3)),
       x = "Variant Type",
       y = "CADD (PHRED) Score",
       fill = "Variant Classification") +
  theme_minimal(base_size = 14) +  # Apply clean theme
  theme(legend.position = "none",  # Remove redundant legend
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_blank())  # Remove minor grid lines

