# Load required libraries
library(ggplot2)
library(dplyr)

# Create a data frame manually based on the table image
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 4)),
  Subcategory = c("Very Strong", "Strong", "Moderate", "Supporting",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong", "Very Strong"),
  Count = c(494, 102, 62, 63, 
            82, 
            98, 98, 362, 0)  # 0 replaces "—" for plotting
)

# Convert factors to maintain order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Very Strong", "Strong", "Moderate", "Supporting", "Indeterminate"))

# Define color scheme
subcat_colors <- c("Very Strong" = "#99CCFF",  # Light Blue
                   "Strong" = "#66B2FF",
                   "Moderate" = "#3399FF",
                   "Supporting" = "#ADD8E6",
                   "Indeterminate" = "#B2B2B2",  # Gray for Indeterminate
                   "Very Strong" = "#FF3333",  # Dark Red for Pathogenic
                   "Strong" = "#FF6666",
                   "Moderate" = "#FF9999",
                   "Supporting" = "#FFCCCC")

# Create the stacked bar plot
ggplot(classification_counts, aes(x = Category, y = Count, fill = Subcategory)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  scale_fill_manual(values = subcat_colors) +  # Apply color scheme
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Variant Count",
       fill = "Subcategory") +
  theme_minimal(base_size = 14) +  # Clean theme
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")  # Legend positioned to the side

###

# Load required libraries

# Load required libraries
library(ggplot2)
library(dplyr)

# Create a data frame manually based on the table image
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 4)),
  Subcategory = c("Very Strong", "Strong", "Moderate", "Supporting",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong", "Very Strong"),
  Count = c(494, 102, 62, 63, 
            82, 
            98, 98, 362, 0)  # 0 replaces "—" for plotting
)

# Convert factors to maintain order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Very Strong", "Strong", "Moderate", "Supporting", "Indeterminate"))

# Calculate percentages
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100)

# Define color scheme
subcat_colors <- c("Very Strong" = "#99CCFF",  # Light Blue (Benign)
                   "Strong" = "#66B2FF",
                   "Moderate" = "#3399FF",
                   "Supporting" = "#ADD8E6",
                   "Indeterminate" = "#B2B2B2",  # Gray for Indeterminate
                   "Very Strong" = "#FF3333",  # Dark Red (Pathogenic)
                   "Strong" = "#FF6666",
                   "Moderate" = "#FF9999",
                   "Supporting" = "#FFCCCC")

# Create the horizontal stacked bar plot with percentages
ggplot(classification_counts, aes(y = Category, x = Count, fill = Subcategory)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +  # "fill" makes it 100% stacked
  geom_text(aes(label = paste0(round(Percent, 1), "%")),  # Add percentage labels
            position = position_fill(vjust = 0.5), size = 5, color = "black") +
  scale_fill_manual(values = subcat_colors) +  # Apply color scheme
  scale_x_continuous(labels = scales::percent_format(scale = 1)) +  # Format x-axis as percentages
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Percentage of Variants",
       y = "Classification Category",
       fill = "Subcategory") +
  theme_minimal(base_size = 14) +  # Clean theme
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")  # Legend positioned to the side

#########

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

# Create the vertical 100% stacked bar plot with percentages inside
ggplot(classification_counts, aes(x = Category, y = Count, fill = Combined)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            position = position_fill(vjust = 0.5),
            size = 4, color = "black") +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  scale_fill_manual(values = custom_palette) +
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Percentage of Variants",
       fill = "Subcategory") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Optionally, save the plot
ggsave("Variant_Classification_StackedBar.png", width = 8, height = 6, dpi = 300)


#######

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
  "Benign (BP4) - Very Strong" = "#3399FF",  # Darkest Blue (strongest Benign)
  "Benign (BP4) - Strong"      = "#66B2FF",  # Medium Blue
  "Benign (BP4) - Moderate"    = "#99CCFF",  # Light Blue
  "Benign (BP4) - Supporting"  = "#ADD8E6",  # Lightest Blue (weakest Benign)
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  "Pathogenic (PP3) - Supporting" = "#FFCCCC",  # Lightest Red (weakest Pathogenic)
  "Pathogenic (PP3) - Moderate"   = "#FF9999",  # Medium-Light Red
  "Pathogenic (PP3) - Strong"     = "#FF6666",  # Medium-Dark Red
  "Pathogenic (PP3) - Very Strong"= "#FF3333"   # Darkest Red (strongest Pathogenic)
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
