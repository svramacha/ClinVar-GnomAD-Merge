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

####

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

# Create a corrected data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 3)),  # Removed PP3 - Very Strong (0 count)
  Subcategory = c("Very Strong", "Strong", "Moderate", "Supporting",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong"),
  Count = c(494, 102, 62, 63, 
            82, 
            98, 98, 362)  # Removed 0 count for "Very Strong" Pathogenic
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# Explicitly order Subcategories so "Very Strong" is at the bottom
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Supporting", "Moderate", "Strong", "Very Strong", "Indeterminate"),
                                            ordered = TRUE)  # Ensure ordered stacking

# Calculate percentages within each Category
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# Define improved color scheme with stronger shades for stronger classifications
custom_palette <- c(
  "Benign (BP4) - Very Strong" = "#0057B7",  # Dark Blue (Strongest Benign)
  "Benign (BP4) - Strong"      = "#0073E6",  # Medium-Strong Blue
  "Benign (BP4) - Moderate"    = "#66B2FF",  # Medium-Light Blue
  "Benign (BP4) - Supporting"  = "#99CCFF",  # Lightest Blue
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  "Pathogenic (PP3) - Supporting" = "#FFCCCC",  # Lightest Red
  "Pathogenic (PP3) - Moderate"   = "#FF6666",  # Medium Red
  "Pathogenic (PP3) - Strong"     = "#CC0000"   # Dark Red (Strongest Pathogenic)
)

# Ensure the correct order of stacking using Subcategory as fill
ggplot(classification_counts, aes(x = Category, y = Count, fill = Subcategory)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars with correct ordering
  geom_text(aes(label = ifelse(Percent > 5, paste0(Count, " (", round(Percent, 1), "%)"), "")),  # Hide small labels
            position = position_stack(vjust = 0.5), size = 4, color = "black") +  
  scale_y_continuous(labels = comma) +  # Use numeric format for y-axis
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

# Save the corrected plot
ggsave("Variant_Classification_StackedBar_Corrected.png", width = 8, height = 6, dpi = 300)


# Load required libraries
library(ggplot2)
library(dplyr)

# Define a logical order for classification categories
phred_uncertain_df$Category <- factor(phred_uncertain_df$Category, 
                                      levels = c("Very Strong Benign", "Strong Benign", "Moderate Benign",
                                                 "Indeterminate", "Supporting Pathogenic", "Moderate Pathogenic",
                                                 "Strong Pathogenic"))

# Define color scheme
category_colors <- c("Very Strong Benign" = "#99CCFF",  # Light Blue
                     "Strong Benign" = "#66B2FF",
                     "Moderate Benign" = "#3399FF",
                     "Indeterminate" = "#B2B2B2",  # Gray for uncertainty
                     "Supporting Pathogenic" = "#FF9999",
                     "Moderate Pathogenic" = "#FF6666",
                     "Strong Pathogenic" = "#FF3333")  # Dark Red

# Create stacked bar plot
ggplot(phred_uncertain_df, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +  # Stacked bar
  scale_fill_manual(values = category_colors) +  # Apply color scheme
  labs(title = "Distribution of Uncertain Variants by Classification",
       subtitle = "Predicted pathogenicity based on logistic regression",
       x = NULL, y = "Count",
       fill = "Classification") +
  theme_minimal(base_size = 14) +  # Clean theme
  theme(axis.text.x = element_blank(),  # Hide x-axis label since it's a single bar
        axis.ticks.x = element_blank(),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(size = 13),
        legend.position = "right")  # Move legend to the side for clarity


###

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

# Create a corrected data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 3)),  
  Subcategory = c("Supporting", "Moderate", "Strong", "Very Strong",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong"),
  Count = c(63, 62, 102, 494,  # Benign counts (Supporting smallest, Very Strong largest)
            82, 
            98, 98, 362)  # Pathogenic counts (Supporting smallest, Strong largest)
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# **Reverse the Subcategory order so "Supporting" is at the bottom and "Very Strong" is at the top**
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = rev(c("Supporting", "Moderate", "Strong", "Very Strong", "Indeterminate")),
                                            ordered = TRUE)  

# Calculate percentages within each Category
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# Define correct color mapping for Benign (Blue), Indeterminate (Gray), and Pathogenic (Red)
custom_palette <- c(
  "Supporting"  = "#99CCFF",  # Lightest Blue (smallest count)
  "Moderate"    = "#66B2FF",  # Medium-Light Blue
  "Strong"      = "#0057B7",  # Medium-Strong Blue
  "Very Strong" = "#003366",  # Darkest Blue (largest count)
  "Indeterminate" = "#B2B2B2",  # Neutral Gray
  "Supporting"  = "#FFCCCC",  # Lightest Red (smallest count)
  "Moderate"    = "#FF6666",  # Medium Red
  "Strong"      = "#CC0000"   # Dark Red (largest count)
)

# Ensure the correct order of stacking using Subcategory
ggplot(classification_counts, aes(x = Category, y = Count, fill = Subcategory)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars with correct ordering
  scale_fill_manual(values = custom_palette, guide = guide_legend(reverse = TRUE)) +  # Reverse legend to match stacking
  scale_y_continuous(labels = comma) +  # Use numeric format for y-axis
  geom_text(aes(label = ifelse(Percent > 5, paste0(Count, " (", round(Percent, 1), "%)"), "")),  # Hide small labels
            position = position_stack(vjust = 0.5), size = 4, color = "black") +  
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Variant Count",
       fill = "Subcategory") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Save the corrected plot
ggsave("Variant_Classification_StackedBar_Corrected.png", width = 8, height = 6, dpi = 300)

######

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

# Create a corrected data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 3)),  
  Subcategory = c("Supporting", "Moderate", "Strong", "Very Strong",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong"),
  Count = c(63, 62, 102, 494,  # Benign counts
            82, 
            98, 98, 362)  # Pathogenic counts
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# Ensure Subcategory levels are in the correct stacking order (smallest at bottom)
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = rev(c("Supporting", "Moderate", "Strong", "Very Strong", "Indeterminate")),
                                            ordered = TRUE)

# Calculate percentages within each Category
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# **Create a proper fill label combining Category and Subcategory for correct color mapping**
classification_counts$Fill_Label <- paste(classification_counts$Category, classification_counts$Subcategory, sep = " - ")

# Define **correct color mapping** for both Benign (Blue) and Pathogenic (Red)
custom_palette <- c(
  "Benign (BP4) - Supporting"  = "#99CCFF",  # Lightest Blue (smallest count)
  "Benign (BP4) - Moderate"    = "#66B2FF",  # Medium-Light Blue
  "Benign (BP4) - Strong"      = "#0057B7",  # Medium-Strong Blue
  "Benign (BP4) - Very Strong" = "#003366",  # Darkest Blue (largest count)
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  "Pathogenic (PP3) - Supporting" = "#FFCCCC",  # Lightest Red (smallest count)
  "Pathogenic (PP3) - Moderate"   = "#FF6666",  # Medium Red
  "Pathogenic (PP3) - Strong"     = "#CC0000"   # Dark Red (largest count)
)

# **Ensure colors are assigned using Fill_Label instead of just Subcategory**
ggplot(classification_counts, aes(x = Category, y = Count, fill = Fill_Label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars with correct ordering
  scale_fill_manual(values = custom_palette, guide = guide_legend(reverse = TRUE)) +  # Reverse legend to match stacking
  scale_y_continuous(labels = comma) +  # Use numeric format for y-axis
  geom_text(aes(label = ifelse(Percent > 5, paste0(Count, " (", round(Percent, 1), "%)"), "")),  # Hide small labels
            position = position_stack(vjust = 0.5), size = 4, color = "black") +  
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Variant Count",
       fill = "Subcategory") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Save the corrected plot
ggsave("Variant_Classification_StackedBar_Corrected.png", width = 8, height = 6, dpi = 300)

####

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

# Create a corrected data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 3)),  
  Subcategory = c("Supporting", "Moderate", "Strong", "Very Strong",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong"),
  Count = c(63, 62, 102, 494,  # Benign counts
            82, 
            98, 98, 362)  # Pathogenic counts
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# **Reverse Subcategory order so "Very Strong" is at the top in the legend**
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Very Strong", "Strong", "Moderate", "Supporting", "Indeterminate"),
                                            ordered = TRUE)  

# Calculate percentages within each Category
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# **Create a proper fill label combining Category and Subcategory for correct color mapping**
classification_counts$Fill_Label <- paste(classification_counts$Category, classification_counts$Subcategory, sep = " - ")

# Define **correct color mapping** for both Benign (Blue) and Pathogenic (Red)
custom_palette <- c(
  "Benign (BP4) - Very Strong" = "#003366",  # Darkest Blue (largest count)
  "Benign (BP4) - Strong"      = "#0057B7",  # Medium-Strong Blue
  "Benign (BP4) - Moderate"    = "#66B2FF",  # Medium-Light Blue
  "Benign (BP4) - Supporting"  = "#99CCFF",  # Lightest Blue (smallest count)
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  "Pathogenic (PP3) - Strong"     = "#CC0000",  # Dark Red (largest count)
  "Pathogenic (PP3) - Moderate"   = "#FF6666",  # Medium Red
  "Pathogenic (PP3) - Supporting" = "#FFCCCC"   # Lightest Red (smallest count)
)

# **Ensure colors are assigned using Fill_Label instead of just Subcategory**
ggplot(classification_counts, aes(x = Category, y = Count, fill = Fill_Label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars with correct ordering
  scale_fill_manual(values = custom_palette, 
                    labels = c("Very Strong", "Strong", "Moderate", "Supporting", "Indeterminate"),
                    guide = guide_legend(reverse = FALSE)) +  # **Reverse legend order**
  scale_y_continuous(labels = comma) +  # Use numeric format for y-axis
  geom_text(aes(label = ifelse(Percent > 5, paste0(Count, " (", round(Percent, 1), "%)"), "")),  # Hide small labels
            position = position_stack(vjust = 0.5), size = 4, color = "black") +  
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Variant Count",
       fill = "Subcategory") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Save the corrected plot
ggsave("Variant_Classification_StackedBar_Corrected.png", width = 8, height = 6, dpi = 300)


#####

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

# Create a corrected data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 3)),  
  Subcategory = c("Supporting", "Moderate", "Strong", "Very Strong",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong"),
  Count = c(63, 62, 102, 494,  # Benign counts
            82, 
            98, 98, 362)  # Pathogenic counts
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# Ensure Subcategory levels are in correct stacking order
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Very Strong", "Strong", "Moderate", "Supporting", "Indeterminate"),
                                            ordered = TRUE)  

# Calculate percentages within each Category
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# **Split the dataset into Benign (BP4) and Pathogenic (PP3)**
benign_data <- classification_counts %>% filter(Category == "Benign (BP4)")
pathogenic_data <- classification_counts %>% filter(Category == "Pathogenic (PP3)")
indeterminate_data <- classification_counts %>% filter(Category == "Indeterminate")

# Define **separate color mappings**
benign_colors <- c("Very Strong" = "#003366",  # Darkest Blue
                   "Strong"      = "#0057B7",  # Medium-Strong Blue
                   "Moderate"    = "#66B2FF",  # Medium-Light Blue
                   "Supporting"  = "#99CCFF")  # Lightest Blue

pathogenic_colors <- c("Very Strong" = "#800000",  # Darkest Red
                       "Strong"      = "#CC0000",  # Medium-Strong Red
                       "Moderate"    = "#FF6666",  # Medium-Light Red
                       "Supporting"  = "#FFCCCC")  # Lightest Red

# Plot both using two separate `geom_bar()` calls
ggplot() +
  # **Benign (BP4) stacked bar plot**
  geom_bar(data = benign_data, aes(x = Category, y = Count, fill = Subcategory), 
           stat = "identity", position = "stack", width = 0.7) +
  
  # **Pathogenic (PP3) stacked bar plot**
  geom_bar(data = pathogenic_data, aes(x = Category, y = Count, fill = Subcategory), 
           stat = "identity", position = "stack", width = 0.7) +
  
  # **Indeterminate bar (single neutral color)**
  geom_bar(data = indeterminate_data, aes(x = Category, y = Count), 
           stat = "identity", fill = "#B2B2B2", width = 0.7) +  # Gray for Indeterminate
  
  # **Benign color mapping**
  scale_fill_manual(name = "Benign Subcategory",
                    values = benign_colors) +
  
  # **Pathogenic color mapping**
  scale_fill_manual(name = "Pathogenic Subcategory",
                    values = pathogenic_colors) +
  
  # **Reverse legend order (strongest at top)**
  guides(fill = guide_legend(reverse = TRUE)) +
  
  # Add text labels for percentages
  geom_text(data = classification_counts, aes(x = Category, y = Count, label = ifelse(Percent > 5, paste0(Count, " (", round(Percent, 1), "%)"), "")), 
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  
  # Titles and themes
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Variant Count") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Save the corrected plot
ggsave("Variant_Classification_StackedBar_SeparateColors.png", width = 8, height = 6, dpi = 300)


##

# Load required libraries
library(ggplot2)
library(dplyr)
library(scales)  # For percentage formatting

# Create a corrected data frame based on your table
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 4)),  
  Subcategory = c("Supporting", "Moderate", "Strong", "Very Strong",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong", "Very Strong"),
  Count = c(63, 62, 102, 494,  # Benign counts
            82, 
            98, 98, 362, 50)  # Pathogenic counts
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# Ensure Subcategory levels are in correct stacking order (smallest on bottom)
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Supporting", "Moderate", "Strong", "Very Strong", "Indeterminate"),
                                            ordered = TRUE)  

# Calculate percentages within each Category
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# **Create a unique fill label for color mapping**
classification_counts$Fill_Label <- paste(classification_counts$Category, classification_counts$Subcategory, sep = " - ")

# **🔹 Corrected Color Mapping**
custom_palette <- c(
  # **Benign (BP4) in Blue**
  "Benign (BP4) - Very Strong" = "#003366",  # Darkest Blue
  "Benign (BP4) - Strong"      = "#0057B7",  # Medium-Strong Blue
  "Benign (BP4) - Moderate"    = "#66B2FF",  # Medium-Light Blue
  "Benign (BP4) - Supporting"  = "#99CCFF",  # Lightest Blue
  
  # **Indeterminate in Gray**
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  
  # **Pathogenic (PP3) in Red**
  "Pathogenic (PP3) - Very Strong" = "#800000",  # Darkest Red
  "Pathogenic (PP3) - Strong"      = "#CC0000",  # Medium-Strong Red
  "Pathogenic (PP3) - Moderate"    = "#FF6666",  # Medium-Light Red
  "Pathogenic (PP3) - Supporting"  = "#FFCCCC"   # Lightest Red
)

# **Ensure colors are assigned using Fill_Label**
ggplot(classification_counts, aes(x = Category, y = Count, fill = Fill_Label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars with correct ordering
  scale_fill_manual(values = custom_palette, 
                    name = "Subcategory",
                    breaks = rev(unique(classification_counts$Fill_Label)),  # **Ensure legend follows the correct order**
                    labels = c("Very Strong", "Strong", "Moderate", "Supporting", "Indeterminate"),
                    guide = guide_legend(reverse = TRUE)) +  # **Reverse legend order**
  scale_y_continuous(labels = comma) +  # Use numeric format for y-axis
  geom_text(aes(label = ifelse(Percent > 5, paste0(Count, " (", round(Percent, 1), "%)"), "")),  # Hide small labels
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  
  # Titles and themes
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Variant Count") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Save the corrected plot
ggsave("Variant_Classification_StackedBar_FinalFix.png", width = 8, height = 6, dpi = 300)

####

# Load required libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)  # For color palettes

# Create a corrected data frame
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 4)),  
  Subcategory = c("Very Strong", "Strong", "Moderate", "Supporting",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong", "Very Strong"),
  Count = c(494, 102, 62, 63,  # Benign counts
            82,  # Indeterminate
            98, 98, 362, 50)  # Pathogenic counts
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# Ensure Subcategory levels are in correct stacking order (Supporting at the bottom)
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = rev(c("Very Strong", "Strong", "Moderate", "Supporting", "Indeterminate")),
                                            ordered = TRUE)  

# Define separate color palettes for Benign (Blue) and Pathogenic (Red)
benign_colors <- brewer.pal(4, "Blues")  # Shades of Blue for Benign (BP4)
pathogenic_colors <- brewer.pal(4, "Reds")  # Shades of Red for Pathogenic (PP3)
indeterminate_color <- "#B2B2B2"  # Neutral Gray

# Map colors to each category
classification_counts$Fill_Color <- ifelse(classification_counts$Category == "Benign (BP4)",
                                           benign_colors[match(classification_counts$Subcategory, levels(classification_counts$Subcategory))],
                                           ifelse(classification_counts$Category == "Pathogenic (PP3)",
                                                  pathogenic_colors[match(classification_counts$Subcategory, levels(classification_counts$Subcategory))],
                                                  indeterminate_color))

# Plot the stacked bar chart
ggplot(classification_counts, aes(x = Category, y = Count, fill = Fill_Color)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  
  # Manually apply the colors
  scale_fill_identity(guide = "legend", name = "Subcategory",
                      breaks = unique(classification_counts$Fill_Color),
                      labels = rev(levels(classification_counts$Subcategory))) +
  
  # Add text labels inside bars
  geom_text(aes(label = ifelse(Count > 0, paste0(Count), "")), 
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  
  # Titles and theme
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Variant Count") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# Save the corrected plot
ggsave("Variant_Classification_StackedBar_Alternative.png", width = 8, height = 6, dpi = 300)
###

# Load required libraries
library(ggplot2)
library(dplyr)

# Create a corrected data frame
classification_counts <- data.frame(
  Category = c(rep("Benign (BP4)", 4), "Indeterminate", rep("Pathogenic (PP3)", 4)),  
  Subcategory = c("Very Strong", "Strong", "Moderate", "Supporting",
                  "Indeterminate",
                  "Supporting", "Moderate", "Strong", "Very Strong"),
  Count = c(494, 102, 62, 63,  # Benign counts
            82,  # Indeterminate
            98, 98, 362, 50)  # Pathogenic counts
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# **Ensure Subcategory levels are in the correct stacking order (Supporting at the bottom, Very Strong at the top)**
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = rev(c("Very Strong", "Strong", "Moderate", "Supporting", "Indeterminate")),
                                            ordered = TRUE)  

# **Create a unique fill label for correct color mapping**
classification_counts$Fill_Label <- paste(classification_counts$Category, classification_counts$Subcategory, sep = " - ")

# **🔹 Define Correct Color Mapping for Separate Benign & Pathogenic Schemes**
custom_palette <- c(
  # **Benign (BP4) in Blue**
  "Benign (BP4) - Very Strong" = "#003366",  # Darkest Blue
  "Benign (BP4) - Strong"      = "#0057B7",  # Medium-Strong Blue
  "Benign (BP4) - Moderate"    = "#66B2FF",  # Medium-Light Blue
  "Benign (BP4) - Supporting"  = "#99CCFF",  # Lightest Blue
  
  # **Indeterminate in Gray**
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  
  # **Pathogenic (PP3) in Red**
  "Pathogenic (PP3) - Very Strong" = "#800000",  # Darkest Red
  "Pathogenic (PP3) - Strong"      = "#CC0000",  # Medium-Strong Red
  "Pathogenic (PP3) - Moderate"    = "#FF6666",  # Medium-Light Red
  "Pathogenic (PP3) - Supporting"  = "#FFCCCC"   # Lightest Red
)

# **Ensure colors are assigned using Fill_Label**
ggplot(classification_counts, aes(x = Category, y = Count, fill = Fill_Label)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  
  # **Manually assign colors using scale_fill_manual()**
  scale_fill_manual(values = custom_palette, 
                    name = "Subcategory",
                    labels = rev(levels(classification_counts$Subcategory)),  # Ensure legend follows the correct order
                    guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  
  # **Add text labels inside bars**
  geom_text(aes(label = ifelse(Count > 0, paste0(Count), "")), 
            position = position_stack(vjust = 0.5), size = 4, color = "black") +
  
  # **Titles and theme**
  labs(title = "Distribution of Variant Classifications",
       subtitle = "Based on CADD Pathogenicity Predictions",
       x = "Classification Category",
       y = "Variant Count") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        legend.position = "right")

# **Save the corrected plot**
ggsave("Variant_Classification_StackedBar_FinalFix.png", width = 8, height = 6, dpi = 300)
###

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
            98, 98, 362, 0)  # Use 0 instead of "—"
)

# Convert Category to a factor with a specific order
classification_counts$Category <- factor(classification_counts$Category, 
                                         levels = c("Benign (BP4)", "Indeterminate", "Pathogenic (PP3)"))

# **Fix stacking order by ensuring Subcategory is a properly ordered factor**
classification_counts$Subcategory <- factor(classification_counts$Subcategory, 
                                            levels = c("Supporting", "Moderate", "Strong", "Very Strong", "Indeterminate"),
                                            ordered = TRUE)  

# Combine Category and Subcategory for coloring
classification_counts <- classification_counts %>%
  mutate(Combined = ifelse(Category == "Indeterminate",
                           paste(Category, Subcategory, sep = " - "),
                           paste(Category, Subcategory, sep = " - ")))

# Calculate percentages within each Category
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

# **Define Correct Color Mapping (Benign = Blue, Pathogenic = Red, Indeterminate = Gray)**
custom_palette <- c(
  "Benign (BP4) - Supporting"  = "#ADD8E6",  # Light Blue
  "Benign (BP4) - Moderate"    = "#3399FF",  # Medium Blue
  "Benign (BP4) - Strong"      = "#66B2FF",  # Darker Blue
  "Benign (BP4) - Very Strong" = "#003366",  # Darkest Blue
  "Indeterminate - Indeterminate" = "#B2B2B2",  # Neutral Gray
  "Pathogenic (PP3) - Supporting" = "#FFCCCC",  # Lightest Red
  "Pathogenic (PP3) - Moderate"   = "#FF9999",  # Medium Red
  "Pathogenic (PP3) - Strong"     = "#FF6666",  # Darker Red
  "Pathogenic (PP3) - Very Strong"= "#FF3333"   # Darkest Red
)

# **Sort the dataset to ensure correct stacking order**
classification_counts <- classification_counts %>%
  arrange(Category, Subcategory)  # This ensures ggplot2 stacks bars correctly

# **Create the vertical stacked bar plot**
ggplot(classification_counts, aes(x = Category, y = Count, fill = Combined)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  scale_y_continuous(labels = comma) +  # Use normal numeric format for y-axis
  scale_fill_manual(values = custom_palette, guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
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

# **Save the corrected plot**
ggsave("Variant_Classification_StackedBar_Corrected.png", width = 8, height = 6, dpi = 300)

### THIS ONE FOR NOW 

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

# **Calculate percentages within each Category**
classification_counts <- classification_counts %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  ungroup()

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

# **Create the vertical stacked bar plot**
ggplot(classification_counts, aes(x = Category, y = Count, fill = factor(Combined, levels = rev(names(custom_palette))))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
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
ggsave("Variant_Classification_StackedBar_Fixed.png", width = 8, height = 6, dpi = 300)


####

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

####

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
  arrange(Category, Subcategory)  # Ensures ggplot2 stacks bars correctly

# **Create the vertical stacked bar plot with dynamic text color**
ggplot(classification_counts, aes(x = Category, y = Count, fill = factor(Combined, levels = rev(names(custom_palette))))) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +  # Stacked bars
  geom_text(aes(label = ifelse(Count > 0, Count, ""),  # Label inside bars (only count)
                color = ifelse(Subcategory %in% c("Very Strong", "Strong"), "white", "black")),  # Change text color
            position = position_stack(vjust = 0.5), size = 4) +
  scale_y_continuous(labels = comma) +  # Use normal numeric format for y-axis
  scale_fill_manual(values = custom_palette, name = "Subcategory",
                    breaks = rev(names(custom_palette)), guide = guide_legend(reverse = TRUE)) +  # Reverse legend order
  scale_color_identity() +  # Allows using custom text colors inside bars
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
ggsave("Variant_Classification_StackedBar_CountsOnly_TextColor.png", width = 8, height = 6, dpi = 300)






#########
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
