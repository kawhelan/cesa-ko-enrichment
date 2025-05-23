# Load required packages
library(tidyverse) # Core data manipulation and visualization (includes dplyr, ggplot2, readr, etc.)
library(patchwork) # For combining multiple plots into a combined figure
library(RColorBrewer)
library(viridis)  # Optional for color-blind-friendly plotting
library(rprojroot)    # Ensures consistent paths regardless of where the .Rmd or .R file is run

# Load data
# Define input directory using project root
input_dir <- file.path(find_rstudio_root_file(), "data")  # All input files should be stored in a 'data/' folder in the project root

# Read in the PANTHER analysis tables
# These are output files from the PANTHER classification system, each containing enrichment results
# for a different Gene Ontology category (BP = Biological Process, MF = Molecular Function, CC = Cellular Component)
bp <- read.delim(file.path(input_dir, "PANTHER_bp_analysis_table.txt"), 
                 skip = 11, header = TRUE, sep = "\t")

mf <- read.delim(file.path(input_dir, "PANTHER_mf_analysis_table.txt"), 
                 skip = 11, header = TRUE, sep = "\t")

cc <- read.delim(file.path(input_dir, "PANTHER_cc_analysis_table.txt"), 
                 skip = 11, header = TRUE, sep = "\t")

# Data cleaning
# Define cleaning function
clean_panther <- function(df) {
  colnames(df)[1:8] <- c("GO_term", "REF_count", "DEG_count", "Expected", "OverUnder", "Fold_Enrichment", "P_value", "FDR") # Rename columns for clarity
  df <- df %>% filter(!is.na(FDR) & FDR < 0.05) # Keep only significant terms (FDR < 0.05) and remove NAs
  # Make sure columns are numeric
  df$Fold_Enrichment <- as.numeric(df$Fold_Enrichment)
  df$FDR <- as.numeric(df$FDR)
  return(df)
}

# Apply cleaning
bp <- clean_panther(bp)
mf <- clean_panther(mf)
cc <- clean_panther(cc)

# Remove GO IDs from labels
clean_labels <- function(term) {
  sub(" \\(GO:\\d+\\)", "", term)
}

bp$GO_term <- clean_labels(bp$GO_term)
mf$GO_term <- clean_labels(mf$GO_term)
cc$GO_term <- clean_labels(cc$GO_term)

# Top 15 bar plot
# Sort BP GO terms by ascending FDR (most significant first) and select the top 15 most significant
bp_top <- bp %>% arrange(FDR) %>% slice(1:15)

plot_bp <- ggplot(bp_top, aes(x = reorder(GO_term, -Fold_Enrichment), y = Fold_Enrichment, fill = Fold_Enrichment)) +  # Reorder GO terms by Fold Enrichment (descending order)
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_viridis(option = "C", direction = -1) + # Apply the Viridis color scale (option "C", reversed direction)
  theme_minimal(base_size = 12) +
  labs(title = "A) GO classification", x = NULL, y = "Fold Enrichment") +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 12, color = "black"))  # Adjust size as needed

plot_bp

# Gene Ontology pie charts
# Define pie chart function
make_pie <- function(df, title) {
  df <- df %>%
    mutate(Label = clean_labels(GO_term),
           Count = as.numeric(DEG_count)) %>%
    arrange(desc(Count)) # Sort GO terms by DEG count in descending order
  # Set number of top GO terms to show in the pie chart and select the top seven terms
  top_n <- 7
  df_top <- df[1:top_n, ]
  other_count <- sum(df$Count) - sum(df_top$Count) # Add 'Others' slice
  # Calculate percentage for each slice
  df_pie <- rbind(
    df_top[, c("Label", "Count")],
    data.frame(Label = "Others", Count = other_count)
  )
  
  df_pie <- df_pie %>%
    mutate(Percent = Count / sum(Count) * 100)
  
  ggplot(df_pie, aes(x = "", y = Percent, fill = Label)) +
    geom_bar(stat = "identity", width = 1, color = "white") + # Draw bars with white borders between
    coord_polar("y") + # Convert bar chart to pie chart
    scale_fill_brewer(palette = "Set3") +
    theme_void(base_size = 12) +
    labs(title = title) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12)
    )
}

# Generate pie charts
plot_cc <- make_pie(cc, "B) Cellular Component")
plot_mf <- make_pie(mf, "C) Molecular Function")
plot_bp_pie <- make_pie(bp, "D) Biological Process")

combined_pies <- plot_cc / plot_mf / plot_bp_pie  # Stack them vertically

combined_pies

# Combine the bar plot and pie charts into one figure
# Place the GO bar plot on the left and stack the three pie charts vertically on the right
final_plot <- (plot_bp | (plot_cc / plot_mf / plot_bp_pie)) + 
  plot_layout(widths = c(1.8, 1.8)) # Set relative widths for the bar plot and pie charts to balance spacing

final_plot

# Optional: save figures as png and pdf
ggsave("GO_summary_figure.png", plot = final_plot, width = 14, height = 12, dpi = 300)
ggsave("GO_summary_figure.pdf", plot = final_plot, width = 10, height = 12)

