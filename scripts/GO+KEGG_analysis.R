# Load required packages

library(tidyverse)

# Load data

# Load GO PANTHER tables
bp <- read.delim("data/PANTHER_bp_analysis_table.txt", skip = 11, header = TRUE, sep = "\t")
mf <- read.delim("data/PANTHER_mf_analysis_table.txt", skip = 11, header = TRUE, sep = "\t")
cc <- read.delim("data/PANTHER_cc_analysis_table.txt", skip = 11, header = TRUE, sep = "\t")

# Load KEGG annotation chart (from DAVID)
kegg <- read.delim("data/KEGG_functional_annotation_chart.txt", header = TRUE, sep = "\t", quote = "")

# Clean data

# Define PANTHER cleaning function
clean_panther <- function(df) {
  colnames(df)[1:8] <- c("GO_term", "REF_count", "DEG_count", "Expected", "OverUnder",
                         "Fold_Enrichment", "P_value", "FDR")
  df %>%
    filter(!is.na(FDR) & FDR < 0.05) %>% # Keep only rows where FDR is not NA and less than 0.05 (significant terms)
    mutate(
      # Make sure columns are numeric
      Fold_Enrichment = as.numeric(Fold_Enrichment),
      FDR = as.numeric(FDR),
      DEG_count = as.numeric(DEG_count),
      REF_count = as.numeric(REF_count),
      GeneRatio = DEG_count / REF_count, # Calculate proportion of DEGs among background genes
      neg_log10_FDR = -log10(FDR), # Calculate negative log10 of FDR for plotting
      GO_term = gsub(".*\\/| \\(GO:\\d+\\)", "", GO_term) # Clean GO term labels by removing slashes and GO IDs
    )
}

# Apply to each GO dataset
bp <- clean_panther(bp)
mf <- clean_panther(mf)
cc <- clean_panther(cc)

# Clean KEGG data
kegg_clean <- kegg %>%
  filter(!is.na(Benjamini) & Benjamini < 0.05) %>%
  mutate(
    GeneRatio = Count / Pop.Hits,
    neg_log10_FDR = -log10(Benjamini),
    Count = as.numeric(Count),
    Term_clean = Term %>%
      gsub("^.*[:~]", "", .) %>%
      gsub(" \\(.*\\)", "", .)
  )

# Sort GO terms by ascending FDR (most significant first) and select the top 30 based on FDR
bp_top <- bp %>% arrange(FDR) %>% slice(1:30)

# Create GO dot plot
GO_dotplot <- ggplot(bp_top, aes(x = GeneRatio, y = reorder(GO_term, GeneRatio))) +
  geom_point(aes(size = DEG_count, color = neg_log10_FDR)) +
  scale_color_gradientn(
    colors = c("red", "orange", "gold", "blue", "purple"),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
  theme_minimal(base_size = 13) +
  xlim(0, 0.5) +
  labs(x = "Gene Ratio", y = "GO Term")

GO_dotplot

# KEGG Pathway Enrichment Dot Plot

# Sort KEGG pathways by ascending Benjamini-adjusted p-value and select the top 30 pathways
kegg_top <- kegg_clean %>% arrange(Benjamini) %>% slice(1:30)

# Create KEGG dot plot
kegg_dotplot <- ggplot(kegg_top, aes(x = GeneRatio, y = reorder(Term_clean, GeneRatio))) +
  geom_point(aes(size = Count, color = neg_log10_FDR)) +
  scale_color_gradientn(
    colors = c("red", "orange", "gold", "blue", "purple"),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
  scale_x_continuous(limits = c(0, 0.4), expand = c(0, 0.01)) +
  theme_minimal(base_size = 12) +
  labs(x = "Gene Ratio", y = "KEGG Pathway")

kegg_dotplot


