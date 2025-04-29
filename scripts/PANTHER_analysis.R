# PANTHER analysis
library(tidyverse)

# Read the PANTHER GO biological process enrichment table
panther <- read.delim("data/PANTHER_bp_analysis_table.txt",  # change input file as needed
                      skip = 11,        # Skip the header metadata
                      header = TRUE,
                      sep = "\t",
                      quote = "",
                      stringsAsFactors = FALSE)

# Rename columns for easier access
colnames(panther) <- c("GO_term", "REF_count", "DEG_count", "Expected", "OverUnder",
                       "Fold_Enrichment", "P_value", "FDR")

# Convert relevant columns to numeric
panther$Fold_Enrichment <- as.numeric(panther$Fold_Enrichment)
panther$FDR <- as.numeric(panther$FDR)

# Filter significant GO terms
panther_sig <- subset(panther, FDR < 0.05)

# Top 15 enriched GO terms by FDR
panther_top <- panther_sig[order(panther_sig$FDR), ][1:15, ]


ggplot(panther_top, aes(x = reorder(GO_term, -Fold_Enrichment), y = Fold_Enrichment)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(title = "Top Enriched GO Biological Processes",
       x = "GO Term",
       y = "Fold Enrichment")
