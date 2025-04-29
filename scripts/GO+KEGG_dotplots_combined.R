library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# Load KEGG functional annotation chart from DAVID
kegg <- read.delim("KEGG_functional_annotation_chart.txt", header = TRUE, sep = "\t", quote = "")

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

# Select top 30 pathways
top_n <- 30
kegg_top <- kegg_clean %>% arrange(Benjamini) %>% slice(1:top_n)

# === KEGG dot plot with tighter right-side legend ===
kegg_dotplot <- ggplot(kegg_top, aes(x = GeneRatio, y = reorder(Term_clean, GeneRatio))) +
  geom_point(aes(size = Count, color = neg_log10_FDR)) +
  scale_color_gradientn(
    colors = c("red", "orange", "gold", "blue", "purple"),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
  scale_x_continuous(limits = c(0, 0.4), expand = c(0, 0.01)) + 
  theme_minimal(base_size = 12) +
  labs(
    title = "B) Top Enriched KEGG Pathways",
    x = "Gene Ratio",
    y = "KEGG Pathway"
  ) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  theme(
    legend.position = "right",
    legend.justification = c("center", "bottom"),
    legend.box.margin = margin(t = 0, r = 0, b = 2, l = -5),  # move legend closer
    legend.key.width = unit(1.2, "lines"),
    legend.key.height = unit(0.9, "lines"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(margin = margin(r = 20)),
    plot.title = element_text(size = 16, face = "bold", hjust = 1.2),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2)
  )

# === GO dot plot with tighter right-side legend ===
GO_dotplot <- GO_dotplot +
  theme(
    legend.position = "right",
    legend.justification = c("center", "bottom"),
    legend.box.margin = margin(t = 0, r = 0, b = 2, l = 0),  # move legend closer
    legend.key.width = unit(1.2, "lines"),
    legend.key.height = unit(0.9, "lines"),
    legend.text = element_text(size = 10),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 16, face = "bold", hjust = 1),
    legend.title = element_text(size = 10)
  )

# === Combine both plots side by side ===
combined_dotplots <- GO_dotplot + kegg_dotplot + plot_layout(widths = c(1, 1.2))

# === Save final output ===
ggsave("GO_KEGG_dotplots_combined.png",
       plot = combined_dotplots,
       width = 16, height = 10, dpi = 600)

# Show the result
print(combined_dotplots)



