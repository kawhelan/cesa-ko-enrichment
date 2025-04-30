library(ggplot2)
library(dplyr)
library(readr)
library(viridis)
library(patchwork)

# === Load and clean the Biological Process table ===
bp <- read.delim("PANTHER_bp_analysis_table.txt", skip = 11, header = TRUE, sep = "\t")
mf <- read.delim("PANTHER_mf_analysis_table.txt", skip = 11, header = TRUE, sep = "\t")
cc <- read.delim("PANTHER_cc_analysis_table.txt", skip = 11, header = TRUE, sep = "\t")

# Rename and clean function
clean_panther <- function(df) {
  colnames(df)[1:8] <- c("GO_term", "REF_count", "DEG_count", "Expected", "OverUnder",
                         "Fold_Enrichment", "P_value", "FDR")
  df <- df %>%
    filter(!is.na(FDR) & FDR < 0.05) %>%
    mutate(
      Fold_Enrichment = as.numeric(Fold_Enrichment),
      FDR = as.numeric(FDR),
      DEG_count = as.numeric(DEG_count),
      REF_count = as.numeric(REF_count),
      GeneRatio = DEG_count / REF_count,
      neg_log10_FDR = -log10(FDR)
    )
  df$GO_term <- gsub(".*\\/| \\(GO:\\d+\\)", "", df$GO_term)
  return(df)
}

# Apply to each dataset
bp <- clean_panther(bp)
mf <- clean_panther(mf)
cc <- clean_panther(cc)

# === Dot plot (Top 30 by FDR) ===
bp_top <- bp %>% arrange(FDR) %>% slice(1:30)

GO_dotplot <- ggplot(bp_top, aes(x = GeneRatio, y = reorder(GO_term, GeneRatio))) +
  geom_point(aes(size = DEG_count, color = neg_log10_FDR)) +
  scale_color_gradientn(
    colors = c("red", "orange", "gold", "blue", "purple"),
    name = expression(-log[10](FDR))
  ) +
  scale_size_continuous(name = "Gene Count", range = c(2, 8)) +
  theme_minimal(base_size = 13) +
  xlim(0, 0.5) +
  labs(
    title = "A) Top Enriched GO Biological Processes",
    x = "Gene Ratio",
    y = "GO Term"
  ) +
  theme(
    legend.key.size = unit(1, "lines"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(GO_dotplot)

# Save dot plot
ggsave("GO_dotplot.png", plot = GO_dotplot, width = 10, height = 8, dpi = 600)

# === Pie chart function ===

# Create a unified color palette for pie charts
pie_colors <- c(
  "#1f78b4",  # strong blue
  "#33a02c",  # vibrant green
  "#e31a1c",  # crimson red
  "#ff7f00",  # bold orange
  "#6a3d9a",  # deep violet
  "#cab2d6",  # soft lavender (replacement for brown)
  "#a6cee3",  # pale blue
  "gray60"    # Others
)


make_pie <- function(df, title, pie_colors) {
  df <- df %>%
    mutate(Label = GO_term, Count = as.numeric(DEG_count)) %>%
    arrange(desc(Count))
  
  top_n <- 7
  df_top <- df[1:top_n, ]
  other_count <- sum(df$Count) - sum(df_top$Count)
  
  df_pie <- rbind(
    df_top[, c("Label", "Count")],
    data.frame(Label = "Others", Count = other_count)
  ) %>%
    mutate(
      Label = factor(Label, levels = c(as.character(df_top$Label), "Others")),
      Percent = Count / sum(Count) * 100
    )
  
  ggplot(df_pie, aes(x = "", y = Percent, fill = Label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    theme_void(base_size = 10) +
    labs(title = title) +
    scale_fill_manual(values = setNames(pie_colors[1:nrow(df_pie)], levels(df_pie$Label))) +
    theme(legend.position = "right")
}


# Generate pie charts
plot_cc <- make_pie(cc, "B) Cellular Component", pie_colors)
plot_mf <- make_pie(mf, "C) Molecular Function", pie_colors)
plot_bp_pie <- make_pie(bp, "D) Biological Process", pie_colors)

# === Combine all into one figure ===
combined_plot <- (GO_dotplot | (plot_cc / plot_mf / plot_bp_pie)) + plot_layout(widths = c(2, 1))

print(combined_plot)

# Save combined figure
ggsave("GO_dotplot_and_piecharts.png", plot = combined_plot, width = 14, height = 8, dpi = 600)

