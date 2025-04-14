
# Prepare DEG–Pfam pairs
deg_pfam_pairs <- deg_annotated %>%
  select(GeneID, Pfam) %>%
  filter(!is.na(Pfam)) %>%
  separate_rows(Pfam, sep = ",") %>%
  distinct(GeneID, Pfam)

# Preview
glimpse(deg_pfam_pairs)

# Count domain frequencies
pfam_deg_counts <- deg_pfam_pairs %>%
  count(Pfam, name = "deg_count") %>%
  arrange(desc(deg_count))

# Preview
head(pfam_deg_counts)

# Step 2: Prepare background gene–Pfam pairs
background_pfam_pairs <- anno %>%
  select(locusName, Pfam) %>%
  filter(!is.na(Pfam)) %>%
  separate_rows(Pfam, sep = ",") %>%
  distinct(locusName, Pfam)

# Preview
glimpse(background_pfam_pairs)

# Count domains in background
pfam_bg_counts <- background_pfam_pairs %>%
  count(Pfam, name = "bg_count") %>%
  arrange(desc(bg_count))


# Wrap fisher.test in a tryCatch() block to handle cases with sparse data and avoid errors
safe_fisher <- safely(~ fisher.test(matrix(c(.x, .y, other_deg, other_bg), nrow = 2))$p.value)

pfam_stats <- pfam_bg_counts %>%
  full_join(pfam_deg_counts, by = "Pfam") %>%
  mutate(
    deg_count = replace_na(deg_count, 0),
    bg_count = replace_na(bg_count, 0)
  ) %>%
  filter(deg_count + bg_count > 0) %>%
  mutate(
    total_deg = n_distinct(deg_pfam_pairs$GeneID),
    total_bg = n_distinct(background_pfam_pairs$locusName)
  ) %>%
  rowwise() %>%
  mutate(
    other_deg = total_deg - deg_count,
    other_bg = total_bg - bg_count,
    pval = tryCatch(
      fisher.test(matrix(c(deg_count, bg_count, other_deg, other_bg), nrow = 2))$p.value,
      error = function(e) NA_real_
    )
  ) %>%
  ungroup() %>%
  mutate(
    padj = p.adjust(pval, method = "fdr")
  ) %>%
  arrange(padj)

# View top enriched domains
pfam_stats %>%
  filter(padj < 0.05) %>%
  head()

pfam_stats <- pfam_stats %>%
  mutate(
    fold_enrichment = (deg_count / total_deg) / (bg_count / total_bg)
  )

top_domains <- pfam_stats %>%
  slice_min(padj, n = 15) %>%
  left_join(pfam_annotated %>% select(Pfam, Description) %>% distinct(), by = "Pfam") %>%
  mutate(
    Description = str_replace_na(Description, "Unknown"),
    Description = str_remove(Description, "^.*?;"),
    Description = str_remove(Description, "anon.*"),
    Description = str_remove(Description, ":.*"),
    Description = str_trim(Description),
    label = paste(Pfam, "-", Description),
    label = fct_reorder(label, fold_enrichment),
    log_padj = -log10(padj)
  )


# Plot
ggplot(top_domains, aes(x = fold_enrichment, y = label)) +
  geom_point(aes(size = deg_count, color = log_padj)) +
  scale_color_gradient(low = "steelblue", high = "firebrick") +
  labs(x = "Fold Enrichment",
       y = "Pfam Domain",
       size = "# DEGs",
       color = "-log10(FDR-adjusted p-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

