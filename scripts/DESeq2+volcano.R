
library(DESeq2)

# Load gene count data
count_data <- read.table("KO1-WT_gene_counts.txt", header = TRUE, row.names = 1, sep = "\t")

# Check the first few rows
head(count_data)

# Prepare metadata with conditions
sample_info <- read.table("testCondition_KO1vsWT.txt", header = TRUE, row.names = 1, sep = "\t")

# Check column names
colnames(count_data)

# Check sample_info row names
rownames(sample_info)

# Replace "." with "-" in column names of count_data to match sample_info
colnames(count_data) <- gsub("\\.", "-", colnames(count_data))

# Now check again to confirm
print(colnames(count_data))
print(rownames(sample_info))

# Should return TRUE if now matching and in order:
all(rownames(sample_info) %in% colnames(count_data))

# Reorder count_data columns to match sample_info row names
count_data <- count_data[, rownames(sample_info)]

# Check again for exact match (should return TRUE)
all(rownames(sample_info) == colnames(count_data))

# DESeq2 object WITHOUT batch since it's all NA
dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = sample_info,
  design = ~ condition  # Only include condition
)

# Run DESeq2
dds <- DESeq(dds)

# Extract results: KO1 vs WT
res <- results(dds, contrast = c("condition", "KO1", "WT"))

# View summary and head
summary(res)
head(res)

# Save results to file
write.table(res, file = "DESeq2_results_KO1_vs_WT.txt", sep = "\t", quote = FALSE, col.names = NA)

# Extract significant genes (adjusted p-value < 0.05)
sig_genes <- res[which(res$padj < 0.05), ]
write.table(sig_genes, file = "DESeq2_significant_genes_padj_0.05.txt", sep = "\t", quote = FALSE, col.names = NA)

# Optional: if using 0.1 threshold
# sig_genes_0_1 <- res[which(res$padj < 0.1), ]
# write.table(sig_genes_0_1, file = "DESeq2_significant_genes_padj_0.1.txt", sep = "\t", quote = FALSE, col.names = NA)

plotMA(res, ylim = c(-5, 5), main = "DESeq2 MA Plot")

vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
vsd_mat <- assay(vsd)
dir.create("data", recursive = TRUE, showWarnings = FALSE)
write.csv(vsd_mat, file = "data/expression_data.csv")

plotPCA(vsd, intgroup = "condition")

library(ggplot2)

# Replace 'res' with your DESeq2 results object if named differently
# Make sure to remove NA values for plotting
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)

# Define significance threshold
padj_cutoff <- 0.05
log2FC_cutoff <- 1  # Set your desired fold change threshold (|LFC| > 1)

# Add a new column to classify significance
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < padj_cutoff & res_df$log2FoldChange > log2FC_cutoff] <- "Upregulated"
res_df$Significance[res_df$padj < padj_cutoff & res_df$log2FoldChange < -log2FC_cutoff] <- "Downregulated"

# Extract gene IDs for downstream analysis (PANTHER)
# Subset significant DEGs with padj < 0.05 and abs(log2FC) > 1
sig_gene_ids <- rownames(res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ])

# Save gene IDs to a text file, one per line (suitable for PANTHER)
writeLines(sig_gene_ids, "PANTHER_gene_list.txt")


volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +  # Scatter points
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Volcano Plot: KO1 vs WT",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value (padj)",
    color = "Regulation"
  ) +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "black") +  # LFC thresholds
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +  # p-value threshold
  theme(legend.position = "top")

print(volcano)

ggsave("volcano_plot_KO1_vs_WT.png", plot = volcano, width = 8, height = 6, dpi = 300)




# publication figure with labels
install.packages("ggrepel")
library(ggrepel)

# Convert results to dataframe and remove NA values
res_df <- as.data.frame(res)
res_df <- na.omit(res_df)

# define thresholds and significance
padj_cutoff <- 0.05  # Adjusted p-value threshold
log2FC_cutoff <- 1  # Fold-change threshold

# Categorize genes
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < padj_cutoff & res_df$log2FoldChange > log2FC_cutoff] <- "Upregulated"
res_df$Significance[res_df$padj < padj_cutoff & res_df$log2FoldChange < -log2FC_cutoff] <- "Downregulated"

# select genes for labeling
# Top 10 most significant (lowest padj) genes for labeling
top_genes <- rownames(res_df[order(res_df$padj), ])[1:10]  # Top 10 by significance
res_df$GeneLabel <- ifelse(rownames(res_df) %in% top_genes, rownames(res_df), NA)

# create volcano plot
volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = GeneLabel)) +
  geom_point(alpha = 0.8, size = 2) +  # Main points
  scale_color_manual(values = c("Upregulated" = "firebrick", "Downregulated" = "steelblue", "Not Significant" = "gray70")) +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "black", size = 0.6) +  # LFC cutoffs
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black", size = 0.6) +  # Padj cutoff
  geom_text_repel(  # Labels
    na.rm = TRUE, 
    size = 3.5, 
    max.overlaps = 15,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = 'grey50'
  ) +
  labs(x = expression(Log[2]~Fold~Change),
    y = expression(-Log[10]~Adjusted~P~value),
    color = "Regulation"
  ) +
  theme_minimal(base_size = 16) +  # Nice clean theme with larger font
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

# display volcano
print(volcano)

# save image
ggsave("volcano_plot_KO1_vs_WT_labeled.png", plot = volcano, width = 8, height = 8, dpi = 600)

