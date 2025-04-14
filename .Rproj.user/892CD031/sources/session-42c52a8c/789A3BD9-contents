library(tidyverse)

# Load the DEGs
deg_path <- "/mnt/homes4celsrs/kerrina.whelan/BIO439-BIg_Data/GitHub/CesA_KO_GO_Enrichment/data/Significant-DEGs_WT-vs-KO1.csv"
degs <- read_csv(deg_path)

# View the the data
glimpse(degs)
head(degs)

# Load the annotation file
anno <- read_tsv("/mnt/homes4celsrs/kerrina.whelan/BIO439-BIg_Data/GitHub/CesA_KO_GO_Enrichment/data/phytozome_annotation.txt")

# View the columns
glimpse(anno)

# Join by GeneID (from degs) and locusName (from anno)
deg_annotated <- degs %>%
  left_join(anno, by = c("GeneID" = "locusName"))

# Check result
glimpse(deg_annotated)

# Separate comma-separated Pfams into individual rows
pfam_all <- deg_annotated %>%
  filter(!is.na(Pfam)) %>%
  separate_rows(Pfam, sep = ",") %>%
  count(Pfam, sort = TRUE)

head(pfam_all, 10)

top_pfam <- pfam_all %>%
  slice_max(n, n = 10)

# Read in only the first three columns of the pfam description table: Pfam ID, short name, and description
pfam_df <- read_table("data/pfamA.txt", col_names = FALSE) %>%
  select(X1, X2, X3:X10) %>%  # X1 = ID, X2 = ShortName, X3+ = start of description
  unite("Description", X3:X10, sep = " ", na.rm = TRUE) %>%
  rename(Pfam = X1, ShortName = X2)

glimpse(pfam_df)

# Join with pfam_all df and clean up labels
pfam_annotated <- pfam_all %>%
  left_join(pfam_df, by = "Pfam") %>%
  mutate(label = paste(Pfam, "-", str_remove(Description, "(anon|:|;).*")))

top10 <- pfam_annotated %>%
  slice_max(n, n = 10)

# Create a lookup table of domain functions
pfam_function_map <- tribble(
  ~Pfam,      ~Function,
  "PF07714", "Kinase/Signaling",
  "PF00069", "Kinase/Signaling",
  "PF00005", "Transporter",
  "PF00249", "Transcription Factor",
  "PF00067", "Detox/Metabolism",
  "PF00847", "Transcription Factor",
  "PF00046", "Transcription Factor",
  "PF07690", "Transporter",
  "PF13855", "Defense/Signaling",
  "PF08263", "Defense/Signaling"
)

# Join with your domain counts
top10_annotated <- top10 %>%
  left_join(pfam_function_map, by = "Pfam")

ggplot(top10_annotated, aes(x = reorder(label, n), y = n, fill = Function)) +
  geom_col() +
  geom_text(aes(label = n), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(title = "Top 10 Pfam Domains in All DEGs",
       x = "Pfam Domain", y = "Count") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Build heatmap of Pfam domain frequencies
# Add a new column to deg_annotated for direction
deg_annotated <- deg_annotated %>%
  mutate(regulation = case_when(
    log2FoldChange > 0 ~ "Upregulated",
    log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Unchanged"
  ))

# Extract Pfam and regulation info, then filter and count Pfam domains by regulation
pfam_regulation_counts <- deg_annotated %>%
  filter(!is.na(Pfam)) %>%
  separate_rows(Pfam, sep = ",") %>%  # in case of multiple domains per gene
  count(Pfam, regulation) %>%
  pivot_wider(names_from = regulation, values_from = n, values_fill = 0)

# Join descriptions
pfam_heatmap_data <- pfam_regulation_counts %>%
  left_join(pfam_df %>% select(Pfam, Description), by = "Pfam") %>%
  mutate(label = paste(Pfam, "-", str_remove(Description, "(anon|:|;).*")))

# Plot heatmap
# Convert to long format
pfam_long <- pfam_heatmap_data %>%
  pivot_longer(cols = c("Upregulated", "Downregulated"),
               names_to = "Direction", values_to = "Count")

# Select top 50 domains by total count
top_pfam <- pfam_long %>%
  group_by(label) %>%
  summarise(total = sum(Count)) %>%
  slice_max(total, n = 40)

# Filter to top
pfam_heat <- pfam_long %>%
  filter(label %in% top_pfam$label)

# Heatmap plot
# Convert to long format
pfam_long <- pfam_heatmap_data %>%
  pivot_longer(cols = c("Upregulated", "Downregulated"),
               names_to = "Direction", values_to = "Count")

# Heatmap plot
ggplot(pfam_heat, aes(x = Direction, y = reorder(label, Count), fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "lightblue", high = "goldenrod") +
  labs(title = "Pfam Domain Frequencies in Up/Downregulated DEGs",
       x = "Regulation", y = "Pfam Domain") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )


ggplot(pfam_heat, aes(x = Direction, y = reorder(label, Count), fill = Count)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "#bdd7e7",   # light blue
    high = "#f16913"   # medium orange
  ) +
  labs(
    title = "Pfam Domain Frequencies in Up/Downregulated DEGs",
    x = "Regulation",
    y = "Pfam Domain",
    fill = "Count"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )

# need to pick blue/red or blue/gold




