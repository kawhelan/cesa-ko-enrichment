
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

# this gives the top 10 Pfam domains but without labels
ggplot(top_pfam, aes(x = reorder(Pfam, n), y = n)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 10 Pfam Domains in All DEGs",
       x = "Pfam Domain",
       y = "Count") +
  theme_minimal()

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


ggplot(pfam_annotated %>% slice_max(n, n = 10),
       aes(x = reorder(label, n), y = n)) +
  geom_col(fill = "darkgreen") +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  coord_flip() +
  labs(title = "Top 10 Pfam Domains in All DEGs",
       x = "Pfam Domain", y = "Count") +
  theme_minimal()





