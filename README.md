# CesA_KO_Pfam_Enrichment

This project analyzes Pfam domain enrichment among differentially expressed genes (DEGs) in a CesA knockout line using GO annotations and Fisher's exact test.


## How to Run

1. Open "CesA_KO_Pfam_Enrichment_final.R" in RStudio.
2. Ensure you have the tidyerse package installed.
3. Run the script from the **project root folder** — it assumes files are in "data/".

## Input Files

- "Significant-DEGs_WT-vs-KO1.csv": Differential gene expression results
- "phytozome_annotation.txt": Annotation file with Pfam domains
- "pfamA.txt": Pfam domain descriptions

## Outputs

- Top 10 Pfam domain bar plot
- Heatmap of Pfam domains by regulation
- Dot plot of statistically enriched Pfam domains

## Author

Kerrina Whelan – final project for BIO539 Big Data Analysis


