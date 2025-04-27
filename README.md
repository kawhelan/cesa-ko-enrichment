# CesA_KO__Enrichment

This project analyzes Pfam domain enrichment, Gene Ontology (GO) enrichment, and KEGG pathway analysis among differentially expressed genes (DEGs) in a CesA knockout line of Physcomitrium patens. The analysis integrates Pfam annotations, GO terms, and KEGG pathways to investigate the functional categories and biological pathways responding to CesA disruption.

Enrichment analyses include:

- Pfam domain enrichment using Fisher’s exact test.
- GO enrichment analysis (Biological Process, Molecular Function, Cellular Component) via the PANTHER classification system.
- KEGG pathway analysis using DAVID functional annotation.


## How to Run

1. Open the Rmd or R script in RStudio
   - For Pfam domain enrichment analysis: Open "Pfam_Enrichment.Rmd" or "CesA_KO_Pfam_Enrichment.R" in RStudio.
   - For GO enrichment analysis: Open
   - For GO enrichment and KEGG pathway analysis: Open
2. Ensure you have the tidyverse package installed.
3. Run the script from the **project root folder** — it assumes files are in "data/".

## Input Files

For Pfam domain enrichment analysis:
- "Significant-DEGs_WT-vs-KO1.csv": Differential gene expression results
- "phytozome_annotation.txt": Annotation file with Pfam domains
- "pfamA.txt": Pfam domain descriptions

For GO enrichment analysis:
- "
For GO enrichment and KEGG pathway analysis:
- "

## Outputs

For Pfam domain enrichment analysis:
- Top 10 Pfam domain bar plot
- Heatmap of Pfam domains by regulation
- Dot plot of statistically enriched Pfam domains

For GO enrichment analysis:
- list plots

For GO and KEGG pathway analysis:
- list plots

## Author

Kerrina Whelan – final project for BIO539 Big Data Analysis


