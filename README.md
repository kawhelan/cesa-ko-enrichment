# CESA_KO_Enrichment

This project analyzes Pfam domain enrichment, Gene Ontology (GO) enrichment, and KEGG pathway analysis among differentially expressed genes (DEGs) in a CESA knockout line of *Physcomitrium patens*. The analysis integrates Pfam annotations, GO terms, and KEGG pathways to investigate the functional categories and biological pathways responding to CESA disruption.

## Research Question

What functional categories, protein domains, and biological pathways are enriched among differentially expressed genes in *Physcomitrium patens* following CESA knockout?

Enrichment analyses include:

- Pfam domain enrichment using Fisher’s exact test.
- GO enrichment analysis (Biological Process, Molecular Function, Cellular Component) via the PANTHER classification system.
- KEGG pathway analysis using DAVID functional annotation.


## How to Run

1. Open the Rmd or R script in RStudio
   - For Pfam domain enrichment analysis: Open "Pfam_Enrichment.Rmd" or "Pfam_Enrichment.R" in RStudio.
   - For GO enrichment analysis: Open "GO_Analysis.Rmd" or "GO_Analysis.R" in RStudio.
   - For GO enrichment and KEGG pathway analysis: Open "GO+KEGG_Analysis.Rmd" or "GO+KEGG_Analysis.R" in RStudio.
2. Ensure you have the tidyverse, patchwork, and viridis packages installed.
3. Run the script from the **project root folder** — it assumes files are in "data/".

## Input Files

- "Significant-DEGs_WT-vs-KO1.csv": Differential gene expression results
- "phytozome_annotation.txt": Annotation file with Pfam domains
- "pfamA.txt": Pfam domain descriptions
- "PANTHER_bp_analysis_table.txt": Output table from PANTHER GO enrichment analysis for Biological Process (BP) terms
- "PANTHER_mf_analysis_table.txt": Output table from PANTHER GO enrichment analysis for Molecular Function (MF) terms
- "PANTHER_cc_analysis_table.txt": Output table from PANTHER GO enrichment analysis for Cellular Component (CC) terms
- "KEGG_functional_annotation_chart.txt": Output table from KEGG pathway enrichment analysis with enriched pathways among DEGs
  

## Outputs

- Top 10 Pfam domain bar plot
- Heatmap of Pfam domains by regulation
- Dot plot of statistically enriched Pfam domains
- Lollipop of significantly enriched Pfam domains
- Top 15 GO biological process bar plot
- Pie charts for GO biological processes, molecular functions, and cellular components
- Combined bar plot and pie chart figure
- GO enrichment dot plot
- KEGG enrichment dot plot


## Author

Kerrina Whelan – final project for BIO 539 Big Data Analysis

## Acknowledgments

This project was completed as part of the final project requirements for BIO 539 Big Data Analysis.



