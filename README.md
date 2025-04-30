# CESA_KO_Enrichment

This project analyzes Pfam domain enrichment, Gene Ontology (GO) enrichment, and KEGG pathway analysis among differentially expressed genes (DEGs) in a CESA knockout line of *Physcomitrium patens*. The analysis integrates Pfam annotations, GO terms, and KEGG pathways to investigate the functional categories and biological pathways responding to CESA disruption.

## Research Question

What functional categories, protein domains, and biological pathways are enriched among differentially expressed genes in *Physcomitrium patens* following CESA knockout?

Enrichment analyses include:

- **Pfam domain enrichment** using Fisher’s exact test.
- **GO enrichment analysis** (Biological Process, Molecular Function, Cellular Component) via the PANTHER classification system.
- **KEGG pathway analysis** using DAVID functional annotation.


## How to Run

1. Open the R Markdown (`.Rmd`) or R script (`R`) file in RStudio.
2. Available scripts:
   - Pfam domain enrichment analysis: `Pfam_Enrichment.Rmd` or `Pfam_Enrichment.R`
   - GO enrichment analysis: `GO_Analysis.Rmd` or `GO_Analysis.R`
   - GO enrichment and KEGG pathway analysis: `GO+KEGG_Analysis.Rmd` or `GO+KEGG_Analysis.R`
3. Make sure the following packages are installed:
   - `tidyverse`
   - `patchwork`
   - `viridis`
   - `rprojroot` (required for setting relative paths when knitting `.Rmd` files)
4. Run the script from the **project root folder** — it assumes files are in `data/`.

## Input Files

- `Significant-DEGs_WT-vs-KO1.csv`: Differential gene expression results
- `phytozome_annotation.txt`: Annotation file with Pfam domains
- `pfamA.txt`: Pfam domain descriptions
- `PANTHER_bp_analysis_table.txt`: GO enrichment (Biological Processes)
- `PANTHER_mf_analysis_table.txt`: GO enrichment (Molecular Functions) 
- `PANTHER_cc_analysis_table.txt`: GO enrichment (Cellular Components)
- `KEGG_functional_annotation_chart.txt`: KEGG enrichment results from DAVID
  

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

## Reports

Rendered HTML versions of the analysis (Pfam, GO, and KEGG enrichment) and the final project report are available in `docs/`. These files were generated from R Markdown and include all code, figures, and narrative explanations.


## Author

**Kerrina Whelan**
Final project for BIO 539: Big Data Analysis

## Acknowledgments

This project was completed as part of the final project requirements for BIO 539: Big Data Analysis.



