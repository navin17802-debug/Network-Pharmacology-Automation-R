# Automated Network Pharmacology Pipeline in R 🧬🌿

[![R Version](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An end-to-end, automated computational workflow built in **R** to streamline network pharmacology analyses. This pipeline replaces manual web-tool clicking with a 100% reproducible script covering target identification, intersection, Protein-Protein Interaction (PPI) networking, network topology (hub genes), and functional enrichment.

---

## 🚀 Workflow Overview (Steps 1 to 6B)
1. **Ligand Collection:** Input SMILES strings for phytochemical/synthetic ligands.
2. **Target Prediction & Disease Curation:** Collection of phyto-targets and disease-associated target genes.
3. **Target Intersection (Step 4):** Programmatic Venn diagram generation to identify common candidate targets.
4. **PPI Network Construction (Step 5A):** Querying the STRING database for protein interactions.
5. **Topological Analysis (Step 5B):** Computing degree centrality via `igraph` to isolate core hub genes.
6. **Functional Enrichment (Steps 6A & 6B):** Gene Ontology (GO) Biological Process and KEGG pathway enrichment via `clusterProfiler`.

---

## 📦 Prerequisites & Dependencies
Make sure you have R and RStudio installed. Run the following command in your R console to install the required packages:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("STRINGdb", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("tidyverse", "igraph", "ggplot2"))
