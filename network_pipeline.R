# =====================================================================
# Project: Automated Network Pharmacology Pipeline (Steps 1 to 6B)
# Description: End-to-end R script for multi-target drug mechanism analysis
# =====================================================================

# 1. Load Required Libraries
library(tidyverse)
library(STRINGdb)
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


# =====================================================================
# STEP 1 & 2: Ligands Dataset (SMILES Strings)
# =====================================================================
smiles_list <- c(
  "CC1C(C(C(O1)NC2C(C(C(C(O2)N)N)O)O)O)NC3C(C(C(C(O3)O)N)O)O",
  "CC1C(C(C(O1)NC2C(C(C(C(O2)O)NC(=O)CCN)N)O)O)N",
  "CNCCc1ccc(O)c(O)c1",
  "CC(C)[N+](C)(C)C(=O)Oc1ccc([N+](C)(C)C)cc1",
  "CN1C(=O)CN=C(c2ccccc2Cl)c3ccccc13",
  "O=C1NC(=S)NC(=O)C1(CC(C)C)C2=CC=CC=C2",
  "C1=CC2=C3C(=C1)C=CC4=C3C=CC5=C4C=CC5=C2",
  "C1=CC=C2C=C(C=CC2=C1)N",
  "CCNCCOC1=CC=CC2=C1CCCC2",
  "CN(C)C1=CC=C(C=C1)C(C#N)C2=CC=CC=C2O"
)

ligands_df <- data.frame(
  Ligand_ID = paste0("Ligand_", seq_along(smiles_list)),
  SMILES = smiles_list,
  stringsAsFactors = FALSE
)


# =====================================================================
# STEP 3 & 4: Target Lists & Venn Diagram (Common Targets)
# =====================================================================
core_24_genes <- c(
  "APP", "BACE1", "PSEN1", "PSEN2", "MAPT", "GSK3B",
  "CDK5", "PIN1", "UBE2N", "TREM2", "APOE", "IL1B",
  "TNF", "CX3CR1", "SOD1", "PINK1", "PARK7", "NFE2L2",
  "GRIN2B", "BDNF", "SNCA", "LRP1", "ABCB1", "CLDN5"
)

phyto_targets   <- c(core_24_genes, "EGFR", "VEGFA", "STAT3", "MMP9", "CASP3", "AKT1")
disease_targets <- c(core_24_genes, "TP53", "MYC", "HIF1A", "PTEN", "BCL2", "ESR1", "MTOR", "JUN")

common_targets <- intersect(phyto_targets, disease_targets)

# Custom Venn Diagram using ggplot2
theta <- seq(0, 2 * pi, length.out = 200)
circle_phyto   <- data.frame(x = cos(theta) - 0.55, y = sin(theta), Set = "Phyto Targets")
circle_disease <- data.frame(x = cos(theta) + 0.55, y = sin(theta), Set = "Disease Targets")
venn_circles   <- rbind(circle_phyto, circle_disease)

p_venn <- ggplot() +
  geom_polygon(data = venn_circles, aes(x = x, y = y, fill = Set), alpha = 0.45, color = "black", linewidth = 0.8) +
  scale_fill_manual(values = c("Phyto Targets" = "#4E79A7", "Disease Targets" = "#F28E2B")) +
  annotate("text", x = -0.95, y = 0, label = paste0("Phyto Only\n", length(setdiff(phyto_targets, disease_targets))), size = 5, fontface = "bold") +
  annotate("text", x = 0,     y = 0, label = paste0("Common\n", length(common_targets)), size = 5.5, fontface = "bold") +
  annotate("text", x = 0.95,  y = 0, label = paste0("Disease Only\n", length(setdiff(disease_targets, phyto_targets))), size = 5, fontface = "bold") +
  coord_fixed() + theme_void() +
  labs(title = "Step 4: Venn Diagram of Common Targets", fill = "Target Group") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), legend.position = "bottom")

print(p_venn)


# =====================================================================
# STEP 5A & 5B: PPI Network & Topological Analysis
# =====================================================================
input_df <- data.frame(gene = common_targets, stringsAsFactors = FALSE)

string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = "")
mapped_genes <- string_db$map(input_df, "gene", removeUnmappedRows = TRUE)
interactions <- string_db$get_interactions(mapped_genes$STRING_id)

igraph_net <- graph_from_data_frame(d = interactions[, c("from", "to")], directed = FALSE)
igraph_net <- simplify(igraph_net)

deg_centrality <- degree(igraph_net)
symbol_map     <- setNames(mapped_genes$gene, mapped_genes$STRING_id)

hub_summary <- data.frame(
  STRING_id  = names(deg_centrality),
  GeneSymbol = unname(symbol_map[names(deg_centrality)]),
  Degree     = as.numeric(deg_centrality),
  stringsAsFactors = FALSE
) %>% arrange(desc(Degree))

top10_hubs <- head(hub_summary, 10)
print(top10_hubs)

# Hub genes bar chart
p_hubs <- ggplot(top10_hubs, aes(x = reorder(GeneSymbol, Degree), y = Degree, fill = Degree)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_gradient(low = "#56B1F7", high = "#132B43") +
  theme_minimal() +
  labs(title = "Step 5B: Top 10 Hub Genes by Degree Centrality", x = "Hub Gene", y = "Degree (Connections)")
print(p_hubs)


# =====================================================================
# STEP 6A & 6B: GO & KEGG Pathway Enrichment Analysis
# =====================================================================
entrez_ids <- bitr(common_targets, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO Biological Process
ego <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
if (nrow(as.data.frame(ego)) > 0) {
  print(dotplot(ego, showCategory = 10, title = "Step 6A: GO Biological Process Enrichment") + theme_minimal())
}

# KEGG Pathway
ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
if (nrow(as.data.frame(ekegg)) > 0) {
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  print(dotplot(ekegg, showCategory = 10, title = "Step 6B: KEGG Pathway Enrichment") + theme_minimal())
}

cat("\n=== PIPELINE EXECUTION SUCCESSFUL! ===\n")