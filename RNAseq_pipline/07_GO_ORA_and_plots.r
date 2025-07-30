######################################################
##### Gene Ontology Over-Representation Analysis #####
######################################################

# Load libraries
library(DESeq2)
library(dplyr)
library(readr)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ggplot2)
library(data.table)

#--------------------------------------
#------------- Analysis ---------------
#--------------------------------------

# Set root directory
root_dir <- "/path/to/the/project"

# Load data
res_df_p8_sig <- read.csv(paste0(root_dir, "/06_fc/ctrl_vs_p8_DESeq2_significant.csv"), header = TRUE, stringsAsFactors = TRUE)
res_df_p15_sig <- read.csv(paste0(root_dir, "/06_fc/ctrl_vs_p15_DESeq2_significant.csv"), header = TRUE, stringsAsFactors = TRUE)
res_p8_vs_p15_sig <- read.csv(paste0(root_dir, "/06_fc/p8_vs_p15_DESeq2_significant.csv"), header = TRUE, stringsAsFactors = TRUE)

# Split selected genes into overexpressed (upregulated) and underexpressed (downregulated) based on log2FoldChange sign
selected_genes_p8_over <- res_df_p8_sig %>% filter(log2FoldChange > 0)
selected_genes_p8_under <- res_df_p8_sig %>% filter(log2FoldChange < 0)

selected_genes_p15_over <- res_df_p15_sig %>% filter(log2FoldChange > 0)
selected_genes_p15_under <- res_df_p15_sig %>% filter(log2FoldChange < 0)

selected_genes_p8_vs_p15_over <- res_p8_vs_p15_sig %>% filter(log2FoldChange > 0)
selected_genes_p8_vs_p15_under <- res_p8_vs_p15_sig %>% filter(log2FoldChange < 0)

# Compute GO for BP
GO_results1_ctrl_vs_p8_over_BP <- enrichGO(gene = selected_genes_p8_over$ENSEMBL_short,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP")

GO_results1_ctrl_vs_p8_under_BP <- enrichGO(gene = selected_genes_p8_under$ENSEMBL_short,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP")

GO_results2_ctrl_vs_p15_over_BP <- enrichGO(gene = selected_genes_p15_over$ENSEMBL_short,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENSEMBL",
                        ont = "BP")
GO_results2_ctrl_vs_p15_under_BP <- enrichGO(gene = selected_genes_p15_under$ENSEMBL_short,
                       OrgDb = "org.Hs.eg.db",
                       keyType = "ENSEMBL",
                       ont = "BP")                        

GO_results3_p8_vs_p15_over_BP <- enrichGO(gene = selected_genes_p8_vs_p15_over$ENSEMBL_short,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENSEMBL",
                        ont = "BP")
GO_results3_p8_vs_p15_under_BP <- enrichGO(gene = selected_genes_p8_vs_p15_under$ENSEMBL_short,
                        OrgDb = "org.Hs.eg.db",
                        keyType = "ENSEMBL",
                        ont = "BP")

# Save results
# write.csv(GO_results1_ctrl_vs_p8_over_BP, file = paste0(root_dir, "/08_ontologies/ctrl_vs_p8_GO_up_res1.csv"), row.names = FALSE)
# write.csv(GO_results1_ctrl_vs_p8_under_BP, file = paste0(root_dir, "/08_ontologies/ctrl_vs_p8_GO_down_res1.csv"), row.names = FALSE)
# write.csv(GO_results2_ctrl_vs_p15_over_BP, file = paste0(root_dir, "/08_ontologies/ctrl_vs_p15_GO_up_res1.csv"), row.names = FALSE)
# write.csv(GO_results2_ctrl_vs_p15_under_BP, file = paste0(root_dir, "/08_ontologies/ctrl_vs_p15_GO_down_res1.csv"), row.names = FALSE)
# write.csv(GO_results3_p8_vs_p15_over_BP, file = paste0(root_dir, "/08_ontologies/p15_vs_p8_GO_up_res1.csv"), row.names = FALSE)
# write.csv(GO_results3_p8_vs_p15_under_BP, file = paste0(root_dir, "/08_ontologies/p15_vs_p8_GO_down_res1.csv"), row.names = FALSE)

# ----------------------
# ------- Plots --------
# ----------------------

# Illustrate top 20 ontologies
pdf(paste0(root_dir, "/07_plots/GO_BP_ctrl_vs_p8_up_barplot.pdf"), height = 12, width = 8)
barplot(GO_results1_ctrl_vs_p8_over_BP, showCategory = 20) +
        labs(title = "GO Enrichment Analysis: Control vs p8",
             subtitle = "Overexpressed Genes") +
        theme(plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(face = "italic"))
dev.off()

pdf(paste0(root_dir, "/07_plots/GO_BP_ctrl_vs_p8_down_barplot.pdf"), height = 12, width = 8)
barplot(GO_results1_ctrl_vs_p8_under_BP, showCategory = 20) +
        labs(title = "GO Enrichment Analysis: Control vs p8",
             subtitle = "Underexpressed Genes") +
        theme(plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(face = "italic"))
dev.off()

pdf(paste0(root_dir, "/07_plots/GO_BP_ctrl_vs_p15_up_barplot.pdf"), height = 12, width = 8)
barplot(GO_results2_ctrl_vs_p15_over_BP, showCategory = 20) +
        labs(title = "GO Enrichment Analysis: Control vs p15",
             subtitle = "Overexpressed Genes") +
        theme(plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(face = "italic"))
dev.off()

pdf(paste0(root_dir, "/07_plots/GO_BP_ctrl_vs_p15_down_barplot.pdf"), height = 12, width = 8)
barplot(GO_results2_ctrl_vs_p15_under_BP, showCategory = 20) +
        labs(title = "GO Enrichment Analysis: Control vs p15",
             subtitle = "Underexpressed Genes") +
        theme(plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(face = "italic"))
dev.off()

pdf(paste0(root_dir, "/07_plots/GO_BP_p15_vs_p8_up_barplot.pdf"), height = 12, width = 8)
barplot(GO_results3_p8_vs_p15_over_BP, showCategory = 20) +
        labs(title = "GO Enrichment Analysis: p15 vs p8",
             subtitle = "Overexpressed Genes") +
        theme(plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(face = "italic"))
dev.off()

pdf(paste0(root_dir, "/07_plots/GO_BP_p15_vs_p8_down_barplot.pdf"), height = 12, width = 8)
barplot(GO_results3_p8_vs_p15_under_BP, showCategory = 20) +
        labs(title = "GO Enrichment Analysis: p15 vs p8",
             subtitle = "Underexpressed Genes") +
        theme(plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(face = "italic"))
dev.off()


