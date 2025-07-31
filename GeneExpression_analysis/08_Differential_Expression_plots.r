##########################################
##### Plots from RNA-seq DE analysis #####
##########################################

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

#--------------------------------------------------------------
#--------------- Differentially expressed genes ---------------
#--------------------------------------------------------------

# Set root directory
root_dir <- "/users/genomics/jmartinez/a_primary_cilia_project"

# Load data
res_p8_vs_ctrl <- read.csv(paste0(root_dir, "/06_fc/ctrl_vs_p8_DESeq2_results1.csv"), header = TRUE, stringsAsFactors = TRUE)
res_p15_vs_ctrl <- read.csv(paste0(root_dir, "/06_fc/ctrl_vs_p15_DESeq2_results1.csv"), header = TRUE, stringsAsFactors = TRUE)
res_p8_vs_p15 <- read.csv(paste0(root_dir, "/06_fc/p8_vs_p15_DESeq2_results1.csv"), header = TRUE, stringsAsFactors = TRUE)

# Generate enhanced volcano plots
create_volcano_plot_enhanced <- function(res_df, title, padj_threshold = 0.05, log2fc_threshold = 0.5) {
  if (!all(c("GeneSymbol", "log2FoldChange", "padj") %in% names(res_df))) 
    stop("Data frame must contain GeneSymbol, log2FoldChange, and padj.")

  res_df_filtered <- res_df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      Expression = factor(case_when(
        log2FoldChange > log2fc_threshold & padj < padj_threshold ~ "Up-regulated",
        log2FoldChange < -log2fc_threshold & padj < padj_threshold ~ "Down-regulated",
        TRUE ~ "Not significant"
      ), levels = c("Up-regulated", "Down-regulated", "Not significant"))
    )
  
  genes_to_label_df <- res_df_filtered %>%
  filter(Expression != "Not significant") %>%
  arrange(padj) %>%
  slice_head(n = 23)

  message(sprintf("For plot '%s': %d genes pass thresholds (padj < %.3f, |log2FC| > %.2f) and will be labeled.",
                  title, nrow(genes_to_label_df), padj_threshold, log2fc_threshold))

  ggplot(res_df_filtered, aes(log2FoldChange, -log10(padj), color = Expression, alpha = Expression != "Not significant")) +
    geom_point(size = 2) +
    scale_color_manual(values = c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8", "Not significant" = "grey70"),
                       name = "Gene Regulation") +
    scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none") +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black", linewidth = 0.5) +
    ggrepel::geom_text_repel(data = genes_to_label_df, aes(label = GeneSymbol),
                             size = 3.2, color = "black", box.padding = 0.4, point.padding = 0.3,
                             segment.color = "grey50", segment.size = 0.3, max.overlaps = Inf, force = 1) +
    labs(title = title,
         x = bquote(~Log[2]~ "Fold Change"),
         y = bquote(-~Log[10]~ "(Adjusted P-value)")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 11, color = "black"),
          legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))
}

# Generate and save volcano plots
pdf(paste0(root_dir, "/07_plots/volcano_ctrl_vs_p8.pdf"), height = 12, width = 8)
print(create_volcano_plot_enhanced(res_p8_vs_ctrl, "Differential Gene Expression: p8 vs Control"))
dev.off()

pdf(paste0(root_dir, "/07_plots/volcano_ctrl_vs_p15.pdf"), height = 12, width = 8)
print(create_volcano_plot_enhanced(res_p15_vs_ctrl, "Differential Gene Expression: p15 vs Control"))
dev.off()

pdf(paste0(root_dir, "/07_plots/volcano_p15_vs_p8.pdf"), height = 12, width = 8)
print(create_volcano_plot_enhanced(res_p8_vs_p15, "Differential Gene Expression: p8 vs p15"))
dev.off()
