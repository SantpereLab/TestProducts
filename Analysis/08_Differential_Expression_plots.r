############################################
##### Plots from RNA-seq DE analysis 1 #####
############################################

# Load libraries
library(dplyr)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(ggplot2)
library(data.table)
library(ggrepel)
library(colorspace)
library(RColorBrewer)
library(scales)

#-------------------------------------------------------------------------
#--------------- Volcano of Differentially Expressed Genes ---------------
#-------------------------------------------------------------------------

# Set root directory
root_dir <- "/path/to/the/project"

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

#------------------------------------------------------------------------
#--------------- Volcano Segregating by Ciliary Function  ---------------
#------------------------------------------------------------------------

# Load functional categories
gene_info <- read_tsv(paste0(root_dir, "/99_general/Reiter2017_BigCategories.txt"))
head(gene_info)

# Merge DESeq2 results with your gene_info
annotated_genes_p8_vs_ctrl <- res_p8_vs_ctrl %>%
  inner_join(gene_info, by = c("GeneSymbol" = "Gene"))
annotated_genes_p15_vs_ctrl <- res_p15_vs_ctrl %>%
  inner_join(gene_info, by = c("GeneSymbol" = "Gene"))
annotated_genes_p8_vs_p15 <- res_p8_vs_p15 %>%
  inner_join(gene_info, by = c("GeneSymbol" = "Gene"))

# Generat volcano function
create_volcano_plot_selected_genes_grouped <- function(
    annotated_df,
    title,
    padj_threshold = 0.05,
    log2fc_threshold = 0.5,
    custom_group_colors = NULL,
    point_size = 3.5,
    label_size = 3.2,
    label_max_overlaps = 20,
    max_labels = NULL
) {
    required_cols <- c("GeneSymbol", "log2FoldChange", "padj", "Group")
    if (!all(required_cols %in% names(annotated_df))) {
        stop("Missing required columns: ", paste(setdiff(required_cols, names(annotated_df)), collapse = ", "))
    }

    plot_df <- annotated_df %>%
        filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
        mutate(
            is_significant = padj < padj_threshold,
            Group = factor(Group)
        )

    if (nrow(plot_df) == 0) {
        return(ggplot() + labs(title = title, subtitle = "No data to plot") + theme_void())
    }

    # Prepare labels data frame
    labels_df <- plot_df %>%
        arrange(padj, desc(abs(log2FoldChange)))  # Order by significance and effect size

    if (!is.null(max_labels)) {
        labels_df <- labels_df %>% head(max_labels)  # Limit to max_labels
    }

    unique_groups <- levels(plot_df$Group)

    base_colors <- custom_group_colors
    if (is.null(base_colors)) {
        palette_name <- "Set2"
        base_colors <- setNames(
            if (length(unique_groups) <= brewer.pal.info[palette_name, "maxcolors"]) {
                brewer.pal(length(unique_groups), palette_name)
            } else {
                scales::hue_pal()(length(unique_groups))
            },
            unique_groups
        )
    } else {
        missing_groups <- setdiff(unique_groups, names(base_colors))
        if (length(missing_groups) > 0) {
            base_colors <- c(base_colors, setNames(rep("grey50", length(missing_groups)), missing_groups))
        }
    }

    plot_df <- plot_df %>%
        mutate(
            color = base_colors[as.character(Group)],
            alpha = ifelse(is_significant, 1, 0.3)
        )

    p <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(
            aes(fill = color, alpha = alpha, shape = Group),
            size = point_size
        ) +
        scale_fill_identity() +
        scale_alpha_identity() +
        scale_shape_manual(name = "Group", values = rep(21:25, length.out = length(unique_groups))) +
        geom_text_repel(
            data = labels_df,
            aes(label = GeneSymbol),
            size = label_size, max.overlaps = label_max_overlaps,
            box.padding = 0.35, point.padding = 0.3,
            segment.color = 'grey50', segment.size = 0.3,
            force = 1.5
        ) +
        geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey30", linewidth = 0.6) +
        geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "grey30", linewidth = 0.6) +
        labs(
            title = title,
            x = bquote(~Log[2]~ "Fold Change"),
            y = bquote(-~Log[10]~ "(Adjusted P-value)")
        ) +
        theme_bw(base_size = 12) +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            axis.title = element_text(size = 14, face = "bold"),
            axis.text = element_text(size = 11),
            legend.position = "bottom",
            legend.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 10),
            panel.grid.major = element_line(color = "grey90"),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7)
        ) +
        guides(shape = guide_legend(
            override.aes = list(
                fill = base_colors[levels(plot_df$Group)],
                alpha = 1, size = 4
            )
        ))

    return(p)
}

########### Execution of Volcano ###############

# Define custom group colors
my_group_colors <- c(
  "Ciliary trafficking" = "#E41A1C",
  "Ciliogenesis" = "#4DAF4A",
  "Motile cilium structure" = "#377EB8",
  "Non-motile cilium structure" = "#984EA3"
)

# Create and save plot
pdf(paste0(root_dir, "/07_plots/volcano_ctrl_vs_p8_specific_enhanced_final_prova.pdf"), width = 8, height = 12)
print(create_volcano_plot_selected_genes_grouped(
  annotated_df = annotated_genes_p8_vs_ctrl,
  title = "Differential Expression of Cilia-Related Genes: Control vs p8",
  custom_group_colors = my_group_colors
))
dev.off()

# Create and save plot
pdf(paste0(root_dir, "/07_plots/volcano_ctrl_vs_p15_specific_enhanced_final.pdf"), width = 8, height = 12)
print(create_volcano_plot_selected_genes_grouped(
  annotated_df = annotated_genes_p15_vs_ctrl,
  title = "Differential Expression of Cilia-Related Genes: Control vs p15",
  custom_group_colors = my_group_colors,
  max_labels = 30
))
dev.off()

# Create and save plot
pdf(paste0(root_dir, "/07_plots/volcano_p8_vs_p15_specific_enhanced_final.pdf"), width = 8, height = 12)
print(create_volcano_plot_selected_genes_grouped(
  annotated_df = annotated_genes_p8_vs_p15,
  title = "Differential Expression of Cilia-Related Genes: p15 vs p8",
  custom_group_colors = my_group_colors
))
dev.off()
