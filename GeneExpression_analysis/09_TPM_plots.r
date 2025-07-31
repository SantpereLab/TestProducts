############################
##### TPM violin plots #####
############################

# Load libraries
library(dplyr)
library(purrr)
library(readr)
library(ggplot2)
library(data.table)
library(tidyr)
library(ggh4x)
library(RColorBrewer)

# ----------------------------------------
# ------------ Violin Plots --------------
# ----------------------------------------

# Set root directory
root_dir <- "/users/genomics/jmartinez/a_primary_cilia_project"

# Load the annotated TPM data
tpm_path <- paste0(root_dir, "/06_fc/tpm_annotated.csv")
tpm_annotated <- read.csv(tpm_path)
head(tpm_annotated)

# Calculate the log2(TPM + 1) and pivot_longer for visualization
tpm_long <- tpm_annotated %>%
  pivot_longer(cols = -c(geneid, length, SYMBOL),
               names_to = "Sample",
               values_to = "TPM") %>%
  mutate(log2_TPM = log2(TPM + 1))
head(tpm_long)

# Generate an extra column with condition information
tpm_long$condition <- gsub("RPE\\.([^.]+)\\..*", "\\1", tpm_long$Sample)

# Filer genes of interest
genes_of_interest_symbol <- c("NQO1", "ADGRG1",
                              "IL32", "IGFBP5", "RPL13AP25", "SBF2-AS1", "FST",
                              "PLK1", "CDCA3", "AURKA", "KIF20A")
genes_of_interest_ensembl <- "ENSG00000294512"
                   
tpm_long <- tpm_long %>%
  filter(SYMBOL %in% genes_of_interest_symbol | geneid %in% genes_of_interest_ensembl)


# Define the order of your categories and the genes within them
category_order <- c(
  "Upregulated in p8 and p15",
  "Downregulated in p8 and p15",
  "Inverse effect in p8 and p15"
)

gene_order <- c(
  # Upregulated
  "NQO1", "ADGRG1",
  # Downregulated
  "IL32", "IGFBP5", "RPL13AP25", "SBF2-AS1", "FST", "ENSG00000294512",
  # Inverse effect
  "PLK1", "CDCA3", "AURKA", "KIF20A"
)

# Prepare the data
tpm_plot_data <- tpm_long %>%

  # Create a label that uses SYMBOL or falls back to geneid
  mutate(facet_label = coalesce(SYMBOL, geneid)) %>%
  
  # Assign each gene to its category
  mutate(category = case_when(
    facet_label %in% c("NQO1", "ADGRG1") ~ "Upregulated in p8 and p15",
    facet_label %in% c("IL32", "IGFBP5", "RPL13AP25", "SBF2-AS1", "FST", "ENSG00000294512") ~ "Downregulated in p8 and p15",
    facet_label %in% c("PLK1", "CDCA3", "AURKA", "KIF20A") ~ "Inverse effect in p8 and p15",
    TRUE ~ "Other" # A fallback for any genes not in your lists
  )) %>%
  
  # Apply the desired ordering by converting to factors
  mutate(
    category = factor(category, levels = category_order),
    facet_label = factor(facet_label, levels = gene_order),
    condition = factor(condition, levels = c("DMSO", "p8", "p15"))
  )

# Plot
p <- ggplot(tpm_plot_data, aes(x = condition, y = log2_TPM, fill = condition)) +
  
  # Geoms for visualization
  geom_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.7, color = "black") +
  
  # Use facet_nested to create two levels of headers
  ggh4x::facet_nested(~ category + facet_label, scales = "free_y") +
  
  # Color palette
  scale_fill_brewer(palette = "Set2") +
  
  # Labels and Titles
  labs(
    title = "Gene Expression Patterns Across Conditions",
    x = "Experimental Condition",
    y = "log2(TPM + 1)",
    fill = "Condition"
  ) +
  
  # Theming
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom",
    
    # Custom styling for the nested facet strips
    ggh4x.facet.nestline = element_line(colour = "gray50"), # Line between nested strips
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold")
  )

pdf(paste0(root_dir, "/07_plots/tpms_shared_genes_p8_p15.pdf"), width = 20, height = 6)
print(p)
dev.off()