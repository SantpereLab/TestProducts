############################################
##### Plots from RNA-seq DE analysis 2 #####
############################################

# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(data.table)
library(ggrepel)
library(colorspace)
library(RColorBrewer)
library(scales)
library(stringr)s
library(tidyverse)
library(GO.db) 

#---------------------------------------------
#--------------- Load DEG data ---------------
#---------------------------------------------

# Set root directory  
root_dir <- "/path/to/the/project"

# Load data
res_p8_vs_ctrl <- read.csv(paste0(root_dir, "/06_fc/ctrl_vs_p8_DESeq2_results1.csv"), header = TRUE, stringsAsFactors = TRUE)
res_p15_vs_ctrl <- read.csv(paste0(root_dir, "/06_fc/ctrl_vs_p15_DESeq2_results1.csv"), header = TRUE, stringsAsFactors = TRUE)

# Load functional categories
gene_info <- read_tsv(paste0(root_dir, "/99_general/Reiter2017_BigCategories.txt"))

#-----------------------------------------------------
#------- Add Ciliary and Ontology Function  ----------
#-----------------------------------------------------

# Get a unique list of all ENSEMBL IDs from the results
all_gene_ids <- as.character(unique(res_p8_vs_ctrl$ENSEMBL_short))

# Get the gene-to-GO mappings
gene_to_go <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_gene_ids,
  columns = c("GO", "ONTOLOGY"),
  keytype = "ENSEMBL"
)

# Remove any rows where the GO annotation is missing
gene_to_go <- na.omit(gene_to_go)
head(gene_to_go)

# Get all the unique GO IDs from our previous result
all_go_ids <- unique(gene_to_go$GO)

# Use GO.db to get the term descriptions
go_id_to_term <- AnnotationDbi::select(
  GO.db,
  keys = all_go_ids,
  columns = c("GOID", "TERM"),
  keytype = "GOID"
)
head(go_id_to_term)

# Join our two mapping tables to create a complete annotation set
full_go_annotation <- left_join(gene_to_go, go_id_to_term, by = c("GO" = "GOID"))

# This aggregates all GO terms for a single gene into one string
summarized_go <- full_go_annotation %>%
  group_by(ENSEMBL) %>%
  summarize(
    GO_IDs = paste(GO, collapse = "; "),
    GO_Terms = paste(TERM, collapse = "; ")
  )

# Join the summarized GO data with our DESeq2 results and custom annotations
annotated_genes_p8_vs_ctrl <- res_p8_vs_ctrl %>%
  left_join(gene_info, by = c("GeneSymbol" = "Gene")) %>%
  left_join(summarized_go, by = c("ENSEMBL_short" = "ENSEMBL"))
annotated_genes_p15_vs_ctrl <- res_p15_vs_ctrl  %>%
  left_join(gene_info, by = c("GeneSymbol" = "Gene")) %>%
  left_join(summarized_go, by = c("ENSEMBL_short" = "ENSEMBL"))

# View the final dataframe
glimpse(annotated_genes_p8_vs_ctrl)
glimpse(annotated_genes_p15_vs_ctrl)

# Save annotated results
write.csv(annotated_genes_p8_vs_ctrl,
          file = paste0(root_dir, "/ctrl_vs_p8_DESeq2_annotated.csv"),
          row.names = FALSE)
write.csv(annotated_genes_p15_vs_ctrl,
          file = paste0(root_dir, "/ctrl_vs_p15_DESeq2_annotated.csv"),
          row.names = FALSE)

# -------------------------------------
# ------------- Plotting --------------
# -------------------------------------

# Filter the p8 dataset
filtered_p8_data <- annotated_genes_p8_vs_ctrl %>%
  filter((!is.na(padj) & (padj < 0.05 | abs(log2FoldChange) > 0.5)))

# Filter the p15 dataset
filtered_p15_data <- annotated_genes_p15_vs_ctrl %>%
  filter((!is.na(padj) & (padj < 0.05 | abs(log2FoldChange) > 0.5)))

# Make a list of every individual category to plot
categories_to_plot <- c(
  "GO:0008017",
  "GO:0003777",
  "Ciliary trafficking",
  "Ciliogenesis",
  "Non-motile cilium structure",
  "Motile cilium structure"
)

# Set grid dimensions for plotting
GRID_ROWS <- 15
GRID_COLS <- 4
MAX_GENES <- GRID_ROWS * GRID_COLS

# Utility function
generate_filtered_tile_plot <- function(dataset, filter_term, condition_name) {

  # --- Part A: Filter for the specific term using your column names ---
  if (str_starts(filter_term, "GO:")) {
    # Search in the single GO_IDs column
    selected_genes <- dataset %>%
      filter(str_detect(replace_na(GO_IDs, ""), filter_term))
  } else {
    # Search in both Function and Group columns for text terms
    selected_genes <- dataset %>%
      filter(str_detect(replace_na(Function, ""), fixed(filter_term)) |
             str_detect(replace_na(Group, ""), fixed(filter_term)))
  }
  
  # --- Part B: Check if any genes were found ---
  if (nrow(selected_genes) == 0) {
    cat(paste0("No genes found for '", filter_term, "' in the '", condition_name, "' dataset.\n"))
    return(NULL)
  }
  
  # --- Part C: Prepare the data for ggplot (create the grid) ---
  if (nrow(selected_genes) > MAX_GENES) {
    cat(paste0("Warning: ", nrow(selected_genes), " genes found for '", filter_term, "'. Truncating to the top ", MAX_GENES, " by p-value.\n"))
    selected_genes <- selected_genes %>% arrange(padj) %>% slice_head(n = MAX_GENES)
  }
    
  # --- Part D: Create the plot ---
  # 1. Create the gene data with grid positions
  gene_positions <- selected_genes %>%
    distinct(GeneSymbol, .keep_all = TRUE) %>%
    filter(!is.na(GeneSymbol)) %>%
    arrange(GeneSymbol) %>%
    mutate(
      x_pos = ((row_number() - 1) %% GRID_COLS) + 1,
      y_pos = floor((row_number() - 1) / GRID_COLS) + 1
    )
    
  # 2. Create the full, empty grid canvas
  master_grid <- expand.grid(x_pos = 1:GRID_COLS, y_pos = 1:GRID_ROWS)
  
  # 3. Join the genes onto the master canvas
  plot_data <- left_join(master_grid, gene_positions, by = c("x_pos", "y_pos"))

  # --- Part E: Create the plot from the master grid data ---
  plot_title <- paste(condition_name, "\n(", filter_term, ")", sep = "")
  
  tile_plot <- ggplot(plot_data, aes(x = x_pos, y = y_pos, fill = log2FoldChange)) +
    # Draw all tiles; genes will be colored, empty cells will not be filled
    geom_tile(color = "white", linewidth = 1.5) + 
    geom_text(aes(label = GeneSymbol), color = "black", size = 3, na.rm = TRUE) +
    # Use blue-red color scale
    scale_fill_gradient2(
      low = "#2166AC",  # Blue
      mid = "#F7F7F7",  # Off-white/light grey
      high = "#B2182B", # Red
      midpoint = 0, 
      name = "log2 Fold Change",
      na.value = NA,
      limits = c(-2, 2) # Set fixed limits for color consistency
    ) +
    labs(title = plot_title) +
    scale_y_reverse() +
    theme_void() +
    theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      legend.position = "right"
    )
    
  return(tile_plot)
}


# Create a named list of your datasets for easy iteration
datasets_to_process <- list(
  "p8_vs_Ctrl" = filtered_p8_data,
  "p15_vs_Ctrl" = filtered_p15_data
)

# --- Main Loop to Generate All Plots ---
for (category in categories_to_plot) {
  for (condition_name in names(datasets_to_process)) {
    
    # Get the current dataset
    current_dataset <- datasets_to_process[[condition_name]]
    
    # Generate the plot for the current combination
    my_plot <- generate_filtered_tile_plot(
      dataset = current_dataset,
      filter_term = category,
      condition_name = condition_name
    )
    
    # Check if the plot was created (i.e., genes were found) and then print it
    if (!is.null(my_plot)) {
      print(my_plot)

      # Save each plot to a file
      safe_filename <- paste0("Plot_", condition_name, "_", str_replace_all(category, ":| ", "_"), ".pdf")
      ggsave(safe_filename, plot = my_plot, width = 8, height = 6, bg = "white")
      cat(paste0("Saved plot to '", safe_filename, "'\n"))
    }
  }
}
