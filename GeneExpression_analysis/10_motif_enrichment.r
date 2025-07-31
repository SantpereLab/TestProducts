#############################################
######### Motif Enrichment Analysis #########
#############################################

# Load libraries
library(RcisTarget)
library(htmlwidgets)
library(data.table)
library(dplyr)
library(VennDiagram)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(DT)

# Set root directory
root_dir <- "/path/to/the/project"

# ---------------------------------
# ------ RcisTarget pipeline ------
# ---------------------------------

# Load annotation
data(motifAnnotations_hgnc)
motifAnnotations[199:202,]

# Load motif rankings
motifRankings <- importRankings(paste0(root_dir, "/99_general/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))

# Load gene sets
gene_id <- "GeneSymbol" # ENSEMBL_short can be used as an alternative 

# P8 gene set
p8_data <- read.csv(paste0(root_dir,"/06_fc/ctrl_vs_p8_DESeq2_significant.csv"))
nrow(p8_data)
p8_up <- p8_data[p8_data$log2FoldChange > 0,]
length(p8_up)
p8_up <- p8_up[, gene_id]
p8_up <- p8_up[!is.na(p8_up)]

p8_down <- p8_data[p8_data$log2FoldChange < 0,]
length(p8_down)
p8_down <- p8_down[, gene_id]
p8_down <- p8_down[!is.na(p8_down)]

# P15 gene set
p15_data <- read.csv(paste0(root_dir, "/06_fc/ctrl_vs_p15_DESeq2_significant.csv"))
nrow(p15_data)
p15_up <- p15_data[p15_data$log2FoldChange > 0 ,]
length(p15_up)
p15_up <- p15_up[, gene_id]
p15_up <- p15_up[!is.na(p15_up)] 

p15_down<- p15_data[p15_data$log2FoldChange < 0 ,]
length(p15_down)
p15_down <- p15_down[, gene_id]
p15_down <- p15_down[!is.na(p15_down)] 


# Check for overlaps between p15 and p8 gene sets
overlap_genes_up <- intersect(p8_up, p15_up)
overlap_genes_down <- intersect(p8_down, p15_down)
inversed_genes <- c(intersect(p8_down, p15_up), intersect(p8_up, p15_down))

length(overlap_genes_up)
length(overlap_genes_down)
length(inversed_genes)

overlap_genes_up
overlap_genes_down
inversed_genes

# Take labels to make venn diagram
intersect(union(p8_up, p8_down), union(p15_up, p15_down))
length(intersect(union(p8_up, p8_down), union(p15_up, p15_down)))
length(union(p8_up, p8_down))
length(union(p15_up, p15_down))

# Create gene sets list
geneSets <- list(
  p8_up = p8_up,
  p15_up = p15_up,
  p8_down = p8_down,
  p15_down = p15_down
) 

# Run motif enrichment analysis
motifEnrichment <- cisTarget(
  geneSets = geneSets,
  motifRankings = motifRankings,
  motifAnnot = motifAnnotations,
  nCores = 8,
  verbose = TRUE
)

############## Save results ##############

# Add logos to the motif enrichment results
motifEnrichment_wLogo <- addLogo(motifEnrichment)

# Export datatable as HTML
saveWidget(
  datatable(motifEnrichment_wLogo[,-c("TF_lowConf"), with=FALSE], 
            escape = FALSE, 
            filter="top", options=list(pageLength=5)),
  file = paste0(root_dir, "/09_motifEnrich/files/motif_enrichment_table_full_updown.html"),
  selfcontained = FALSE
)

# --------------------------------------------------------
# -------- Further Analysis: Search of Regulons ----------
# --------------------------------------------------------


# Extract the unique list of high-confidence TFs for each enriched set.
get_enriched_tfs <- function(setName, enrichment_df) {
  # Filter for the specific gene set
  subset_df <- enrichment_df[enrichment_df$geneSet == setName, ]
  # Get all TFs from the high-confidence column
  tfs_raw_strings <- subset_df$TF_highConf
  # Split strings that have multiple TFs (e.g., "SOX9;SOX10")
  tfs_no_annot <- gsub(" \\(.*", "", tfs_raw_strings)
  tfs_list <- strsplit(tfs_no_annot, ";\\s*")
  all_tfs <- unlist(tfs_list)
  # Return a single, unique vector of TF names
  return(unique(all_tfs[all_tfs != ""]))
}

enriched_tfs_in_p8_up <- get_enriched_tfs("p8_up", motifEnrichment)
enriched_tfs_in_p8_down <- get_enriched_tfs("p8_down", motifEnrichment)
enriched_tfs_in_p15_up <- get_enriched_tfs("p15_up", motifEnrichment)
enriched_tfs_in_p15_down <- get_enriched_tfs("p15_down", motifEnrichment)


# Identify candidate activators and repressor regulons by cross-referencing
p8_motif_up.p8_up <- intersect(enriched_tfs_in_p8_up, p8_up)
p8_motif_down.p8_up<- intersect(enriched_tfs_in_p8_down, p8_up)
p8_motif_up.p8_down <- intersect(enriched_tfs_in_p8_up, p8_down)
p8_motif_down.p8_down <- intersect(enriched_tfs_in_p8_down, p8_down)

p8_motif_up.p8_up
p8_motif_down.p8_up
p8_motif_up.p8_down
p8_motif_down.p8_down

p15_motif_up.p15_up <- intersect(enriched_tfs_in_p15_up, p15_up)
p15_motif_down.p15_up<- intersect(enriched_tfs_in_p15_down, p15_up)
p15_motif_up.p15_down <- intersect(enriched_tfs_in_p15_up, p15_down)
p15_motif_down.p15_down <- intersect(enriched_tfs_in_p15_down, p15_down)

p15_motif_up.p15_up
p15_motif_down.p15_up
p15_motif_up.p15_down
p15_motif_down.p15_down

# ------------------------ Check all the genes regulated by candidate TF ---------------------- #

extract_enriched_genes <- function(tf_name, gene_set_name, cisTarget_results) {
  
  # Create a pattern to match the whole word for the TF to avoid partial matches (e.g., "NFE2" vs "NFE2L1")
  tf_pattern <- paste0("\\b", tf_name, "\\b")
  
  # Filter rows for the specified gene set and TF
  filtered_results <- cisTarget_results[
    cisTarget_results$geneSet == gene_set_name & grepl(tf_pattern, cisTarget_results$TF_highConf),
  ]
  
  # Get the strings from the 'enrichedGenes' column (e.g., "GENE1;GENE2")
  gene_strings <- filtered_results$enrichedGenes
  
  # Split each string by the semicolon into a list of character vectors
  list_of_gene_vectors <- strsplit(gene_strings, ";")
  
  # Collapse the list into a single vector, find unique names, and sort them
  all_genes <- unlist(list_of_gene_vectors)
  unique_sorted_genes <- sort(unique(all_genes))
  
  return(unique_sorted_genes)
}

# Call the function to get the list of genes for NFE2, FOSL2 and JUN from p15_up
nfe2_regulated_genes <- extract_enriched_genes(
  tf_name = "NFE2",
  gene_set_name = "p15_up",
  cisTarget_results = motifEnrichment
)

fosl2_regulated_genes <- extract_enriched_genes(
  tf_name = "FOSL2",
  gene_set_name = "p15_up",
  cisTarget_results = motifEnrichment
)

jun_regulated_genes <- extract_enriched_genes(
  tf_name = "JUN",
  gene_set_name = "p15_up",
  cisTarget_results = motifEnrichment
)

nfe2_regulated_genes
fosl2_regulated_genes
jun_regulated_genes