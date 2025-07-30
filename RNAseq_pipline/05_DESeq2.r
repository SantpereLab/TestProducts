############################################
##### Differential Expression Analysis #####
############################################

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

############## Data preparation ##############

# Define file paths
phenodata_path <- paste0(root_dir, "/phenodata.csv")
counts_dir <- paste0(root_dir, "/05_counts/strand2")

# Read the phenotype table
phenodata <- read.csv(phenodata_path, header = TRUE, stringsAsFactors = TRUE)

# Ensure replicate and condition is a factor
phenodata$Replicate <- factor(phenodata$Replicate)
phenodata$Condition <- factor(phenodata$Condition)

# Extract IDs 
ids <- phenodata$name

# Initialize an empty list to store data
count_data_list <- list()

# Iterate over ids and find corresponding count files
for (id in ids) {
    print(paste("Processing sample:", id))
    
    count_file <- list.files(counts_dir, pattern = paste0(id, ".*txt$"), full.names = TRUE)
    
    if (length(count_file) == 1) {  # Ensure a single match
        print(paste("Found count file:", count_file)) 
        
        count_table <- read.delim(count_file, header = TRUE, stringsAsFactors = FALSE, skip = 1)
        
        # Extract relevant columns: Geneid and counts
        count_data <- count_table[, c(1, 7)]
        colnames(count_data) <- c("Geneid", id)
        
        # Store data
        count_data_list[[id]] <- count_data
    } else {
        print(paste("Count file not found or multiple files found for SRR code:", id))
    }
}

# Merge all count tables by Geneid
counts <- Reduce(function(x, y) full_join(x, y, by = "Geneid"), count_data_list)

# Set Geneid as row names
rownames(counts) <- counts$Geneid
counts <- counts[, -1]  # Remove Geneid column

# Generate deseq object with counts and condition information
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = phenodata,
                              design = ~Condition)

############## Data visualization ##############

vsd <- vst(dds, blind = TRUE)

# Export PCA plot to PDF
pdf("PCA_plot.pdf")
plotPCA(vsd, intgroup = c("Condition"))
dev.off()


############## Deseq Analysis ##############

# Filter genes with less than 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds_modified <- dds[keep,]

# Set uninduced as reference level
dds_modified$Condition <- relevel(dds_modified$Condition, ref = "control")

# Run DESeq and save results
dds_modified <- DESeq(dds_modified)
resultsNames(dds_modified)
res_p8 <- results(dds_modified, name = "Condition_p8_vs_control")
res_p15 <- results(dds_modified, name = "Condition_p15_vs_control")
res_p8_vs_p15 <- results(dds_modified, contrast = c("Condition", "p8", "p15"))

# Extract log2FoldChange and ENSEMBL identifier
res_df_p8 <- as.data.frame(res_p8)
res_df_p15 <- as.data.frame(res_p15)
res_p8_vs_p15 <- as.data.frame(res_p8_vs_p15)

# Remove version number from gene names (to extract geneNames)
res_df_p8$ENSEMBL <- rownames(res_df_p8)
res_df_p15$ENSEMBL <- rownames(res_df_p15)
res_p8_vs_p15$ENSEMBL <- rownames(res_p8_vs_p15)
res_df_p8$ENSEMBL_short <- gsub("\\..*", "",row.names(res_df_p8))
res_df_p15$ENSEMBL_short <- gsub("\\..*", "",row.names(res_df_p15))
res_p8_vs_p15$ENSEMBL_short <- gsub("\\..*", "",row.names(res_p8_vs_p15))

# Map ENSEMBL IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = res_df_p8$ENSEMBL_short, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Add gene symbols to the results data frames
res_df_p8$GeneSymbol <- gene_symbols
res_df_p15$GeneSymbol <- gene_symbols
res_p8_vs_p15$GeneSymbol <- gene_symbols

# Reorder columns to place GeneSymbol first
res_df_p8 <- res_df_p8[, c("ENSEMBL", "ENSEMBL_short", "GeneSymbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")]
res_df_p15 <- res_df_p15[, c("ENSEMBL", "ENSEMBL_short", "GeneSymbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")]
res_p8_vs_p15 <- res_p8_vs_p15[, c("ENSEMBL", "ENSEMBL_short", "GeneSymbol", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")]

# Remove version number from gene symbols
rownames(res_df_p8) <- NULL
rownames(res_df_p15) <- NULL
rownames(res_p8_vs_p15) <- NULL

# Sort files by padj value desc
res_df_p8 <- res_df_p8 %>% arrange(padj)
res_df_p15 <- res_df_p15 %>% arrange(padj)
res_p8_vs_p15 <- res_p8_vs_p15 %>% arrange(padj)

write.csv(res_df_p8, file = paste0(root_dir, "/06_fc/ctrl_vs_p8_DESeq2_results1.csv"), row.names = FALSE)
write.csv(res_df_p15, file = paste0(root_dir, "/06_fc/ctrl_vs_p15_DESeq2_results1.csv"), row.names = FALSE)
write.csv(res_p8_vs_p15, file = paste0(root_dir, "/06_fc/p8_vs_p15_DESeq2_results1.csv"), row.names = FALSE)

sum(is.na(res_df_p8$padj))
sum(is.na(res_df_p8$pvalue))
dim(res_df_p8)

head(res_df_p8, 20)

# Take only signiifcant genes
res_df_p8_sig <- res_df_p8 %>% filter(padj < 0.05)
res_df_p15_sig <- res_df_p15 %>% filter(padj < 0.05)
res_p8_vs_p15_sig <- res_p8_vs_p15 %>% filter(padj < 0.05)

# Save results
write.csv(res_df_p8_sig, file = paste0(root_dir, "/06_fc/ctrl_vs_p8_DESeq2_significant.csv"), row.names = FALSE)
write.csv(res_df_p15_sig, file = paste0(root_dir, "/06_fc/ctrl_vs_p15_DESeq2_significant.csv"), row.names = FALSE)
write.csv(res_p8_vs_p15_sig, file = paste0(root_dir, "/06_fc/p8_vs_p15_DESeq2_significant.csv"), row.names = FALSE)