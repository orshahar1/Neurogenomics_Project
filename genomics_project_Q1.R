
### Genomics project part 1

################################################################
#Or Shahar 206582017
#Daniel Brooker 315015594
################################################################

# This R script performs differential gene expression analysis on RNA-seq data
# The data was previously processed using Kallisto for transcript quantification
# Kallisto was run in Python to generate abundance estimates for each sample
# The Kallisto output files are used as input for this analysis

# Install  packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tximport")
BiocManager::install("readr")
BiocManager::install("biomaRt")

# Load  libraries
library(tximport)
library(readr)
library(GenomicFeatures)
library(biomaRt)

# Load the list of Kallisto output files
files <- c("C:/Users/danie/OneDrive/מסמכים/R_programing/C1_output.tsv",
           "C:/Users/danie/OneDrive/מסמכים/R_programing/C2_output.tsv",
           "C:/Users/danie/OneDrive/מסמכים/R_programing/C3_output.tsv",
           "C:/Users/danie/OneDrive/מסמכים/R_programing/KO1_output.tsv",
           "C:/Users/danie/OneDrive/מסמכים/R_programing/KO2_output.tsv",
           "C:/Users/danie/OneDrive/מסמכים/R_programing/KO3_output.tsv")

# Path to the GTF file
gtf_file <- "C:/Users/danie/OneDrive/מסמכים/R_programing/Mus_musculuss.GRCm39.112.gtf"

# Create a TxDb object from the GTF file
# This step creates a transcript database from the GTF file
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

# Create a mapping between transcripts and genes
# Purpose: To establish a relationship between individual transcripts and their corresponding genes
# Importance: Essential for aggregating transcript-level data to gene-level
tx2gene <- select(txdb, keys(txdb, "TXNAME"), "GENEID", "TXNAME")

# Import the Kallisto data using tximport
# Purpose: To convert transcript-level abundance estimates from Kallisto to gene-level count data
# Importance: 
#   1. Allows for gene-level analysis
#   2. Prepares data in a format suitable for downstream differential expression analysis  in DESeq2
#   3. Accounts for uncertainty in transcript quantification when summarizing to gene-level
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, txOut = FALSE, ignoreTxVersion = TRUE)

# Edit column names 
colnames(txi$counts) <- c("C1", "C2", "C3", "KO1", "KO2", "KO3")

# Check and display the first few rows of the count data
head(txi$counts)

# Load DESeq2 library for differential expression analysis
library(DESeq2)

# Create a data frame with condition information
# This matches the experimental design: 3 control samples and 3 knockout samples
coldata <- data.frame (condition = factor(c("control", "control", "control", "knockout", "knockout", "knockout")))

# Round the count values to integers
# DESeq2 requires integer count data
txi$counts <- round(txi$counts)

# Create a DESeqDataSet object
# This combines the count data with the experimental design
dds <- DESeqDataSetFromMatrix(countData = txi$counts, colData = coldata, design = ~ condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

# Extract the results of differential expression analysis
res <- results(dds)
head(res)

# Sort results by log2FoldChange
res_ordered <- res[order(res$log2FoldChange, decreasing = TRUE), ]

# Extract Ensemble IDs after sorting
ensembl_ids <- rownames(res_ordered)

# Connect to Ensemble database for mouse
# This allows us to map Ensemble IDs to gene symbols for easier interpretation
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Map Ensemble IDs to gene symbols
gene_mapping <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = 'ensembl_gene_id', values = ensembl_ids, mart = ensembl)

# Display the mapping
head(gene_mapping)

# Add gene names to the results table
res_ordered$gene_name <- gene_mapping$external_gene_name[match(rownames(res_ordered), gene_mapping$ensembl_gene_id)]

# Display updated results
head(res_ordered)

# Filter for significantly differential expressed genes (p-adj < 0.05)
significant_genes <- subset(res_ordered, padj < 0.05)

# Save the list of significant genes to a file
write.table(significant_genes, file="significant_genes_with_names.txt", sep="\t", row.names=FALSE)

# Filter genes higher expressed in control
higher_in_control <- subset(significant_genes, log2FoldChange > 0)

# Filter genes higher expressed in knockout
higher_in_knockout <- subset(significant_genes, log2FoldChange < 0)

# Save genes higher expressed in control to a file
write.table(higher_in_control, file="higher_in_control.txt", sep="\t", row.names=FALSE)

# Save genes higher expressed in knockout to a file
write.table(higher_in_knockout, file="higher_in_knockout.txt", sep="\t", row.names=FALSE)