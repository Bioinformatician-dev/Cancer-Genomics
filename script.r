# Load necessary libraries
library(TCGAbiolinks)
library(DESeq2)
library(ggplot2)
library(maftools)
library(dplyr)

# Step 1: Download TCGA Data
# Query for RNA-Seq data (gene expression)
query_rna <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "HTSeq - Counts")

# Query for mutation data
query_mutation <- GDCquery(project = "TCGA-BRCA",
                            data.category = "Simple Nucleotide Variation",
                            data.type = "Masked Somatic Mutation")

# Download the data
GDCdownload(query_rna)
GDCdownload(query_mutation)

# Prepare the data
rna_data <- GDCprepare(query_rna)
mutation_data <- GDCprepare(query_mutation)

# Step 2: Preprocess RNA-Seq Data
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = assay(rna_data),
                              colData = colData(rna_data),
                              design = ~ 1)  # Adjust design formula as needed

# Normalize the data
dds <- DESeq(dds)
normalized_counts <- assay(dds)

# Convert to a data frame for further analysis
normalized_counts_df <- as.data.frame(normalized_counts)

# Step 3: Analyze Mutations
# Convert mutation data to MAF format
maf_data <- read.maf(maf = mutation_data)

# Summary of MAF data
print(maf_data)

# Plot the summary of mutations
plotmafSummary(maf = maf_data, main = "Mutational Landscape")

# OncoPlot to visualize mutations across samples
oncoplot(maf = maf_data, top = 10)  # Show top 10 mutated genes

# Step 4: Identify Differentially Expressed Genes (DEGs)
# Assuming you have a column in colData indicating tumor vs normal
dds$condition <- factor(ifelse(grepl("01", dds$sample), "Tumor", "Normal"))

# Run DESeq2 analysis
dds <- DESeq(dds)
results <- results(dds)

# Extract significant results
sig_genes <- subset(results, padj < 0.05)

# View significant genes
head(sig_genes[order(sig_genes$pvalue), ])

# Step 5: Visualization of DEGs
# Create a volcano plot
ggplot(data = as.data.frame(results), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  xlim(c(-10, 10)) +
  ylim(c(0, 30)) +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")

# Step 6: Integrate Results
# Merge significant genes with mutation data
# Assuming 'gene' is a common column in both datasets
biomarkers <- merge(as.data.frame(sig_genes), maf_data@data, by = "gene", all.x = TRUE)

# Save biomarkers to a CSV file
write.csv(biomarkers, file = "potential_biomarkers.csv")

# Optional: Save normalized RNA-Seq data
write.csv(normalized_counts_df, file = "normalized_rna_counts.csv")
