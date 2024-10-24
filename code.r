# Load necessary libraries
library(TCGAbiolinks)
library(DESeq2)
library(ggplot2)

# Download TCGA data
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Gene expression", 
                  data.type = "Gene expression quantification", 
                  workflow.type = "HTSeq - Counts")
GDCdownload(query)
data <- GDCprepare(query)

# Preprocess data
dds <- DESeqDataSetFromMatrix(countData = data$counts, 
                              colData = data$colData, 
                              design = ~ condition)
dds <- DESeq(dds)

# Extract results
results <- results(dds)

# Identify significant mutations
sig_results <- results[which(results$padj < 0.05), ]

# Visualize expression changes
ggplot(data = as.data.frame(sig_results), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  theme_minimal() +
  labs(title = "Significant Expression Changes in TCGA-BRCA",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value")

# Save significant results
write.csv(as.data.frame(sig_results), file = "significant_mutations.csv")
