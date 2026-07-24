# DESeq2 analysis: ALIVE vs PARALYZED at 18 nM ivermectin
# Species: Anopheles dirus
# Output: significant tolerance genes, top 40 genes, volcano plot

# ------------------------------------------------------------
# 1) Setup
# ------------------------------------------------------------

library(DESeq2)
library(ggplot2)

folder <- "C:/Users/19198826/Downloads/RotationProject2/Sequencing/DSEQ Data"
setwd(folder)

# Optional: confirm location and files
getwd()
list.files()

# ------------------------------------------------------------
# 2) Import data
# ------------------------------------------------------------

# featureCounts output
fc <- read.delim("Adirus_gene_counts.txt", comment.char = "#", check.names = FALSE)

# Metadata table
metadata <- read.csv("an_dirus_deseq_metadata_with_group.csv", row.names = 1)

# ------------------------------------------------------------
# 3) Prepare metadata
# ------------------------------------------------------------

metadata$batch <- as.factor(metadata$batch)
metadata$phenotype <- as.factor(metadata$phenotype)
metadata$concentration_nM <- as.numeric(as.character(metadata$concentration_nM))

# Keep only Batch2 samples exposed to 18 nM ivermectin
keep <- metadata$batch == "Batch2" & metadata$concentration_nM == 18
metadata_sub <- metadata[keep, , drop = FALSE]

# ------------------------------------------------------------
# 4) Prepare count matrix
# ------------------------------------------------------------

# featureCounts stores annotation columns in the first 6 columns
counts <- fc[, 7:ncol(fc)]
rownames(counts) <- fc$Geneid
counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

# Make sample names cleaner if needed
colnames(counts) <- basename(colnames(counts))
colnames(counts) <- sub("\\.sorted\\.bam$", "", colnames(counts))
colnames(counts) <- sub("\\.bam$", "", colnames(counts))

# Keep only the selected samples
counts_sub <- counts[, rownames(metadata_sub), drop = FALSE]

# Confirm sample order matches
stopifnot(all(colnames(counts_sub) == rownames(metadata_sub)))

# ------------------------------------------------------------
# 5) Differential expression analysis
# ------------------------------------------------------------

# Set PARALYZED as the reference so results are interpreted as ALIVE vs PARALYZED
metadata_sub$phenotype <- relevel(factor(metadata_sub$phenotype), ref = "PARALYZED")

# Build DESeq2 object
DDS <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = metadata_sub,
  design = ~ phenotype
)

# Remove low-count genes (total counts across selected samples < 10)
DDS <- DDS[rowSums(counts(DDS)) >= 10, ]

# Run DESeq2
DDS <- DESeq(DDS)

# Extract results
res <- results(DDS)
res <- res[order(res$padj), ]

# Save all results
write.csv(as.data.frame(res), "DESeq2_All_Results_Alive_vs_Paralyzed.csv", row.names = TRUE)

# ------------------------------------------------------------
# 6) Tolerance genes (higher in ALIVE)
# ------------------------------------------------------------

# Significant genes with higher expression in ALIVE
# padj < 0.05 and log2FoldChange > 1

tolerance_genes <- subset(
  as.data.frame(res),
  padj < 0.05 & log2FoldChange > 1
)

# Order by significance
if (nrow(tolerance_genes) > 0) {
  tolerance_genes <- tolerance_genes[order(tolerance_genes$padj), ]
}

paralyzed_genes <- subset(
  as.data.frame(res),
  padj < 0.05 & log2FoldChange < -1)
# Save all tolerance genes
write.csv(
  tolerance_genes,
  "Alive_vs_Paralyzed_Tolerance_Genes.csv",
  row.names = TRUE
)

# Save all Paralysed genes
write.csv(paralyzed_genes, "Paralyzed_genes.csv", row.names = TRUE)

# Save top 40 tolerance genes
top40_tolerance <- head(tolerance_genes, 40)
write.csv(
  top40_tolerance,
  "Top40_Tolerance_Genes_Alive_vs_Paralyzed.csv",
  row.names = TRUE
)

# Save top 40 Paralysed genes
top40_paralyzed <- head(paralyzed_genes[order(paralyzed_genes$padj), ], 40)

write.csv(
  top40_paralyzed,
  "Top40_Paralyzed.csv",
  row.names = TRUE
)
# ------------------------------------------------------------
# 7) Volcano plot
# ------------------------------------------------------------

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

res_df$significant <- "Not significant"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "Higher in ALIVE"
res_df$significant[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "Higher in PARALYZED"

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(
    "Not significant" = "grey70",
    "Higher in ALIVE" = "red",
    "Higher in PARALYZED" = "blue"
  )) +
  labs(
    title = "Differentially Expressed Genes (IVM 18 nM): Alive vs Paralysed An. dirus",
    x = "log2 fold change (ALIVE vs PARALYZED)",
    y = "-log10(adjusted p-value)",
    color = NULL
  ) +
  theme_minimal()

print(volcano_plot)

# Save plot files

ggsave(
  filename = "Differentially_Expressed_Genes_IVM_18nM_Alive_vs_Paralysed_An_dirus.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "Differentially_Expressed_Genes_IVM_18nM_Alive_vs_Paralysed_An_dirus.pdf",
  plot = volcano_plot,
  width = 8,
  height = 6
)
write.csv(as.data.frame(res),
          "DESeq2_All_Results.csv")
save.image("DESeq2_Alive_vs_Paralyzed_Workspace.RData")
