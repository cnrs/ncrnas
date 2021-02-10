# @ charles.k.w
# charles.k.w@msn.com
# V2021.01.10

# Usage: Rscript diffExp.R currentWorkDir samplesheet genecount
# options(echo=TRUE) # if you want see commands in output file

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3){
	print("Usage: Rscript diffExp.R currentWorkDir samplesheet genecount")
	exit(1)
}

currentWorkDir  <- "/home/wangk/lab/shenlan/result/DAY7_vs_WT5"
genecountsFile  <- "CIRC.GENECOUNT.txt"
samplesheetFile <- "samplesheet.csv"
genecountsAnno  <- "GENECOUNTS.DAY7_vs_WT5.txt"

lfc_cutoff  <- log2 (1.5)  #0.584962501
pval_cutoff <- - log10(0.05) #1.301029996
padj_cutoff <- 0.05

setwd(currentWorkDir)

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("refGenome")

#library(refGenome)
#library(Rsubread)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(biomaRt)
library(ReactomePA)
library(DOSE)
library(KEGG.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
library(GO.db)
library(topGO)
library(dplyr)
library(gage)
library(ggsci)
library(AnnotationDbi)

organismDb <- org.Mm.eg.db

# - - - - - - - - - - - -
genecountsData <- read.table(genecountsFile, header = TRUE, sep = "\t", quote = "", row.names = 1)

head(genecountsData)

rownum <- nrow(genecountsData)
colnum <- ncol(genecountsData)

metasheetData <- read.table(samplesheetFile, header = TRUE, sep = ",", quote = "", row.names = NULL)
metasheetData <- metasheetData[match(colnames(genecountsData), metasheetData$Sample), ]

head (metasheetData)

# An example of the sample sheet file:
###  Sample,Group,Replicate,Condition
###  A_M6,M6,Rep1,treatment
###  B_M6,M6,Rep2,treatment
###  A_ME1,ME1,Rep1,control
###  B_ME1,ME1,Rep2,control


datasetMatrix <- DESeqDataSetFromMatrix(countData = genecountsData, colData = metasheetData, design = ~ Condition)

# Find differential expressed genes using DESEq2
datasetMatrix <- DESeq(datasetMatrix)

# - - - - - - - - - - - - - 
# Differential expression analysis
# - - - - - - - - - - - - - 

# comparison using treatment and control, and dort the result acooding to pvalues
comparison_res <- results(datasetMatrix, contrast=c("Condition", "Day7", "WT5"))
deseq = comparison_res[order(comparison_res$pvalue),]
deseq.df <- as.data.frame(deseq)


# Subset for only significant genes (q < padj_cutoff & FC > lfc_cutoff)
results_sig <- subset(deseq.df, padj <= padj_cutoff & abs(log2FoldChange) >= lfc_cutoff)
head(results_sig)

# - - - - - - - - - - - - -
# Export table result data
# - - - - - - - - - - - - -

# Write normalized gene counts for all samples to a .txt file
# write.table(x = as.data.frame(counts(datasetMatrix), normalized = T), file = 'gene_count_normalized.txt', sep = '\t', quote = F, col.names = NA)
# Write the annotated results table to a .txt file
# write.table(x = as.data.frame(deseq.df), file = "all_filtered_gene_annotated.txt", sep = '\t', quote = F, col.names = NA)

# Write significant normalized gene counts to a .txt file
# write.table(x = counts(datasetMatrix[row.names(results_sig)], normalized = T), file = 'gene_count_normalized_significant.txt', sep = '\t', quote = F, col.names = NA)
# Write significant annotated results table to a .txt file
# write.table(x = as.data.frame(results_sig), file = "all_filtered_gene_annotated_significant.txt", sep = '\t', quote = F, col.names = NA)

results_gene_count <- as.data.frame(counts(datasetMatrix, normalized = T))
results_gene_count$GENEID <- row.names(results_gene_count)

results_annotation <- deseq.df
results_annotation$GENEID <- row.names(results_annotation)

file_merged <- merge(results_gene_count, results_annotation, by = "GENEID", sort = F)
write.table(file_merged, genecountsAnno, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)


if (file.exists("./figures")){
} else {
	dir.create("./figures")
}

# -------------------------
# Convert all samples to rlog
# perform the rlog transformation and used as gene expression matrix for heatmap and PCA application
# -------------------------

datasetMatrix_rlog <- rlog(datasetMatrix, blind=FALSE)
head(assay(datasetMatrix_rlog), 3)

# -------------------------
# Sample distance heatmap
# ------------------------
sampleDists <- dist( t( assay(datasetMatrix_rlog) ) )
sampleDists

# We visualize the distances in a heatmap in a Figure below, using the function pheatmap from the pheatmap package.
# Heatmap of sample-to-sample distances using the rlog-transformed values

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- metasheetData$Sample #rownames(metasheetData)
colnames(sampleDistMatrix) <- metasheetData$Sample #rownames(metasheetData)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
p <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         angle_col = 90, col=colors, fontsize = 8)

ggsave(p, filename="figures/Heatmap_sample_to_sample_distances.pdf", width=9, height=8, units=c("cm"), colormodel="srgb")


# - - - - - - - - - - - - - 
# PCA plot
# - - - - - - - - - - - - - 
# PCA plot using the rlog-transformed values. Each sample is color-coded according to their affiliations.

# Plot PCA by column variable
p <- plotPCA(datasetMatrix_rlog, intgroup = "Group", ntop = 1000)       # + theme_bw() + ggsave('figures/pca_plot.pdf')
ggsave(p, filename="figures/pca_plot.pdf", width=9, height=5, units=c("cm"), colormodel="srgb")


# - - - - - - - - - - - - - 
# MA plot
# - - - - - - - - - - - - -
pdf(file='figures/MA_plot.pdf', width=5, height=4, colormodel="srgb")
plotMA(deseq, ylim = c(-5, 5))
dev.off()


# - - - - - - - - - - - - -
# Volcano plot
# - - - - - - - - - - - - -
# Load libraries
library(ggplot2)
library(RColorBrewer)

# Gather Log-fold change and FDR-corrected pvalues from deseq2 results
data <- data.frame(pval = -log10(deseq$padj), lfc = deseq$log2FoldChange, row.names = row.names(deseq))

# Remove any rows that have NA as an entry
data <- na.omit(data)

# Color the points which are up or down
## If fold-change > 1.5 and pvalue < 0.05 (Increased significant)
## If fold-change < -1.5 and pvalue < 0.05 (Decreased significant)
data <- mutate(data, color = case_when(data$lfc >= lfc_cutoff & data$pval >= pval_cutoff   ~ "Increased",
                                       data$lfc <= - lfc_cutoff & data$pval >= pval_cutoff ~ "Decreased",
                                       data$lfc > - lfc_cutoff & data$lfc < lfc_cutoff     ~ "Nonsignificant",
                                       data$pval < pval_cutoff                             ~ "Nonsignificant"))

# Make a basic ggplot2 object with x-y values


p <- ggplot(data, aes(x = lfc, y = pval, color = color)) +
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 1.0, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality", values = c(Increased = "#CD4F39", Decreased = "#008B00", Nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + # change overall theme
  theme(legend.position = "right") + # change the legend
  xlab(expression(log[2]("Treatment" / "Control"))) + # Change X-Axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # Change Y-Axis label
  geom_hline(yintercept = 1.301029996, colour = "darkgrey") + # Add p-adj value cutoff line
  scale_y_continuous(trans = "log1p") # Scale yaxis due to large p-value

ggsave(p, filename="figures/Volcano_plot.pdf", width=12, height=8, units=c("cm"), colormodel="srgb")
