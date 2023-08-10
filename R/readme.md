## Differential expression (DE) analysis in R
Count data from kallisto pseudoalignment are imported into R, further analysis was done using RStudios. The steps below are repeated for each set of data that were produced with the core, soft-core, pangenome and PAO1 reference. 
- Define function to save pheatmaps to be used later 
- Installation  any required packages, most packages can be found with BiocManager

```R

#function to save pheatmaps
save_pheatmap_pdf <- function(x, filename, width=11, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

##### Package installations

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.17")

BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("tximportData")
BiocManager::install("rhdf5")
BiocManager::install("vsn")
BiocManager::install("apeglm")

##### Load R packages 
library(tximportData)
library(tximport)
library(readr)
library(dplyr)
library(rhdf5)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(vsn)

##### Import data into R

# Set the directory, folders "CORE", "SOFTCORE", "PAN", "REF" are stored under the kallisto folder
dir <- "~/kallisto"
ref <- "REF"
list.files(dir)

# Import metadata for samples (e.g SRA, sample_name, group)
samples <- read.table(file.path(dir, "samples.tsv"), header = TRUE)
samples

# Import kallisto abundance files
files <-file.path(dir, ref, samples$sra_accession, "abundance.h5")
all(file.exists(files))
names(files) <- samples$sample

# Import Gene Id and Gene name
tx_id <- read.table(file.path(dir, "tx_id.txt"), sep = ',', header = TRUE )
colnames(tx_id) <- "tx_id"
tx_id$gene_id <- tx_id$PGD.Gene.ID

# Import kallisto file content into 
txi <- tximport(files, type = "kallisto", tx2gene = tx_id)
names(txi)

tpm <- txi$abundance
colSums(txi$counts)

# Some TPM stats and plot
tpm_stats <- transform(tpm, SD=rowSds(tpm), AVG=rowMeans(tpm), MED=rowMedians(tpm))

ggplot(tpm_stats) +
  aes(x = SD, y = MED) +
  geom_hex(shape=16, size=2) +
  labs(y="Medians", x="Standard Deviations", 
       title="Transcripts per million (TPM)", 
       subtitle= "unfiltered, non-normalised data",
       caption="here's an interesting caption") +
  theme_bw()
ggsave(paste(dir, ref, "TPM_raw.pdf", sep = '/' ))


#### DESeq2

sample_info <- samples
rownames(sample_info) <- sample_info$sample
colData <- sample_info[,c("sra_accession","group")]
colData$sra_accession <- factor(colData$sra_accession)
colData$group <- factor(colData$group)

# Create DESeq2 object
dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ group)
head(dds)

# Counts per sample and counts per gene
write.csv(colSums(counts(dds)), paste(dir,ref,"counts_per_sample.csv", sep = '/'), row.names = TRUE)
write.csv(rowSums(counts(dds)), paste(dir,ref,"counts_per_gene.csv", sep = '/'), row.names = TRUE)

# Filter low count reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq method
dds <- DESeq(dds)
res <- results(dds)

# Order by lowest padj values 
head(res[order(res$padj),])

write.csv(res, paste(dir,ref,"result_deseq2_invivo_vs_invitro.csv", sep = '/'), row.names = TRUE)

# LFC shrinkage/data transformation for visualisation ----

library(apeglm)

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="group_sputum_vs_in_vitro", type="apeglm")
resLFC

resOrdered <- res[order(res$pvalue),]

# exploring genes with padj < 0.05 and log fold change >1 
res05 <- results(dds, alpha = 0.05)
summary(res05)

resSig <- (res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >1),])
summary(resSig)

# Upregulated DEGS in in vivo samples
upreg <- resSig[which(resSig$log2FoldChange >1),]
up_30 <- upreg[order(-upreg$log2FoldChange),] 
write.table(up_30, paste(dir,ref,'UP_reg_invivo_all.txt', sep = '/'), col.names = T, row.names = T, sep = '\t', quote = F)

# Downregulated DEGS in in vivo samples
downreg <- resSig[which(resSig$log2FoldChange <1),]
down_30 <- downreg[order(downreg$log2FoldChange),]
write.table(down_30, paste(dir,ref,'DOWN_reg_invivo_all.txt', sep = '/'), col.names = T, row.names = T, sep = '\t', quote = F)


# Significant DEGs stored into tables
write.table(resSig, paste(dir,ref,'DEG_padj_0.05_logFC_1.txt', sep = '/'), col.names = T, row.names = T, sep = '\t', quote = F)
write.table(rownames(resSig), paste(dir,ref, 'DEG_padj_0.05_logFC_1_gene_names.txt', sep = '/'), col.names = F, row.names = F, sep = '\t', quote = F)
write.table(rownames(resSig), paste(dir,ref,'All_tested_genes_gene_names.txt', sep = '/'), col.names = F, row.names = F, sep = '\t', quote = F)


# MA plot
plotMA(res, ylim=c(-2,2))
resLFC <- lfcShrink(dds_all, coef="group_sputum_vs_in_vitro", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))


# Dispersion plot
plotDispEsts(dds)


# Variance and PCA
library(gridExtra)

# SD mean plot
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)
msd <- meanSdPlot(assay(dds))

pdf(paste(dir, ref, 'PCA.pdf', sep = '/'))
plotPCA(rld, intgroup="group")
dev.off()

pdf(paste(dir, ref, 'msd_all.pdf', sep = '/'), width=11, height=8)
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))
meanSdPlot(assay(ntd))
dev.off()

# Heatmap 
library("RColorBrewer")

# Sample to sample distance 
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Plot and save heatmap
xx <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
save_pheatmap_pdf(xx, paste(dir,ref,'sample_dist_heatmap.pdf', sep='/'))

# Volcano plots?
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

pdf(paste(dir,ref,'volcano_plot.pdf', sep = '/'), width = 10, height=11 )
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'In vivo versus In vitro',
                pCutoff = 0.001,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 5.0,
                xlim = c(-6, 6),
                ylim = c(0, 20))
dev.off()

```
