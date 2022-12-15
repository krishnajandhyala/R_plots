library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(EnhancedVolcano)
library(devtools)
library(pcaExplorer)
library(tidyverse)
# Read the counts and design files
counts <- as.matrix(read.csv('LO_vs_LOensembl_counts.tsv', sep='\t', row.names='Geneid'))
counts

col_data <- read.csv('LO_vs_LO_metadata.txt', sep='\t', row.names=1)
col_data
# Check that the rownames of design and column names of counts match
all(rownames(col_data) %in% colnames(counts))
all(rownames(col_data) == colnames(counts))
# Create a DESeq Dataset from counts matrix
dds <- DESeqDataSetFromMatrix(countData = counts, colData = col_data, design = ~condition)
dds
# See counts from the dds object
counts(dds)
# Filter transcripts with low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# Set the reference condition
dds$condition <- relevel(dds$condition, ref='LO')
dds$condition
# Differential Expression analysis
des <- DESeq(dds)
res_1 <- results(des, alpha = 0.05, contrast = c("condition","LO_ensembl","LO"),)
summary(res_1)
res_2 <- results(des, alpha = 0.05, contrast = c("condition","hiPSC","LO"),)
summary(res_2)

write.table(res_1, file = 'results_1.txt', sep = "\t", row.names = T)
write.table(res_2, file = 'results_2.txt', sep = "\t", row.names = T)

#PCA
vsd <- vst(des, blind = F)
rld <- rlog(dds, blind=FALSE)
rld
#PCA <- plotPCA(rld, intgroup = c("condition")) #+ geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black")
PCA <- pcaplot(rld, intgroup = c("condition"),ntop = 1000, pcX = 1, pcY = 3, title = "LO vs LO_ensembl PCA", ellipse = TRUE)
PCA

#Volcano plot_1
keyvals <- ifelse(
  res_1$log2FoldChange < -1 & res_1$pvalue < 0.05, 'red',
  ifelse(res_1$log2FoldChange > 1 & res_1$pvalue < 0.05, 'royalblue',
         'grey'))

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'red'] <- 'low'

#selectLab_italics = paste0(
#  "italic('",
#  c('AGER','CXCR4','FOXA2','KRT5','MUC5AC','NKX2-1','SCGB3A2','SOX17','SOX2','SPC24','SPC25'),
#  "')")
volcano1 <- EnhancedVolcano(res_1,lab = rownames(res_1), x = 'log2FoldChange',y = 'pvalue',selectLab = c('AGER','CXCR4','FOXA2','KRT5','MUC5AC','NKX2-1','SCGB3A2','SOX17','SOX2','SPC24','SPC25'),title = 'LO vs LO ensembl dataset',subtitle = "Volcano Plot", colCustom = keyvals,pCutoff = 0.05, FCcutoff = 2)
volcano1

#Volcano plot_2
keyvals <- ifelse(
  res_2$log2FoldChange < -1 & res_2$pvalue < 0.05, 'red',
  ifelse(res_2$log2FoldChange > 1 & res_2$pvalue < 0.05, 'royalblue',
         'grey'))

keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'royalblue'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'red'] <- 'low'

volcano2 <- EnhancedVolcano(res_2,lab = rownames(res_2), x = 'log2FoldChange',y = 'pvalue',selectLab = c('AGER','CXCR4','FOXA2','KRT5','MUC5AC','NKX2-1','SCGB3A2','SOX17','SOX2','SPC24','SPC25'),title = 'LO vs hiPSC ensembl dataset',subtitle = "Volcano Plot", colCustom = keyvals,pCutoff = 0.05, FCcutoff = 2)
volcano2
#volcano <- EnhancedVolcano(res, lab = rownames(res), selectLab = c('AGER','CXCR4','FOXA2','KRT5','MUC5AC','NKX2-1','SCGB3A2','SOX17','SOX2','SPC24','SPC25'), x = 'log2FoldChange',y = 'pvalue',title = 'Lung Organoids vs hiPSC',pCutoff = 0.05, FCcutoff = 2)

#filtered <- subset(res_1, (padj < 0.05) & (abs(log2FoldChange > 1 | (abs(log2FoldChange < -1)))))
#summary(filtered)

#write.table(filtered, file = "filtered.tsv", sep = "\t", row.names = T)
#NHEK_res <- res_1[(res_1$padj < 0.05) & (res_1$lfcSE > 2) & !(is.na(res_1$padj)),]

dds
min(as.data.frame(assay(dds))$LB01)
class(keep)
max(keep$LB01)
max(res_1$log2FoldChange)

pal <- colorRampPalette(c("#67001F", 'white', "#053061"))
selected <- order(res_1$padj, decreasing = F)[1:100]
assay(rld)[selected,]
#heatmap <- pheatmap(assay(rld)[selected,], scale =  cluster_rows = T, show_rownames = TRUE, cluster_cols = TRUE, annotation_col = col_data, breaks = c(seq(1, 10000, by = 100)))
heatmap_1 <- pheatmap(assay(rld)[selected,], cluster_rows = T, show_rownames = T, cluster_cols = T, border_color = F, main = "LO vs LO ensembl Dataset", 
         scale = "row", color = pal(100), breaks = c(seq(-2, 2, by = 0.04)),fontsize = 12,legend_breaks = c(-2,-1,0,1,2), legend_labels = c(-200,-100,0,100,200))
heatmap_1

pal <- colorRampPalette(c("#67001F", 'white', "#053061"))
selected <- order(res_2$padj, decreasing = F)[1:10]
assay(rld)[selected,]
#heatmap <- pheatmap(assay(rld)[selected,], scale =  cluster_rows = T, show_rownames = TRUE, cluster_cols = TRUE, annotation_col = col_data, breaks = c(seq(1, 10000, by = 100)))
heatmap_2 <- pheatmap(assay(rld)[selected,], cluster_rows = T, show_rownames = T, cluster_cols = T, border_color = F, main = "LO vs hiPSC Ensembl dataset",
                    scale = "row", color = pal(100), breaks = c(seq(-2, 2, by = 0.04)),fontsize = 12,legend_breaks = c(-2,0,2))
heatmap_2

# Get genes upregulated in condition 1 using adjusted pvalue and log2FoldChange
up_res_1 <- res_1[(res_1$padj < 0.05) & (res_1$log2FoldChange > 1) & !(is.na(res_1$padj)),]
 # Save the details of DE genes into files
write.table(as.data.frame(up_res_1), file = 'upregulated_1.tsv', row.names = T, quote = F, sep = "\t")
# Get genes downregulated in sample 2 using adjusted pvalue and log2FoldChange
down_res_1 <- res_1[(res_1$padj < 0.05) & (res_1$log2FoldChange < -1) & !(is.na(res_1$padj)),]
 # Save the details of DE genes into files
write.table(as.data.frame(down_res_1), file = 'downregulated_1.tsv', row.names = T, quote = F, sep = "\t")

# Get genes upregulated in condition 2 using adjusted pvalue and log2FoldChange
up_res_2 <- res_2[(res_2$padj < 0.05) & (res_2$log2FoldChange > 1) & !(is.na(res_2$padj)),]
# Save the details of DE genes into files
write.table(as.data.frame(up_res_2), file = 'upregulated_2.tsv', row.names = T, quote = F, sep = "\t")
# Get genes downregulated in sample 2 using adjusted pvalue and log2FoldChange
down_res_2 <- res_2[(res_2$padj < 0.05) & (res_2$log2FoldChange < -1) & !(is.na(res_2$padj)),]
# Save the details of DE genes into files
write.table(as.data.frame(down_res_2), file = 'downregulated_2.tsv', row.names = T, quote = F, sep = "\t")

#Load R package
library(clusterProfiler)
library(org.Hs.eg.db)
 # Get all supported key types
keytypes(org.Hs.eg.db)
# Get names of upregulated genes in Sample2
up_genes_1 <- rownames(up_res_1)
# Run Gene Ontology analysis
up_GOBP_1 <- enrichGO(up_genes_1, OrgDb=org.Hs.eg.db, ont = "BP", keyType = 'SYMBOL', qvalueCutoff=0.05)
# See significant categories as dotplot
dotplot1_1 <- dotplot(up_GOBP_1, showCategory = 20)
down_genes_1 <- rownames(down_res_1)
down_GOBP_1 <- enrichGO(down_genes_1, OrgDb=org.Hs.eg.db, ont = "BP", keyType = 'SYMBOL', qvalueCutoff=0.05)
dotplot2_1 <- dotplot(down_GOBP_1, showCategory = 20)
#write.table(as.data.frame(down_GOBP), file = 'down_GO_analysis.tsv', row.names = F, quote = F, sep = "\t")
#write.table(as.data.frame(up_GOBP), file = 'up_GO_analysis.tsv', row.names = F, quote = F, sep = "\t")
# Use compareCluster function to do the analysis for both lists simultaneously
gene_list_1 <- list(up = paste(up_genes_1), down = paste(down_genes_1))
list_GOBP_1 <- compareCluster(gene_list_1, fun = "enrichGO", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", qvalueCutoff=0.05)
dotplotBP_1 <- dotplot(list_GOBP_1, showCategory=10, title = "GO Biological Process (LO vs LO ensembl dataset)")
dotplotBP_1

# Get names of upregulated genes in Sample3
up_genes_2 <- rownames(up_res_2)
# Run Gene Ontology analysis
up_GOBP_2 <- enrichGO(up_genes_2, OrgDb=org.Hs.eg.db, ont = "BP", keyType = 'SYMBOL', qvalueCutoff=0.05)
# See significant categories as dotplot
dotplot1_2 <- dotplot(up_GOBP_2, showCategory = 20)
down_genes_2 <- rownames(down_res_2)
down_GOBP_2 <- enrichGO(down_genes_2, OrgDb=org.Hs.eg.db, ont = "BP", keyType = 'SYMBOL', qvalueCutoff=0.05)
dotplot2_2 <- dotplot(down_GOBP_2, showCategory = 20)
#write.table(as.data.frame(down_GOBP), file = 'down_GO_analysis.tsv', row.names = F, quote = F, sep = "\t")
#write.table(as.data.frame(up_GOBP), file = 'up_GO_analysis.tsv', row.names = F, quote = F, sep = "\t")
# Use compareCluster function to do the analysis for both lists simultaneously
gene_list_2 <- list(up = paste(up_genes_2), down = paste(down_genes_2))
list_GOBP_2 <- compareCluster(gene_list_2, fun = "enrichGO", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP", qvalueCutoff=0.05)
dotplotBP_2 <- dotplot(list_GOBP_2, showCategory=10, title = "GO Biological Process (LO vs hiPSC ensembl dataset)")
dotplotBP_2

#CC_1
up_GOCC_1 <- enrichGO(up_genes_1, OrgDb=org.Hs.eg.db, ont = "CC", keyType = 'SYMBOL', qvalueCutoff=0.05)
# See significant categories as dotplot
dotplot4_1 <- dotplot(up_GOCC_1, showCategory = 20)
down_GOCC_1 <- enrichGO(down_genes_1, OrgDb=org.Hs.eg.db, ont = "CC", keyType = 'SYMBOL', qvalueCutoff=0.05)
dotplot5_1 <- dotplot(down_GOCC_1, showCategory = 20)
# Use compareCluster function to do the analysis for both lists simultaneously
list_GOCC_1 <- compareCluster(gene_list_1, fun = "enrichGO", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC", qvalueCutoff=0.05)
dotplotCC_1 <- dotplot(list_GOCC_1, showCategory=10, title = "GO Cellular Component (LO vs LO ensembl dataset)")
dotplotCC_1

#CC_2
up_GOCC_2 <- enrichGO(up_genes_2, OrgDb=org.Hs.eg.db, ont = "CC", keyType = 'SYMBOL', qvalueCutoff=0.05)
# See significant categories as dotplot
dotplot4_2 <- dotplot(up_GOCC_2, showCategory = 20)
down_GOCC_2 <- enrichGO(down_genes_2, OrgDb=org.Hs.eg.db, ont = "CC", keyType = 'SYMBOL', qvalueCutoff=0.05)
dotplot5_2 <- dotplot(down_GOCC_2, showCategory = 20)
# Use compareCluster function to do the analysis for both lists simultaneously
list_GOCC_2 <- compareCluster(gene_list_2, fun = "enrichGO", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "CC", qvalueCutoff=0.05)
dotplotCC_2 <- dotplot(list_GOCC_2, showCategory=10, title = "GO Cellular Componentb (LO vs hiPSC ensembl dataset)")
dotplotCC_2

#MF_1
up_GOMF_1 <- enrichGO(up_genes_1, OrgDb=org.Hs.eg.db, ont = "MF", keyType = 'SYMBOL', qvalueCutoff=0.05)
# See significant categories as dotplot
dotplot7_1 <- dotplot(up_GOMF_1, showCategory = 20)
down_GOMF_1 <- enrichGO(down_genes_1, OrgDb=org.Hs.eg.db, ont = "MF", keyType = 'SYMBOL', qvalueCutoff=0.05)
dotplot8_1 <- dotplot(down_GOMF_1, showCategory = 20)
# Use compareCluster function to do the analysis for both lists simultaneously
list_GOMF_1 <- compareCluster(gene_list_1, fun = "enrichGO", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF", qvalueCutoff=0.05)
dotplotMF_1 <- dotplot(list_GOMF_1, showCategory=10, title = "GO Molecular Function (LO vs LO ensembl dataset)")
dotplotMF_1

#MF_2
up_GOMF_2 <- enrichGO(up_genes_2, OrgDb=org.Hs.eg.db, ont = "MF", keyType = 'SYMBOL', qvalueCutoff=0.05)
# See significant categories as dotplot
dotplot7_2 <- dotplot(up_GOMF_2, showCategory = 20)
down_GOMF_2 <- enrichGO(down_genes_2, OrgDb=org.Hs.eg.db, ont = "MF", keyType = 'SYMBOL', qvalueCutoff=0.05)
dotplot8_2 <- dotplot(down_GOMF_2, showCategory = 20)
# Use compareCluster function to do the analysis for both lists simultaneously
list_GOMF_2 <- compareCluster(gene_list_2, fun = "enrichGO", OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "MF", qvalueCutoff=0.05)
dotplotMF_2 <- dotplot(list_GOMF_2, showCategory=10, title = "GO Molecular Function (LO vs hiPSC ensembl dataset)")
dotplotMF_2


#KEGG_1
up_KEGG_1 <- enrichKEGG(gene = bitr(up_genes_1, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
barplot9_1 <- barplot(up_KEGG_1, showCategory = 20)
barplot9_1
down_KEGG_1 <- enrichKEGG(gene = bitr(down_genes_1, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
barplot10_1 <- barplot(down_KEGG_1, showCategory = 20)
barplot10_1
#KEGG_2
up_KEGG_2 <- enrichKEGG(gene = bitr(up_genes_2, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
dotplot11_1 <- dotplot(down_KEGG_1, showCategory = 20)
dotplot11_1
down_KEGG_2 <- enrichKEGG(gene = bitr(down_genes_2, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
dotplot12_1 <- dotplot(down_KEGG_1, showCategory = 20)
dotplot12_1

#res$label[which(rownames(res) %in% c('AGER','SPC25','FOXA2','SPC24','NKX2-1','KRT5','SOX17','SCGB3A2','SOX2','CXCR4'), arr.ind = T)] = c('AGER','SPC25','FOXA2','SPC24','NKX2-1','KRT5','SOX17','SCGB3A2','SOX2','CXCR4')
#volcano <- EnhancedVolcano(res,lab = NA, x = 'log2FoldChange',y = 'pvalue',title = 'Lung organoids difference with ensembl dataset',colCustom = keyvals,pCutoff = 0.05, FCcutoff = 2) + geom_text_repel(aes(label = label))