## 2026.07.20
## Processing H3K27me3 counts with DiffBind
#Goal :  Getting the differential peaks of H3K27me3 with LFAP spikein for all D2 comparaison
suppressMessages(library(DESeq2))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
setwd("/Users/maxim/Documents/Marlet/Results/H3K27me3/")
load("MMF_H3K27me3_spikein_featureCounts.RData")

## Explore count matrix
dim(counts$counts)
colnames(counts$counts)

## Create metadata
condition <- c("D0_TKO_H3","D0_TKO_H3","D0_TKO_K4M","D0_TKO_K4M",
               "D0_WT_H3", "D0_WT_H3","D0_WT_K4M","D0_WT_K4M",
               "D2_TKO_H3","D2_TKO_H3","D2_TKO_K4M","D2_TKO_K4M",
               "D2_WT_H3", "D2_WT_H3","D2_WT_K4M","D2_WT_K4M")
rep <- rep(c(1,2), times = 8)
metaData <- data.frame(Sample = colnames(counts$counts), condition = condition, replicate = rep )
metaData
rm(condition, rep, gtf_file)

## Create Granges object for DESeq
dds <- DESeqDataSetFromMatrix(countData = counts$counts,
                              colData = metaData,
                              design = ~ condition)
dds

## Add normalization factors from bigWigs; cap at six digits
normFactors <- c(5.424656, 2.118055, 2.751415, 3.300762,
                 1.818131, 2.815352, 3.566436, 3.307109,
                 1.627687, 3.939744, 1.566919, 12.68346,
                 1.377817, 4.568384, 1.383347, 4.459780)
## Invert normalization (bigWig multiplies, DESeq2 divides)
normFactors <- 1/normFactors
normMatrix <- matrix(normFactors,
                     nrow  = nrow(dds),      
                     ncol  = ncol(dds),
                     byrow = TRUE)           
rownames(normMatrix) <- rownames(dds)
colnames(normMatrix) <- colnames(dds)

## Normalize dds object
normalizationFactors(dds) <- normMatrix
## Preview normalization matrix
normalizationFactors(dds)[1:5, 1:4]

## -- What is the count distribution before and after normalization?
Counts_norm <- counts(dds, normalized=TRUE)
par(mfrow = c(1,2))
hist(log2(counts$counts), main = "Original counts")
hist(log2(Counts_norm), main = "Normalized counts")
par(mfrow = c(1,1))
## -- Adding all reads in the quantified peaks
cat("Original reads:\n")
sort(colSums(counts$counts))
cat("Normalized reads:\n")
sort(colSums(Counts_norm))
## -- How many fold between min and max number of reads?
## Original counts: 8.381734
cat("Original reads max/min:  ", max(colSums(counts$counts))/min(colSums(counts$counts)), "\n")
## Normalized counts: 2.30931 (with non-inversed normalization matrix, it used to be 67.85)
cat("Normalized reads max/min:  ", max(colSums(Counts_norm))/min(colSums(Counts_norm)), "\n")

rm(Counts_norm, normMatrix)


## ------------------ ##
## --  Clustering  -- ##
## ------------------ ##

# The rlog tends to work well on small datasets (n < 30), potentially outperforming the VST
# when there is a wide range of sequencing depth across samples.

## Cluster samples
rld <- rlog(dds, blind = TRUE)
pcaData <- plotPCA(rld, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

## Add factors
pcaData$condition <- factor(pcaData$condition, levels = c("D0_WT_H3", "D0_WT_K4M", "D0_TKO_H3", "D0_TKO_K4M", "D2_WT_H3", "D2_WT_K4M", "D2_TKO_H3", "D2_TKO_K4M"))
pcaData$replicate <- factor(pcaData$replicate, levels = c(1, 2))

## Plot PCA
p <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(size=5) +
  scale_color_manual(values = c("#ff5d5e", "#ffab28", "#56a3fd", "#00d987",
                                "#D10000", "#D18100", "#0262CF", "#007548")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()
p
#ggsave(filename = "../../Plots/H3K27me3/PCA_rlog_H3K27me3_spikeIn_counts.pdf", plot = p, width = 16, height = 12, units = "cm")

## Plot hierarchical clustering
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
## Add labels
rownames(sampleDistMatrix) <- paste(rld$condition, rld$replicate, sep = ".")
colnames(sampleDistMatrix) <- paste(rld$condition, rld$replicate, sep = ".")
## Plot
p <- pheatmap(sampleDistMatrix,
         color = colorRampPalette(rev(RColorBrewer::brewer.pal(9,"Blues")))(255),
         border_color = "gray70",
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists, 
         treeheight_row = 40,
         treeheight_col = 40,
         angle_col = 45)
p
#ggsave(filename = "../../Plots/H3K27me3/HierClustering_rlog_H3K27me3_spikeIn.pdf", plot = p, width = 16, height = 14, units = "cm")

## Clear environment
rm(p, pcaData, rld, sampleDistMatrix, percentVar, sampleDists)


## -------------------------------------------
## -------------------------------------------

## -- Process with DESeq directly (and not through DiffBind)
dds <- DESeq(dds)

## Specify contrasts
G1 <- c("D0_WT_H3",  "D0_TKO_H3",  "D0_WT_H3",  "D0_WT_K4M",  "D2_WT_H3",  "D2_TKO_H3",  "D2_WT_H3",  "D2_WT_K4M")
G2 <- c("D0_WT_K4M", "D0_TKO_K4M", "D0_TKO_H3", "D0_TKO_K4M", "D2_WT_K4M", "D2_TKO_K4M", "D2_TKO_H3", "D2_TKO_K4M")
G1_G2 <- data.frame("Group_1" = G1, "Group_2" = G2)

## Iterate over every contrast combination
for(i in 1:nrow(G1_G2)){
  ## Select comparisons
  group_1 <- G1_G2$Group_1[i]
  group_2 <- G1_G2$Group_2[i]
  cat("Comparing ",group_1," vs ",group_2,"\n",sep = "")
  
  ## Specify log file and Plot_folder
  Name_log <- paste("log_files/",group_1,"_vs_",group_2,"_H3K27me3_spikeIn_DE.txt", sep = "")
  Plot_folder <- "../../Plots/H3K27me3/"
  cat("Comparing ",group_1," vs ",group_2,"\n\n",file = Name_log, sep = "", append = FALSE)
  
  ## Get normalization factors
  tmp <- normalizationFactors(dds)[1,]
  tmp_1 <- tmp[grep(pattern = group_1, x = names(tmp))]
  tmp_2 <- tmp[grep(pattern = group_2, x = names(tmp))]
  
  cat("Normalization factors ",group_1,": ",file = Name_log, sep = "", append = TRUE)
  cat(tmp_1,"\n",file = Name_log, sep = "  ", append = TRUE)
  cat("Normalization factors ",group_2,": ",file = Name_log, sep = "", append = TRUE)
  cat(tmp_2,"\n\n",file = Name_log, sep = "  ", append = TRUE)
  rm(tmp, tmp_1, tmp_2)
  
  ## -- Differential Enrichment -- ##
  res <- results(dds, alpha = 0.05, contrast = c("condition", group_2, group_1))
  
  tmp_up <- sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = TRUE)
  tmp_down <- sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = TRUE)
  tmp_low <- sum(is.na(res$padj))
  
  cat("\n** -- Differential Enrichment -- **\n\n", file = Name_log, sep = "", append = TRUE)
  cat("Total domains:\t\t",res@nrows,"\n",file = Name_log, sep = "", append = TRUE)
  cat("LFC > 0 (up):\t\t",tmp_up," (",round(tmp_up*100/res@nrows, 2),"%)\n",file = Name_log, sep = "", append = TRUE)
  cat("LFC < 0 (down):\t\t",tmp_down," (",round(tmp_down*100/res@nrows, 2),"%)\n",file = Name_log, sep = "", append = TRUE)
  cat("Domains with low counts:\t",tmp_low," (",round(tmp_low*100/res@nrows, 2),"%)\n\n",file = Name_log, sep = "", append = TRUE)
  
  ## Assumming logFC = 1
  tmp_up <- sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE)
  tmp_down <- sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)
  cat("LFC > 1 (up):\t\t",tmp_up," (",round(tmp_up*100/res@nrows, 2),"%)\n",file = Name_log, sep = "", append = TRUE)
  cat("LFC < -1 (down):\t\t",tmp_down," (",round(tmp_down*100/res@nrows, 2),"%)\n\n",file = Name_log, sep = "", append = TRUE)

  ## Export res in nice format
  tmp <- as.data.frame(res)
  tmp[,c("baseMean", "log2FoldChange", "lfcSE", "stat")] <- round(tmp[,c("baseMean", "log2FoldChange", "lfcSE", "stat")], digits = 4)
  tmp[,c("pvalue", "padj")] <- signif(tmp[,c("pvalue","padj")], digits = 5)
  tmp <- tmp |> tibble::rownames_to_column(var = "coord")
  Name_res <- paste(group_1,"_vs_",group_2,"_H3K27me3_spikeIn.txt", sep = "")
  write.table(x = tmp, file = Name_res, append = FALSE, quote = FALSE, sep = '\t', eol = '\n', na = "", row.names = FALSE, col.names = TRUE)
  rm(tmp, Name_res)
  
  ## Stats of domains with counts: remove rows with padj = NA (they were not used for Bonferroni correction either way)
  res_df <- as.data.frame(res)
  res_df <- res_df |> filter(!is.na(padj))
  cat("Remaining domains:\t",nrow(res_df),"\n\n",file = Name_log, sep = "", append = TRUE)
  
  ## Mark differentially-enriched (logFC > 1)
  res_df$direction <- "none"
  res_df$direction[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "up"
  res_df$direction[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "down"
  res_df$direction <- factor(res_df$direction, levels = c("down", "none", "up"))
  res_df <- res_df |> tibble::rownames_to_column(var = "name")
  lim_fc <- max(abs(res_df$log2FoldChange))
  
  tmp_none <- sum(res_df$direction == "none")
  tmp_up <- sum(res_df$direction == "up")
  tmp_down <- sum(res_df$direction == "down")
  cat("Not Diff Enriched:\t\t",tmp_none," (",round(tmp_none*100/nrow(res_df), 2),"%)\n",file = Name_log, sep = "", append = TRUE)
  cat("Up in ",group_2,":\t\t",tmp_up," (",round(tmp_up*100/nrow(res_df), 2),"%)\n",file = Name_log, sep = "", append = TRUE)
  cat("Down in ",group_2,":\t\t",tmp_down," (",round(tmp_down*100/nrow(res_df), 2),"%)\n\n",file = Name_log, sep = "", append = TRUE)
  rm(tmp_none, tmp_up, tmp_down, tmp_low)
  
  ## Define color groups
  color_groups <- c("#ff5d5e", "#ffab28", "#56a3fd", "#00d987", "#ff5d5e", "#ffab28", "#56a3fd", "#00d987")
  names(color_groups) <-  c("D0_WT_H3", "D0_WT_K4M", "D0_TKO_H3", "D0_TKO_K4M", "D2_WT_H3", "D2_WT_K4M", "D2_TKO_H3", "D2_TKO_K4M") 
    
  ## -- Plot volcano
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
    geom_point(size = 1, alpha = 0.7, na.rm = TRUE) +
    scale_color_manual(values = as.character(c(color_groups[group_1], "gray50", color_groups[group_2]))) +
    geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed", color = 'gray30') +
    geom_vline(xintercept = c(-1,1), linewidth = 0.3, linetype = "dashed", color = 'gray30') +
    coord_cartesian(xlim = c(-lim_fc, lim_fc)) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
  p
  Name_plot <- paste("Volcano_plot/",group_1,"_vs_",group_2,"_logFC_padj_volcano.pdf", sep = "")
  ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = p, width = 12, height = 12, units = "cm")
  
  ## -- Plot MA
  ## Re-store res to keep low_count genes
  res_full <- as.data.frame(res)
  res_full$direction <- "none"
  res_full$direction[res_full$padj < 0.05 & res_full$log2FoldChange > 1] <- "up"
  res_full$direction[res_full$padj < 0.05 & res_full$log2FoldChange < -1] <- "down"
  res_full$direction <- factor(res_full$direction, levels = c("down", "none", "up"))
  res_full <- res_full |> tibble::rownames_to_column(var = "name")
  lim_fc <- max(abs(res_full$log2FoldChange))
  
  p <- ggplot(res_full, aes(x = log10(baseMean), y = log2FoldChange, color = direction)) +
    geom_point(size = 0.5, alpha = 1) +
    scale_color_manual(values = as.character(c(color_groups[group_1], "gray50", color_groups[group_2]))) +
    #geom_hline(yintercept = c(-1,1), linewidth = 0.3, linetype = "dashed", color = 'blue') +
    geom_hline(yintercept = 0, linewidth = 0.5, linetype = "solid", color = 'black') +
    #coord_cartesian(ylim = c(-lim_fc, lim_fc)) +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")
  p
  Name_plot <- paste("MA_plot/",group_1,"_vs_",group_2,"_meanExp_logFC_MA.pdf", sep = "")
  ggsave(filename = Name_plot, device = "pdf", path = Plot_folder, plot = p, width = 12, height = 12, units = "cm")
}



# ## ---- D0 ---- ##
# 
# ## -- WT D0 H3 vs WT D0 K4M 
# D0_db_H3_WT_vs_K4M_WT <- results(dds, contrast=c("condition","D0_WT_K4M", "D0_WT_H3"), alpha=0.05)
# D0_db_H3_WT_vs_K4M_WT
# summary(D0_db_H3_WT_vs_K4M_WT)
# plotMA(D0_db_H3_WT_vs_K4M_WT)
# ## The plot also looks slightly different; before all observarions were below logFC=0
# ## Inverse: Now we get DE (mostly down).
# 
# ## -- TKO D0 H3 vs TKO D0 K4M 
# D0_db_H3_TKO_vs_K4M_TKO <- results(dds, contrast=c("condition","D0_TKO_K4M", "D0_TKO_H3"), alpha=0.05)
# D0_db_H3_TKO_vs_K4M_TKO
# summary(D0_db_H3_TKO_vs_K4M_TKO)
# plotMA(D0_db_H3_TKO_vs_K4M_TKO)
# 
# ## -- WT D0 H3 vs TKO D0 H3 
# D0_db_H3_WT_vs_H3_TKO <- results(dds, contrast=c("condition","D0_TKO_H3", "D0_WT_H3"), alpha=0.05)
# D0_db_H3_WT_vs_H3_TKO
# summary(D0_db_H3_WT_vs_H3_TKO)
# plotMA(D0_db_H3_WT_vs_H3_TKO)
# 
# ## -- WT D0 K4M vs TKO D0 K4M 
# D0_db_K4M_WT_vs_K4M_TKO <- results(dds, contrast=c("condition","D0_TKO_K4M", "D0_WT_K4M"), alpha=0.05)
# D0_db_K4M_WT_vs_K4M_TKO
# summary(D0_db_K4M_WT_vs_K4M_TKO)
# plotMA(D0_db_K4M_WT_vs_K4M_TKO)
# 
# ## ---- D2 ---- ##
# 
# ## -- WT D2 H3 vs WT D2 K4M 
# D2_db_H3_WT_vs_K4M_WT <- results(dds, contrast=c("condition","D2_WT_K4M", "D2_WT_H3"), alpha=0.05)
# D2_db_H3_WT_vs_K4M_WT
# summary(D2_db_H3_WT_vs_K4M_WT)
# plotMA(D2_db_H3_WT_vs_K4M_WT)
# 
# ## -- TKO D2 H3 vs TKO D2 K4M 
# D2_db_H3_TKO_vs_K4M_TKO <- results(dds, contrast=c("condition","D2_TKO_K4M", "D2_TKO_H3"), alpha=0.05)
# D2_db_H3_TKO_vs_K4M_TKO
# summary(D2_db_H3_TKO_vs_K4M_TKO)
# plotMA(D2_db_H3_TKO_vs_K4M_TKO)
# ## Inverse: big gain
# 
# ## -- WT D2 H3 vs TKO D2 H3 
# D2_db_H3_WT_vs_H3_TKO <- results(dds, contrast=c("condition","D2_TKO_H3", "D2_WT_H3"), alpha=0.05)
# D2_db_H3_WT_vs_H3_TKO
# summary(D2_db_H3_WT_vs_H3_TKO)
# plotMA(D2_db_H3_WT_vs_H3_TKO)
# 
# ## -- WT D2 K4M vs TKO D2 K4M 
# D2_db_K4M_WT_vs_K4M_TKO <- results(dds, contrast=c("condition","D2_TKO_K4M", "D2_WT_K4M"), alpha=0.05)
# D2_db_K4M_WT_vs_K4M_TKO
# summary(D2_db_K4M_WT_vs_K4M_TKO)
# plotMA(D2_db_K4M_WT_vs_K4M_TKO)
 

