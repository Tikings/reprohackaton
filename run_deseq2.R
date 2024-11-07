if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages
BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("apeglm")
BiocManager::install("biomaRt")

library(ggplot2)
library(DESeq2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(httr)


library(tximport)
library(apeglm)
library(biomaRt)



count_data <- as.matrix(read.table("counts_final.txt", 
                                   header = TRUE, 
                                   row.names = 1, 
                                   sep = "\t", 
                                   skip = 1)[, -c(1:5)])

#print(count_data)


sample_names <- colnames(count_data)
sample_conditions <- c(
  "SRR10379721.bam" = "persister",
  "SRR10379722.bam" = "persister",
  "SSR10379723.bam" = "persister", 
  #erreur potentielle du T mais ca va être corrigé le SRR est tapé à la main
  "SRR10379723.bam" = "persister",
  "SRR10379724.bam" = "control",
  "SRR10379725.bam" = "control",
  "SRR10379726.bam" = "control" )

col_data <- data.frame(
  row.names = sample_names,
  condition = sample_conditions[sample_names] )

#print(col_data)


dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
summary(res)


res <- lfcShrink(dds, coef="condition_persister_vs_control", type="apeglm") 

plotMA(res, ylim=c(-5,5))





res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05, "Significant", "Non-significant")

print(res_df)

# Create the customized MA plot
ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 1.0, alpha=0.5) + # Adjust point size as needed
  scale_color_manual(values = c("Significant" = "red","Non-significant" = "black")) +
  scale_x_log10() +  # Log scale for base mean
  geom_hline(yintercept = 0, linetype = "dashed") +  # Add a dashed line at y = 0
  theme_minimal() +  # Use a minimal theme for a clean look
  labs(x = "Mean of normalized counts", y = "Log2 fold change") +
  theme(
    legend.position = "none",  # Remove legend if not needed
    panel.grid.minor = element_blank()  # Clean up minor grid lines
  ) +
  ylim(c(-4, 4)) # Set the y-axis limits

##### à suivre 



#  biomaRt  fait le pont entre l'ID Ensembl et l'annotation des gènes
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "aureowiki_gene")
gene_annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),mart = ensembl)

res <- as.data.frame(res)
res <- left_join(res, gene_annotations, by = c("row.names" = "ensembl_gene_id"))


# KEGG Translation genes
# est ce que je dois prendre tous les gènes de translation ?
kegg_translation <- GET("http://rest.kegg.jp/link/aureowiki/translation")
kegg_translation_genes <- content(kegg_translation, as = "text")
kegg_translation_genes <- read.table(text = kegg_translation_genes, sep="\t", header = FALSE)

translation_genes <- unique(kegg_translation_genes$V2)

# ajout d'une colonne avec vrai faux si gène à encercler
res$translation_gene <- ifelse(res$Row.names %in% translation_genes, TRUE, FALSE)


summary(res)





# modification potentielle des seuils optionelle,
#il faut chercher quelles combinaisons ont été utilisées dans l'article 
#seuil pour lfcThreshold
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

#seuil pour padj
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)


resultsNames(dds)


#coef names à changer
# quel type utiliser : apeglm, ashr, normal ?
res <- lfcShrink(dds, coef="condition_persister_vs_control", type="apeglm") 

plotMA(res, ylim=c(-5,5))
#ggsave("MA_plot.png")


ggplot(as.data.frame(res), aes(x = log2(baseMean), y = log2FoldChange, color = padj < 0.05)) +
  geom_point(aes(shape = translation_gene), size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  scale_shape_manual(values = c(19, 1)) +  
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal() +
  geom_text_repel(aes(label = ifelse(Row.names %in% c("frr", "infA", "infC", "infB", "tsf", "pth"), external_gene_name, "")),
                  color = "black", size = 3)
# récupéreration des noms à la main on peut faire autrement ? 
ggsave("MA_plot_enhanced.png")



