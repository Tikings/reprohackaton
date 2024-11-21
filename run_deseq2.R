if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}




# Install Bioconductor packages
BiocManager::install("DESeq2")
BiocManager::install("tximport")
BiocManager::install("apeglm")

library(ggplot2)
library(DESeq2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(httr)


library(tximport)
library(apeglm)
library(biomaRt)



count_data <- as.matrix(read.table("./counts_final.txt", 
                                   header = TRUE, 
                                   row.names = 1, 
                                   sep = "\t", 
                                   skip = 1)[, -c(1:5)])


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
  condition = factor(sample_conditions[sample_names]) )

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data
                              , design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res
#plotMA(res, ylim=c(-4,4))


res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05,
                             "Significant",
                             "Non-significant")


ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 1.0, alpha=0.5) + 
  scale_color_manual(values = c("Significant" = "red","Non-significant" = "black")) +
  scale_x_log10() +  
  geom_hline(yintercept = 0, linetype = "dashed") +  
  theme_minimal() +  
  labs(x = "Mean of normalized counts", y = "Log2 fold change") +
  theme(
    legend.position = "none",  
  ) + ylim(c(-4, 4)) 


ggsave("MA_plot_1.png")





# pas de lfc shrink du tout ? l'image sans semble bcp plus proche de celle du papier que celle avec

#res <- lfcShrink(dds, coef="condition_persister_vs_control", type="apeglm") 
# quel type utiliser : apeglm, ashr, normal ?


##### à suivre 

library(readr)
library(tidyverse)

gene_annotations <- read_csv("NCTC8325.csv", show_col_types = FALSE)

gene_annotations <- gene_annotations %>%
  separate(
    col = `locus tag;pan gene symbol;symbol;synonym;Gene ID` ,
    into = c("locus_tag", "pan_gene_symbol", "symbol", "synonym", 
             "Gene_ID"),
    sep = ";",
    fill = "right"
  )





res <- as.data.frame(res)

res <- res %>% rownames_to_column(var = "row.names")

res$row.names <- sub("^gene-", "", res$row.names)

res

res <- left_join(res, gene_annotations, by = c("row.names" = "locus_tag"))

res

## plus que KEGG à faire 


