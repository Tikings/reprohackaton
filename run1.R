library(ggplot2)
library(DESeq2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(httr)
library(readr)
library(tidyverse)
library(EnrichmentBrowser)


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
  geom_point(aes(color = significant), size = 2.0, alpha=0.5) + 
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

res <- lfcShrink(dds, coef="condition_persister_vs_control", type="apeglm") 
# quel type utiliser : apeglm, ashr, normal ?


##### à suivre 



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




kegg_pathway_ids <- c("sao03010", "sao00970", "sao03013", "sao03015", "sao03008")
kegg_gs <- getGenesets(org = "sao", db = "kegg", cache = FALSE)

selected_gene_sets <- kegg_gs[grep(paste(kegg_pathway_ids, collapse = "|"), names(kegg_gs))]

additional_genes <- c("SAOUHSC_01203")
translation_gene_list <- unique(c(unlist(selected_gene_sets), additional_genes))
aatRNA_gene_list <- unique(unlist(selected_gene_sets["sao00970_Aminoacyl-tRNA_biosynthesis"]))


selected_gene_sets_df_AtRNA <- as.data.frame(selected_gene_sets$`sao00970_Aminoacyl-tRNA_biosynthesis`)
selected_gene_sets_df_ribo <- as.data.frame(selected_gene_sets$sao03010_Ribosome)


translation_genes_df <- res  %>%  filter(Gene_ID %in% translation_gene_list)



# Define thresholds
padj_threshold <- 0.05
logfc_threshold <- 1




# Filter significant genes based on thresholds
significant_translation_genes <- translation_genes_df %>%
  filter(padj < padj_threshold & abs(log2FoldChange) > logfc_threshold) %>%
  pull(row.names)

# Update res to focus on translation-related genes only
res <- translation_genes_df %>%
  mutate(
    Significant = ifelse(row.names %in% significant_translation_genes, "Significant", "Non-significant"),
    Translation_label = ifelse(Gene_ID %in% translation_gene_list, symbol, NA),
    AAtRNA_label = ifelse(Gene_ID %in% aatRNA_gene_list, "AAtRNA synthetase", NA)
  )

getOption("max.print")
res
getOption("max.print")
print(translation_genes_df$pan_gene_symbol)

print(res$pan_gene_symbol)










p <- ggplot(res, aes(x = log2(baseMean), y = log2FoldChange, color = Significant)) +
  geom_point(size = 2) +  
  scale_color_manual(values = c("Non-significant" = "grey", "Significant" = "red")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = expression(Log[2] * " Base Mean"),
    y = expression(Log[2] * " Fold Change"),
    title = "MA-plot of genes related to translation"
  ) +
  
  theme_minimal() +
  theme(
    text = element_text(size = 14),  
    legend.position = "bottom" ) +
  
  ylim(c(-6, 6)) +  
  xlim(c(0, 20))  



p <- p +
  geom_text_repel(
    data = filter(res, !is.na(Translation_label) & Significant == "Significant"),
    aes(label = Translation_label),
    size = 3, max.overlaps = 10
  )

p <- p +
  geom_point(
    data = filter(res, !is.na(AAtRNA_label)),
    shape = 21, color = "black", fill = NA, size = 4, stroke = 1.2
  )


print(p)

ggsave(p,"MA_plot_enhanced.png")
