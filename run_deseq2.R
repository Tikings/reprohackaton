



library(DESeq2)
library(tximport)
library(apeglm)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(biomaRt)
library(httr)
library(ggrepel)



count_data <- as.matrix(read.csv("counts.csv", row.names = 1))
print(count_data)

# Charger les métadonnées des échantillons
col_data <- read.csv("coldata.csv", row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
summary(res)




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
res <- lfcShrink(dds, coef="condition_treated_vs_control", type="apeglm") 
# quel type utiliser : apeglm, ashr, normal ?

resultsNames(dds)




plotMA(res, ylim=c(-5,5))
ggsave("MA_plot.png")


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



