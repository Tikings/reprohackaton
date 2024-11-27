


library(ggplot2)
library(DESeq2)
library(dplyr)
library(ggrepel)
library(pheatmap)
library(httr)

library(readr)
library(tidyverse)

library(EnrichmentBrowser)






count_data_last <- as.matrix(read.table("./counts-last.txt", 
                                        header = TRUE, 
                                        row.names = 1, 
                                        sep = "\t", 
                                        skip = 1)[, -c(1:5)])


sample_names_last<- colnames(count_data_last)
new_names <- sub(".*?(SRR[0-9]+).*", "\\1", sample_names_last)
colnames(count_data_last) <- new_names
count_data<- count_data_last




sample_names <- colnames(count_data)
sample_conditions <- c(
  "SRR10379721" = "persister",
  "SRR10379722" = "persister",
  "SRR10379723" = "persister",
  "SRR10379724" = "control",
  "SRR10379725" = "control",
  "SRR10379726" = "control" )

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


p<-ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 2.0, alpha=0.5) + 
  scale_color_manual(values = c("Significant" = "red","Non-significant" = "black")) +
  scale_x_log10() +  
  geom_hline(yintercept = 0, linetype = "dashed") +  
  theme_minimal() +  
  labs(x = "Mean of normalized counts", y = "Log2 fold change") +
  theme(
    legend.position = "none",  
  ) + ylim(c(-4, 4))


ggsave("MA_plot_1.pdf",plot=p,device='pdf', width=10,height=8)





# pas de lfc shrink du tout ? l'image sans semble bcp plus proche de celle du papier que celle avec

res <- lfcShrink(dds, coef="condition_persister_vs_control", type="apeglm") 
# quel type utiliser : apeglm, ashr, normal ?


##### à suivre 

gene_annotations <- read_csv("NCTC8325.csv")

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


res <- left_join(res, gene_annotations, by = c("row.names" = "locus_tag"))

# faire la fusion 



## plus que KEGG à faire 


pwys <- download.kegg.pathways("sao")
gs <- get.kegg.genesets(pwys)

sao00970<-gs$`sao00970_Aminoacyl-tRNA_biosynthesis`
sao03010<-gs$sao03010_Ribosome


sao03012 <- c("SAOUHSC_01246", "SAOUHSC_01786", "SAOUHSC_02274",
                "SAOUHSC_01234", "SAOUHSC_00530", "SAOUHSC_00529",
                "SAOUHSC_02359", "SAOUHSC_01236", "SAOUHSC_02358",
                "SAOUHSC_00475", "SAOUHSC_02489", "SAOUHSC_01625",
                "SAOUHSC_00771", "SAOUHSC_00956" )
  






res <- as.data.frame(res)

res <- res %>%
  rename( gene_id= row.names)

res <- res %>%
  mutate(
    GeneType = case_when(
      gene_id %in% sao03012 ~ "Typical Member",
      gene_id %in% sao00970 ~ "AA-tRNA Synthetase",
      gene_id %in% sao03010 ~ "Ribosome",
      TRUE ~ "Other"
    ),
    Significant = ifelse(padj < 0.05, "Significant", "Non-significant")
  )

res_filtered <- res %>%
  filter(GeneType != "Other")


p2 <- ggplot(res_filtered, aes(x = log2(baseMean), y = log2FoldChange, color = Significant)) +
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


p2 <- p2 +
  geom_text_repel(
    data = filter(res_filtered, GeneType == "Typical Member"& Significant == "Significant"), 
    aes(label = pan_gene_symbol),  
    size = 6, 
    max.overlaps = 20,
    color = "black",              
    box.padding = 0.8,            
    point.padding = 0.5,          
    segment.color = "black",      
    segment.size = 0.5            
  
  )

p2 <- p2 +
  geom_point(
    data = filter(res_filtered, GeneType == "AA-tRNA Synthetase"),  
    shape = 21, color = "black", fill = NA, size = 3, stroke = 1.2
  )


ggsave("MA_plot_enhanced.pdf",plot=p2,device='pdf', width=10,height=8)














