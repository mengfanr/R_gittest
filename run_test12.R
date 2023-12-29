####TCGA差异分析热图####
setwd("xena")
library(tidyverse)
exp <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

logFC_cutoff <- 1
type1 = (DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)
type2 = (DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
table(DEG$change)
library(pheatmap)
cg = rownames(DEG)[DEG$change !="NOT"]
exp_diff <- exp[cg,]
group_list=factor(ifelse(substr(colnames(exp),14,16) == "01A","T","N"),levels = c("N","T"))
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(exp_diff)
pheatmap(exp_diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()