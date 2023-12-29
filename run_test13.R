####TCGA差异分析火山图####
# setwd("xena")
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


install.packages("ggpubr")
install.packages("ggthemes")
library(ggpubr)
library(ggthemes)

DEG$logP <- -log10(DEG$padj)
ggscatter(DEG,
          x = "log2FoldChange", y = "logP") +
  theme_base()

#增加基因上下调信息
ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base()

#添加分界线
ggscatter(DEG, x = "log2FoldChange", y = "logP", xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1) +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")
dev.off()

#添加基因标签信息
DEG$Label = ""   #新加一列label
DEG <- DEG[order(DEG$padj), ]   #对差异基因的p值进行从小到大的排序
DEG$Gene <- rownames(DEG)
#高表达的基因中，选择fdr值最小的5个
up.genes <- head(DEG$Gene[which(DEG$change == "UP")], 5)
#低表达的基因中，选择fdr值最小的5个
down.genes <- head(DEG$Gene[which(DEG$change == "DOWN")], 5)
#将up.genes和down.genes合并，并加入到Label中
DEG.top5.genes <- c(as.character(up.genes), as.character(down.genes))
DEG$Label[match(DEG.top5.genes, DEG$Gene)] <- DEG.top5.genes

ggscatter(DEG, x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("blue", "black", "red"),
          size = 1,
          label = DEG$Label,
          font.label = 8,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(Adjust P-value)") +
  theme_base() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

dev.off()