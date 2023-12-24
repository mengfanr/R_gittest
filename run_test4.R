setwd("xena")
counts_01A = read.table("LIHC_counts_mRNA_01A.txt", sep = "\t", check.names = 0,
                    stringsAsFactors = TRUE, header = 1, row.names = 1)

counts = read.table("LIHC_counts_mRNA_all.txt", sep = "\t", check.names = 0,
                    stringsAsFactors = TRUE, header = 1, row.names = 1)
library(tidyverse)
#安装BiocManager
if(!require(DESeq2))BiocManager::install('DESeq2')
library(DESeq2)
counts = counts[apply(counts, 1, function(x) sum(x > 1) > 32), ]# stupid function fuck, for each row, count the number of elements that is larger than 1, if the number is greater than 32, get it.
counts_tmp = apply(counts, 1, function(x) sum(x > 1) > 32)
conditions=data.frame(sample=colnames(counts),
                      group=factor(ifelse(substr(colnames(counts),14,16) == "01A","T","N"))) %>% 
  column_to_rownames("sample")
dds = DESeqDataSetFromMatrix(countData = counts, colData = conditions, design = ~ group)
dds <- DESeq(dds)
resultsNames(dds)
res = results(dds)
save(res,file = "LIHC_DEG.rda")
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 3, padj < 0.05)#根据自己需要 pvalue statiscal meaning and abs(log2FoldChange)-> comparing with normal expression, how many times this gene has been expressed
