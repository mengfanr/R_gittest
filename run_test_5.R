library(tidyverse)
setwd("xena")
fpkm1 = read.table(file = 'TCGA-LIHC.htseq_fpkm.tsv', sep = '\t', header = TRUE) 
rownames(fpkm1) <- fpkm1[,1]  
fpkm1 = fpkm1[,-1]
table(substr(colnames(fpkm1),14,16))
fpkm1 <- fpkm1[,substr(colnames(fpkm1),14,16)%in% c("01A","11A")]
table(substr(colnames(fpkm1),14,16))
rownames(fpkm1) <- substr(rownames(fpkm1),1,15)
fpkm <- fpkm1