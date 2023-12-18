# # day_number = 2
# # day_number_base = as.character(day_number)
# # day_base = "DAY"
# # day_dir = paste(day_base,day_number_base,sep='')
# # if (!dir.exists(day_dir))
# #   dir.create(day_dir)
# # install.packages(tidyverse)
# setwd("xena")
# library(tidyverse)
# counts = read.table(file='TCGA-LIHC.htseq_counts.tsv', sep='\t', header=TRUE)
# rownames(counts)  = counts[,1] 
# counts = counts[,  -1]
# counts = counts[, substr(colnames(counts), 14, 16) %in% c("01A","11A")]
# rownames(counts) = substr(rownames(counts), 1, 15)
# counts = ceiling(2^(counts)-1)
# write.table(counts, "counts1.txt", sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)
# # str_tmp = table(substr(colnames(counts), 14, 16), exclude=c('01B', '02A','02B'))
# # # str_tmp = table(substr(colnames(str_tmp), 14, 16), exclude='02A')
# # # str_tmp = table(substr(colnames(str_tmp), 14, 16), exclude='02B')
# # # table()

#########################################
setwd("xena")
counts = read.table("counts1.txt", sep = "\t", check.names = 0,
                    stringsAsFactors = TRUE, header = 1, row.names = 1)
G_info0 = read.table("gene_length_table.txt", sep = "\t", check.names = 0,
                     stringsAsFactors = TRUE, header = 1, row.names = 1)
G_info = G_info0[which(G_info0$genetype=="protein_coding"), ]
comgene = intersect(rownames(G_info), rownames(counts))
counts = counts[comgene, ]
counts = ceiling(2^counts-1)
counts$Gene = as.character(G_info$genename)
counts = counts[!duplicated(counts$Gene),]
rownames(counts) = counts$Gene
counts = counts[,-ncol(counts)]
write.table(counts, file = "LIHC_counts_mRNA_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
tumor = colnames(counts)[substr(colnames(counts),14,16) == "01A"]
counts_01A = counts[, tumor]
write.table(counts, file = "LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
