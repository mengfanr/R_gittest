# day_number = 2
# day_number_base = as.character(day_number)
# day_base = "DAY"
# day_dir = paste(day_base,day_number_base,sep='')
# if (!dir.exists(day_dir))
#   dir.create(day_dir)
# install.packages(tidyverse)
library(tidyverse)
counts = read.table(file='TCGA-LIHC.htseq_counts.tsv', sep='\t', header=TRUE)
rownames(counts)  = counts[,1] 
counts = counts[,  -1]


