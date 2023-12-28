####计算患者免疫评分与肿瘤纯度#####
setwd("TCGA_ESTIMATE")  #设置工作目录
#安装包
library(utils) #这个包应该不用下载，自带的
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)
#读取肿瘤患者01A表达谱
expr <- read.table("LIHC_fpkm_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)


#计算免疫评分
filterCommonGenes(input.f = "LIHC_fpkm_mRNA_01A.txt",   #输入文件名
                  output.f = "LIHC_fpkm_mRNA_01A.gct",   #输出文件名
                  id = "GeneSymbol")   #行名为gene symbol
estimateScore("LIHC_fpkm_mRNA_01A.gct",   #刚才的输出文件名
              "LIHC_fpkm_mRNA_01A_estimate_score.txt",   #新的输出文件名（即估计的结果文件）
              platform="affymetrix")   #默认平台

#3. 输出每个样品的打分
result <- read.table("LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result <- result[,-1]   
colnames(result) <- result[1,]   
result <- as.data.frame(t(result[-1,]))

rownames(result) <- colnames(expr)
write.table(result, file = "LIHC_fpkm_mRNA_01A_estimate_score.txt",sep = "\t",row.names = T,col.names = NA,quote = F) # 保存并覆盖得分

