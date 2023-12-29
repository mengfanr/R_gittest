####ROC####
#读取生存信息tsv文件
chooseBioCmirror()
setwd("C:/Users/Administrator/Desktop/R_test/R_gittest/ROC")
library(tidyverse)
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
#保存整理好的生存信息
write.table(surv, file = "survival.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#提取上次cox作图的10个基因
gene <- c("SPP1","PAGE1","G6PD","MAGEA4",'CDCA8',
          'TRIM54','KIF2C','KIF20A','ANLN',"SLC7A11")
exp10 <- expr[gene,] %>% t() %>% as.data.frame()
#整合表达谱与生存信息
exp_sur <- cbind(exp10,surv)
write.table(exp_sur, file = "exp_sur.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#准备R包
install.packages("ROCR")
install.packages("rms")
library(ROCR)
library(rms)
#构建ROC预测模型
ROC1 <- prediction(exp_sur$SPP1,exp_sur$OS)   #构建ROC预测模型 SPP1 changable
ROC2 <- performance(ROC1,"tpr","fpr")   #计算预测模型的TPR/FPR值
AUC <- performance(ROC1,"auc")   #计算曲线下面积(AUC)值

AUC<- 0.5604839 #改 根据结果对AUC进行赋值

#1.4 绘制ROC曲线
plot(ROC2,
     col="red",   #曲线的颜色
     xlab="False positive rate", ylab="True positive rate",   #x轴和y轴的名称
     lty=1,lwd=3,
     main=paste("AUC=",AUC))
abline(0, 1, lty=2, lwd=3)   #绘制对角线
dev.off()
