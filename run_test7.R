setwd("C:/Users/Administrator/Desktop/R_test/R_gittest/GSE84402")
###加载R包
library(tidyverse)
library(GEOquery)
###下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
###提取子集
gset[[1]]
#通过pData函数获取分组信息
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
#设置参考水平
group_list <- ifelse(str_detect(pdata$source_name_ch1, "hepatocellular carcinoma"), "tumor",
                     "normal")
#因子型
group_list = factor(group_list,
                    levels = c("normal","tumor"))

##读取上节课整理好的表达数据exp##
exp <- read.table("exp.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
#差异分析
library(limma)
design=model.matrix(~group_list)
##### 直接跑 就是个线性最小二乘 ####
fit=lmFit(exp,design)#在给定一系列阵列的情况下，拟合每个基因的线性模型
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)#eBayes函数指在微阵列线性模型拟合下，通过经验Bayes对标准误差向一个共同值的缓和，计算缓和t-统计、缓和f-统计和微分表达式的对数概率。 用法举列：EBAYES（配合，比例=0.01，标准偏差，系数，极限值=C（0.1,4） 
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
##标记上下调基因
logFC=1
P.Value = 0.05
k1 = (deg$P.Value < P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)

##热图无脑运行 走你##
cg = rownames(deg)[deg$change !="stable"]
diff=exp[cg,]
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(diff) 
pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=3,
         fontsize_col=3)
dev.off()
