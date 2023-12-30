# 学生的脑子要是跟代码一样听话就好了 还是机器对我最好 我只需要承认我错了就行
####cibersort####
setwd("cibersort")   
install.packages('e1071')
install.packages('parallel')
#install.packages("BiocManager")
BiocManager::install("preprocessCore", version = "3.13")
library(e1071)
library(parallel)
library(preprocessCore)
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'LIHC_fpkm_mRNA_01A.txt'   #肿瘤患者表达谱
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersort.Rdata")   #保存中间文件

load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞
#可视化（阿琛老师）
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框写
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板