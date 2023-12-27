####cox回归分析####
#设置工作目录
setwd("C:/Users/Administrator/Desktop/R_test/R_gittest/cox")
#安装加载R包
# install.packages("survival")
# install.packages("forestplot")
library(survival)
library(forestplot)
library(tidyverse)
#下载生存信息
#xena官网：https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
#读取生存信息tsv文件
surv = read.table(file = 'TCGA-LIHC.survival.tsv', sep = '\t', header = TRUE) 
#整理生存信息数据
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
#读取表达数据
expr <- read.table("LIHC_fpkm_mRNA_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
comgene <- intersect(colnames(expr),rownames(surv))
table(substr(comgene,14,16))
expr <- expr[,comgene]
surv <- surv[comgene,]
#表达数据整理完毕
#读取tcga差异分析结果
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 2, padj < 0.05)
#整合
deg_expr <- expr[rownames(res_deseq2),] %>% t() %>% as.data.frame()
surv.expr <- cbind(surv,deg_expr)
Coxoutput <- NULL 
for(i in 3:ncol(surv.expr)){
  g <- colnames(surv.expr)[i]
  cox <- coxph(Surv(OS.time,OS) ~ surv.expr[,i], data = surv.expr) # 单变量cox模型
  coxSummary = summary(cox)
  
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1],
                                           lower = as.numeric(coxSummary$conf.int[,3][1]),
                                           upper = as.numeric(coxSummary$conf.int[,4][1]),
                                           stringsAsFactors = F),
                                stringsAsFactors = F)
}


write.table(Coxoutput, file = "cox results.txt",sep = "\t",row.names = F,col.names = T,quote = F)
###筛选top基因
pcutoff <- 0.001
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),] # 取出p值小于阈值的基因
topgene <- topgene[1:10,]

#3. 绘制森林图
##3.1 输入表格的制作
tabletext <- cbind(c("Gene",topgene$gene),
                   c("HR",format(round(as.numeric(topgene$HR),3),nsmall = 3)),
                   c("lower 95%CI",format(round(as.numeric(topgene$lower),3),nsmall = 3)),
                   c("upper 95%CI",format(round(as.numeric(topgene$upper),3),nsmall = 3)),
                   c("pvalue",format(round(as.numeric(topgene$p),3),nsmall = 3)))
##3.2 绘制森林图
forestplot(labeltext=tabletext,
           mean=c(NA,as.numeric(topgene$HR)),
           lower=c(NA,as.numeric(topgene$lower)), 
           upper=c(NA,as.numeric(topgene$upper)),
           graph.pos=5,# 图在表中的列位置
           graphwidth = unit(.25,"npc"),# 图在表中的宽度比
           fn.ci_norm="fpDrawDiamondCI",# box类型选择钻石
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),# box颜色
           
           boxsize=0.4,# box大小固定
           lwd.ci=1,
           ci.vertices.height = 0.1,ci.vertices=T,# 显示区间
           zero=1,# zero线横坐标
           lwd.zero=1.5,# zero线宽
           xticks = c(0.5,1,1.5),# 横坐标刻度根据需要可随意设置
           lwd.xaxis=2,
           xlab="Hazard ratios",
           txt_gp=fpTxtGp(label=gpar(cex=1.2),# 各种字体大小设置
                          ticks=gpar(cex=0.85),
                          xlab=gpar(cex=1),
                          title=gpar(cex=1.5)),
           hrzl_lines=list("1" = gpar(lwd=2, col="black"), # 在第一行上面画黑色实线
                           "2" = gpar(lwd=1.5, col="black"), # 在第一行标题行下画黑色实线
                           "12" = gpar(lwd=2, col="black")), # 在最后一行上画黑色实线
           lineheight = unit(.75,"cm"),# 固定行高
           colgap = unit(0.3,"cm"),
           mar=unit(rep(1.5, times = 4), "cm"),
           new_page = F
)
dev.off()