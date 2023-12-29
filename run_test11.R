setwd("C:/Users/Administrator/Desktop/R_test/R_gittest/timeROC")
#R包
install.packages("timeROC")
install.packages("survival")
library(timeROC)
library(survival)
library(tidyverse)

#2.2 数据的整理与载入
exp_sur <- read.table("exp_sur.txt", header=T,sep="\t", check.names=F, row.names=1)
exp_sur$OS.time <- exp_sur$OS.time/365
exp_sur_01A <- exp_sur[substr(rownames(exp_sur),14,16) == "01A",]
write.table(exp_sur_01A, file = "exp_sur_01A.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#2.3 构建ROC曲线函数
ROC3 <- timeROC(T=exp_sur_01A$OS.time,   #结局时间
                delta=exp_sur_01A$OS,   #结局指标
                marker=exp_sur_01A$SPP1,   #预测变量
                cause=1,   #阳性结局指标数值
                weighting="marginal",   #计算方法，默认为marginal
                times=c(1, 3, 5),   #时间点，选取1年，3年和5年的生存率
                iid=TRUE)
ROC3   #查看模型变量信息

#2.4 绘制ROC曲线
plot(ROC3,
     time=1, col="red")   #time是时间点，col是线条颜色
plot(ROC3,
     time=3, col="green", add=TRUE)   #add指是否添加在上一张图中
plot(ROC3,
     time=5, col="blue", add=TRUE)
legend("bottomright",
       c("Year-1", "Year-3", "Year-5"),
       col=c("red", "green", "blue"),
       lty=1, lwd=2)   #添加标签信息

dev.off()