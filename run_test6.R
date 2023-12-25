####GSE84402####
setwd("GSE84402")
###加载R包
library(tidyverse)
chooseBioCmirror()
BiocManager::install('GEOquery')
library(GEOquery)
###下载数据，如果文件夹中有会直接读入
# chooseBioCmirror()
gset = getGEO('GSE84402', destdir=".", AnnotGPL = F, getGPL = F)
class(gset)
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
##2.2 通过exprs函数获取表达矩阵并校正
exp <- exprs(gset[[1]])
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()
###数据校正
library(limma) 
exp=normalizeBetweenArrays(exp)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=2)
range(exp)
exp <- log2(exp+1)
range(exp)
dev.off()
#使用R包转换id
index = gset[[1]]@annotation
if(!require("hgu133plus2.db"))#平台对应平台对应的包的名称
  BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
ls("package:hgu133plus2.db")
ids <- toTable(hgu133plus2SYMBOL)
head(ids)
#length(unique(ids$symbol))
#table(sort(table(ids$symbol)))
#id转换
library(tidyverse)
exp <- as.data.frame(exp)
exp <- exp %>% mutate(probe_id=rownames(exp))
exp <- exp %>% inner_join(ids,by="probe_id") 
exp <- exp[!duplicated(exp$symbol),]
rownames(exp) <- exp$symbol
exp <- exp[,-(29:30)]
write.table(exp, file = "exp.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####GEO手动注释####
####GSE31056####
setwd("GSE31056")
###加载R包
library(tidyverse)
BiocManager::install('GEOquery')
library(GEOquery)
###下载数据，如果文件夹中有会直接读入
chooseBioCmirror()
gset = getGEO('GSE31056', destdir=".", AnnotGPL = F, getGPL = F)
#有时会报错  Increase it by setting `Sys.setenv("VROOM_CONNECTION_SIZE")`
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
class(gset)
###提取子集
gset[[1]]
#读取表达谱
exp <- exprs(gset[[1]])
#把表达谱转为数据框格式
exp <- as.data.frame(exp)
##转换id
#读取GPL文件
comname <- intersect(rownames(exp),rownames(GPL))
exp <- exp[comname,]
GPL <- GPL[comname,]
exp1 <- cbind(GPL,exp)
exp1 <- exp1[!duplicated(exp1$SYMBOL),]
rownames(exp1) <- exp1$SYMBOL
exp1 <- exp1[,-(1:5)]
write.table(exp1, file = "exp1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)