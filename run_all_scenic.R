setwd("/public/home/zhangzb/test/work/SCENIC/all_raw_sum")
library(SCENIC)
library(SCopeLoomR)
library(doParallel)
# loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")
# loom <- open_loom(loomPath, mode = 'r+')
# exprMat <- get_dgem(loom)

# cellInfo <- get_cell_annotation(loom)

exprMat <- as.matrix(read.csv('/public/home/zhangzb/test/work/SCENIC/raw_data/all_raw_sum.csv',
                               header = T, row.names = 1)) #不能有全为0的列
# dim(exprMat)
# exprMat <- exprMat[, -c(grep("E175B", colnames(exprMat)))]
# exprMat <- exprMat[, -c(grep("meninges", colnames(exprMat)))]

scenicOptions <- initializeScenic(org="mgi", dbDir = "//public//home//zhangzb//test//work//cisTarget_databases//V1" , nCores=10)
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- exprMat_filtered
exprMat_filtered_log <- log2(exprMat_filtered + 1)
exprMat_filtered_log[is.na(exprMat_filtered_log)] <- 0 #去掉nan值
runGenie3(exprMat_filtered_log, scenicOptions)
#SCENIC流程其实是包装好的4个函数，
exprMat_log <- exprMat
exprMat_log <- log2(exprMat + 1)
exprMat_log[is.na(exprMat_log)] <- 0
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
scenicOptions <-runSCENIC_1_coexNetwork2modules(scenicOptions) #GENIE3/GRNBoost：基于共表达情况鉴定每个TF的潜在靶点

library(BiocParallel)
scenicOptions <-runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #cisTarget：基于DNA-motif 分析选择潜在的直接结合靶点；

library(foreach)
scenicOptions <- initializeScenic(org="mgi", dbDir = "//public//home//zhangzb//test//work//cisTarget_databases//V1" , nCores=1)
scenicOptions <-runSCENIC_3_scoreCells(scenicOptions, exprMat_log)#AUCell：分析每个细胞的regulons活性 

f1 <- readRDS("int/3.4_regulonAUC.Rds")
sum <- f1@assays@data@listData$AUC
write.csv(sum,"int/all_raw_sum.csv")

scenicOptions <-runSCENIC_4_aucell_binarize(scenicOptions)#细胞聚类：基于regulons的活性鉴定稳定的细胞状态并对结果进行探索
tsneAUC(scenicOptions, aucType = "AUC")
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file = "int/scenicOptions.Rds")

