# /public/home/tangy/miniconda3/envs/r413/bin/R
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(ggplot2)
cluster <- "/public/home/zhangzb/test/work/monocle3/DEG/Hypothalamus.csv"
DGE <- read.csv(cluster,
                header = T,
                row.names = 1)
gene_list <- row.names(DGE) 
f <- read.table("/public/home/zhangzb/test/work/zzb/Finalize/GO/development-down.txt",
                sep = "\n")
gene_list <- f$V1      
gene_symbol <- bitr(geneID = gene_list,
                fromType = "SYMBOL",
                toType=c("SYMBOL", "ENTREZID"),
                OrgDb="org.Mm.eg.db")

gene <- gene_symbol[, 2]
CC <- enrichGO(gene = gene,  #基因列表(转换的ID)
        keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
        OrgDb=org.Mm.eg.db,  #物种对应的org包
        ont = "all",   #CC细胞组件，MF分子功能，BF生物学过程
        pvalueCutoff = 0.01,  #p值阈值
        pAdjustMethod = "fdr",  #多重假设检验校正方式
        minGSSize = 1,   #注释的最小基因集，默认为10
        maxGSSize = 500,  #注释的最大基因集，默认为500
        qvalueCutoff = 0.01,  #p值阈值
        readable = TRUE)  #基因ID转换为基因名
kegg <- enrichKEGG(gene = gene, organism = "mmu", pvalueCutoff =0.05, qvalueCutoff =0.2)
# setwd("/public/home/zhangzb/monkey-heart/seurat/goenrich")
# a <- paste0(unlist(strsplit(i, "[.]"))[1], ".pdf")

b <- barplot(CC,
        #GO富集分析结果
        # x = "Count",  #横坐标,默认Count,也可以为GeneRation
        #color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = 10,  #展示前20个，默认为10个
        #size = NULL,  
        #title = "CC_barplot"  #设置图片的标题
        )
b <- barplot(kegg,
        #GO富集分析结果
        # x = "Count",  #横坐标,默认Count,也可以为GeneRation
        #color = "p.adjust",  #右纵坐标,默认p.adjust,也可以为pvalue和qvalue
        showCategory = 10,  #展示前20个，默认为10个
        #size = NULL,  
        #title = "CC_barplot"  #设置图片的标题
        )
ggsave(b, filename = paste0("/public/home/zhangzb/test/work/zzb/Finalize/GO/development-up_GO.pdf"))

        
