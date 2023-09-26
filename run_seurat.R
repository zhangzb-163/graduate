library(Seurat)
library(stringr)
library(ggplot2)
library(data.table)
library(SeuratObject)
library(pheatmap)
library(dplyr)
library(RColorBrewer)

se <- Load10X_Spatial("/public/home/hongyz/workspace/mouse-brain-full/spaceranger/E165A/outs")
brain <- SCTransform(se, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:50)
brain <- FindClusters(brain, verbose = FALSE, resolution = 5)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:50)
SpatialDimPlot(brain, label = TRUE, label.size = 3, alpha = 0.6)
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(15, 26)), 
                                        facet.highlight = TRUE, ncol = 3)

DefaultAssay(brain) <- "SCT"
DEG <- FindMarkers(brain, ident.1 = c(15), group.by = "seurat_clusters")
DEG <- DEG[which((DEG$avg_log2FC >= 1 | DEG$avg_log2FC <= -1)
                    & DEG$p_val_adj <= 0.05), ]
DEG[order(DEG$avg_log2FC, decreasing = T), ]
VlnPlot(brain, features = c("Tac1"), group.by = "seurat_clusters") #Tac1 Eomes Gpr88 Magel2

#亚聚类
brain@meta.data$seurat_clusters <- as.vector(brain@meta.data$seurat_clusters)
brain <- FindSubCluster(brain, cluster = 6, graph.name = "SCT_snn", resolution = 0.5)
brain@meta.data[grep('_',brain$sub.cluster),'seurat_clusters'] <- brain@meta.data[grep('_',brain$sub.cluster),'sub.cluster']

write.csv(brain@meta.data[,c(6,7)], "/public/home/zhangzb/test/work/tem-duc/test-seurat.csv")
write.csv(brain@meta.data[,c(6,7)], "/public/home/zhangzb/test/work/raw_data/seurat-cluster/E155A-seurat.csv")


a$seurat_clusters <- as.vector(a$seurat_clusters)
a[grep('_',a$sub.cluster),'seurat_clusters'] <- a[grep('_',a$sub.cluster),'sub.cluster']



#整合数据
sample <- c("E135A", "E135B", "E155A", "E155B", "E165A", "E165B",
            "E175A1", "E175A2", "P0A1", "P0A2")
for(i in sample){
    f <- read.csv(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                          i, "_rawcount.csv"),header = T, row.names = 1)
    f <- CreateSeuratObject(f, project = i)
    if(i == "E135A"){
        seurat.list <- f
    }
    else {
        seurat.list <- merge(seurat.list, f)
    }
}
seurat.list <- SplitObject(seurat.list, split.by = "orig.ident")
seurat.list2 <- lapply(X = seurat.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = seurat.list2, nfeatures = 3000)
seurat.list2 <- PrepSCTIntegration(object.list = seurat.list2, anchor.features = features)
seurat.list2 <- lapply(X = seurat.list2, FUN = RunPCA, features = features)
seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list2, normalization.method = "SCT",
anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
seurat_combined <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", dims = 1:30)
seurat.data <- ScaleData(seurat_combined, verbose = FALSE)
seurat.data <- RunPCA(seurat.data, npcs = 30, verbose = FALSE,  features = features)
seurat.data <- RunUMAP(seurat.data, reduction = "pca", dims = 1:30)
seurat.data <- FindNeighbors(seurat.data, reduction = "pca", dims = 1:30)
seurat.data <- FindClusters(seurat.data, resolution = 0.5)
# p1 <- DimPlot(seurat.data, reduction = "umap", group.by = "orig.ident")
# p2 <- DimPlot(seurat.data, reduction = "umap", group.by = "seurat_clusters")
# p1 + p2
saveRDS(seurat.data, "/public/home/zhangzb/test/work/seurat/all_SCT_integrate.rds")

# seurat.data@meta.data$development <- substr(seurat.data@meta.data$orig.ident, start = 1, stop = 3)
# seurat.data@meta.data$sample <- seurat.data@meta.data$orig.ident 
seurat.data <- readRDS("/public/home/zhangzb/test/work/seurat/all_SCT_integrate.rds")
DefaultAssay(seurat.data) <- "integrated"
seurat.data <- PrepSCTFindMarkers(seurat.data)


FeaturePlot(seurat.data)
# DEG <- DEG[which((DEG$avg_log2FC >= 1 | DEG$avg_log2FC <= -1)
#                     & DEG$p_val_adj <= 0.01), ]
# DEG <- DEG[order(DEG$avg_log2FC, decreasing = T), ]

regions <- c("Cortex", "Thalamus", "Hypothalamus", "Striatum", "SVZ", "Meninges", "Amygdalar",
            "Olfactory", "Hippocampus", "GE", "CXTsp", "Habenular", "Choroid", "3V", "Noregion")
for(i in regions){
    DEG <- FindMarkers(seurat.data, ident.1 = i, group.by = "development")
    write.csv(DEG, paste0("/public/home/zhangzb/test/work/zzb/integrated_DEG/",i ,"_DEG.csv")) 
}
#绘制region DEG 热图
regions <- c("Cortex", "Thalamus", "Hypothalamus", "Striatum", "SVZ", "Meninges", "Amygdalar",
            "Olfactory", "Hippocampus", "GE", "CXTsp", "Habenular", "Choroid", "3V")
up_gene <- c()
for(i in regions){
    DEG <- read.csv(paste0("/public/home/zhangzb/test/work/zzb/integrated_DEG/",i ,"_DEG.csv"), header = T, row.names = 1)
    DEG <- DEG[which((DEG$avg_log2FC >= 1 | DEG$avg_log2FC <= -1) & DEG$p_val_adj <= 0.01), ]
    # DEG <- DEG[which(DEG$pct.2 <= 0.3),]
    DEG <- DEG[order(DEG$avg_log2FC, decreasing = T),]
    if(dim(DEG)[1] > 15){
        DEG <- DEG[c(1:15),]
    }
    up_gene <- c(up_gene, row.names(DEG))

    
}
up_gene <- unique(up_gene)
count <- seurat.data@assays$integrated@data[up_gene,]
# count <- log10(count + 1)
annotation_col <- seurat.data@meta.data[,c(8, 10)]
annotation_col <- annotation_col[which(annotation_col$region != "Noregion"),]
library(forcats)
annotation_col <- annotation_col %>% mutate(x=fct_relevel(region, regions))%>% 
                arrange(x)
annotation_col <- annotation_col[,-3]
count <- count[,row.names(annotation_col)]
colors <- colorRampPalette(colors = c("blue","white","red"))(50)
p <- pheatmap(count, cluster_rows = FALSE, cluster_cols = FALSE, color = colors,
                            annotation_col=annotation_col, show_colnames=FALSE,
                            show_rownames = FALSE, fontsize_row = 2,scale = "row",
                            )
ggsave(p, filename = "/public/home/zhangzb/test/work/zzb/Finalize/integrate_region_DEG_heatmap2.pdf")

pheatmap(f, cluster_cols = FALSE, scale = "row", color = colors)
# 绘制发育时间点DEG 热图
time <- c("E135", "E155", "E165", "E175", "P0")
up_gene <- c()
for(i in time){
    DEG <- read.csv(paste0("/public/home/zhangzb/test/work/zzb/integrated_DEG/",i ,"_DEG.csv"), header = T, row.names = 1)
    DEG <- DEG[which(DEG$avg_log2FC >= 0.5 & DEG$p_val_adj <= 0.01), ]
    # DEG <- DEG[which(DEG$pct.2 <= 0.3),]
    DEG <- DEG[order(DEG$avg_log2FC, decreasing = T),]
    if(dim(DEG)[1] > 15){
        DEG <- DEG[c(1:15),]
    }
    up_gene <- c(up_gene, row.names(DEG))
}
up_gene <- unique(up_gene)
count <- seurat.data@assays$integrated@data[up_gene,]
annotation_col <- seurat.data@meta.data[,c(8), drop = F]
count <- count[,row.names(annotation_col)]
colors <- colorRampPalette(colors = c("blue","white","red"))(50)
p <- pheatmap(count, cluster_rows = FALSE, cluster_cols = FALSE, color = colors,
                            annotation_col=annotation_col, show_colnames=FALSE,
                            show_rownames = TRUE, fontsize_row = 2,scale = "row",
                            )
ggsave(p, filename = "/public/home/zhangzb/test/work/zzb/test.pdf")


#计算细胞周期

sample <- c("E135A", "E135B", "E155A", "E155B", "E165A", "E165B",
            "E175A1", "E175A2", "P0A1", "P0A2")
for (i in sample){
    f <- read.csv(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                          i, "_rawcount.csv"), header = T, row.names = 1, check.names = FALSE)
    se <- CreateSeuratObject(f)
    pbmc <- NormalizeData(se)
    g2m_genes <- cc.genes$g2m.genes ## 获取G2M期marker基因
    g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(pbmc)) #提取pbmc矩阵中的G2M期marker基因
    s_genes <- cc.genes$s.genes   #获取S期marker基因 
    s_genes <- CaseMatch(search=s_genes, match=rownames(pbmc)) #提取pbmc矩阵中的S期marker基因
    pbmc <- CellCycleScoring(se, g2m.features=g2m_genes, s.features=s_genes)
    pbmc <- ScaleData(pbmc, features = rownames(pbmc))
    pbmc <- RunPCA(pbmc, features = c(s_genes, g2m_genes))
    write.csv(pbmc@meta.data, paste0("/public/home/zhangzb/test/work/seurat/cell_cycle/",
                                     i, "_cellcycle.csv"))
}

for(i in sample){
    f <- read.csv(paste0("/public/home/zhangzb/test/work/seurat/cell_cycle/",
                  i, "_cellcycle.csv"), header = T, row.names = 1, check.names = FALSE)
    if(i == "E135A"){
        cycle <- f
    }
    else{
        cycle <- rbind(cycle, f)
    }
}
cycle <- cycle[, 6, drop = F]

seurat.data@meta.data <- cbind(seurat.data@meta.data,
                              cycle[row.names(seurat.data@meta.data),,drop = F])
seurat.data@meta.data[which(seurat.data@meta.data$Phase == "G1"), "Phase"] <- "None"
colors <- c("blue", "grey", "red")
p <- DimPlot(seurat.data, reduction = "umap", group.by = "Phase", cols = colors)
colors <- if (blend) {     c("lightgrey", "#ff0000", "#00ff00") } else {    
    c("lightyellow", "#ff0000") }
color2 <- c("grey", "lightgrey", "blue")
p <- FeaturePlot(seurat.data, reduction = "umap", features = "nCount_RNA")
ggsave(p, filename = "/public/home/zhangzb/test/work/seurat/cell_cycle/all-nCount_RNA.pdf")
seurat.data <- RunTSNE(seurat.data)

library(showtext)
t <- c("E13.5", "E15.5", "E16.5", "E17.5", "P0")
ratio <- as.data.frame(t)
ratio1 <- ratio
ratio2 <- ratio
k <- 0
for (i in t){
    k <- k+1
    a <- subset(seurat.data@meta.data, development == i)
    print(paste0(i, ":", dim(a)[1]))
    b <- subset(a, Phase == "None")
    c <- (1 - dim(b)[1]/dim(a)[1])*100
    ratio1[k,"tio"] <- c
    ratio1[k, "num"] <- dim(a)[1] - dim(b)[1]
    ratio1$type <- "CellCycle"
    
    ratio2[k,"tio"] <- 100 - c
    ratio2[k, "num"] <- dim(b)[1]
    ratio2$type <- "None"
    
}
ratio3 <- rbind(ratio1, ratio2)
ratio3$type <- factor(ratio3$type, levels=c("CellCycle", "None"))
p <- ggplot(ratio3,aes(t,num, fill=type))+
        geom_bar(stat="identity",position="stack")+
        theme_bw()+
        theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20, face = "bold"), #坐标轴标题
        # axis.line = element_line(size = 5) #坐标轴
        text = element_text(size=20, face = "bold") #坐标轴刻度字体大小 family = "Times New Roman",
        # plot.margin=unit(rep(6,10),'cm') #画布的大小
        )+
        labs(x='Development',y='Number')+
        scale_fill_manual(values = c("purple", "grey"))+
        coord_flip()
        
        
        
ggsave(p, filename = "/public/home/zhangzb/test/work/seurat/cell_cycle/cellcycle_ratio.pdf", width = 20)
# DimPlot(pbmc,group.by = 'Phase')
# SpatialDimPlot(pbmc, group.by = "Phase", label = T, label.size = 3, alpha = 0.6)

cyc <- subset(seurat.data@meta.data, Phase != "None")
none <- subset(seurat.data@meta.data, Phase == "None")

##绘制时空可变基因折线图








#按时期SCT integrate 数据
sample <- c("E175", "E135", "E155", "E165", "P0")
integrate <- function(a, b, c = NULL){
    f1 <- as.data.frame(fread(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                      a, "_rawcount.csv"),
                    header = T, data.table = T))
    row.names(f1) <- f1[, 1]
    f1 <- f1[, -1]
    f1 <- CreateSeuratObject(f1, project = a)
    f2 <- as.data.frame(fread(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                      b, "_rawcount.csv"),
                    header = T, data.table = T))
    row.names(f2) <- f2[, 1]
    f2 <- f2[, -1]
    f2 <- CreateSeuratObject(f2, project = b)
    if (!is.null(c)){
            f3 <- as.data.frame(fread(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                      c, "_rawcount.csv"),
                    header = T, data.table = T))
            row.names(f3) <- f3[, 1]
            f3 <- f3[, -1]
            f3 <- CreateSeuratObject(f3, project = c)
    }
    seurat.list1 <- list(f1, f2)
    if (!is.null(c)){
        seurat.list1 <- list(f1, f2, f3)
    }
    # seurat.list1 <- SplitObject(seurat.list1, split.by = "orig.ident")
    seurat.list2 <- lapply(X = seurat.list1, FUN = SCTransform)
    # if (length(seurat.list2) == 2){
    #     features = intersect(row.names(seurat.list2[[1]]@assays$SCT), row.names(seurat.list2[[2]]@assays$SCT))
    # }
    # if(length(seurat.list2) == 3){
    #    features = intersect(row.names(seurat.list2[[1]]@assays$SCT), row.names(seurat.list2[[2]]@assays$SCT))
    #    features = intersect(features, row.names(seurat.list2[[3]]@assays$SCT))
    # }
    features <- SelectIntegrationFeatures(object.list = seurat.list2, nfeatures = 3000)
    seurat.list2 <- PrepSCTIntegration(object.list = seurat.list2, anchor.features = features)
    seurat.list2 <- lapply(X = seurat.list2, FUN = RunPCA, features = features)
    seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list2, normalization.method = "SCT",
    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
    seurat_combined <- IntegrateData(anchorset = seurat.anchors, normalization.method = "SCT", dims = 1:30)
    op <- list(seurat_combined, features)
    return(op)
}

for (i in sample){
    a <- NULL
    b <- NULL
    c <- NULL
    k <- NULL
    if (i == "E175" ){
        a <- paste0(i, "A1")
        b <- paste0(i, "A2")
        # c <- paste0(i, "B")
        k <- integrate(a, b, c)
        seurat.data <- k[[1]]
        features <- k[[2]]
    }
    if (i == "P0" ){
        a <- paste0(i, "A1")
        b <- paste0(i, "A2")
        # c <- paste0(i, "B")
        k <- integrate(a, b, c)
        seurat.data2 <- k[[1]]
        features2 <- k[[2]]
    }
    if (i != "E175" & i != "P0") {
        a <- paste0(i, "A")
        b <- paste0(i, "B")
        k <- integrate(a, b)
        seurat.data2 <- k[[1]]
        features2 <- k[[2]]
    }
    if (i != "E175"){
        seurat.data <- merge(seurat.data, seurat.data2)
        features <- c(features, features2)
    }
}
features <- unique(features)
seurat.data <- ScaleData(seurat.data, verbose = FALSE)
seurat.data <- RunPCA(seurat.data, npcs = 30, verbose = FALSE,  features = features)
seurat.data <- RunUMAP(seurat.data, reduction = "pca", dims = 1:30)
seurat.data <- FindNeighbors(seurat.data, reduction = "pca", dims = 1:30)
seurat.data <- FindClusters(seurat.data, resolution = 0.5)
# p1 <- DimPlot(seurat.data, reduction = "umap", group.by = "orig.ident")
# p2 <- DimPlot(seurat.data, reduction = "umap", group.by = "seurat_clusters")
# p1 + p2
saveRDS(seurat.data, "/public/home/zhangzb/test/work/seurat/all_SCT_integrate.rds")
# data <- t(as.data.frame(seurat.data@assays$integrated@data))
# write.csv(data, "/public/home/zhangzb/test/work/pyscenic/all_SCT_integrate/all_SCT_all.csv")

#nor
sample <- c("E175", "E135", "E155", "E165", "P0")
integrate <- function(a, b, c = NULL){
    f1 <- as.data.frame(fread(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                      a, "_rawcount.csv"),
                    header = T, data.table = T))
    row.names(f1) <- f1[, 1]
    f1 <- f1[, -1]
    f1 <- CreateSeuratObject(f1, project = a)
    f2 <- as.data.frame(fread(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                      b, "_rawcount.csv"),
                    header = T, data.table = T))
    row.names(f2) <- f2[, 1]
    f2 <- f2[, -1]
    f2 <- CreateSeuratObject(f2, project = b)
    if (!is.null(c)){
            f3 <- as.data.frame(fread(paste0("/public/home/zhangzb/test/work/raw_data/raw_count/",
                                      c, "_rawcount.csv"),
                    header = T, data.table = T))
            row.names(f3) <- f3[, 1]
            f3 <- f3[, -1]
            f3 <- CreateSeuratObject(f3, project = c)
    }
    seurat.list1 <- list(f1, f2)
    if (!is.null(c)){
        seurat.list1 <- list(f1, f2, f3)
    }
    # seurat.list1 <- SplitObject(seurat.list1, split.by = "orig.ident")
    seurat.list2 <- lapply(X = seurat.list1, FUN = function(x){
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    })
    features <- SelectIntegrationFeatures(object.list = seurat.list2)
    anchors <- FindIntegrationAnchors(
    object.list = seurat.list2,
    anchor.features = 32285,
    verbose = FALSE
    )
    seurat_combined <- IntegrateData(anchorset = anchors, verbose = FALSE)
    op <- list(seurat_combined, features)
    return(op)
}

for (i in sample){
    a <- NULL
    b <- NULL
    c <- NULL
    k <- NULL
    if (i == "E175"){
        a <- paste0(i, "A1")
        b <- paste0(i, "A2")
        # c <- paste0(i, "B")
        k <- integrate(a, b, c)
        seurat.data <- k[[1]]
        features <- k[[2]]
    }
    if (i == "P0"){
        a <- paste0(i, "A1")
        b <- paste0(i, "A2")
        # c <- paste0(i, "B")
        k <- integrate(a, b, c)
        seurat.data2 <- k[[1]]
        features2 <- k[[2]]
    }
    if (i != "E175" & i != "P0") {
        a <- paste0(i, "A")
        b <- paste0(i, "B")
        k <- integrate(a, b)
        seurat.data2 <- k[[1]]
        features2 <- k[[2]]
    }
    if (i != "E175"){
        seurat.data <- merge(seurat.data, seurat.data2)
        features <- c(features, features2)
    }
}
features <- unique(features)
seurat.data <- ScaleData(seurat.data, verbose = FALSE)
seurat.data <- RunPCA(seurat.data, npcs = 30, verbose = FALSE,  features = features)
seurat.data <- RunUMAP(seurat.data, reduction = "pca", dims = 1:30)
seurat.data <- FindNeighbors(seurat.data, reduction = "pca", dims = 1:30)
seurat.data <- FindClusters(seurat.data, resolution = 0.5)
# p1 <- DimPlot(seurat.data, reduction = "umap", group.by = "orig.ident")
# p2 <- DimPlot(seurat.data, reduction = "umap", group.by = "seurat_clusters")
# p1 + p2
saveRDS(seurat.data, "/public/home/zhangzb/test/work/seurat/all_nor_integrated.rds")
# data <- t(as.data.frame(seurat.data@assays$integrated@data))
# write.csv(data, "/public/home/zhangzb/test/work/pyscenic/all_nor_integrate/all_nor_all.csv")

