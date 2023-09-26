# /public/home/dengchx/anaconda3/envs/RR/bin/R
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)

#10X P0的数据
pdf_path <- "/public/home/zhangzb/test/work/cellchat/"
data <- readRDS("/public/home/zhangzb/test/work/seurat/all_SCT_integrate.rds")
data <- subset(data, subset = region != "Noregion")
sub_data <- subset(data, subset = development == "P0")
count <- sub_data@assays$RNA@data
count <- CreateSeuratObject(count)
count <- NormalizeData(count)
count <- count@assays$RNA@data
meta <- sub_data@meta.data
count <- count[,row.names(meta)]
#stereo的数据
pdf_path <- "/public/home/zhangzb/test/work/stereo/bin110/cellchat/"
data <- readRDS("/public/home/zhangzb/test/work/stereo/bin110/SS200000977TL_A1_bin110_sub.rds")
data <- subset(data, subset = region != "Noregion")
meta <- data@meta.data
count <- data@assays$Spatial@data
count <- CreateSeuratObject(count)
count <- NormalizeData(count)
count <- count@assays$RNA@data
count <- count[,row.names(meta)]

cellchat <- createCellChat(object = count, meta = meta, group.by = "region")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 1)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)



sample <- "stereo"
if(sample == "10X"){
  cellchat <- readRDS("/public/home/zhangzb/test/work/cellchat/10X-P0_cellchat.rds")
}else {
   cellchat <- readRDS("/public/home/zhangzb/test/work/cellchat/stereo-110_cellchat.rds")
}
pdf_path <- paste0("/public/home/zhangzb/test/work/cellchat/",
                  sample,"/")
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste0(pdf_path, "脑区之间交互数量和强度.pdf"))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights")
dev.off()
pdf(paste0(pdf_path, "每个脑区之间交互强度.pdf"))
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pathway <- cellchat@netP$pathways  #获得所有信号通路
pdf(paste0(pdf_path, "test3.pdf"))
pathways.show <- c("MK")
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
dev.off()

tf_import <- read.csv("/public/home/zhangzb/test/work/pyscenic/P0_integrate/P0_all_adj.csv",
                      header = T) 
chat_10x <- readRDS("/public/home/zhangzb/test/work/cellchat/10X-P0_cellchat.rds")
chat_stereo <- readRDS("/public/home/zhangzb/test/work/cellchat/stereo-110_cellchat.rds")
pathway_10x <- chat_10x@netP$pathways
pathway_stereo <- chat_stereo@netP$pathways
dif_pathway <- setdiff(pathway_10x, pathway_stereo)
nr4a2_tf <- tf_import[which(tf_import$TF == "Nr4a2"),]
row.names(nr4a2_tf) <- nr4a2_tf$target
# nr4a2_tf <- nr4a2_tf[which(nr4a2_tf$importance >= 0.1),]
d <- list()
for (i in dif_pathway){
  lr_10x <- extractEnrichedLR(chat_10x, signaling = i, geneLR.return = TRUE)[2]
  lr_10x <- intersect(lr_10x$geneLR, row.names(nr4a2_tf))
  if(!is.na(lr_10x[1])){
    o <- nr4a2_tf[lr_10x,]
    d[[i]] <- o
  }
}
a <- CellChatDB$interaction
b <- a[a$pathway_name=="EGF",]

pdf_path <- "/public/home/zhangzb/test/work/cellchat/10X/"
pathway <- names(d)[9]
pdf(paste0(pdf_path,pathway, "-LR-贡献图.pdf"))
netAnalysis_contribution(chat_10x, signaling = pathway)
dev.off()

#信号通路热图
pdf_path <- "/public/home/zhangzb/test/work/cellchat/10X/"
pathway <-"EGF"
pdf(paste0(pdf_path,pathway, "-脑区热图.pdf"))
netVisual_heatmap(chat_10x, signaling = pathway, color.heatmap = "Reds")
dev.off()
#对给定的信号通路和相关信号基因
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = TRUE)

plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE)








options(browser = 'firefox')
library(tidyverse)
library(networkD3)
MisLinks = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/Network/link.csv",sep = ",")
MisNodes = read.delim("https://www.bioladder.cn/shiny/zyp/bioladder2/demoData/Network/node.csv",sep = ",")

MisLinks <- read.csv("/public/home/zhangzb/test/work/cellchat/10X/10X-Nr4a2-TF.csv",
                     header = T, row.names = 1)
MisNodes <- read.csv("/public/home/zhangzb/test/work/cellchat/10X/10X-Nr4a2-TF-node.csv",
                      header = T, sep = "\t")
# MisNodes <- MisNodes[1,] 
Node2index = list()
Node2index[MisNodes$name] = 0:length(MisNodes$name)

MisLinks = MisLinks %>%
  mutate(source2 = unlist(Node2index[source])) %>%
  mutate(target2 = unlist(Node2index[target]))

# 定义颜色
group2project = paste(unique(MisNodes$group),collapse = '","')
color2project = paste(unique(MisNodes$color),collapse = '","')
my_color <- paste0('d3.scaleOrdinal().domain(["',group2project,'"]).range(["',color2project,'"])')
# MisLinks$source <- 1
# MisLinks$target2 <- 2
# MisNodes$TF <- "Nr4a2"
my_color <- "d3.scaleOrdinal().domain([\"TF\",\"target\"]).range([\"green\",\"black\"])"
a <- forceNetwork(Links = MisLinks, 
             Nodes = MisNodes,
             Source = "source", 
             Target = "target2",
             Value ="value",
             NodeID = "name",
             Group = "group",
             opacity= 0.9,        # 透明度
             Nodesize="size",
             zoom = TRUE,       # 是否可以缩放
             opacityNoHover=1,  # 鼠标没有悬浮在节点上时，文字的透明度(0-1)
             colourScale = JS(my_color),   # 节点颜色，JavaScript
             radiusCalculation = JS("d.nodesize"),
             legend=T, 
             fontSize = 3,
             linkColour= MisLinks$color,
             )
library(magrittr)
a %>% saveNetwork(file = "/public/home/zhangzb/test/work/cellchat/10X/10X-配受体基因.html")

library(igraph)
library(tidyverse)
library(ggplot2)

link <- read.csv("/public/home/zhangzb/test/work/cellchat/10X/10X-Nr4a2-TF.csv",
                     header = T, row.names = 1)
node <- read.csv("/public/home/zhangzb/test/work/cellchat/10X/10X-Nr4a2-TF-node.csv",
                      header = T, sep = "\t")
link <- link[1:31,]
node$size <- if_else(node$group == "TF", 16, 11)
node$color <- if_else(node$group == "TF", "#80b1d3", "#fb8072")

g <- graph_from_data_frame(
          d = link, vertices=node
)
# E(g)$width <- link$value
pdf("/public/home/zhangzb/test/work/cellchat/10X/10X-配受体基因.pdf")
plot(g)
dev.off()

