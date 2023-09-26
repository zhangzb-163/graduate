library(monocle3)
library(stringr)
library(dplyr)
library(ggplot2)
library(data.table)
options(warn=-1)
svg <- read.table("/public/home/zhangzb/test/work/zzb/Finalize/SVG_gene.txt")
all <- as.data.frame(fread("/public/home/zhangzb/test/work/combat/all_combat.csv",
                header = T, data.table = T))
row.names(all) <- all[, 1]
all <- all[,-1]
all <- all[,-c(grep("Nregion", colnames(all)))]
all <- all[,-c(grep("Meninges", colnames(all)))]
all <- all[svg$V1,]
all <- t(all)
# SVG <- read.csv("/public/home/zhangzb/test/work/raw_data/cluster-by-hyz/top_1000_union.csv",
#                col.names =FALSE)
# all1 <- all[, SVG[,1]]
regions <- c("3V", "Hypothalamus")  #  , "Striatum" "Choroid", "SVZ", 

part <- data.frame()
for (i in regions){
part_1 <- all[c(grep(i, row.names(all))),]
part <- rbind(part, part_1)
}
part <- t(part)
a <- read.csv("/public/home/zhangzb/test/work/monocle3/test/sub_Hypothalamus.csv", 
               header = T, row.names = 1)
b <- setdiff(colnames(part), row.names(a))
part <- part[,b]
# part <- t(all)
# time <- c("E135B", "E155A", "E165A", "E175A1", "P0A1", "E175B") #"E135B", "E155B", "E165B", "E175A1", "P0A1", "E135A", "Nregions", "meninges" 
# file_colnames_1 <- unique(str_split_fixed(colnames(part), pattern = "_", 2)[,1])
# file_colnames_2 <- unique(str_split_fixed(file_colnames_1, pattern = "\\-", 2)[,2])
# file_colnames_3 <- unique(str_split_fixed(file_colnames_1, pattern = "\\-", 2)[,1])
# file_colnames <- c(file_colnames_1, file_colnames_2, file_colnames_3)
# for (i in time){
#        if(i %in% file_colnames){
#               part <- part[, -c(grep(i, colnames(part)))]
#               }
# }

#构建cell_metadata

cell_metadata <- data.frame(Type = str_split_fixed(colnames(part), pattern = "_", 2)[,1])
cell_metadata$nUMI <- apply(part, 2, sum)
cell_metadata$batch <- substr(str_split_fixed(cell_metadata$Type, pattern = '\\-',2)[,2],
                              start = 1, stop = 3)
rownames(cell_metadata) <- colnames(part)
cell_metadata.1 <- subset(cell_metadata, batch != "E17" & batch != "P0A")
cell_metadata.1$Type2 <- unlist(lapply(cell_metadata.1$Type, 
                                function(x)substr(x, start = 1, stop = nchar(x)-1)))
cell_metadata.2 <- subset(cell_metadata, batch == "E17" | batch == "P0A")
cell_metadata.2$Type2 <- unlist(lapply(cell_metadata.2$Type, 
                                   function(x)substr(x, start = 1, stop = nchar(x)-2)))
cell_metadata <- rbind(cell_metadata.1, cell_metadata.2)
cell_metadata$batch <- str_split_fixed(cell_metadata$Type2, pattern = '\\-',2)[,2]
cell_metadata <- cell_metadata[colnames(part),]
cell_metadata[which(cell_metadata$batch == "E135"), "batch"] <- "E13.5"
cell_metadata[which(cell_metadata$batch == "E155"), "batch"] <- "E15.5"
cell_metadata[which(cell_metadata$batch == "E165"), "batch"] <- "E16.5"
cell_metadata[which(cell_metadata$batch == "E175"), "batch"] <- "E17.5"
cell_metadata$Type3 <- paste0(str_split_fixed(row.names(cell_metadata), pattern = '\\-',2)[,1],
                               "-", cell_metadata$batch)
#构建gene_annotation
gene_annotation <- data.frame(gene_short_name = row.names(part))
gene_annotation$num_cells_expressed <- apply(part, 1, function(x)sum(x > 0 ))
row.names(gene_annotation) <- row.names(part)

part <- as.matrix(part)
cds <- new_cell_data_set(part,
                     cell_metadata = cell_metadata,
                     gene_metadata = gene_annotation)

cds <- preprocess_cds(cds, num_dim = 62)
# plot_pc_variance_explained(cds)
## Step 2: Remove batch effects with cell alignment
cds <- align_cds(cds, alignment_group = "batch")

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
# plot_cells(cds)

#替换为seurat的umap图
# umap.cds <- cds@int_colData$reducedDims$UMAP
# umap.seurat <- Embeddings(brain.combined, reduction = "umap")
# umap.seurat <- umap.seurat[row.names(umap.cds),]
# cds@int_colData$reducedDims$UMAP <- umap.seurat
## Step 4: Cluster the cells
cds <- cluster_cells(cds)
# plot_cells(cds, reduction_method="UMAP", color_cells_by="Type2", group_label_size= 3) 

## Step 5: Learn a graph
cds <- learn_graph(cds)
p <- plot_cells(cds, color_cells_by = "Type3", label_groups_by_cluster=FALSE,
       label_leaves=FALSE, label_branch_points=FALSE, group_label_size= 6)
ggsave(p, filename = "/public/home/zhangzb/test/work/zzb/Finalize/发育轨迹/下丘脑_发育轨迹.pdf")
## Step 6: Order cells
cds <- order_cells(cds)

# 自动选择根节点
# get_earliest_principal_node <- function(cds, time_bin="SVZ-E135B"){
#   cell_ids <- which(colData(cds)[, "Type"] == time_bin)
  
#   closest_vertex <-
#   cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <-
#   igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#   (which.max(table(closest_vertex[cell_ids,]))))]
  
#   root_pr_nodes
# }
# cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
 
p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, 
               label_leaves = FALSE,  label_branch_points = FALSE)
## Step 7: analyzing subset
cds_subset <- choose_cells(cds)
cds_sub <- choose_graph_segments(cds)
sub_col <- cds_subset@colData
write.csv(sub_col,"/public/home/zhangzb/test/work/monocle3/test/sub_Hypothalamus.csv")

#DEG
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.005))
result <- subset(ciliated_cds_pr_test_res, q_value < 0.005)
result <- result[order(result$q_value, decreasing = T),]
result <- result[order(result$morans_I),]
Track_genes_sig <- result %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
a <- c("Nkx2-1", "Zeb2", "Dlx1", "Dlx2", "Dlx5", "Dlx6")  #"Unc5b", 
Cortex_tf <- c("Cdon", "Celsr1", "Dbi", "E2f5", "Eomes",
       "Hmgn2", "Neurog2", "Notch1", "Sox3", "Ssrp1", "Tead2", "Tgif2")  #"Pcnt", 
Thalamus_tf <- c("Gbx2", "Lhx9", "Zic4", "Etv1", "Esrrg", "Nr4a2", "Sox1", "Adarb1",
                 "Prox1", "Lef1", "Sp9", "Mef2a", "Six3", "Zeb2")
Hypothalamus_tf <- c("Tead2", "Tcf12", "Tgif2", "Hbp1", "Rfx4", "Sox21", "Sox9", "Zfhx3",
                     "E2f5", "Nfia", "Sox1", "Supt20", "Luzp2", "Msi1", "Nanos1","Ctbp1",
                     "Zfp62", "Zscan18", "Zkscan3")
a <- read.table("/public/home/zhangzb/test/work/cisTarget_databases/allTFs_mm.txt")
a <- intersect(a$V1, svg$V1)
b <- Hypothalamus_tf 
b <- c("Slit1", "Slit2", "Robo1", "Robo2")
b <- c("Zfp62")
p <- plot_genes_in_pseudotime(cds[b,], color_cells_by="batch",
                              min_expr=0.3, ncol = 2)
ggsave(paste0("/public/home/zhangzb/test/work/zzb/Finalize/发育轨迹/monocle3/丘脑-Nr4a2.pdf"),
       plot = p, width = 10, height = 10)
for (i in a){
       if(!(i %in% row.names(result))){
              print(i)
       }
}
print(p)
ggsave("/public/home/zhangzb/test/work/monocle3/part_result/Genes_Jitterplot.jpg",
       plot = p, width = 10, height = 20)

p <- plot_cells(cds, genes=Track_genes_sig, show_trajectory_graph=FALSE,
                label_cell_groups=FALSE,  label_leaves=FALSE)
ggsave("/public/home/zhangzb/test/work/monocle3/part_result/Genes_Featureplot.jpg",
       plot = p, width = 20, height = 15)
write.csv(result, "/public/home/zhangzb/test/work/monocle3/DEG/Hypothalamus.csv")

