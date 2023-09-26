#%%
from operator import index
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
from enum import unique
import numpy as np
import math

idx_full = {
    "E135A": "V10M17-100-E135A",
    "E135B": "V10M17-085-E135B",
    "E155A": "V10M17-100-E155A",
    "E155B": "V10M17-085-E155B",
    "E165A": "V10M17-100-E165A",
    "E165B": "V10M17-085-E165B",
    "E175A1": "V10M17-101-E175A1",
    "E175A2": "V10M17-101-E175A2",
    "E175B": "V10M17-085-E175B",
    "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
test = "E175A1"

cluster = pd.read_csv("/public/home/zhangzb/test/work/tem-duc/test-seurat.csv",
                     header=0, index_col=0)
cluster.index = [i.replace("-", ".") for i in cluster.index]
cluster.index = [test + "_" + i for i in cluster.index ]
if test == "E135B" or test == "E135A":
        flie = pd.read_csv(f"/public/home/zhangzb/test/work/raw_data/raw_count/{test}_rawcount.csv",header=0,index_col=0).T
        # flie.index = [i.split("_")[1] for i in flie.index]
        flie.index = [i.replace("-", ".") for i in flie.index]
        cluster = cluster.T
        cluster = cluster[[i for i in flie.index]]
        cluster = cluster.T
        if type(cluster.iloc[0, 1]) == str:
                for i in cluster.index:
                        a = cluster.loc[i, "seurat_clusters"]
                        if "_" in a:
                                b = a.replace("_", ".")
                                cluster.loc[i, "seurat_clusters"] = b 
                cluster["seurat_clusters"] = pd.to_numeric(cluster["seurat_clusters"])   
        sort = sorted(list(set(cluster.iloc[:,1].tolist())))
        t_dic = {}
        for i in range(0, len(sort)):
                t_dic[sort[i]] = i
        for j in t_dic:
                cluster.replace(j, t_dic[j], inplace=True)



position = pd.read_csv(f"/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/{test}-coor.csv",
                     header=0, index_col=0).T
# position.index = [i.split("_")[1] for i in position.index]
position.columns = [i.replace("-", ".") for i in position.columns]
position = position[[i for i in cluster.index]].T
HE_path = f"/public/home/hongyz/Data/SpatialTranscriptomics/20210225/HE/{idx_full[test]}.tif"
HE_image = Image.open(HE_path)
fig, ax = plt.subplots(figsize=(13, 10), dpi=400) 
ax.set_xticks([]) 
ax.set_yticks([])
ax.imshow(HE_image)
position = position.reindex(cluster.index)
#%%同时画多个类
sub_num = [6]
sub_cluster = pd.DataFrame()
for i in sub_num:
        a = cluster[cluster["seurat_clusters"] == i]
        sub_cluster = pd.concat([sub_cluster, a], axis=0)
position_T = position.T
sub_position = position_T[[i for i in sub_cluster.index ]].T
sub_position = sub_position.reindex(sub_cluster.index)
scatter = ax.scatter(
                        sub_position["X"],
                        sub_position["Y"],
                        c = sub_cluster["seurat_clusters"],
                        cmap='tab20',
                        alpha=0.5
    )
clusters = np.max(np.array(cluster["seurat_clusters"]))
legend1 = ax.legend(*scatter.legend_elements(num= clusters + 1),
                    loc="best")
ax.add_artist(legend1)
plt.show()

# %% 绘制sub——cluster
num_cluster = cluster.iloc[:,1].tolist()
rows = math.ceil((len(set(num_cluster)) + 1)/5)

fig, axes = plt.subplots(rows, 5, figsize = (50, rows * 10))
[ax.axis("off") for ax in axes.flatten()]
for n, ax in zip(set(num_cluster), axes.flatten()):
        draw_series = cluster[cluster["seurat_clusters"] == n]
        draw_coor = position.reindex(index = draw_series.index)
        [axis.set_visible(False) for axis in ax.spines.values()]
        ax.imshow(HE_image)
        ax.scatter(draw_coor["X"], draw_coor["Y"],
                s = 16)
        ax.set_title(f"{test} cluster {n}"
        )
plt.savefig("/public/home/zhangzb/test/work/tem-duc/E175A1-cluster.pdf")
plt.show()


# %%绘制华大数据聚类结果
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
from enum import unique
import numpy as np
bin = "bin110"
cluster = pd.read_csv(f"/public/home/zhangzb/test/work/stereo/{bin}/position_cluster.csv",
                      header=0, index_col=0)
# cluster = cluster[cluster["seurat_clusters"] == 5]
fig, ax = plt.subplots(figsize=(10, 14))
ax.set_xticks([]) 
ax.set_yticks([])
sc = ax.scatter(cluster["x"],
             cluster["y"],
             c = cluster["seurat_clusters"],
             s = 5,
             alpha=1,
             cmap='tab20')
clusters = np.max(np.array(cluster["seurat_clusters"]))
legend1 = ax.legend(*sc.legend_elements(num= clusters + 1),
                    loc="best")
ax.add_artist(legend1)
# fig.savefig(f"/public/home/zhangzb/test/work/stereo/{bin}/cluster.pdf")
plt.show()

a = [14]
b = pd.DataFrame()
c = cluster.T
for i in a:
        b_1 = cluster[cluster["seurat_clusters"] == i]
        b = pd.concat([b, b_1], axis = 0)
b["test"] = 0
c = c[[i for i in c.columns if i not in b.index]].T
c["test"] = 1
test = pd.concat([b, c], axis=0)
color = ["red", "blue"]
fig, ax = plt.subplots(figsize=(10, 14))
ax.set_xticks([]) 
ax.set_yticks([])
for i in [0, 1]:
       ax.scatter(
        test.loc[test["test"] == i, "x"],
        test.loc[test["test"] == i, "y"],
        c = color[i]
       ) 
plt.show()
deg = pd.read_csv("/public/home/zhangzb/test/work/stereo/bin110/DEG.csv",
                  header = 0, index_col=0)
sub_deg = deg[deg["cluster"] == a[0]]
sub_deg = sub_deg[(sub_deg["avg_log2FC"] >= 1) | (sub_deg["avg_log2FC"] <= -1)]
sub_deg = sub_deg[sub_deg["p_val_adj"] <= 0.05]
sub_deg = sub_deg.sort_values(by = ["avg_log2FC"], ascending=False)

#%%
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
from enum import unique
import numpy as np
bin = "bin110"
cluster = pd.read_csv(f"/public/home/zhangzb/test/work/stereo/{bin}/position_cluster.csv",
                      header=0, index_col=0)
region = list(set(cluster["region"].to_list()))
region.sort(reverse=False)
colors = ["yellow","darkorange","tan","skyblue","darkviolet",
"darkgreen","khaki","wheat","red","cyan","peru","olive",
"yellow","cyan","mediumaquamarine","skyblue","purple","fuchsia",
"indigo","khaki"]

fig, ax = plt.subplots(figsize=(10, 14))
ax.set_xticks([]) 
ax.set_yticks([])       
for i in range(len(region)):
        plt.scatter(
                cluster.loc[cluster["region"] == region[i], "x"],
                cluster.loc[cluster["region"] == region[i], "y"],
                s = 15,
                alpha=1,
                c = colors[i],
                label = region[i]
        )
font1 = {'size':8}
plt.legend(loc = 'best', prop=font1)
fig.savefig("/public/home/zhangzb/test/work/zzb/Finalize/region-stereo.pdf")
# plt.show()
# %%绘制华大基因表达图
from matplotlib import pyplot as plt
import pandas as pd
bin = "bin110"
raw_count = pd.read_csv(f"/public/home/zhangzb/test/work/stereo/{bin}/SS200000977TL_A1_{bin}_raw_count.csv",
                        header=0,
                        index_col=0).T
raw_count.index = [int(i.split("X")[1]) for i in raw_count.index]
position = pd.read_csv(f"/public/home/zhangzb/test/work/stereo/{bin}/position_cluster.csv",
                        header=0,
                        index_col=0).T
gene = ["Dlk1", "Ptpru", "Klhl1", "Vip"]
for i in gene:
        draw_count = raw_count[[i]]
        # draw_count = draw_count[draw_count[i] > 0]
        draw_position = position[[j for j in draw_count.index]].T
        draw_position.reindex(index = draw_count.index)
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.set_xticks([]) 
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False) 
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        sc = ax.scatter(draw_position["x"],
                        draw_position["y"],
                        c = draw_count[i],
                        cmap = "OrRd",
                        s =1,
                        vmin=0)
        ax.set_title(f"{bin}_{i}")
        fig.colorbar(sc)
        fig.savefig(f"/public/home/zhangzb/test/work/stereo/{bin}/draw_gene/{i}.pdf")
wt_count = pd.read_csv("/public/home/zhangzb/test/work/raw_data/raw_count/P0A2_rawcount.csv",
                       header =0 ,index_col=0).T
wt_coor = pd.read_csv("/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/P0A2-coor.csv",
                               header = 0, index_col=0).T
for i in gene:
        draw_count = wt_count[[i]]
        # draw_count = draw_count[draw_count[i] > 0]
        draw_position = wt_coor[[j for j in draw_count.index]].T
        draw_position.reindex(index = draw_count.index)
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.set_xticks([]) 
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False) 
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        sc = ax.scatter(draw_position["X"],
                        draw_position["Y"],
                        c = draw_count[i],
                        cmap = "OrRd",
                        s =1,
                        vmin=0)
        ax.set_title(f"WT_{i}")
        fig.colorbar(sc)
        fig.savefig(f"/public/home/zhangzb/test/work/stereo/{bin}/draw_gene/WT-{i}.pdf")
# %%绘制10X与stereo 共同聚类的结果
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
type = "10X"
cluster = pd.read_csv("/public/home/zhangzb/test/work/stereo/test_cluster.csv",
                      header = 0, index_col=0)
fig, ax = plt.subplots(figsize=(5, 5))
ax.set_xticks([]) 
ax.set_yticks([])
if type == "stereo":
        position = pd.read_csv("/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/P0A2-coor.csv",
                               header = 0, index_col=0)
        position.index = [i.split("_")[1] for i in position.index]
        cluster2 = cluster[cluster["orig.ident"] == "10X"]
        position.reindex(index=cluster2.index)
        sc = ax.scatter(position["X"],
                        position["Y"],
                        c = cluster2["seurat_clusters"],
                        s =3,
                        cmap = "tab20")
        clusters = np.max(np.array(cluster2["seurat_clusters"]))
        legend1 = ax.legend(*sc.legend_elements(num= clusters + 1),
                        loc="best")
        ax.add_artist(legend1)
        plt.show()
else:
        cluster2 = cluster[cluster["orig.ident"] == "stereo"]
        cluster2.index = [int(i.split("X")[1]) for i in cluster2.index]
        position = pd.read_csv("/public/home/zhangzb/test/work/stereo/bin160/position_cluster.csv",
                               header = 0, index_col=0)
        position.reindex(index=cluster2.index)
        sc = ax.scatter(position["x"],
                        position["y"],
                        c = cluster2["seurat_clusters"],
                        s=3,
                        cmap = "tab20")
        clusters = np.max(np.array(cluster2["seurat_clusters"]))
        legend1 = ax.legend(*sc.legend_elements(num= clusters + 1),
                        loc="best")
        ax.add_artist(legend1)
        plt.show()      
# %% 绘制细胞周期评分
from operator import index
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
from enum import unique
import numpy as np
import math
from matplotlib.backends.backend_pdf import PdfPages
idx_full = {
    "E135A": "V10M17-100-E135A",
    "E135B": "V10M17-085-E135B",
    "E155A": "V10M17-100-E155A",
    "E155B": "V10M17-085-E155B",
    "E165A": "V10M17-100-E165A",
    "E165B": "V10M17-085-E165B",
    "E175A1": "V10M17-101-E175A1",
    "E175A2": "V10M17-101-E175A2",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
colors = ["red", "blue", "grey"]
ti = ["S", "G2M", "None"]

for i in idx_full:
        cycle = pd.read_csv(f"/public/home/zhangzb/test/work/seurat/cell_cycle/{i}_cellcycle.csv",
                        header = 0, index_col=0)
        cycle = cycle.replace("G1", "None")
        position = pd.read_csv(f"/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/{i}-coor.csv",
                        header = 0, index_col=0)
        position = position.reindex(cycle.index)
        HE_path = f"/public/home/hongyz/Data/SpatialTranscriptomics/20210225/HE/{idx_full[i]}.tif"
        HE_image = Image.open(HE_path)
        with PdfPages(f'/public/home/zhangzb/test/work/seurat/cell_cycle/{i}.pdf') as pdf:
                fig, ax = plt.subplots(figsize=(13, 10), dpi=400) 
                ax.set_xticks([]) 
                ax.set_yticks([])
                ax.imshow(HE_image)
                for j in range(len(ti)):
                        ax.scatter(
                                position.loc[cycle[cycle.Phase == ti[j]].index.tolist(), "X"],
                                position.loc[cycle[cycle.Phase == ti[j]].index.tolist(), "Y"],
                                s = 16,
                                c = colors[j],
                                label = ti[j],
                                alpha=0.8
                        )
                plt.legend(loc = 'best')
                plt.title(f"{i}_ cellcycle phase")
                pdf.savefig()
                plt.close()
                fig, ax = plt.subplots(figsize=(13, 10), dpi=400) 
                ax.set_xticks([]) 
                ax.set_yticks([])
                ax.imshow(HE_image)
                sc = ax.scatter(
                        position["X"],
                        position["Y"],
                        c = cycle["nFeature_RNA"],
                        cmap = "viridis_r"
                )
                fig.colorbar(sc)
                plt.title(f"{i}_Feature")
                pdf.savefig()
                plt.close()
                fig, ax = plt.subplots(figsize=(13, 10), dpi=400) 
                ax.set_xticks([]) 
                ax.set_yticks([])
                ax.imshow(HE_image)
                sc2 = ax.scatter(
                        position["X"],
                        position["Y"],
                        c = cycle["nCount_RNA"],
                        cmap = "OrRd"
                )
                fig.colorbar(sc2)
                plt.title(f"{i}_Count")
                pdf.savefig()
                plt.close()
# %%
from matplotlib import pyplot as plt
import pandas as pd
inf = pd.read_csv("/public/home/zhangzb/test/work/seurat/cell_cycle/all_meta.csv",
                  header = 0, index_col=0)
umap = pd.read_csv("/public/home/zhangzb/test/work/seurat/cell_cycle/seurat.umap.csv",
                   header = 0, index_col=0)
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_xlabel("UMAP_1")
ax.set_ylabel("UMAP_2")
sc = ax.scatter(
        umap["UMAP_1"],
        umap["UMAP_2"],
        s=5,
        c = inf["nCount_SCT"],
        cmap = "OrRd"
)
fig.colorbar(sc)
# plt.title("ngene_count")
plt.show()

#%% 绘制单基因定义脑区
from operator import index
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
from enum import unique
import numpy as np
import math

idx_full = {
    "E135A": "V10M17-100-E135A",
    "E135B": "V10M17-085-E135B",
    "E155A": "V10M17-100-E155A",
    "E155B": "V10M17-085-E155B",
    "E165A": "V10M17-100-E165A",
    "E165B": "V10M17-085-E165B",
    "E175A1": "V10M17-101-E175A1",
    "E175A2": "V10M17-101-E175A2",
    "E175B": "V10M17-085-E175B",
    "P0B": "V10M17-100-P0B",
    "P0A1": "V10M17-101-P0A1",
    "P0A2": "V10M17-101-P0A2",
}
count = pd.read_csv("/public/home/zhangzb/test/work/raw_data/raw_count/all_rawcount.csv",
                    header = 0, index_col=0).T

sample = "E165A"
HE_path = f"/public/home/hongyz/Data/SpatialTranscriptomics/20210225/HE/{idx_full[sample]}.tif"
coor = pd.read_csv(f"/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/{sample}-coor.csv",
                   header = 0, index_col=0).T
gene = "Mpped1"
sub_count = count[[gene]]
sub_count = sub_count.T
sub_count = sub_count[[i for i in sub_count.columns if sample in i]].T
precent = sub_count.quantile(0.90)[0]
sub_count = sub_count[sub_count[gene] >= precent]
sub_coor = coor [[i for i in coor.columns if i in sub_count.index]].T

HE_image = Image.open(HE_path)
fig, ax = plt.subplots(figsize=(13, 10), dpi=400) 
ax.set_xticks([]) 
ax.set_yticks([])
ax.imshow(HE_image)
sc = ax.scatter(
        sub_coor["X"],
        sub_coor["Y"],
        s = 16
)
ax.set_title(f"{sample}_{gene}_95%")
plt.show()
# %% 绘制seurat E175A2手动注释的结果
from operator import index
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
from enum import unique
import numpy as np
import math
import matplotlib as mpl
cluster = pd.read_csv("/public/home/zhangzb/test/work/tem-duc/E175A1_manu.csv",
                      header = 0, index_col= 0)
cluster.index = ["E175A1_" + i for i in cluster.index]
coor = pd.read_csv("/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/E175A1-coor.csv",
                   header = 0, index_col= 0)
manual_cluster	= list(set(cluster["manual_cluster"].tolist()))
manual_cluster = sorted(manual_cluster)
precision_cluster = list(set(cluster["precision_cluster"].tolist()))
precision_cluster = sorted(precision_cluster)
colors_list =  [mpl.colors.rgb2hex(i) for i in plt.cm.tab20.colors]

HE_path = "/public/home/hongyz/Data/SpatialTranscriptomics/20210225/HE/V10M17-101-E175A1.tif"
HE_image = Image.open(HE_path)
fig, ax = plt.subplots(figsize=(13, 10)) 
ax.set_xticks([]) 
ax.set_yticks([])
ax.imshow(HE_image)
for i in range(len(manual_cluster)):
        ax.scatter(
                coor.loc[cluster[cluster.manual_cluster == manual_cluster[i]].index.tolist(), "X"],
                coor.loc[cluster[cluster.manual_cluster == manual_cluster[i]].index.tolist(), "Y"],
                c = colors_list[i],
                s = 16,
                label = manual_cluster[i]
        )
plt.legend(loc = 'best')
plt.savefig("/public/home/zhangzb/test/work/tem-duc/E175A1_粗略注释.pdf")

HE_path = "/public/home/hongyz/Data/SpatialTranscriptomics/20210225/HE/V10M17-101-E175A1.tif"
HE_image = Image.open(HE_path)
fig, ax = plt.subplots(figsize=(13, 10)) 
ax.set_xticks([]) 
ax.set_yticks([])
ax.imshow(HE_image)
for i in range(len(precision_cluster)):
        ax.scatter(
                coor.loc[cluster[cluster.precision_cluster == precision_cluster[i]].index.tolist(), "X"],
                coor.loc[cluster[cluster.precision_cluster == precision_cluster[i]].index.tolist(), "Y"],
                c = colors_list[i],
                s = 16,
                label = precision_cluster[i]
        )
plt.legend(loc = 'best')
plt.savefig("/public/home/zhangzb/test/work/tem-duc/E175A1_精细注释.pdf")
# %%
