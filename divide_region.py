#%%  拆分spot 手动将每个spot注释为region+时期
from enum import unique
import json
import pandas as pd
import numpy as np
import scanpy as sc
import loompy as lp

#从loom提取regulor矩阵
f_pyscenic_output = "/public/home/zhangzb/test/work/pyscenic/combat_all_sum_integrate/combat_all_sum_all_pyscenic_output.loom"
with lp.connect( f_pyscenic_output, mode='r+', validate=False ) as lf:
    auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID).T
    # regulons = lf.ra.Regulons
    all_count = auc_mtx
    all_count = all_count[[i for i in all_count.columns if "Meninges" not in i]]
##

with open (r"/public/home/zhangzb/test/work/raw_data/seurat-cluster/regions.json") as f: #clusters 所在的区域
    load_dict = json.load(f)
regions = load_dict["regions"]
keys = list(regions.keys())
# values = list(regions.values())
# values = sum(values,[])
cluster = pd.read_csv("/public/home/zhangzb/test/work/raw_data/seurat-cluster/seurat_all.csv",
                      header= 0,index_col= 0)
# cluster.index = [i.replace("-", ".") for i in cluster.index]
# cluster = cluster.T
# cluster = cluster[[i for i in cluster.columns if "P0B" not in i]].T
all_count = pd.read_csv("/public/home/zhangzb/test/work/combat/all_combat.csv", 
                        header=0,
                        index_col=0) #除P0B以外所有数据combat之后的结果
# all_count.columns = [i.replace("-", ".") for i in all_count.columns]
# all_count = all_count[[i for i in all_count.columns if "P0B" not in i]]

#在spot前面加上区域的名称
def function1(value, dictionay):
    """
    对于字典,输入values值返回对应的key值
    """
    a = list(dictionay.keys())
    b = sum(list(dictionay.values()), [])
    for i in a:
        c = list(dictionay[i])
        if value not in b:
            return "Nregions"
            break
        if value in c:
            return i
            break

for i in all_count.columns:
    a = cluster.loc[i, "seurat_clusters"]
    # a = i.split("_")[0] +"_"+str(a)
    b = function1(a, regions)
    all_count.rename(columns={i: b +"-"+ i }, inplace=True)
all_count.to_csv("/public/home/zhangzb/test/work/combat/all_combat.csv")

# 求regulor按照region-time 的平均值或和
# re_time = all_count.columns.to_list()
# re_time = set([i.split("_")[0] for i in re_time])
# re_time_frame = pd.DataFrame()
# for i in re_time:
#     a = all_count[[j for j in all_count.columns if i in j]]
#     b = pd.DataFrame(a.sum(axis=1), columns = [i])
#     re_time_frame = pd.concat([re_time_frame, b], axis=1)

#将各个区域按照时间点进行求和或求平均值
all_time = list(set([i.split("_")[0].split("-")[1] for i in all_count.columns]))
all_time.sort()
all_regions = list(set([i.split("_")[0].split("-")[0] for i in all_count.columns]))
all_regions.append("Nregions")

re_time = []
for i in regions:
    for j in regions[i]:
        j = j.split("_")[0]
        re_time.append(i + "-" + j)
re_time = list(set(re_time))
for i in re_time:
    if "P0B" in i:
        re_time.remove(i)

sum_regions = pd.DataFrame()
mean_regions = pd.DataFrame()
for i in all_regions:
    for j in all_time:
        if i + "-" + j in re_time:
            a = all_count[[idx for idx in all_count.columns if i + "-" + j  in idx]]
            b = pd.DataFrame(a.sum(axis=1), columns= [i + "-" + j])
            c = pd.DataFrame(a.mean(axis=1), columns= [i + "-" + j])
            sum_regions = pd.concat([sum_regions, b], axis=1)
            mean_regions = pd.concat([mean_regions, c], axis=1)

sum_regions.to_csv("/public/home/zhangzb/test/work/pyscenic/combat_all_sum_integrate/combat_sum.csv")
mean_regions.to_csv("/public/home/zhangzb/test/work/pyscenic/combat_all_sum_integrate/combat_mean.csv")

#run /public/home/zhangzb/test/work/src/R/run_all_scenic.R 之后接着往下

#%% 将得到的AUC评分矩阵，进行PCA
import matplotlib.pyplot as plt
from sklearn import decomposition
sum = all_count.T
p = "-"
method = "test"
# sum = pd.read_csv(f"/public/home/zhangzb/test/work/tem-duc/seurat_norinte_sum/int/sum_norintegrate.csv", 
#                   header=0, index_col=0).T
# sum = sum[[i for i in sum.columns if "E175B" not in i]]
# sum = sum.T
pca = decomposition.PCA(n_components=2)
pca_sum = pca.fit(sum)
de_sum = pca_sum.transform(sum)
label = sum.index.tolist()
labels = [i.split(p)[1] for i in label]
position_sum = pd.DataFrame(de_sum, columns=["X","Y"])

re_name = []
for i in label:
    if "3V" in i or "GE" in i:
        a = i.split(p)[0][0:2]
    else:    
        a = i.split(p)[0][0:3]
        # a = a.capitalize()
    # if "SVZ" in i or "CXT" in i:
    #     a = i.split(p)[0][0:3]
    b = i.split(p)[1]
    c = a + p + b
    re_name.append(c)
color = ['darkred', 'red', 'chocolate', 'sandybrown', 'navy', 'darkblue', 'limegreen',
    'lime', 'purple', 'darkmagenta', '#6495ED']
position_sum["labels"] = labels
period = list(set(labels))
period.sort()
fig, ax = plt.subplots(figsize=(12,12))
for i in range(len(period)):
    print(color[i])
    plt.scatter(position_sum.loc[period[i] == position_sum["labels"], "X"],#传入数据x
                position_sum.loc[period[i] == position_sum["labels"], "Y"],#传入数据y
                s = 80,#散点图形（marker）的大小
                c = color[i],#marker颜色
                marker = '.',#marker形状
                #marker=matplotlib.markers.MarkerStyle(marker = markers[i],fillstyle='full'),#设置marker的填充
                alpha=1,#marker透明度，范围为0-1
                # facecolors='r',#marker的填充颜色，当上面c参数设置了颜色，优先c
                edgecolors='none',#marker的边缘线色
                linewidths=1,#marker边缘线宽度，edgecolors不设置时，该参数不起作用
                label = period[i]
                )#后面图例的名称取自label
for i in range(len(label)):
    x = position_sum["X"][i]
    y = position_sum["Y"][i]
    plt.text(x, y, re_name[i], fontsize = 8)
plt.xlabel("pc1")
plt.ylabel("pc2")
font1 = {'size':5}
plt.legend(loc = 'best', prop=font1)
# plt.show()
plt.savefig(f'/public/home/zhangzb/test/work/zzb/PCA_发育趋势2.pdf')
            #   dpi=900)
plt.close()
#%%
#z-score normalization
# auc_mtx = auc_mtx[[i for i in auc_mtx.columns if "Noregion" not in i]]
# auc_mtx = auc_mtx[[i for i in auc_mtx.columns if "meninges" not in i]]

# re_col = []
# for i in regulonAUC.columns:
#     if "3V" not in i and "GE" not in  i and "SVZ" not in i and "CXT" not in i:
#         b = i.split("-")[0]
#         c = i.split("-")[1]
#         b = b.capitalize()
#         a = b+"-"+c
#     else:
#         a = i
#     re_col.append(a)
# regulonAUC .columns = re_col
regulonAUC = auc_mtx.loc[(auc_mtx != 0).any(axis = 1)]
nor_regulonAUC = (regulonAUC - np.mean(regulonAUC)) / np.std(regulonAUC)
# nor_regulonAUC = nor_regulonAUC.T
full_regions = set(nor_regulonAUC.columns.tolist())
  
cor_dict = {}
for i in full_regions:
    count_mean_1 = nor_regulonAUC[[i]].T.mean()
    a = []
    for j in full_regions:
        count_mean_2 = nor_regulonAUC[[j]].T.mean()
        p = np.corrcoef(count_mean_1,count_mean_2)[0,1]
        a.append(p)
    cor_dict[i] = a
cor_metric = pd.DataFrame(cor_dict)
cor_metric.index = cor_metric.columns
cor_metric.to_csv("/public/home/zhangzb/test/work/zzb/cor_metric.csv")


# %%
