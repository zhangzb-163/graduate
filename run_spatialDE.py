import pandas as pd
import SpatialDE
import NaiveDE

def SVG(sample1, sample2):
    data = pd.DataFrame()
    for s in [sample1, sample2]: 
        drct = {}
        count = pd.read_csv(f"/public/home/zhangzb/test/work/raw_data/raw_count/{s}_rawcount.csv",
                            header=0, index_col=0).T
        position = pd.read_csv(f"/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/{s}-coor.csv",
                            header=0, index_col=0)
        count = count.loc[position.index]
        a = pd.DataFrame(count.sum(axis=1), columns = ["total_counts"])
        position = pd.concat([position, a], axis=1)
        norm_expr = NaiveDE.stabilize(count.T).T
        resid_expr = NaiveDE.regress_out(position, norm_expr.T, 'np.log(total_counts)').T
        # sample_resid_expr = resid_expr.sample(n=1000, axis=1, random_state=1)
        X = position[['X', 'Y']]
        results = SpatialDE.run(X, resid_expr)
        sub_results = results[['g', 'l', 'qval']]
        # sub_results = sub_results[sub_results['qval'] <= 0.01]
        sub_results.columns = [s +"_"+ l for l in sub_results.columns]
        data = pd.concat([data, sub_results], axis = 1)
    return data
sample = ["E135", "E155", "E165", "E175", "P0"]
svg_res = pd.DataFrame()
for i in sample:
    if i == "E175" or i == "P0":
        sample1 = i + "A1"
        sample2 = i + "A2"
        res = SVG(sample1, sample2)
    else:
        sample1 = i + "A"
        sample2 = i + "B"
        res = SVG(sample1, sample2)
    svg_res = pd.concat([svg_res, res], axis = 1)

svg_res.to_csv("/public/home/zhangzb/test/work/spatialDE/all_spatialDE.csv")  


#%%整合spark和spatialDE的数据
import pandas as pd
import numpy as np
from venn import venn
from matplotlib import pyplot as plt
svg_spatialDE = pd.read_csv("/public/home/zhangzb/test/work/spatialDE/all_spatialDE.csv",
                      header = 0, index_col=0)
svg_spatialDE = svg_spatialDE[[i for i in svg_spatialDE.columns if "_l" not in i]]
svg_spark = pd.read_csv("/public/home/zhangzb/test/work/spark/all_spark.csv",
                        header = 0, index_col=0)
qval = 0.05
def SVG_res(sample1, sample2):
    svg_data = pd.DataFrame()
    for i in [sample1, sample2]:
        spatialDE_res = svg_spatialDE[[j for j in svg_spatialDE.columns if i in j]]
        spatialDE_res = spatialDE_res[spatialDE_res[i + "_qval"] <= qval].sort_values(by = i + "_qval")
        spatialDE_res = spatialDE_res.iloc[:,[0]]
        spatialDE_res.index = range(0, len(spatialDE_res))
        spark_res = svg_spark[[j for j in svg_spark.columns if i in j]]
        spark_res = spark_res[spark_res[i + "_adjusted_pvalue"] <= qval].sort_values(by = i + "_adjusted_pvalue")
        spark_res = spark_res.iloc[:,[0]]
        spark_res.index = range(0, len(spark_res))
        a = pd.concat([spatialDE_res, spark_res], axis=1)
        a.columns = [i+"_spatialDE", i+"_spark"]
        svg_data = pd.concat([svg_data, a], axis=1)
    return svg_data
        

sample = ["E135", "E155", "E165", "E175", "P0"]
svg_res_data = pd.DataFrame() #筛选spark和spatialDE的结果
for i in sample:
    if i == "E175" or i == "P0":
        sample1 = i + "A1"
        sample2 = i + "A2"
        res = SVG_res(sample1, sample2)
    else:
        sample1 = i + "A"
        sample2 = i + "B"
        res = SVG_res(sample1, sample2)
    svg_res_data = pd.concat([svg_res_data, res], axis = 1)

def Intersect(t_data):
    a = t_data
    b = a.iloc[:,[0]].dropna().iloc[:,0].to_list()
    c = a.iloc[:,[1]].dropna().iloc[:,0].to_list()
    d = [k for k in b if k in c]
    return(d)

all_res = pd.DataFrame() #每个样品取spatialDE和spark交集的结果
results = pd.DataFrame() #取每个时期的SVG的并集
for i in sample:
    if i == "E175" or i == "P0":
        sample1 = i + "A1"
        sample2 = i + "A2"
        a1 = svg_res_data[[j for j in svg_res_data.columns if sample1 in j]]
        b1 = Intersect(a1)
        a2 = svg_res_data[[j for j in svg_res_data.columns if sample2 in j]]
        b2 = Intersect(a2)
        c = list(set(b1) | set(b2))
        c= pd.DataFrame(c ,columns = [i])
        results = pd.concat([results, c], axis = 1)
        b1 = pd.DataFrame(b1, columns = [sample1])
        b2 = pd.DataFrame(b2, columns = [sample2])
        all_res = pd.concat([all_res, b1, b2], axis=1)
    else:
        sample1 = i + "A"
        sample2 = i + "B"
        a1 = svg_res_data[[j for j in svg_res_data.columns if sample1 in j]]
        b1 = Intersect(a1)
        a2 = svg_res_data[[j for j in svg_res_data.columns if sample2 in j]]
        b2 = Intersect(a2)
        c = list(set(b1) | set(b2))
        c= pd.DataFrame(c ,columns = [i])
        results = pd.concat([results, c], axis = 1)
        b1 = pd.DataFrame(b1, columns = [sample1])
        b2 = pd.DataFrame(b2, columns = [sample2])
        all_res = pd.concat([all_res, b1, b2], axis=1)
results_1 = []  #所有样品加起来一个有多少个SVG  
for i in sample:
    a = results[[i]].dropna()
    for j in a.iloc[:,0]:
        if j not in results_1:
            results_1.append(j)
results_2 = pd.DataFrame() #查看每个时期独特的SVG
for i in sample:
    a = results[[i]].dropna()
    b = results[[k for k in results.columns if i not in k]]
    c = []
    for j in a.iloc[:,0]:
        if j not in b.values:
            c.append(j)
    results_2 = pd.concat([results_2, pd.DataFrame(c, columns = [i])], axis = 1)

def p_val(method, stage):
    """
    查看每个时期独特的SVG在原来的p_val值
    method:指定spark或spatialDE的结果
    stage:为样品的发育时间点
    """
    svg_d = method
    unique_svg = results_2[[i for i in results_2.columns if stage in i]].dropna()[stage].to_list()
    if stage == "E175" or stage == "P0":
        a1 = stage + "A1"
        a2 = stage +"A2"
    if stage != "E175" and stage != "P0":
        a1 = stage + "A"
        a2 = stage + "B"
    b1 = svg_d[[i for i in svg_d.columns if a1 in i]].dropna().T
    b1.columns = b1.iloc[0,:]
    b1 = b1[[i for i in unique_svg if i in b1.columns]].T
    
    b2 = svg_d[[i for i in svg_d.columns if a2 in i]].dropna().T
    b2.columns = b2.iloc[0,:]
    b2 = b2[[i for i in unique_svg if i in b2.columns]].T
    c = [i for i in b1.index if i in b2.index]
    data1 = pd.DataFrame()
    for j in c:
        if b1.loc[j, b1.columns[1]] <= qval and b2.loc[j, b2.columns[1]] <= qval:
            b1_1 = pd.DataFrame(b1.loc[j]).T 
            b2_2 = pd.DataFrame(b2.loc[j]).T
            b3 = pd.concat([b1_1, b2_2], axis = 1)
            data1 = pd.concat([data1, b3], axis = 0)
    b1.index = range(0, len(b1))
    b2.index = range(0, len(b2))
    data2 = pd.concat([b1, b2], axis = 1)
    return(data1, data2)
sample_t = "P0"
results_3 = p_val(svg_spatialDE, sample_t) #查看svg的p_val
results_4 = [] #查看在results_3中是否有spark和spatialDE的交集 
for i in results_3[0].index:
    if i in results[[sample_t]].dropna().values:
        results_4.append(i)


# %%
allen = pd.read_table("/public/home/zhangzb/test/work/raw_data/Allen-development.txt",
                      header = None, index_col = None, sep = "\t")
allen = sum(allen.values.tolist(), [])
sample = ["E135", "E155", "E165", "E175", "P0"]
results_5 = [] #查看在每个软件每个样品的并集
for i in sample:
    a = svg_spatialDE[[j for j in svg_spatialDE.columns if i in j]]
    a1 = a.iloc[:,0:2].dropna()
    a1 = a1[a1.iloc[:,1] <= qval].iloc[:,0].to_list()
    a2 = a.iloc[:,2:4].dropna()
    a2 = a2[a2.iloc[:,1] <= qval].iloc[:,0].to_list()
    a3 = list(set(a1) | set(a2))
    for h in a3:
        if h not in results_5:
            results_5.append(h)
hvg = pd.read_csv("/public/home/zhangzb/test/work/seurat/seurat_HVG.csv", 
                  header = 0, index_col = 0)
results_6 = [] #seurat HVG，先在同一时期取交集，再到所有样品取并集
for i in sample:
    a = hvg[[j for j in hvg.columns if i in j]]
    b = list(set(a.iloc[:,0].to_list()).intersection(set(a.iloc[:,1])))
    for j in b :
        if j not in results_6:
            results_6.append(j)
results_6_2 = [] #每个脑区上调的DEG合集
results_6_3 = {} #每个脑区上调的DEG和SVG交集
regions = ["Cortex", "Hippocampus", "Thalamus", "Hypothalamus", "Olfactory",
             "Striatum", "Amygdalar", "Habenular", "3V", "Choroid", "Meninges",
             "CXTsp", "GE", "SVZ"]
for i in regions:
    f = pd.read_csv(f"/public/home/zhangzb/test/work/zzb/integrated_DEG/{i}_DEG.csv",
                    header = 0, index_col = 0)
    f = f[f["p_val_adj"] <= 0.01]
    f = f[(f["avg_log2FC"] >= 2) | (f["avg_log2FC"] <= -2)] #| (f["avg_log2FC"] <= -1)
    results_6_2.append(f.index.to_list())
    print(f"{i}:{len(f)}")
    a = [j for j in f.index.to_list() if j in results_1]
    results_6_3[i] = a
results_6_2 = list(set(sum(results_6_2, [])))
#绘制韦恩图
all_data1 = {"SVG":set(results_1), "region-DEG": set(results_6_2)}
all_data2 = {"SVG": set(results_1), "Allen":set(allen)}
pdf = venn(all_data1)
plt.savefig("/public/home/zhangzb/test/work/zzb/SVG-DEG.pdf")
#%% 从表达量上来筛序时空可变基因
import scipy.stats as stats
all_count = pd.read_csv("/public/home/zhangzb/test/work/raw_data/raw_count/all_rawcount.csv",
                        header = 0, index_col = 0)
test = ["E135", "E155", "E165", "E175", "P0"]
sample1 = all_count[[j for j in all_count.columns if test[0] in j]].T
sample2 = all_count[[j for j in all_count.columns if test[1] in j]].T
sample3 = all_count[[j for j in all_count.columns if test[2] in j]].T
sample4 = all_count[[j for j in all_count.columns if test[3] in j]].T
sample5 = all_count[[j for j in all_count.columns if test[4] in j]].T
results_7 = [] #查看是否有随着发育表达量减少的gene
results_8 = [] #查看是否有随着发育表达量增加的gene

for i in results_1:
    a = sample1[[i]]
    a1 = a[a[i] > 0]
    a1 = len(a1)/len(a)
    # a1 = (a.sum(axis = 0) / len(a)).values[0]
    # a1 = np.median(a)
    b = sample2[[i]]
    b1 = b[b[i] > 0]
    b1 = len(b1)/len(b)
    # b1 = (b.sum(axis = 0) / len(b)).values[0]
    # b1 = np.median(b)
    c = sample3[[i]]
    c1 = c[c[i] > 0]
    c1 = len(c1)/len(c)
    # c1 = (c.sum(axis = 0) / len(c)).values[0]
    # c = np.median(c)
    d = sample4[[i]]
    d1 = d[d[i] >0]
    d1 = len(d1)/len(d)
    # d1 = (d.sum(axis = 0) / len(d)).values[0]
    e = sample5[[i]]
    e1= e[e[i]>0]
    e1 = len(e1)/len(e)
    # e1 = (e.sum(axis = 0) / len(e)).values[0]
    if d1 == 0 or e1 ==0:
        continue
    if (a1 - b1)/a1 >=0.15 and d1 > e1 or abs(d1-e1)/d1 <=0.02 and stats.mannwhitneyu(a[i], b[i])[1] <= 0.01:
        results_7.append(i)
    elif (b1 - c1)/b1 >= 0.15 and d1 > e1 or abs(d1-e1)/d1 <=0.02 and stats.mannwhitneyu(b[i], c[i])[1] <= 0.01:
        results_7.append(i)
           
    elif (b1 - a1)/b1 >= 0.5 and d1 < e1 or abs(e1-d1)/e1 <=0.02 and stats.mannwhitneyu(b[i], a[i])[1] <= 0.01:
        results_8.append(i)
    elif (c1 - b1)/c1 >= 0.5 and d1 < e1 or abs(e1-d1)/e1 <=0.02 and stats.mannwhitneyu(c[i], b[i])[1] <= 0.01:
        results_8.append(i)
print(len(results_7))
print(len(results_8))

for i in sample:
    svg = results[[i]].dropna().iloc[:,0].tolist()
    hvg = pd.read_csv(f"/public/home/zhangzb/test/work/zzb/integrated_DEG/{i}_DEG.csv",
                      header = 0, index_col = 0)
    hvg = hvg[(hvg.avg_log2FC >= 0.5)&(hvg.p_val_adj <= 0.01)]
    a = [i for i in hvg.index if i in svg]
    results_7.append(a)
results_7 = sum(results_7, [])
# %% 绘制表达趋势折线图
time = ["E135", "E155", "E165", "E175", "P0"]
y_all = [] #有表达的spot占的比例
gene ="Mef2c"
draw_gene = pd.DataFrame()
for i in time:
    a = all_count[[j for j in all_count.columns if i in j]].T
    b = a[[gene]]
    b.columns = [i]
    b2 = b[b[i] >0 ]
    # print(f"{i}:{len(b2)/len(b)}")
    y_all.append((len(b2)/len(b))*100)
    draw_gene = pd.concat([draw_gene, b], axis=1)
plt.plot(time, y_all)
plt.scatter(time,y_all, c=  "red")
plt.grid(linestyle="--", alpha=0.3)
plt.xlabel("Development", fontdict={"size":14})
plt.ylabel("expression spot ratio(%)", fontdict={"size":14})
for i in range(0,4):
    c = i
    d = i + 1
    test = stats.mannwhitneyu(draw_gene[[time[c]]].dropna(), draw_gene[[time[d]]].dropna())
    print(f"{time[c]}-{time[d]}:{test[1]}")
# %%
