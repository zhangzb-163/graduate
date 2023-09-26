from operator import index
from matplotlib import pyplot as plt
from PIL import Image
import pandas as pd
from enum import unique
import numpy as np
import math
import matplotlib as mpl
#%%
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

all_meta = pd.read_csv("/public/home/zhangzb/test/work/zzb/all_meta.csv",
                       header = 0 ,index_col =0)
all_region = list(set(all_meta["region"].to_list()))
all_region = sorted(all_region)
colors_list = [mpl.colors.rgb2hex(i) for i in plt.cm.tab20.colors]
colors_list = ['#1f77b4','#aec7e8','#ff7f0e','#ffbb78','#2ca02c','#98df8a','#d62728','#ff9896',
               '#9467bd','#c5b0d5','#8c564b','#c49c94','#e377c2','#f7b6d2','#bcbd22']
region_color = dict(zip(all_region, colors_list))

for i in idx_full:
    draw_sample = all_meta[all_meta["orig.ident"]==i]
    coor = pd.read_csv(f"/public/home/hongyz/workspace/mouse-brain-full/Data/coor_df/{i}-coor.csv",
                   header = 0, index_col=0).T
    coor = coor[[j for j in draw_sample.index]].T
    HE_path = f"/public/home/hongyz/Data/SpatialTranscriptomics/20210225/HE/{idx_full[i]}.tif"
    HE_image = Image.open(HE_path)
    fig, ax = plt.subplots(figsize=(13, 10)) 
    ax.set_xticks([]) 
    ax.set_yticks([])
    ax.imshow(HE_image)
    for d in region_color:
        if d in draw_sample["region"].to_list():
            ax.scatter(
                    coor.loc[draw_sample[draw_sample.region == d].index.tolist(), "X"],
                    coor.loc[draw_sample[draw_sample.region == d].index.tolist(), "Y"],
                    c = region_color[d],
                    s = 16,
                    label = d
            )
    plt.legend(loc = 'best')
    plt.savefig(f"/public/home/zhangzb/test/work/zzb/Finalize/all_region/{i}_region.pdf")
    plt.close()