
import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
from pathlib import Path
import sys
import os
import warnings
warnings.filterwarnings("ignore")
dir_root = os.getcwd()
TILE_PATH = Path(dir_root+"/data/result")
TILE_PATH.mkdir(parents=True, exist_ok=True)

parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
#ratio= sys.argv[1]
ratio = float(parameter_settings_f['parameter'][0])

ratio = float(ratio)
subcelltype_color = {}
def bubbleSort(arr,serial):
    for i in range(1, len(arr)):
        for j in range(0, len(arr)-i):
            if arr[j] > arr[j+1]:
                serial[j],serial[j + 1] = serial[j + 1],serial[j]
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
    return arr , serial
two_time_anndata_file = dir_root+"/data/intermediate_result/tow_time_adata.h5ad"
pseudo_file = dir_root+"/data/pseudotime.csv"
autoencode_distance_comp_file = dir_root+"/data/intermediate_result/tow_time_auto_encode_distances.csv"
pic_distance_comp_file = dir_root+"/data/intermediate_result/tow_time_pic_distances.csv"
choose_hvg_distanc_comp_file = dir_root+"/data/intermediate_result/tow_time_choose_hvg_distances.csv"
TL_edges_save_file = dir_root+"/data/result/tow_time_TL_edges.csv"
samll_TL_edges_save_file = dir_root+"/data/result/tow_time_TL_edges_little_message.csv"
single_cell_adata = sc.read(two_time_anndata_file)
before_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "before"]
after_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "after"]
all_single_cell_csv = before_single_cell_csv.append(after_single_cell_csv)
before_num = len(before_single_cell_csv.index)
after_num = len(after_single_cell_csv.index)
two_time_autoencode_csv = pd.read_csv(autoencode_distance_comp_file,index_col=0)
two_time_autoencode_np = two_time_autoencode_csv.values
pic_distance_csv = pd.read_csv(pic_distance_comp_file,index_col=0)
pic_distance_csv_np = pic_distance_csv.values
two_time_choose_hvg_csv = pd.read_csv(choose_hvg_distanc_comp_file,index_col=0)
two_time_choose_hvg_np  = two_time_choose_hvg_csv.values
pseudo_csv = pd.read_csv(pseudo_file,index_col=0)
all_index = pic_distance_csv.index.append(pic_distance_csv.columns)
new_pseudo_time = pd.DataFrame(index=all_single_cell_csv.index,columns=["pseudotime"])
for i in new_pseudo_time.index:
    cell_name = all_single_cell_csv.loc[i,"single_cell_name"].split(" ")[0]
    new_pseudo_time.loc[i,"pseudotime"] =  pseudo_csv.loc[cell_name,"pseudotime"]
two_distence = (1 - ratio) * two_time_autoencode_np  + ratio * pic_distance_csv_np
similar_threshold = 25
edges = []
pesudo_all = 0
distance_num = 0
break_num = 0
for i in range(after_num):
    after_time_pesuda = new_pseudo_time.loc[two_time_autoencode_csv.columns[i],"pseudotime"]
    in_sort = np.argsort(two_distence[:,i] )
    min_similar = two_distence[in_sort[0],i]
    maybe_point = in_sort[0:similar_threshold]

    cell_gene_similar_value = []
    for pic_cell in maybe_point:

        cell_gene_similar_value.append(two_time_choose_hvg_np[pic_cell,i])

    cell_gene_similar_value,maybe_point = bubbleSort(cell_gene_similar_value,maybe_point)

    j = 0

    min_pesudo = 1
    min_pesudo_loc = maybe_point[0]
    while(True):
        before_time_pesuda = new_pseudo_time.loc[two_time_autoencode_csv.index[maybe_point[j]],"pseudotime"]
        pic_min_loc = maybe_point[j]
        in_pesudo = abs(after_time_pesuda - before_time_pesuda)
        text_presudo = abs(after_time_pesuda - before_time_pesuda)
        if(min_pesudo > in_pesudo):
            min_pesudo = in_pesudo
            min_pesudo_loc = pic_min_loc
        j = j + 1
        if((after_time_pesuda > before_time_pesuda or text_presudo < 0.05) or j >= similar_threshold):
            break
    if(j == similar_threshold):
        distance_num = distance_num + pic_distance_csv_np[pic_min_loc][i]

        break_num =  break_num  + 1
    else:
        distance_num = distance_num + pic_distance_csv_np[pic_min_loc][i]
        edges.append((pic_distance_csv.index[maybe_point[j]],pic_distance_csv.columns[i],1))

    pesudo_all = pesudo_all + abs(after_time_pesuda - before_time_pesuda)
    edges_set = set()
    for ed in edges:
        edges_set.add(ed[0])



nodes  = all_single_cell_csv.index
G = nx.Graph()
G.add_nodes_from(nodes)
G.add_weighted_edges_from(edges)
all_num = len(all_single_cell_csv.index)
try:
    column = ["before_cell","before_cell_type","before_sub_cell","before_pseudotime","before_row","before_col","before_pic_row","before_pic_col",
            "similar",
            "after_cell","after_cell_type","after_sub_cell","after_pseudotime","after_row","after_col","after_pic_row","after_pic_col"]
    edges_csv = pd.DataFrame(columns=column)
    for i in range(len(edges)) :
        before_cell_series = all_single_cell_csv.loc[edges[i][0]]
        before_time_pesuda = new_pseudo_time.loc[edges[i][0],"pseudotime"]
        after_cell_series = all_single_cell_csv.loc[edges[i][1]]
        after_time_pesuda = new_pseudo_time.loc[edges[i][1],"pseudotime"]
        series = pd.Series({"before_cell":edges[i][0],
                            "before_cell_type":before_cell_series["celltype"],
                            "before_sub_cell":before_cell_series["subtype"],
                            "before_pseudotime" : before_time_pesuda,
                            "before_row":before_cell_series["imagerow"],
                            "before_col":before_cell_series["imagecol"],
                            "before_pic_row":before_cell_series["pic_row"],
                            "before_pic_col":before_cell_series["pic_col"],
                            "similar":edges[i][2],
                            "after_cell":edges[i][1],
                            "after_cell_type":after_cell_series["celltype"],
                            "after_sub_cell":after_cell_series["subtype"],
                            "after_pseudotime":after_time_pesuda,
                            "after_row":after_cell_series["imagerow"],
                            "after_col":after_cell_series["imagecol"],
                            "after_pic_row":after_cell_series["pic_row"],
                            "after_pic_col":after_cell_series["pic_col"],
                            },name = i)
        edges_csv = edges_csv.append(series)
    edges_csv.to_csv(TL_edges_save_file)
except KeyError as e:
    print("Cell does not have this property:" + e.args[0])
    column = ["before_cell",
            "similar",
            "after_cell"]
    edges_csv = pd.DataFrame(columns=column)
    for i in range(len(edges)) :
        before_cell_series = all_single_cell_csv.loc[edges[i][0]]
        before_time_pesuda = new_pseudo_time.loc[edges[i][0],"pseudotime"]
        after_cell_series = all_single_cell_csv.loc[edges[i][1]]
        after_time_pesuda = new_pseudo_time.loc[edges[i][1],"pseudotime"]
        series = pd.Series({"before_cell":edges[i][0],
                            "before_cell_type":before_cell_series["celltype"],
                            "similar":edges[i][2],
                            "after_cell":edges[i][1],
                            "after_cell_type":after_cell_series["celltype"],
                            },name = i)
        edges_csv = edges_csv.append(series)
    edges_csv.to_csv(samll_TL_edges_save_file)
print("TL_get is over")