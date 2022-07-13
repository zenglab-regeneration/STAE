import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
from sklearn.decomposition import PCA
import copy
import scipy.stats
import copy
import sys
import os
dir_root=os.getcwd()

two_time_anndata_file = dir_root+"/data/tow_time_adata.h5ad"
#导入文件每个两个时间点的mapping结果和这个两个时间点的sc_gene表达
before_sc_file = dir_root+"/data/before_sc_data.csv"
after_sc_file = dir_root+"/data/after_sc_data.csv"
single_cell_adata = sc.read(two_time_anndata_file)
before_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "before"]
after_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "after"]
before_num = len(before_single_cell_csv.index)
after_num = len(after_single_cell_csv.index)
before_choose_sc_csv = pd.read_csv(before_sc_file,index_col=0)
after_choose_sc_csv = pd.read_csv(after_sc_file,index_col=0)
choose_before_sc_data_np = np.zeros([len(before_choose_sc_csv.index),len(before_single_cell_csv.index)])
choose_before_sc_data_csv = pd.DataFrame(data=choose_before_sc_data_np,index= before_choose_sc_csv.index,columns=before_single_cell_csv.index)
#single_celll_name_list = []
for i in choose_before_sc_data_csv.columns:
    single_cell_name = before_single_cell_csv.loc[i ,"single_cell_name"]
    #single_celll_name_list.append(single_cell_name)
    choose_before_sc_data_csv.loc[:,i] = before_choose_sc_csv.loc[:,single_cell_name]
choose_after_sc_data_np = np.zeros([len(after_choose_sc_csv.index),len(after_single_cell_csv.index)])
choose_after_sc_data_csv = pd.DataFrame(data=choose_after_sc_data_np,index= after_choose_sc_csv.index,columns=after_single_cell_csv.index)
for i in choose_after_sc_data_csv.columns:
    single_cell_name = after_single_cell_csv.loc[i ,"single_cell_name"]
    #single_celll_name_list.append(single_cell_name)
    choose_after_sc_data_csv.loc[:,i] = after_choose_sc_csv.loc[:,single_cell_name]
#谱系特异性gene marker
gene_list_file = dir_root+"/data/marker_gene.csv"
choose_hvg_distance_file = dir_root+"/data/tow_time_choose_hvg_distances.csv"

before_sc_csv = choose_before_sc_data_csv
after_sc_csv = choose_after_sc_data_csv
gene_list_csv = pd.read_csv(gene_list_file,index_col=0)
before_gene_choose_csv = before_sc_csv.loc[gene_list_csv.index]
after_gene_choose_csv = after_sc_csv.loc[gene_list_csv.index]
result_np = np.zeros([len(before_gene_choose_csv.columns),len(after_gene_choose_csv.columns)])
result_csv = pd.DataFrame(data=result_np , index = before_gene_choose_csv.columns,columns=after_gene_choose_csv.columns)
before_gene_choose_np = before_gene_choose_csv.values
after_gene_choose_np = after_gene_choose_csv.values
for i in range(len(before_gene_choose_csv.columns)):
    for j in range(len(after_gene_choose_csv.columns)):
        result_np[i][j] = np.linalg.norm(before_gene_choose_np[:,i] - after_gene_choose_np[:,j])
result_csv.to_csv(choose_hvg_distance_file) 
print("hvg_gene_distance is over")