import numpy as np
import pandas as pd
import scanpy as sc
from rich.progress import track
import os
dir_root = os.getcwd()
two_time_anndata_file = dir_root+"/data/intermediate_result/tow_time_adata.h5ad"
single_cell_adata = sc.read(two_time_anndata_file)
before_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "before"]
after_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "after"]
before_num = len(before_single_cell_csv.index)
after_num = len(after_single_cell_csv.index)
choose_before_sc_data_csv = pd.DataFrame(data=single_cell_adata.X.T[:,0:before_single_cell_csv.shape[0]],index= single_cell_adata.var.index,columns=before_single_cell_csv.index)
choose_after_sc_data_csv = pd.DataFrame(data=single_cell_adata.X.T[:,before_single_cell_csv.shape[0]:single_cell_adata.X.T.shape[1]],index= single_cell_adata.var.index,columns=after_single_cell_csv.index)
gene_list_file = dir_root+"/data/marker_gene.csv"
choose_hvg_distance_file = dir_root+"/data/intermediate_result/tow_time_choose_hvg_distances.csv"
before_sc_csv = choose_before_sc_data_csv
after_sc_csv = choose_after_sc_data_csv
gene_list_csv = pd.read_csv(gene_list_file,index_col=0)
before_gene_choose_csv = before_sc_csv.loc[gene_list_csv.index]
after_gene_choose_csv = after_sc_csv.loc[gene_list_csv.index]
result_np = np.zeros([len(before_gene_choose_csv.columns),len(after_gene_choose_csv.columns)])
result_csv = pd.DataFrame(data=result_np , index = before_gene_choose_csv.columns,columns=after_gene_choose_csv.columns)
before_gene_choose_np = before_gene_choose_csv.values
after_gene_choose_np = after_gene_choose_csv.values
for i in track(sequence =  range(len(before_gene_choose_csv.columns)),transient = True):
    for j in range(len(after_gene_choose_csv.columns)):
        result_np[i][j] = np.linalg.norm(before_gene_choose_np[:,i] - after_gene_choose_np[:,j])
result_csv.to_csv(choose_hvg_distance_file) 
print("hvg_gene_distance is over")