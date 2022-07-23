

import numpy as np
import pandas as pd
import scanpy as sc
from rich.progress import track
import warnings
import os
warnings.filterwarnings("ignore")
dir_root = os.getcwd()
two_time_anndata_file =  dir_root+"/data/intermediate_result/tow_time_adata.h5ad"
loc_distance_file =dir_root+"/data/intermediate_result/tow_time_pic_distances.csv"
choose_hvg_distance_file = dir_root+"/data/intermediate_result/tow_time_choose_hvg_distances.csv"
autoencode_reduce_file = dir_root+"/data/intermediate_result/tow_time_auto_encode_dim_red_hvg_32.csv"
autoencode_distance_comp_file = dir_root+"/data/intermediate_result/tow_time_auto_encode_distances.csv"
single_cell_adata = sc.read(two_time_anndata_file)
before_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "before"]
after_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "after"]
before_num = len(before_single_cell_csv.index)
after_num = len(after_single_cell_csv.index)
pic_distance_csv = pd.read_csv(loc_distance_file,index_col=0)
pic_distance_csv_np = pic_distance_csv.values
pic_distance_csv_np  = (pic_distance_csv_np - np.min(pic_distance_csv_np)) * (1 / (np.max(pic_distance_csv_np) - np.min(pic_distance_csv_np)))
pic_distance_csv= pd.DataFrame(data = pic_distance_csv_np,index = before_single_cell_csv.index,columns = after_single_cell_csv.index)
pic_distance_csv.to_csv(loc_distance_file)
two_time_choose_hvg_csv = pd.read_csv(choose_hvg_distance_file,index_col=0)
two_time_choose_hvg_np  = two_time_choose_hvg_csv.values
two_time_choose_hvg_np  = (two_time_choose_hvg_np - np.min(two_time_choose_hvg_np)) * (1 / (np.max(two_time_choose_hvg_np) - np.min(two_time_choose_hvg_np)))
two_time_choose_hvg_csv = pd.DataFrame(data=two_time_choose_hvg_np,index=before_single_cell_csv.index,columns=after_single_cell_csv.index)
two_time_choose_hvg_csv.to_csv(choose_hvg_distance_file)
autoencode_reduce_dr = pd.read_csv(autoencode_reduce_file,index_col=0)
autoencode_reduce_dr_np = autoencode_reduce_dr.values
two_time_autoencode_np = np.zeros([before_num,after_num])
for i in track(sequence =  range(before_num),transient = True):
    for j in range(after_num):
        two_time_autoencode_np[i][j] = np.linalg.norm(autoencode_reduce_dr_np[:,i] - autoencode_reduce_dr_np[:,j+before_num])
two_time_autoencode_np  = (two_time_autoencode_np - np.min(two_time_autoencode_np)) * (1 / (np.max(two_time_autoencode_np) - np.min(two_time_autoencode_np)))
two_time_autoencode_csv = pd.DataFrame(data = two_time_autoencode_np,index = before_single_cell_csv.index,columns = after_single_cell_csv.index)
two_time_autoencode_csv.to_csv(autoencode_distance_comp_file)
print("prepare_comp is over")