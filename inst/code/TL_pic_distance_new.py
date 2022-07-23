import numpy as np
import pandas as pd
from rich.progress import track
from pathlib import Path
import warnings
import os
warnings.filterwarnings("ignore")
dir_root = os.getcwd()
TILE_PATH = Path(dir_root+"/data/intermediate_result")
TILE_PATH.mkdir(parents=True, exist_ok=True)
begin_file_name = dir_root+"/data/before_iterative_mapping_result.csv"
before_single_cell_csv = pd.read_csv(begin_file_name,index_col = 0)
after_mapping_file =dir_root+"/data/after_iterative_mapping_result.csv"
after_single_cell_csv = pd.read_csv(after_mapping_file,index_col = 0)
loc_distance_file = dir_root+"/data/intermediate_result/tow_time_pic_distances.csv"
before_num = len(before_single_cell_csv.index)
after_num = len(after_single_cell_csv.index)
pic_distance_csv_np = np.zeros([before_num,after_num])
pic_distance_csv = pd.DataFrame(data=pic_distance_csv_np,index=before_single_cell_csv.index,columns=after_single_cell_csv.index)
for i in track(sequence =  range(before_num),transient = True):
    before_loc_np = np.array([before_single_cell_csv.loc[before_single_cell_csv.index[i],"pic_row"],before_single_cell_csv.loc[before_single_cell_csv.index[i],"pic_col"]])
    for j in range(after_num):
        after_loc_np = np.array([after_single_cell_csv.loc[after_single_cell_csv.index[j],"pic_row"],after_single_cell_csv.loc[after_single_cell_csv.index[j],"pic_col"]])
        pic_distance_csv_np[i][j] = np.linalg.norm(before_loc_np - after_loc_np)
pic_distance_csv.to_csv(loc_distance_file)
print("TL_pic_distance_new is over")