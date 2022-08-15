import numpy as np
import pandas as pd
from rich.progress import track
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")
import os
dir_root = os.getcwd()
TILE_PATH = Path(dir_root+"/data/intermediate_result")
TILE_PATH.mkdir(parents=True, exist_ok=True)
begin_file_name = dir_root+"/data/before_iterative_mapping_result.csv"
before_single_cell_csv = pd.read_csv(begin_file_name,index_col = 0)
after_file_name =dir_root+"/data/after_iterative_mapping_result.csv"
after_single_cell_csv = pd.read_csv(after_file_name,index_col = 0)
loc_distance_file = dir_root+"/data/intermediate_result/tow_time_pic_distances.csv"
before_num = len(before_single_cell_csv.index)
after_num = len(after_single_cell_csv.index)
all_single_cell_csv = before_single_cell_csv.append(after_single_cell_csv)
cell_loc = all_single_cell_csv.loc[:,("pic_col","pic_row")].values
from scipy.spatial.distance import pdist, squareform
distance_np = pdist(cell_loc, metric = "euclidean")
distance_np_X =squareform(distance_np)
pic_distance_csv_np = distance_np_X[0:before_num,before_num:after_num + before_num]
pic_distance_csv_np  = (pic_distance_csv_np - np.min(pic_distance_csv_np)) * (1 / (np.max(pic_distance_csv_np) - np.min(pic_distance_csv_np)))
pic_distance_csv= pd.DataFrame(data = pic_distance_csv_np,index = before_single_cell_csv.index,columns = after_single_cell_csv.index)

pic_distance_csv.to_csv(loc_distance_file)
print("TL_pic_distance_new is over")
