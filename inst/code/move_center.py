import pandas as pd
import os
dir_root = os.getcwd()

#parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')

begin_file_name = dir_root+"/data/before_iterative_mapping_result.csv"
begin_mapping_csv = pd.read_csv(begin_file_name,index_col = 0)
begin_length_col = max(begin_mapping_csv.imagecol) - min(begin_mapping_csv.imagecol)
begin_length_row = max(begin_mapping_csv.imagerow) - min(begin_mapping_csv.imagerow)
begin_center_col = min(begin_mapping_csv.imagecol) + begin_length_col / 2
begin_center_row = min(begin_mapping_csv.imagerow) + begin_length_row / 2
begin_mapping_csv["pic_col"] = begin_mapping_csv["imagecol"] - begin_center_col
begin_mapping_csv["pic_row"] = begin_mapping_csv["imagerow"] - begin_center_row
begin_mapping_csv.to_csv(begin_file_name)
file_name = dir_root+"/data/after_iterative_mapping_result.csv"
mapping_csv = pd.read_csv(file_name,index_col = 0)
length_col = max(mapping_csv.imagecol) - min(mapping_csv.imagecol)
length_row = max(mapping_csv.imagerow) - min(mapping_csv.imagerow)
center_col = min(mapping_csv.imagecol) + length_col / 2
center_row = min(mapping_csv.imagerow) + length_row / 2
mapping_csv["pic_col"] = mapping_csv["imagecol"] - center_col
mapping_csv["pic_row"] = mapping_csv["imagerow"] - center_row
mapping_csv["pic_col"] = mapping_csv["pic_col"] * begin_length_col / length_col * 1.2
mapping_csv["pic_row"] = mapping_csv["pic_row"] * begin_length_row / length_row * 1.2
mapping_csv.to_csv(file_name)
print("move_center is over")