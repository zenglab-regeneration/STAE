#计算两个面的细胞之间gene表达谱之间的距离性和连接性
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import plotly.express as px
import plotly.graph_objects as go
import sys
import os
dir_root=os.getcwd()
def mkdir(path):
    # 去除首位空格
    path=path.strip()
    # 去除尾部 \ 符号
    path=path.rstrip("\\")
    # 判断路径是否存在
    # 存在     True
    # 不存在   False
    isExists=os.path.exists(path)
    # 判断结果
    if not isExists:
        # 如果不存在则创建目录
        # 创建目录操作函数
        os.makedirs(path) 
        print (path+' 创建成功')
        return True
    else:
        # 如果目录存在则不创建，并提示目录已存在
        print(path+' 目录已存在')
        return False
file_name = dir_root+"/data/result/"
mkdir(file_name)
two_time_anndata_file = dir_root+"/data/tow_time_adata.h5ad"
#导入文件每个两个时间点的mapping结果和这个两个时间点的sc_gene表达
before_sc_file = dir_root+"/data/before_sc_data.csv"
after_sc_file = dir_root+"/data/after_sc_data.csv"
#"/home/sunhang/data/pl_new/4/" + after_time + "/mapping_result/" + after_time + "_iterative_mapping_result_celltype_" + method + ".csv"
before_mapping_file = dir_root+"/data/before_iterative_mapping_result.csv"
after_mapping_file =dir_root+"/data/after_iterative_mapping_result.csv"
before_sc_csv = pd.read_csv(before_sc_file,index_col = 0)
after_sc_csv = pd.read_csv(after_sc_file,index_col = 0)
before_mapping_csv = pd.read_csv(before_mapping_file,index_col = 0 )
before_mapping_csv["time"] = "before"
after_mapping_csv = pd.read_csv(after_mapping_file,index_col = 0)
after_mapping_csv["time"] = "after"
comb_mapping_result_csv = before_mapping_csv.append(after_mapping_csv)
comb_sc_csv = pd.merge(before_sc_csv,after_sc_csv,left_index=True, right_index=True)
comb_sc_csv = comb_sc_csv.fillna(0)
in_sc_csv = pd.DataFrame(index=comb_sc_csv.index)
#生成single_cell_adata数据，将gene表达和mapping结果一起放进去
for sc_name in comb_mapping_result_csv.index :  
    #print(sc_name)
    celltype = comb_mapping_result_csv.loc[sc_name]["celltype_x"]
    single_cell_name = comb_mapping_result_csv.loc[sc_name]["single_cell_name"]
    in_sc_csv[sc_name] = comb_sc_csv.loc[:,single_cell_name]
single_cell_adata = ad.AnnData(in_sc_csv.values.T)
single_cell_adata.obs = comb_mapping_result_csv
single_cell_adata.write(two_time_anndata_file)
print("TL_gene_simliar_comp is over")