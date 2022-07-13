#生成准备文件
#对计算出来的矩阵，进行归一化操作
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
#from sklearn.decomposition import PCA
import copy
#import scipy.stats
from sklearn import preprocessing
import sys
import os
#import scipy.stats
#time_points = ["S1.5d"]
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


two_time_anndata_file =  dir_root+"/data/tow_time_adata.h5ad"
#pseudo_file = "/home/sunhang/data/planarian/new_connectivities/cell_pseudotime_n_tip_5000.csv"
loc_distance_file =dir_root+"/data/tow_time_pic_distances.csv"
#pseudo_time_file = "/home/sunhang/data/planarian/change/new_connectivities/"+ before_time +"_"+ after_time + "_pseudo_time.csv"
choose_hvg_distance_file = dir_root+"/data/tow_time_choose_hvg_distances.csv"
#autoencode_reduce_file = "/home/sunhang/data/planarian/new2_connectivities/"+ before_time +"_"+ after_time + "_auto_encode_dim_red_256.csv"
autoencode_reduce_file = dir_root+"/data/tow_time_auto_encode_dim_red_hvg_32.csv"
#autoencode_distance_comp_file = "/home/sunhang/data/planarian/intergrate_connectivities/"+ before_time +"_"+ after_time + "_auto_encode_dim_red_256_distance.csv"

hvg_distance_comp_file = dir_root+"/data/tow_time_hvg_distance.csv"

single_cell_adata = sc.read(two_time_anndata_file)
before_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "before"]
after_single_cell_csv = single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "after"]
before_num = len(before_single_cell_csv.index)
after_num = len(after_single_cell_csv.index)

#先算autoencode的相似性
autoencode_reduce_dr = pd.read_csv(autoencode_reduce_file,index_col=0)
autoencode_reduce_dr_np = autoencode_reduce_dr.values
two_time_autoencode_np = np.zeros([before_num,after_num])
#two_time_autoencode_csv = pd.DataFrame(data = two_time_autoencode_np,index = before_single_cell_csv.index,columns = after_single_cell_csv.index)
for i in range(before_num):
    #print(i)
    for j in range(after_num):
        two_time_autoencode_np[i][j] = np.linalg.norm(autoencode_reduce_dr_np[:,i] - autoencode_reduce_dr_np[:,j+before_num])
#two_time_autoencode_np  = two_time_autoencode_csv.values
#归一化
#print(np.max(pca_distance_csv_np) - np.min(pca_distance_csv_np))
two_time_autoencode_np  = (two_time_autoencode_np - np.min(two_time_autoencode_np)) * (1 / (np.max(two_time_autoencode_np) - np.min(two_time_autoencode_np)))
#print(two_time_autoencode_csv)
#print(two_time_autoencode_np)
two_time_autoencode_csv = pd.DataFrame(data = two_time_autoencode_np,index = before_single_cell_csv.index,columns = after_single_cell_csv.index)
#print(two_time_autoencode_csv)
two_time_autoencode_csv.to_csv(autoencode_reduce_file)
#算pic上的相似性

pic_distance_csv = pd.read_csv(loc_distance_file,index_col=0)
pic_distance_csv_np = pic_distance_csv.values
#归一化
#print(np.max(pca_distance_csv_np) - np.min(pca_distance_csv_np))
pic_distance_csv_np  = (pic_distance_csv_np - np.min(pic_distance_csv_np)) * (1 / (np.max(pic_distance_csv_np) - np.min(pic_distance_csv_np)))
pic_distance_csv= pd.DataFrame(data = pic_distance_csv_np,index = before_single_cell_csv.index,columns = after_single_cell_csv.index)
pic_distance_csv.to_csv(loc_distance_file)
#算hvg的相似性
sc.pp.highly_variable_genes(single_cell_adata, flavor="seurat", n_top_genes=3000)
two_times_single_cell_gene = pd.DataFrame(data = single_cell_adata.X.T , index= single_cell_adata.var.index,columns=single_cell_adata.obs.index)
cc = []
num = 0
for i in single_cell_adata.var.highly_variable:
    #print(i)
    if i == True:
        num = num +1
    cc.append(i)
    #print(num)|
two_times_single_cell_gene_hvg = two_times_single_cell_gene.iloc[cc]
two_times_single_cell_gene_hvg_np = two_times_single_cell_gene_hvg.values
two_times_hvg_dis_np = np.zeros([before_num,after_num])
#two_times_hvg_dis = pd.DataFrame(data=two_times_hvg_dis_np,index=before_single_cell_csv.index,columns=after_single_cell_csv.index)
for i in range(before_num):
    #print(i)
    for j in range(after_num):
        two_times_hvg_dis_np[i][j] = np.linalg.norm(two_times_single_cell_gene_hvg_np[:,i] - two_times_single_cell_gene_hvg_np[:,j+before_num])
#two_times_hvg_dis_np = two_times_hvg_dis.values
#归一化
#print(np.max(pca_distance_csv_np) - np.min(pca_distance_csv_np))
two_times_hvg_dis_np  = (two_times_hvg_dis_np - np.min(two_times_hvg_dis_np)) * (1 / (np.max(two_times_hvg_dis_np) - np.min(two_times_hvg_dis_np)))
two_times_hvg_dis = pd.DataFrame(data=two_times_hvg_dis_np,index=before_single_cell_csv.index,columns=after_single_cell_csv.index)
two_times_hvg_dis.to_csv(hvg_distance_comp_file)
#算choose_hvg的相似性
two_time_choose_hvg_csv = pd.read_csv(choose_hvg_distance_file,index_col=0)
two_time_choose_hvg_np  = two_time_choose_hvg_csv.values
#归一化
#print(np.max(pca_distance_csv_np) - np.min(pca_distance_csv_np))
two_time_choose_hvg_np  = (two_time_choose_hvg_np - np.min(two_time_choose_hvg_np)) * (1 / (np.max(two_time_choose_hvg_np) - np.min(two_time_choose_hvg_np)))
two_time_choose_hvg_csv = pd.DataFrame(data=two_time_choose_hvg_np,index=before_single_cell_csv.index,columns=after_single_cell_csv.index)
two_time_choose_hvg_csv.to_csv(choose_hvg_distance_file)
print("yes")