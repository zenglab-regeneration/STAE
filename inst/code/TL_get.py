import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import networkx as nx
import scanpy as sc
import sys
import os
import random
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
result_file = dir_root+"/data/pl_new/4/TL/" + "result"
mkdir(result_file)
subcelltype_color = {}

method = "high_it_Hvg"
must_fix = [
    [1, 0, 0, 25],
    [1, 0, 0, 25],
    [1, 0, 0, 25],
    [1, 0, 0, 25],
]
def bubbleSort(arr,serial):
    for i in range(1, len(arr)):
        for j in range(0, len(arr)-i):
            if arr[j] > arr[j+1]:
                serial[j],serial[j + 1] = serial[j + 1],serial[j]
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
    return arr , serial
#for ratio in np.arange(0,1.1,0.1):
ratio = 0.1

resolution = 4



two_time_anndata_file = dir_root+"/data/tow_time_adata.h5ad"
pseudo_file = dir_root+"/data/pseudotime.csv"
#new_file_name = "/home/sunhang/data/pl_new/4/TL/" + before_time +"_"+ after_time + "/intergrate_connectivities"
#autoencode_distance_comp_file = "/home/sunhang/data/planarian/intergrate_connectivities/"+ before_time +"_"+ after_time + "_auto_encode_dim_red_256_distance.csv"
autoencode_distance_comp_file = dir_root+"/data/tow_time_auto_encode_dim_red_hvg_32.csv"
before_time_mapping_file = dir_root+"/data/before_iterative_mapping_result.csv"
after_time_mapping_file = dir_root+"/data/after_iterative_mapping_result.csv"
hvg_distance_comp_file = dir_root+"/data/tow_time_hvg_distance.csv"
pic_distance_comp_file = dir_root+"/data/tow_time_pic_distances.csv"
choose_hvg_distanc_comp_file = dir_root+"/data/tow_time_choose_hvg_distances.csv"
#TL_edges_save_file = "/home/sunhang/data/planarian/intergrate_connectivities/"+ before_time +"_"+ after_time + "_TL_edges_256.csv"
TL_edges_save_file = dir_root+"/data/tow_time_TL_edges_32.csv"
#single_cell_adata = sc.read(two_time_anndata_file)
before_single_cell_csv = pd.read_csv(before_time_mapping_file,index_col= 0)
after_single_cell_csv = pd.read_csv(after_time_mapping_file,index_col= 0)
all_single_cell_csv = before_single_cell_csv.append(after_single_cell_csv)
before_num = len(before_single_cell_csv.index)
after_num = len(after_single_cell_csv.index)
two_time_autoencode_csv = pd.read_csv(autoencode_distance_comp_file,index_col=0)
two_time_autoencode_np = two_time_autoencode_csv.values
#print(two_time_autoencode_csv)
pic_distance_csv = pd.read_csv(pic_distance_comp_file,index_col=0)
pic_distance_csv_np = pic_distance_csv.values
#print(pic_distance_csv)
two_times_hvg_dis  = pd.read_csv(choose_hvg_distanc_comp_file, index_col=0)
two_times_hvg_dis_np = two_times_hvg_dis.values
#print(two_times_hvg_dis)
two_time_choose_hvg_csv = pd.read_csv(choose_hvg_distanc_comp_file,index_col=0)
two_time_choose_hvg_np  = two_time_choose_hvg_csv.values
#print(two_time_choose_hvg_csv)
pseudo_csv = pd.read_csv(pseudo_file,index_col=0)
all_index = pic_distance_csv.index.append(pic_distance_csv.columns)
new_pseudo_time = pd.DataFrame(index=all_single_cell_csv.index,columns=["pseudotime"])
#print(new_pseudo_time.index)
for i in new_pseudo_time.index:
    #print(pseudo_csv.loc[i.split(".")[0],"pseudotime"])
    cell_name = all_single_cell_csv.loc[i,"single_cell_name"].split(" ")[0]
    new_pseudo_time.loc[i,"pseudotime"] =  pseudo_csv.loc[cell_name,"pseudotime"]
#must_fix[i]
wuwu = 1
haha = 0
lala = 0
huhu = 25
two_distence = wuwu * two_time_autoencode_np + haha * two_time_choose_hvg_np + lala * two_times_hvg_dis_np + ratio * pic_distance_csv_np
num = huhu
break_num = 0
distance_num = 0
#print(huhu)
similar_threshold = huhu
edges = []
maybe_point = []
num = 0
#pca_distance_csv_np = pca_distance_csv.valueswo
#pic_distance_csv_np = pic_distance_csv.values
pesudo_all = 0
for i in range(after_num):
    after_time_pesuda = new_pseudo_time.loc[two_time_autoencode_csv.columns[i],"pseudotime"]
    in_sort = np.argsort(two_distence[:,i] )
    min_similar = two_distence[in_sort[0],i]
    maybe_point = in_sort[0:similar_threshold]
    #print(maybe_point)
    #after_time_pesuda = before_time_pesuda - 1
    cell_gene_similar_value = []
    for pic_cell in maybe_point:
        #print(pic_distance_csv_np[pic_cell,i])
        cell_gene_similar_value.append( 0 * pic_distance_csv_np[pic_cell,i] +
                                                                    1 * two_time_choose_hvg_np[pic_cell,i] + 
                                                                    0 * two_times_hvg_dis_np[pic_cell,i])
    #print(cell_gene_similar_value[0])
    #print(pic_cell_gene_similar_value)
    cell_gene_similar_value,maybe_point = bubbleSort(cell_gene_similar_value,maybe_point)
    #print(maybe_point)
    #after_time_pesuda = before_time_pesuda - 1
    j = 0
    #print(before_time_pesuda)
    #print(after_time_pesuda)
    #print(after_time_pesuda > before_time_pesuda)
    #print(abs(after_time_pesuda - before_time_pesuda) < 0.1)
    min_pesudo = 1
    min_pesudo_loc = maybe_point[0]
    while(True):
        #print(1)
        #pic_min_loc = maybe_point[j]
        before_time_pesuda = new_pseudo_time.loc[two_time_autoencode_csv.index[maybe_point[j]],"pseudotime"]
        #print(j)
        #print(before_time_pesuda)
        pic_min_loc = maybe_point[j]
        in_pesudo = abs(after_time_pesuda - before_time_pesuda)
        text_presudo = abs(after_time_pesuda - before_time_pesuda)
        if(min_pesudo > in_pesudo):
            min_pesudo = in_pesudo
            min_pesudo_loc = pic_min_loc

        """
        pic_min_num = pic_distance_csv_np[maybe_point[0],i]
        pic_min_loc = maybe_point[0]
        for j in maybe_point:
            if( pic_distance_csv_np[maybe_point[0],i] < pic_min_num ):
                pic_min_loc = j
                pic_min_num = pic_distance_csv_np[maybe_point[j],i]
        print(pic_min_loc)
        """

        j = j + 1
        #text_presudo = abs(after_time_pesuda - before_time_pesuda)
        #after_time_pesuda = new_pseudo_time.loc[pca_distance_csv_new.index[pic_min_loc],"pseudotime"]
        #edges.append((pic_distance_csv.index[maybe_point[j]],pic_distance_csv.columns[i],1))
        #print(before_time_pesuda)
        #print(after_time_pesuda)
        #print("#############")
        if((after_time_pesuda > before_time_pesuda or text_presudo < 0.05) or j >= similar_threshold):
            break
    #print(j)
    if(j == similar_threshold):
        distance_num = distance_num + pic_distance_csv_np[pic_min_loc][i]
        #edges.append((pic_distance_csv.index[pic_min_loc],pic_distance_csv.columns[i],1))
        break_num =  break_num  + 1
    else:
        distance_num = distance_num + pic_distance_csv_np[pic_min_loc][i]
        #print(pic_distance_csv_np[pic_min_loc][i])
        edges.append((pic_distance_csv.index[maybe_point[j]],pic_distance_csv.columns[i],1))

    pesudo_all = pesudo_all + abs(after_time_pesuda - before_time_pesuda)
    edges_set = set()
    for ed in edges:
        edges_set.add(ed[0])


#print("所有伪时间的差值和"+str(pesudo_all))
#print("阈值之外的细胞数目"+str(break_num))
#print(distance_num)
print(len(edges_set))
print("##########")
        #loc_similar = 
        #print(maybe_point)
        #print(pca_distance_csv_new_np[in_sort[0],i])
nodes  = all_single_cell_csv.index
#定义graph
G = nx.Graph()
G.add_nodes_from(nodes)
G.add_weighted_edges_from(edges)
all_num = len(all_single_cell_csv.index)

column = ["before_cell","before_cell_type","before_sub_cell","before_pseudotime","before_row","before_col","before_pic_row","before_pic_col",
        "similar",
        "after_cell","after_cell_type","after_sub_cell","after_pseudotime","after_row","after_col","after_pic_row","after_pic_col"]
edges_csv = pd.DataFrame(columns=column)
for i in range(len(edges)) :
    before_cell_series = all_single_cell_csv.loc[edges[i][0]]
    #print(before_cell_series)

    before_time_pesuda = new_pseudo_time.loc[edges[i][0],"pseudotime"]
    after_cell_series = all_single_cell_csv.loc[edges[i][1]]
    after_time_pesuda = new_pseudo_time.loc[edges[i][1],"pseudotime"]
    series = pd.Series({"before_cell":edges[i][0],
                        "before_cell_type":before_cell_series["celltype_x"],
                        "before_sub_cell":before_cell_series["final_subtype"],
                        "before_pseudotime" : before_time_pesuda,
                        "before_row":before_cell_series["imagerow"],
                        "before_col":before_cell_series["imagecol"],
                        "before_pic_row":before_cell_series["pic_row"],
                        "before_pic_col":before_cell_series["pic_col"],
                        "similar":edges[i][2],
                        "after_cell":edges[i][1],
                        "after_cell_type":after_cell_series["celltype_x"],
                        "after_sub_cell":after_cell_series["final_subtype"],
                        "after_pseudotime":after_time_pesuda,
                        "after_row":after_cell_series["imagerow"],
                        "after_col":after_cell_series["imagecol"],
                        "after_pic_row":after_cell_series["pic_row"],
                        "after_pic_col":after_cell_series["pic_col"],
                        },name = i)
    edges_csv = edges_csv.append(series)
num = 0
all_num = 0
for i in edges_csv.index:
    if(edges_csv.loc[i,"before_cell_type"] == edges_csv.loc[i,"after_cell_type"] or edges_csv.loc[i,"before_cell_type"] == "tgs1+_Neoblast" ):
        num = num + 1
    all_num = all_num + 1
edges_csv.to_csv(TL_edges_save_file)