import plotly.graph_objects as go
import pandas as pd
import networkx as nx
import scanpy as sc
import sys
import warnings
import os
warnings.filterwarnings("ignore")
dir_root = os.getcwd()

parameter_settings_f=pd.read_csv(dir_root+'/data/parameter_settings.csv')
# cell_type_choose= sys.argv[1]
cell_type_choose = str(parameter_settings_f['parameter'][0])

two_time_anndata_file = dir_root+"/data/intermediate_result/tow_time_adata.h5ad"
TL_edges_save_file = dir_root+"/data/result/tow_time_TL_edges.csv"
single_cell_adata = sc.read(two_time_anndata_file)
loc_csv = single_cell_adata.obs
loc_csv.loc[single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "before"].index,"y"] = 0
loc_csv.loc[single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "before"].index,"color"] = '#9c1b45'
loc_csv.loc[single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "after"].index,"y"] = 20
loc_csv.loc[single_cell_adata.obs.loc[single_cell_adata.obs["time"] == "after"].index,"color"] = '#f26d44'
G = nx.Graph()
node_traces = []
sub_celltype = []
nodes  = loc_csv.index
G.add_nodes_from(nodes)
for j in loc_csv.index:
    G.nodes[j]["pos"] = [loc_csv.loc[j,"pic_row"],loc_csv.loc[j,"pic_col"],loc_csv.loc[j,"y"],
                        loc_csv.loc[j]["color"]]
node_trace = go.Scatter3d(
    x = loc_csv.pic_row,
    y = loc_csv.y,
    z = loc_csv.pic_col,

    marker=dict(
    size=2,
    color = "#D9D9D9",
    ),
    mode= "markers")
node_traces.append(node_trace)
TL_edges_csv = pd.read_csv(TL_edges_save_file,index_col=0)

edges = []
edges_csv = TL_edges_csv.loc[TL_edges_csv["before_cell_type"] == cell_type_choose]
for j in edges_csv.index:
    edges.append((edges_csv.loc[j,"before_cell"], edges_csv.loc[j,"after_cell"],1) )

edges_csv["before_y"] = 0
edges_csv["after_y"] = 20
node_trace = go.Scatter3d(
x = edges_csv.before_pic_row,
y = edges_csv.before_y,
z = edges_csv.before_pic_col,
marker=dict(
size=2,
color = '#9c1b45',
),
mode= "markers")
node_traces.append(node_trace)
G.add_weighted_edges_from(edges)
edge_traces = []
edge_x = []
edge_y = []
edge_z = []
edge_color = []
line_color_list = []
for edge in G.edges():

    x0, z0, y0 ,edge_color0= G.nodes[edge[0]]['pos'][0:4]
    x1, z1, y1 ,edge_color1= G.nodes[edge[1]]['pos'][0:4]
    edge_x.append(x0)
    edge_x.append(x1)
    edge_x.append(None)
    edge_y.append(y0)
    edge_y.append(y1)
    edge_y.append(None)
    edge_z.append(z0)
    edge_z.append(z1)
    edge_z.append(None)
    edge_color.append(edge_color0)
    edge_color.append(edge_color1)
    edge_color.append('#888')
    edge_trace = go.Scatter3d(
        x=edge_x, y=edge_y, z = edge_z,
        legendgrouptitle_text="Line Title",
        line=dict(width=1, color=edge_color),
        hoverinfo='none',
        mode='lines')
    edge_traces.append(edge_trace)
fig = go.Figure(node_traces + edge_traces)
noaxis=dict(showbackground=False,
            showgrid=False,
            showline=False,
            showticklabels=False,
            ticks='',
            title='',
            zeroline=False)
camera = dict(
    eye=dict(x= 3, y= 0.15, z=0.2)
)
fig.update_layout(
    height=800,
    width=800,
    scene=dict(
        xaxis = noaxis,
        yaxis = noaxis,
        zaxis = noaxis,
        aspectratio=dict(x=1, y=1, z=1)#改变画布空间比例为1：1：1
        ),
    margin=dict(r=0, l=0, b=0, t=0),
    scene_camera = camera
)

fig.write_html(dir_root+"/data/result/3D.html")
fig.write_image(dir_root+"/data/result/3D.pdf")
print("pic is over")
