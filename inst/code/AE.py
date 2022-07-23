import torch
import torch.nn as nn
from torch import optim
import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.preprocessing import MinMaxScaler
import os
dir_root = os.getcwd()
class AutoEncoder(nn.Module):
    def __init__(self,n_feature,n_hidden_1, n_hidden_2,n_output):
        super(AutoEncoder,self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_feature,n_hidden_1),
            nn.ReLU(),
            nn.Linear(n_hidden_1,n_hidden_2),
            nn.ReLU(),
            nn.Linear(n_hidden_2,n_output),
            nn.ReLU(),
        )
        self.decoder = nn.Sequential(
            nn.Linear(n_output, n_hidden_2),
            nn.BatchNorm1d(num_features=n_hidden_2),
            nn.Linear(n_hidden_2, n_hidden_1),
            nn.BatchNorm1d(num_features=n_hidden_1),
            nn.Linear(n_hidden_1, n_feature),
            nn.BatchNorm1d(num_features=n_feature),
        )

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return encoded, decoded
two_time_anndata_file = dir_root+"/data/intermediate_result/tow_time_adata.h5ad"
autoencode_file = dir_root+"/data/intermediate_result/tow_time_auto_encode_dim_red_hvg_32.csv"
before_sc_file = dir_root+"/data/before_sc_data.csv"
after_sc_file = dir_root+"/data/after_sc_data.csv"

before_sc_csv = pd.read_csv(before_sc_file,index_col=0)
after_sc_csv = pd.read_csv(after_sc_file,index_col=0)
two_sc_csv = pd.merge(before_sc_csv,after_sc_csv,left_index=True,right_index=True)
single_cell_adata = sc.read(two_time_anndata_file)
sc.pp.highly_variable_genes(single_cell_adata, flavor="seurat", n_top_genes=3000)
gene_num = []
gene_name = []
num = 0
for i in single_cell_adata.var.highly_variable:
    if i == True:
        gene_num.append(str(num))
        gene_name.append(two_sc_csv.index[num])
    num = num + 1

two_sc_csv_hvg = two_sc_csv.loc[gene_name]
two_sc_np = two_sc_csv_hvg.values.T
two_sc_csv_hvg_T = pd.DataFrame(data=two_sc_csv_hvg.values.T,index=two_sc_csv_hvg.columns,columns=two_sc_csv_hvg.index)
two_sc_csv_hvg_T_np = two_sc_csv_hvg_T.values
N = len(two_sc_csv_hvg_T )
single_cell_feature = len(two_sc_csv_hvg.index)
feature = torch.FloatTensor(two_sc_csv_hvg_T_np)
index_all = np.array(range(0, N), dtype=np.int64)
np.random.shuffle(index_all)

idx_train = index_all[:int(N*0.7)]
idx_test = index_all[int(N*0.7):]

feature_train = feature[idx_train]
feature_test = feature[idx_test]

model = AutoEncoder(n_feature=single_cell_feature,n_hidden_1=512,n_hidden_2=128,n_output=32)
model.cuda()
optimizer = optim.Adam(model.parameters(), lr=0.01)
loss_func = nn.MSELoss()
for epoch in range(500):
    enc,dec = model(feature_train.cuda())
    loss = loss_func(dec, feature_train.cuda())
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
enc,dec = model(feature_test.cuda())
loss = loss_func(dec, feature_test.cuda())
#print(loss.item())
single_cell_adata = sc.read(two_time_anndata_file)
two_times_single_cell_gene = pd.DataFrame(data = single_cell_adata.X.T , index= single_cell_adata.var.index,columns=single_cell_adata.obs.index)

two_times_single_cell_gene = two_times_single_cell_gene.loc[gene_name]
two_times_single_cell_gene_T = pd.DataFrame(data=two_times_single_cell_gene.values.T,index=two_times_single_cell_gene.columns,columns=two_times_single_cell_gene.index)
two_times_single_cell_gene_T_value = two_times_single_cell_gene_T.values
two_times_single_cell_gene_T_value = torch.FloatTensor(two_times_single_cell_gene_T_value)
encoded_data,_ = model(two_times_single_cell_gene_T_value.cuda())
encoded_data = encoded_data.cpu()
encoded_data_np = encoded_data.detach().numpy()
mm = MinMaxScaler()
X_train_standard = mm.fit_transform(encoded_data_np)
sc_dim_reduce_standard = pd.DataFrame(data=X_train_standard.T,columns = single_cell_adata.obs.index)
sc_dim_reduce_standard.to_csv(autoencode_file)
print("AE is over")
