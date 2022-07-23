
import pandas as pd
import anndata as ad
import warnings
import os
warnings.filterwarnings("ignore")
dir_root = os.getcwd()
two_time_anndata_file = dir_root+"/data/intermediate_result/tow_time_adata.h5ad"
before_sc_file = dir_root+"/data/before_sc_data.csv"
after_sc_file = dir_root+"/data/after_sc_data.csv"
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
in_sc_csv = comb_sc_csv.loc[:,comb_mapping_result_csv.loc[:,"single_cell_name"].values]
in_sc_csv.columns = comb_mapping_result_csv.index
single_cell_adata = ad.AnnData(in_sc_csv.values.T)
single_cell_adata.obs = comb_mapping_result_csv
single_cell_adata.var = pd.DataFrame(index = in_sc_csv.index)
single_cell_adata.write(two_time_anndata_file)
print("TL_sample_get_adata is over")