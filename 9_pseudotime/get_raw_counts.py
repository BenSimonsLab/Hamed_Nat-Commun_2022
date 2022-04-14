Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3.6")

library(reticulate)
# use_python("/usr/bin/python3.6")
# py_config()


repl_python()


import numpy as np
import pandas as pd
import scanpy as sc


def get_merged_adata(sampleIDs):
    list_adata = []
    
    for sampleID in sampleIDs:
        # print(sampleID)
        adata_temp = sc.read_10x_h5('data/samples/' + sampleID + '/outs/filtered_feature_bc_matrix.h5')
        adata_temp.var_names_make_unique()
        list_adata.append(adata_temp)
        
    adata_merged = list_adata[0].concatenate(list_adata[1:],
                                             join='outer',
                                             fill_value=0,
                                             batch_key = 'sampleID',
                                             batch_categories = sampleIDs)
    
    adata_merged.obs['cellID'] = [cellID.replace('-1-', '-') for cellID in adata_merged.obs.index.tolist()]
    # adata_merged.obs.index = pd.Index(adata_merged.obs['cellID'])
    
    return(adata_merged)


def get_raw_hvg_counts(fileID):
    
    # file_path_main = "data/forebrain_"+fileID+".h5ad"
    # adata_processed = sc.read_h5ad(file_path_main)
    # 
    # sampleIDs = adata_processed.obs.sampleID.unique()
    # adata = get_merged_adata(sampleIDs)
    # adata.obs.index = pd.Index(adata.obs['cellID'])
    # 
    # adata = adata[adata_processed.obs_names.values, :].copy()
    # 
    # counts = adata.X.to_dense()
    
    adata_path = 'data/forebrain_'+fileID+'.h5ad'
    
    adata = sc.read(adata_path)
    
    sampleIDs = np.unique(adata.obs.sampleID.values)
    # cellIDs = np.concatenate([adata.obs_names.values])
    cellIDs = adata.obs_names.values
    # geneIDs = adata.var_names.values
    adata_raw = get_merged_adata(sampleIDs)
    adata_raw.obs.index = pd.Index(adata_raw.obs['cellID'])
    # adata_raw = adata_raw[cellIDs,geneIDs].copy()
    adata_raw = adata_raw[cellIDs,:].copy()
    sc.pp.filter_genes(adata_raw, min_cells=3)
    
    t = adata_raw.X.toarray()
    counts = pd.DataFrame(data=t, index=adata_raw.obs_names, columns=adata_raw.var_names)
    
    return(counts)





# fileID = "normal_merged_subset_neuronal_glial_subset_RG_OPC_NSC_lineage_cleaned"
# counts = get_raw_hvg_counts(fileID)


fileID = "normal_merged_subset_neuronal_glial"
counts = get_raw_hvg_counts(fileID)
counts.to_csv("data/temp/cmat_normal_merged_subset_neuronal_glial.csv.gz")



adata_path = 'data/forebrain_'+fileID+'.h5ad'
adata = sc.read(adata_path)
adata.obs.to_csv("data/temp/metadata_normal_merged_subset_neuronal_glial.csv.gz")

exit
# in R
counts = t(py$counts)