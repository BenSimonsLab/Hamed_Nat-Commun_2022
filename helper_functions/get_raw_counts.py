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


def get_raw_counts(fileID):
    
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