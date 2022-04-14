import numpy as np
import pandas as pd
import scanpy as sc


# export raw counts imcluding all cells and genes without QC

def get_merged_adata(sampleIDs):
    list_adata = []
    
    for sampleID in sampleIDs:
        print(sampleID)
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


adata_raw = get_merged_adata(sampleIDs=['HAM12604', 'HAM12605',
                                        'HAM13754', 'HAM13755',
                                        'HAM13072', 'HAM13073',
                                        'HAM11199', 'HAM11200',
                                        'HAM10572', 'HAM10573',
                                        'HAM10364', 'HAM10365',
                                        'HAM10288', 'HAM10289',
                                        'HAM11302', 'HAM11303',
                                        'HAM10493', 'HAM10494',
                                        'HAM12214', 'HAM12215',
                                        'HAM13943', 'HAM13944'])


t = adata_raw.X.toarray()
counts = pd.DataFrame(data=t, index=adata_raw.obs_names, columns=adata_raw.var_names)


counts.to_csv('data/cmat_raw_forebrain_atlas.csv')
counts.to_csv('data/cmat_raw_forebrain_atlas.csv.gz', compression='gzip')



####################
## anndata object ##
####################

import numpy as np
import pandas as pd
import scanpy as sc

# add metadata to main h5ad object
adata_forebrain = sc.read_h5ad("data/forebrain_normal.h5ad")

metadata_full = pd.read_csv("data/metadata_forebrain_atlas.csv.gz", low_memory=False, index_col=0)
metadata_full = metadata_full.loc[:,[x not in adata_forebrain.obs.columns.values for x in metadata_full.columns]]


metadata_adata = adata_forebrain.obs
metadata_adata_full = metadata_adata.merge(metadata_full, on='cellID', how='left')
adata_forebrain.obs = metadata_adata_full

adata_forebrain.write('data/forebrain_atlas.h5ad')
