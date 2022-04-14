#########################
## Collect QC Metadata ##
#########################

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
                                             batch_key = 'sampleID',
                                             batch_categories = sampleIDs)
    
    adata_merged.obs['cellID'] = [cellID.replace('-1-', '-') for cellID in adata_merged.obs.index.tolist()]
    adata_merged.obs.index = pd.Index(adata_merged.obs['cellID'])
    
    return(adata_merged)



# import data
metadata_exp = pd.read_csv("data/metadata_experiment.tsv", sep='\t')
filter_set1 = metadata_exp["set"] == "set1"
filter_set = filter_set1
sampleIDs = list(metadata_exp[filter_set].sampleID.sort_values())

adata = get_merged_adata(sampleIDs)

# add metadata
metadata_exp = pd.read_csv("data/metadata_experiment.tsv", sep='\t')

meta_merged = pd.merge(adata.obs[['cellID', 'sampleID']], metadata_exp, on='sampleID', how='left')
scrublet_scores = pd.read_csv("data/scrublet-scores/scrublet_scores.csv.gz")
meta_merged = pd.merge(meta_merged, scrublet_scores, on='cellID', how='left')
meta_merged = meta_merged.set_index('cellID')

adata.obs = meta_merged

adata
# run to calculate n_genes
sc.pp.filter_cells(adata, min_genes=1)
# sc.pp.filter_genes(adata, min_cells=1)

# mt_counts & n_counts
mito_genes = adata.var_names.str.startswith('mt-')
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

adata


adata.obs.to_csv("data/metadata_QC_normal.csv.gz", compression='gzip')