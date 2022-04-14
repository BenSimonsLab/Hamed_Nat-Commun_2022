import numpy as np
import pandas as pd
import scanpy as sc

from collections import Counter


fileID = "normal_subset_neuronal_glial_cleaned"

adata = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")


# add root cells:

root_cellIDs = pd.read_csv("data/URD/"+fileID+"/root_cells.csv").cellID.values
adata.obs["URD_root"] = [x in root_cellIDs for x in adata.obs_names]


root_cellIDs = pd.read_csv("data/URD/"+fileID+"/root_cells_2.csv").cellID.values
adata.obs["URD_root_2"] = [x in root_cellIDs for x in adata.obs_names]


# add tip cells groups
tip_clusters = pd.read_csv("data/URD/"+fileID+"/tips_merged.csv", index_col=0)

tip_clusters_merged = pd.merge(adata.obs["sampleID"], tip_clusters, how='left', left_on="cellID", right_on="cellID", validate="one_to_one")

adata.obs["URD_tips"] = tip_clusters_merged.clusterID.values


###############
## Subsample ##
###############

adata_sub = sc.pp.subsample(adata, n_obs=10000, random_state=2533, copy=True)
adata_sub.write('data/forebrain_'+fileID+'_subsample_10k.h5ad')
# Export metadata
adata_sub.obs.to_csv("data/URD/metadata_normal_subset_neuronal_glial_cleaned_subsample_10k.csv.gz")

adata_sub = sc.pp.subsample(adata, n_obs=20000, random_state=2533, copy=True)
adata_sub.write('data/forebrain_'+fileID+'_subsample_20k.h5ad')
# Export metadata
adata_sub.obs.to_csv("data/URD/metadata_normal_subset_neuronal_glial_cleaned_subsample_20k.csv.gz")

adata_sub = sc.pp.subsample(adata, n_obs=30000, random_state=2533, copy=True)
adata_sub.write('data/forebrain_'+fileID+'_subsample_30k.h5ad')
# Export metadata
adata_sub.obs.to_csv("data/URD/metadata_normal_subset_neuronal_glial_cleaned_subsample_30k.csv.gz")

adata_sub = sc.pp.subsample(adata, n_obs=40000, random_state=2533, copy=True)
adata_sub.write('data/forebrain_'+fileID+'_subsample_40k.h5ad')
# Export metadata
adata_sub.obs.to_csv("data/URD/metadata_normal_subset_neuronal_glial_cleaned_subsample_40k.csv.gz")

adata_sub = sc.pp.subsample(adata, n_obs=50000, random_state=2533, copy=True)
adata_sub.write('data/forebrain_'+fileID+'_subsample_50k.h5ad')
# Export metadata
adata_sub.obs.to_csv("data/URD/metadata_normal_subset_neuronal_glial_cleaned_subsample_50k.csv.gz")


#######################
## Export raw counts ##
#######################

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


def get_raw_hvg_counts(fileID):
    adata_path = 'data/forebrain_'+fileID+'.h5ad'
    
    adata = sc.read(adata_path)
    
    sampleIDs = np.unique(adata.obs.sampleID.values)
    # cellIDs = np.concatenate([adata.obs_names.values])
    cellIDs = adata.obs_names.values
    # geneIDs = adata.var_names.values
    adata_raw = get_merged_adata(sampleIDs)
    adata_raw.obs.index = pd.Index(adata_raw.obs['cellID'])
    # adata_raw = adata_raw[cellIDs,geneIDs].copy()
    adata_raw = adata_raw[cellIDs,adata.var_names].copy()
    
    t = adata_raw.X.toarray()
    counts = pd.DataFrame(data=t, index=adata_raw.obs_names, columns=adata_raw.var_names)
    
    return(counts)


fileID = "normal_subset_embryonic_RG_cleaned"
counts = get_raw_counts(fileID)
counts.to_csv("data/URD/cmat_"+fileID+".csv") # NB: gz compression is very slow
# Export metadata
adata_path = 'data/forebrain_'+fileID+'.h5ad'
adata = sc.read(adata_path)
adata.obs.to_csv("data/URD/metadata_"+fileID+".csv.gz")

fileID = "normal_merged_RG_NSC_lineage"
counts = get_raw_counts(fileID)
counts.to_csv("data/URD/cmat_"+fileID+".csv") # NB: gz compression is very slow
# Export metadata
adata_path = 'data/forebrain_'+fileID+'.h5ad'
adata = sc.read(adata_path)
adata.obs.to_csv("data/URD/metadata_"+fileID+".csv.gz")



fileID = "normal_subset_neuronal_glial_cleaned_subsample_10k"
counts = get_raw_counts(fileID)
counts.to_csv("data/URD/cmat_"+fileID+".csv.gz")


fileID = "normal_subset_neuronal_glial_cleaned_subsample_10k"
counts = get_raw_hvg_counts(fileID)
counts.to_csv("data/URD/cmat_hvg_"+fileID+".csv.gz")

fileID = "normal_subset_neuronal_glial_cleaned_subsample_20k"
counts = get_raw_hvg_counts(fileID)
counts.to_csv("data/URD/cmat_hvg_"+fileID+".csv.gz")

fileID = "normal_subset_neuronal_glial_cleaned_subsample_30k"
counts = get_raw_hvg_counts(fileID)
counts.to_csv("data/URD/cmat_hvg_"+fileID+".csv.gz")

fileID = "normal_subset_neuronal_glial_cleaned_subsample_40k"
counts = get_raw_hvg_counts(fileID)
counts.to_csv("data/URD/cmat_hvg_"+fileID+".csv.gz")

fileID = "normal_subset_neuronal_glial_cleaned_subsample_50k"
counts = get_raw_hvg_counts(fileID)
counts.to_csv("data/URD/cmat_hvg_"+fileID+".csv.gz")


# manually add root cells to subsample
fileID = "normal_subset_neuronal_glial_cleaned_subsample_30k"
fileID = "normal_subset_neuronal_glial_cleaned_subsample_10k"
fileID = "normal_subset_neuronal_glial_cleaned_subsample_40k"
fileID = "normal_subset_neuronal_glial_cleaned_subsample_50k"

main_fileID = "normal_subset_neuronal_glial_cleaned"

adata_path = 'data/forebrain_'+fileID+'.h5ad'
adata = sc.read(adata_path)

root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_2.csv").cellID.values
adata.obs["URD_root_2"] = [x in root_cellIDs for x in adata.obs_names]

root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_3.csv").cellID.values
adata.obs["URD_root_3"] = [x in root_cellIDs for x in adata.obs_names]

root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_4.csv").cellID.values
adata.obs["URD_root_4"] = [x in root_cellIDs for x in adata.obs_names]

root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_5.csv").cellID.values
adata.obs["URD_root_5"] = [x in root_cellIDs for x in adata.obs_names]

root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_6.csv").cellID.values
adata.obs["URD_root_6"] = [x in root_cellIDs for x in adata.obs_names]


root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_em1.csv").cellID.values
adata.obs["URD_root_em1"] = [x in root_cellIDs for x in adata.obs_names]

root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_em2.csv").cellID.values
adata.obs["URD_root_em2"] = [x in root_cellIDs for x in adata.obs_names]

root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_em3.csv").cellID.values
adata.obs["URD_root_em3"] = [x in root_cellIDs for x in adata.obs_names]

root_cellIDs = pd.read_csv("data/URD/"+main_fileID+"/root_cells_em4.csv").cellID.values
adata.obs["URD_root_em4"] = [x in root_cellIDs for x in adata.obs_names]



# Export metadata
adata.obs.to_csv("data/URD/metadata_"+fileID+".csv.gz")
