#!/usr/bin/env python3

#####################
## Scanpy Pipeline ##
#####################

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad



def run_adata_processing(adata, fileID):
    # expects cleaned and scaled adata as input
    
    # HVGs
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var['highly_variable']].copy()
    # regress counts
    sc.pp.regress_out(adata, ['n_counts'])
    sc.pp.scale(adata, max_value=10)
    # PCA
    sc.tl.pca(adata)
    # neighbourhood graph
    sc.pp.neighbors(adata, n_neighbors=15)
    
    # UMAP
    sc.tl.umap(adata)
    
    # Leiden clustering
    list_leiden_res = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    for leiden_res in list_leiden_res:
        sc.tl.leiden(adata, resolution=leiden_res, key_added='leiden_'+str(leiden_res).replace(".", "_"))
        adata.obs['annot_leiden_'+str(leiden_res).replace(".", "_")] = adata.obs['leiden_'+str(leiden_res).replace(".", "_")]
        
    # save results
    adata.write('data/Richards_2020/GBM_Richards_2020_'+fileID+'.h5ad')
    return


def process_raw_Richards_2020():
    cmat_richards = pd.read_csv('data/Richards_2020/cmat_GBM_3_patients_tumour.csv', index_col=0)
    metadata_richards = pd.read_csv('data/Richards_2020/metadata_GBM_3_patients_tumour.csv', index_col=0)
    
    cmat_richards = cmat_richards.loc[:, metadata_richards.index.values]
    
    adata_Richards_2020 = ad.AnnData(cmat_richards.values.T, obs=metadata_richards, var=pd.DataFrame(index=cmat_richards.index.values))
    adata = adata_Richards_2020
    
    sc.pp.filter_cells(adata, min_genes=1) # cells are already filtered
    sc.pp.filter_genes(adata, min_cells=3)
    
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    
    # Normalise and log1p
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    
    run_adata_processing(adata, "GBM_3_patients_tumour")
    return



def process_mm_raw_Richards_2020():
    cmat_richards = pd.read_csv('data/Richards_2020/cmat_mm_GBM_3_patients_tumour.csv', index_col=0)
    metadata_richards = pd.read_csv('data/Richards_2020/metadata_GBM_3_patients_tumour.csv', index_col=0)
    
    cmat_richards = cmat_richards.loc[:, metadata_richards.index.values]
    
    adata_Richards_2020 = ad.AnnData(cmat_richards.values.T, obs=metadata_richards, var=pd.DataFrame(index=cmat_richards.index.values))
    adata = adata_Richards_2020
    
    sc.pp.filter_cells(adata, min_genes=1) # cells are already filtered
    sc.pp.filter_genes(adata, min_cells=3)
    
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    
    # Normalise and log1p
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    
    run_adata_processing(adata, "mm_GBM_3_patients_tumour")
    return


######################################
## Richards 2020 - Process raw data ##
######################################

# process_raw_Richards_2020()
# process_mm_raw_Richards_2020()

########################################
## Richards 2020 - Split into samples ##
########################################


def run_pipe_subset_sampleID(main_fileID, subsetID, sampleIDs_selected=None, sampleIDs_removed=None):
        
        fileID = subsetID
        file_path_main = "data/Richards_2020/GBM_Richards_2020_"+main_fileID+".h5ad"
        
        adata = sc.read_h5ad(file_path_main)
        
        if sampleIDs_selected and not sampleIDs_removed:
            filter_sub = [x in sampleIDs_selected for x in adata.obs["sampleID"]]
        elif sampleIDs_removed and not sampleIDs_selected:
            filter_sub = [x not in sampleIDs_removed for x in adata.obs["sampleID"]]
        else:
            raise Exception('Supply either sampleIDs_selected or sampleIDs_removed')
        
        adata_tmp = adata[filter_sub,:]
        adata_sub = sc.AnnData(adata_tmp.raw.X, obs=adata_tmp.obs, var=adata_tmp.raw.var)
        
        adata_sub.raw = adata_sub
        
        # remove leiden and annot
        filter_meta_cols = [not x.startswith(('annot_', 'leiden_')) for x in adata_sub.obs.columns]
        adata_sub.obs = adata_sub.obs.loc[:, filter_meta_cols]
        
        run_adata_processing(adata_sub, fileID)
        
        return



adata_Richards_2020 = sc.read_h5ad("data/Richards_2020/GBM_Richards_2020_"+"GBM_3_patients_tumour"+".h5ad")
adata_Richards_2020.obs.sampleID

sampleID_list = list(set(adata_Richards_2020.obs.sampleID))

sampleID = 'G910_S1'

for sampleID in sampleID_list:
    run_pipe_subset_sampleID("GBM_3_patients_tumour", sampleID, sampleIDs_selected=[sampleID])

# by patient
run_pipe_subset_sampleID("GBM_3_patients_tumour", "G910_full", sampleIDs_selected=['G910_S1', 'G910_S2', 'G910_S3', 'G910_S4', 'G910_S5'])
run_pipe_subset_sampleID("GBM_3_patients_tumour", "G945_full", sampleIDs_selected=['G945_S1', 'G945_S2', 'G945_S3'])
run_pipe_subset_sampleID("GBM_3_patients_tumour", "G1003_full", sampleIDs_selected=['G1003_S1', 'G1003_S2', 'G1003_S3', 'G1003_S4'])


run_pipe_subset_sampleID("mm_GBM_3_patients_tumour", "mm_G910_full", sampleIDs_selected=['G910_S1', 'G910_S2', 'G910_S3', 'G910_S4', 'G910_S5'])
run_pipe_subset_sampleID("mm_GBM_3_patients_tumour", "mm_G945_full", sampleIDs_selected=['G945_S1', 'G945_S2', 'G945_S3'])
run_pipe_subset_sampleID("mm_GBM_3_patients_tumour", "mm_G1003_full", sampleIDs_selected=['G1003_S1', 'G1003_S2', 'G1003_S3', 'G1003_S4'])



############################################
## Integrate with neuronal & glial subset ##
############################################

import numpy as np
import pandas as pd
import scanpy as sc


def run_adata_processing(adata, fileID):
    # expects cleaned and scaled adata as input
    
    # setup directories
    import os
    if not os.path.exists("output/"+fileID+"/"):
        os.makedirs("output/"+fileID+"/")
        
    sc.settings.figdir = 'output/'+fileID+'/'
    
    # HVGs
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var['highly_variable']].copy()
    # regress counts
    sc.pp.regress_out(adata, ['n_counts'])
    sc.pp.scale(adata, max_value=10)
    # PCA
    sc.tl.pca(adata)
    # neighbourhood graph
    sc.pp.neighbors(adata, n_neighbors=15)
    
    # UMAP
    sc.tl.umap(adata)
    
    # Leiden clustering
    list_leiden_res = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    
    for leiden_res in list_leiden_res:
        sc.tl.leiden(adata, resolution=leiden_res, key_added='leiden_'+str(leiden_res).replace(".", "_"))
        adata.obs['annot_leiden_'+str(leiden_res).replace(".", "_")] = adata.obs['leiden_'+str(leiden_res).replace(".", "_")]
        
    # save results
    adata.write('data/forebrain_'+fileID+'.h5ad')
    return


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


def run_pipe(sampleIDs, fileID):
    adata = get_merged_adata(sampleIDs)
    
    # import metadata    
    metadata_exp = pd.read_csv("data/metadata_experiment.tsv", sep='\t')
    
    meta_merged = pd.merge(adata.obs[['cellID', 'sampleID']], metadata_exp, on='sampleID', how='left')
    scrublet_scores = pd.read_csv("data/scrublet-scores/scrublet_scores.csv.gz")
    meta_merged = pd.merge(meta_merged, scrublet_scores, on='cellID', how='left')
    meta_merged = meta_merged.set_index('cellID')
    
    adata.obs = meta_merged
    
    # filter cells based on QC
    metadata_qc_normal = pd.read_csv('data/metadata_QC_normal_pass.csv.gz')
    metadata_qc = pd.concat([metadata_qc_normal])
    
    filter_QC = [x in metadata_qc.cellID.values for x in adata.obs_names]
    adata = adata[filter_QC, :].copy()
    
    # remove sex specific genes (Zeisel 2018)
    sex_genes = ['Xist', 'Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d']
    filter_sex_genes = [x not in sex_genes for x in adata.var_names]
    adata = adata[:, filter_sex_genes]
    
    # not really needed, but generates n_genes slot
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    mito_genes = adata.var_names.str.startswith('mt-')
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    # add the total counts per cell as observations-annotation to adata
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    
    # Normalise and log1p        
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    
    run_adata_processing(adata, fileID)
    return


# integration by patient

fileID_sample = "GBM_Richards_2020_mm_G910_full"
fileID_sample = "GBM_Richards_2020_mm_G945_full"
fileID_sample = "GBM_Richards_2020_mm_G1003_full"

list_adata = []

adata_normal = sc.read_h5ad("data/forebrain_normal_subset_neuronal_glial_cleaned.h5ad")
adata_richards = sc.read_h5ad("data/Richards_2020/"+fileID_sample+".h5ad")



list_adata = [adata_normal, adata_richards]


adata_merged = list_adata[0].concatenate(list_adata[1:],
                                         batch_key = 'subsetID',
                                         batch_categories = ['normal', 'Richards 2020'],
                                         index_unique=None)


adata = sc.AnnData(adata_merged.raw.X, obs=adata_merged.obs, var=adata_merged.raw.var)

adata.raw = adata

adata_normal = 0
adata_richards = 0
adata_merged = 0

run_adata_processing(adata, "normal_subset_neuronal_glial_cleaned_"+fileID_sample)




# by sample


# by patient



#############
## Harmony ##
#############

import numpy as np
import pandas as pd
import scanpy as sc


def run_harmony_integration(adata):
    
    import rpy2.robjects as robjects 
    
    num_pcs = 20
    pca = adata.obsm['X_pca'][:, :num_pcs]
    batch = adata.obs['subsetID']
    
    robjects.globalenv['pca'] = pca
    robjects.globalenv['batch'] = batch
    
    
    #########################
    robjects.r('''
    library(harmony)
    library(magrittr)
    
    # Harmony determinism!
    set.seed(1)
    
    hem = HarmonyMatrix(pca, batch, theta=4, verbose=FALSE, do_pca=FALSE)
    hem = data.frame(hem)
    ''')
    
    ##################
    
    hem = robjects.globalenv['hem']
    
    
    adata.obsm['X_harmony'] = hem.values
    sc.pp.neighbors(adata, use_rep='X_harmony')
    
    
    sc.tl.umap(adata)
    return(adata)
    


fileID_sample = "GBM_Richards_2020_G910_full"
fileID_sample = "GBM_Richards_2020_G945_full"
fileID_sample = "GBM_Richards_2020_G1003_full"



adata = sc.read_h5ad("data/forebrain_normal_subset_neuronal_glial_cleaned_"+fileID_sample+".h5ad")
adata_integration = run_harmony_integration(adata)



# sc.pl.umap(adata, color=['subsetID'])

# robjects.globalenv['sql_scores'] = sql_scores
# robjects.r("results = fetch(dbSendQuery(mydb, sql_scores))")



# def run_pipe(sampleIDs, fileID):
#     adata = get_merged_adata(sampleIDs)
# 
#     # import metadata    
#     metadata_exp = pd.read_csv("data/metadata_experiment.tsv", sep='\t')
# 
#     meta_merged = pd.merge(adata.obs[['cellID', 'sampleID']], metadata_exp, on='sampleID', how='left')
#     scrublet_scores = pd.read_csv("data/scrublet-scores/scrublet_scores.csv.gz")
#     meta_merged = pd.merge(meta_merged, scrublet_scores, on='cellID', how='left')
#     meta_merged = meta_merged.set_index('cellID')
# 
#     adata.obs = meta_merged
# 
#     # filter cells based on QC
#     metadata_qc_normal = pd.read_csv('data/metadata_QC_normal_pass.csv.gz')
#     metadata_qc = pd.concat([metadata_qc_normal])
# 
#     filter_QC = [x in metadata_qc.cellID.values for x in adata.obs_names]
#     adata = adata[filter_QC, :].copy()
# 
#     # remove sex specific genes (Zeisel 2018)
#     sex_genes = ['Xist', 'Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d']
#     filter_sex_genes = [x not in sex_genes for x in adata.var_names]
#     adata = adata[:, filter_sex_genes]
# 
#     # not really needed, but generates n_genes slot
#     sc.pp.filter_cells(adata, min_genes=200)
#     sc.pp.filter_genes(adata, min_cells=3)
# 
#     mito_genes = adata.var_names.str.startswith('mt-')
#     # for each cell compute fraction of counts in mito genes vs. all genes
#     # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
#     adata.obs['percent_mito'] = np.sum(
#         adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
#     # add the total counts per cell as observations-annotation to adata
#     adata.obs['n_counts'] = adata.X.sum(axis=1).A1
# 
#     # Normalise and log1p        
#     sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
#     sc.pp.log1p(adata)
#     adata.raw = adata
# 
#     run_adata_processing(adata, fileID)
#     return