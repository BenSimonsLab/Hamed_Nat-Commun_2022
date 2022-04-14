#!/usr/bin/env python3

#####################
## Scanpy Pipeline ##
#####################

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


###############################
## Integration - Zeisel 2018 ##
###############################

import numpy as np
import pandas as pd
import scanpy as sc


list_adata = []
# 
# for subsetID in subsetIDs:
#     # print(sampleID)
#     file_path_main = "data/brain-regeneration_"+subsetID+".h5ad"
#     adata_temp = sc.read_h5ad(file_path_main)
#     list_adata.append(adata_temp)


adata_normal = sc.read_h5ad("data/forebrain_normal_subset_neuronal_glial_cleaned.h5ad")
adata_zeisel = adata_zeisel = sc.read_h5ad('data/Zeisel_2018/Zeisel_2018-raw.h5ad')
filter_neuronal_glial = [clusterID in ['Glia', 'Neurons'] for clusterID in adata_zeisel.obs["TaxonomyRank1"]]
adata_zeisel = adata_zeisel[filter_neuronal_glial, :].copy()

################################################################################
# remove sex specific genes (Zeisel 2018)
sex_genes = ['Xist', 'Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d']
filter_sex_genes = [x not in sex_genes for x in adata_lamanno.var_names]
adata_zeisel = adata_zeisel[:, filter_sex_genes]

# not really needed, but generates n_genes slot
sc.pp.filter_cells(adata_zeisel, min_genes=200)
sc.pp.filter_genes(adata_zeisel, min_cells=3)

mito_genes = adata_zeisel.var_names.str.startswith('mt-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata_zeisel.obs['percent_mito'] = np.sum(
adata_zeisel[:, mito_genes].X, axis=1).A1 / np.sum(adata_zeisel.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata_zeisel
adata_zeisel.obs['n_counts'] = adata_zeisel.X.sum(axis=1).A1

# Normalise and log1p        
sc.pp.normalize_per_cell(adata_zeisel, counts_per_cell_after=1e4)
sc.pp.log1p(adata_zeisel)
adata_zeisel.raw = adata_zeisel
################################################################################


list_adata = [adata_normal, adata_zeisel]


adata_merged = list_adata[0].concatenate(list_adata[1:],
                                         batch_key = 'subsetID',
                                         batch_categories = ['normal', 'Zeisel 2018'],
                                         index_unique=None)


adata = sc.AnnData(adata_merged.raw.X, obs=adata_merged.obs, var=adata_merged.raw.var)

adata.raw = adata

adata_normal = 0
adata_kalamakis = 0
adata_merged = 0

run_adata_processing(adata, "normal_subset_neuronal_glial_cleaned_Zeisel_2018")