############
## SCENIC ##
############


# Export from adata
# get raw untransformed data (tranformation has impact on results)

import pandas as pd
import numpy as np
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
    adata_merged.obs.index = pd.Index(adata_merged.obs['cellID'])
    
    return(adata_merged)


def export_cmat_csv(fileID):
    adata_path = 'data/forebrain_'+fileID+'.h5ad'
    
    adata = sc.read(adata_path)
    
    sampleIDs = np.unique(adata.obs.sampleID.values)
    cellIDs = np.concatenate([adata.obs_names.values])
    
    adata_raw = get_merged_adata(sampleIDs)
    adata_raw = adata_raw[cellIDs,:].copy()
    sc.pp.filter_genes(adata_raw, min_cells=3)
    
    t = adata_raw.X.toarray()
    cmat = pd.DataFrame(data=t, index=adata_raw.obs_names, columns=adata_raw.var_names)
    cmat.to_csv('data/scenic/cmat_'+fileID+'.csv.gz')
    return



# untransformed, raw data

fileID = "normal_subset_embryonic_RG_cleaned"
fileID = "normal_subset_ependymal_cleaned"
fileID = "normal_subset_OPCs"
fileID = "normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned"
fileID = "normal_subset_embryonic_RG_cleaned_wo_gliogenic"
fileID = "normal_merged_RG_NSC_lineage"

fileID = "normal"


export_cmat_csv(fileID)


# CLI

# docker aertslab/pyscenic:0.10.0
# docker pull aertslab/pyscenic:0.10.2


fileID="normal_subset_embryonic_RG_cleaned"
fileID="normal_subset_ependymal_cleaned"
fileID="normal_subset_OPCs"
fileID="normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned"
fileID="normal_subset_embryonic_RG_cleaned_wo_gliogenic"
fileID="normal_merged_RG_NSC_lineage"


docker run -it --rm \
    -v /home/ubuntu/data:/scenicdata \
    aertslab/pyscenic:0.10.2 pyscenic grn \
        --num_workers 10 \
        -o /scenicdata/brain-development/data/scenic/cmat_adjacencies_${fileID}.csv \
        /scenicdata/brain-development/data/scenic/cmat_${fileID}.csv.gz \
        /scenicdata/references/scenic/mm_mgi_tfs.txt


docker run -it --rm \
    -v /home/ubuntu/data:/scenicdata \
    aertslab/pyscenic:0.10.2 pyscenic ctx \
        /scenicdata/brain-development/data/scenic/cmat_adjacencies_${fileID}.csv \
        /scenicdata/references/scenic/mm9-500bp-upstream-10species.mc9nr.feather \
        /scenicdata/references/scenic/mm9-tss-centered-5kb-10species.mc9nr.feather \
        /scenicdata/references/scenic/mm9-tss-centered-10kb-10species.mc9nr.feather \
        --annotations_fname /scenicdata/references/scenic/motifs-v9-nr.mgi-m0.001-o0.0.tbl \
        --expression_mtx_fname /scenicdata/brain-development/data/scenic/cmat_${fileID}.csv.gz \
        --mode "dask_multiprocessing" \
        --output /scenicdata/brain-development/data/scenic/regulons_${fileID}.csv.gz \
        --num_workers 10


docker run -it --rm \
    -v /home/ubuntu/data:/scenicdata \
    aertslab/pyscenic:0.10.2 pyscenic aucell \
        /scenicdata/brain-development/data/scenic/cmat_${fileID}.csv.gz \
        /scenicdata/brain-development/data/scenic/regulons_${fileID}.csv.gz \
        -o /scenicdata/brain-development/data/scenic/auc_mtx_${fileID}.csv.gz \
        --num_workers 10
