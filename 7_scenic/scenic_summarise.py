import numpy as np
import pandas as pd
import scanpy as sc

from pyscenic.rss import regulon_specificity_scores


def scenic_enrichment(fileID, leiden_res):
    file_path = "data/forebrain_"+fileID+".h5ad"
    adata = sc.read_h5ad(file_path)
    auc_mtx = pd.read_csv('data/scenic/auc_mtx_'+fileID+'.csv.gz', index_col=0)
    
    rss_cellType = regulon_specificity_scores(auc_mtx, adata.obs[leiden_res][auc_mtx.index.values].values)
    
    rss_cellType.to_csv("data/scenic/rss_"+fileID+"_"+leiden_res+".csv.gz")
    
    return


# scenic_enrichment('normal_subset_embryonic_RG_cleaned', 'leiden_0_33')
# scenic_enrichment('normal_subset_embryonic_RG_cleaned_wo_gliogenic', 'annot_superset_num')
# scenic_enrichment('normal_merged_RG_NSC_lineage', "annot_merged_coarse")
# scenic_enrichment('normal_merged_RG_NSC_lineage', "annot_merged_fine")
# scenic_enrichment('normal_merged_RG_NSC_lineage', "annot_merged_coarse_num")
# scenic_enrichment('normal_merged_RG_NSC_lineage', "annot_merged_fine_num")


scenic_enrichment('normal_subset_ependymal_cleaned', 'leiden_0_15')
scenic_enrichment('normal_subset_OPCs', 'leiden_0_1')

scenic_enrichment('normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned', 'leiden_0_13')
