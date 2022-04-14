import scvelo as scv
import scanpy as sc


def run_velocyto(fileID):
    adata = scv.read('data/forebrain_'+fileID+'.h5ad')
    
    scv.settings.figdir = 'output/'+fileID+'/'
    
    ldata = scv.read('data/set1_normal_splicing.h5ad')
    ldata.var_names_make_unique()
    ldata = ldata[adata.obs_names.values, :].copy()
    
    adata_velo = scv.utils.merge(adata, ldata)
    
    scv.pp.filter_and_normalize(adata_velo)
    scv.pp.moments(adata_velo)
    
    scv.tl.velocity(adata_velo, mode='stochastic')
    scv.tl.velocity_graph(adata_velo)
    scv.tl.velocity_pseudotime(adata_velo)
    
    adata_velo.obs['velocity_pseudotime'].to_csv('output/'+fileID+'/velocity_pseudotime_'+fileID+'.csv.gz')
    scv.pl.velocity_embedding_stream(adata_velo, color='velocity_pseudotime', basis='umap', save='stochastic_'+fileID+'_velocity_pseudotime.png', dpi=300)
    # get list_leiden_res
    list_leiden_res = [x for x in adata.obs_keys() if x.startswith('leiden_')]
    
    for leiden_res in list_leiden_res:
        scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
        
    return
        

# run the stochastic scvelo RNA velocity inference for different subsets and save output

# list_fileIDs = ["normal", "normal_prenatal_RGCs", "normal_RGCs_cleaned", "normal_subset_ependymal_cleaned"]
# 
# for fileID in list_fileIDs:
#     run_velocyto(fileID)



# 2020-11-04
run_velocyto("normal_merged_subset_neuronal_glial_subset_RG_OPC_NSC_lineage_cleaned")