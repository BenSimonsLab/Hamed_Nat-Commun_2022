
import numpy as np
import pandas as pd
import scanpy as sc



def add_leiden_res(fileID, leiden_res_float):
    # expects cleaned and scaled adata as input
    adata = sc.read('data/forebrain_'+fileID+'.h5ad')
    
    leiden_res = 'leiden_'+str(leiden_res_float).replace(".", "_")
    
    # setup directories
    sc.settings.figdir = 'output/'+fileID+'/'
    
    sc.tl.leiden(adata, resolution=leiden_res_float, key_added=leiden_res)
    
    adata.obs['annot_'+leiden_res] = adata.obs[leiden_res]
    
    # save results
    adata.write('data/forebrain_'+fileID+'.h5ad')
    return


def subcluster_leiden_res(fileID, leiden_res, list_clusterIDs, sub_leiden_res_float):
    # expects cleaned and scaled adata as input
    adata = sc.read('data/forebrain_'+fileID+'.h5ad')
    
    leiden_key = leiden_res+'_cl_'+'_'.join(list_clusterIDs)+'_'+'leiden_'+str(sub_leiden_res_float).replace(".", "_")
    sc.tl.leiden(adata, restrict_to=(leiden_res, list_clusterIDs), resolution=sub_leiden_res_float, key_added=leiden_key)
    
    adata.obs['annot_'+leiden_key] = adata.obs[leiden_key]
    
    # save results
    adata.write('data/forebrain_'+fileID+'.h5ad')
    return



add_leiden_res("normal", 0.41)
add_leiden_res("normal", 0.42)
add_leiden_res("normal", 0.43)
add_leiden_res("normal", 0.44)
add_leiden_res("normal", 0.45)
add_leiden_res("normal", 0.46)
add_leiden_res("normal", 0.47)
add_leiden_res("normal", 0.48)
add_leiden_res("normal", 0.49)
add_leiden_res("normal", 0.5)
add_leiden_res("normal", 0.51)
add_leiden_res("normal", 0.52)
add_leiden_res("normal", 0.53)
add_leiden_res("normal", 0.54)
add_leiden_res("normal", 0.55)
add_leiden_res("normal", 0.56)
add_leiden_res("normal", 0.57)
add_leiden_res("normal", 0.58)
add_leiden_res("normal", 0.59)
add_leiden_res("normal", 0.6)



subcluster_leiden_res('normal_subset_neuronal_glial_cleaned', 'leiden_0_4', ['2'], 0.05)
subcluster_leiden_res('normal_subset_neuronal_glial_cleaned', 'leiden_0_4', ['2'], 0.075)
subcluster_leiden_res('normal_subset_neuronal_glial_cleaned', 'leiden_0_4', ['2'], 0.1)
subcluster_leiden_res('normal_subset_neuronal_glial_cleaned', 'leiden_0_4', ['2'], 0.15)
subcluster_leiden_res('normal_subset_neuronal_glial_cleaned', 'leiden_0_4', ['2'], 0.2)
subcluster_leiden_res('normal_subset_neuronal_glial_cleaned', 'leiden_0_4', ['2'], 0.25)



add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.31)
add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.32)
add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.33)
add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.34)
add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.35)
add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.36)
add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.37)
add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.38)
add_leiden_res("normal_subset_embryonic_RG_cleaned", 0.39)


add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.31)
add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.32)
add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.33)
add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.34)
add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.35)
add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.36)
add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.37)
add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.38)
add_leiden_res("normal_subset_juvenile_RG_TAPs", 0.39)



add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.11)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.12)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.13)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.14)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.15)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.16)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.17)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.18)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_ependymal_cleaned", 0.19)

add_leiden_res("normal_subset_ependymal_cleaned", 0.11)
add_leiden_res("normal_subset_ependymal_cleaned", 0.12)
add_leiden_res("normal_subset_ependymal_cleaned", 0.13)
add_leiden_res("normal_subset_ependymal_cleaned", 0.14)
add_leiden_res("normal_subset_ependymal_cleaned", 0.15)
add_leiden_res("normal_subset_ependymal_cleaned", 0.16)
add_leiden_res("normal_subset_ependymal_cleaned", 0.17)
add_leiden_res("normal_subset_ependymal_cleaned", 0.18)
add_leiden_res("normal_subset_ependymal_cleaned", 0.19)


add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.11)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.12)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.13)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.14)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.15)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.16)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.17)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.18)
add_leiden_res("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", 0.19)