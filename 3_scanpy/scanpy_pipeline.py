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


######################
## From raw samples ##
######################

# run on 2020-10-31
# sc.logging.print_header()
## scanpy==1.6.0 anndata==0.7.4 umap==0.3.10 numpy==1.19.2 scipy==1.5.2 pandas==1.1.3 scikit-learn==0.23.2 statsmodels==0.12.1 python-igraph==0.7.1 louvain==0.6.1 leidenalg==0.7.0


# normal
run_pipe(sampleIDs=['HAM12604', 'HAM12605',
                    'HAM13754', 'HAM13755',
                    'HAM13072', 'HAM13073',
                    'HAM11199', 'HAM11200',
                    'HAM10572', 'HAM10573',
                    'HAM10364', 'HAM10365',
                    'HAM10288', 'HAM10289',
                    'HAM11302', 'HAM11303',
                    'HAM10493', 'HAM10494',
                    'HAM12214', 'HAM12215',
                    'HAM13943', 'HAM13944'], fileID='normal')


run_pipe(sampleIDs=['HAM12604', 'HAM12605'], fileID='normal_E12_5')
run_pipe(sampleIDs=['HAM13754', 'HAM13755'], fileID='normal_E14_5')
run_pipe(sampleIDs=['HAM13072', 'HAM13073'], fileID='normal_E16_5')
run_pipe(sampleIDs=['HAM11199', 'HAM11200'], fileID='normal_P0')
run_pipe(sampleIDs=['HAM10572', 'HAM10573'], fileID='normal_P2')
run_pipe(sampleIDs=['HAM10364', 'HAM10365'], fileID='normal_P7')
run_pipe(sampleIDs=['HAM10288', 'HAM10289'], fileID='normal_P13')
run_pipe(sampleIDs=['HAM11302', 'HAM11303'], fileID='normal_P19')
run_pipe(sampleIDs=['HAM10493', 'HAM10494'], fileID='normal_P39')
run_pipe(sampleIDs=['HAM12214', 'HAM12215'], fileID='normal_P111')
run_pipe(sampleIDs=['HAM13943', 'HAM13944'], fileID='normal_P365')







################################################################################
################################################################################

#############
## Subsets ##
#############


def run_pipe_subset(main_fileID, subsetID, subset_leidenID, clusterIDs_selected=None, clusterIDs_removed=None):
        
        fileID = main_fileID+'_'+subsetID
        file_path_main = "data/forebrain_"+main_fileID+".h5ad"
        
        adata = sc.read_h5ad(file_path_main)
        
        if clusterIDs_selected and not clusterIDs_removed:
            filter_sub = [x in clusterIDs_selected for x in adata.obs[subset_leidenID]]
        elif clusterIDs_removed and not clusterIDs_selected:
            filter_sub = [x not in clusterIDs_removed for x in adata.obs[subset_leidenID]]
        else:
            raise Exception('Supply either clusterIDs_selected or clusterIDs_removed')
        
        adata_tmp = adata[filter_sub,:]
        adata_sub = sc.AnnData(adata_tmp.raw.X, obs=adata_tmp.obs, var=adata_tmp.raw.var)
        
        adata_sub.raw = adata_sub
        
        # remove leiden and annot
        filter_meta_cols = [not x.startswith(('annot_', 'leiden_')) for x in adata_sub.obs.columns]
        adata_sub.obs = adata_sub.obs.loc[:, filter_meta_cols]
        
        run_adata_processing(adata_sub, fileID)
        
        return



# 2020-11-01
# # ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
# subset_neuronal_glial
# run_pipe_subset(main_fileID, subsetID, subset_leidenID, clusterIDs_selected=None, clusterIDs_removed=None):

# >>> sc.logging.print_header()
# scanpy==1.6.0 anndata==0.7.4 umap==0.3.10 numpy==1.19.2 scipy==1.5.2 pandas==1.1.3 scikit-learn==0.23.2 statsmodels==0.12.1 python-igraph==0.7.1 louvain==0.6.1 leidenalg==0.7.0

# E12_5, leiden_0_1, remove 11, 8, 4, 13
run_pipe_subset("normal_E12_5", "subset_neuronal_glial", "leiden_0_1",
                clusterIDs_removed=['4', '8', '11', '13'])
# E14_5, leiden_0_1, remove 5, 6, 7, 8
run_pipe_subset("normal_E14_5", "subset_neuronal_glial", "leiden_0_1",
                clusterIDs_removed=['5', '6', '7', '8'])
# E16_5, leiden_0_7, remove 20, 16, 12, 15, 26, 22, 25, 27, 28, 19, 18
run_pipe_subset("normal_E16_5", "subset_neuronal_glial", "leiden_0_7",
                clusterIDs_removed=['12', '15', '16', '18', '19', '20', '22', '25', '26', '27', '28'])
# P0, leiden_0_7, remove 23, 18, 25, 19, 11, 20, 26, 21 (might lose a few relevant cells), 22
run_pipe_subset("normal_P0", "subset_neuronal_glial", "leiden_0_7",
                clusterIDs_removed=['11', '18', '19', '20', '21', '22', '23', '25', '26'])
# P2, leiden_0_3, remove 10, 15, 16, 11, 6, 14, 13
run_pipe_subset("normal_P2", "subset_neuronal_glial", "leiden_0_3",
                clusterIDs_removed=['6', '10', '11', '13', '14', '15', '16'])
# P7, leiden_0_8, remove 20, 21, 14, 13, 11, 19, 0, 6, 22, 15, 18
run_pipe_subset("normal_P7", "subset_neuronal_glial", "leiden_0_8",
                clusterIDs_removed=['0', '6', '11', '13', '14', '15', '18', '19', '20', '21', '22'])
# P13, leiden_0_4, remove 16, 12, 2, 9, 0, 15, 13, 7
run_pipe_subset("normal_P13", "subset_neuronal_glial", "leiden_0_4",
                clusterIDs_removed=['0', '2', '7', '9', '12', '13', '15', '16'])
# P19, leiden_0_7, remove 23, 19, 26, 25, 18, 22, 1, 2, 7, 11, 21, 4, 24, 16
run_pipe_subset("normal_P19", "subset_neuronal_glial", "leiden_0_7",
                clusterIDs_removed=['1', '2', '4', '7', '11', '16', '18', '19', '21', '22', '23', '24', '25', '26'])
# P39, leiden_0_1, remove 6, 0, 8, 7
run_pipe_subset("normal_P39", "subset_neuronal_glial", "leiden_0_1",
                clusterIDs_removed=['0', '6', '7', '8'])
# P111, leiden_0_2, remove 6, 10, 11, 7, 8, 14, 0, 3
run_pipe_subset("normal_P111", "subset_neuronal_glial", "leiden_0_2",
                clusterIDs_removed=['0', '3', '6', '7', '8', '10', '11', '14'])
# P365, leiden_0_3, remove 18, 16, 14, 7, 6, 0, 17, 8, 11, 15, 9, 13, 10
run_pipe_subset("normal_P365", "subset_neuronal_glial", "leiden_0_3",
                clusterIDs_removed=['0', '6', '7', '8', '9', '10', '11', '13', '14', '15', '16', '17', '18'])




# run_pipe_subset("normal_merged_subset_neuronal_glial", "subset_RG_OPC_NSC_lineage", "leiden_0_3",
#                 clusterIDs_selected=['5', '2', '7', '8', '6', '18', '4', '11'])
# run_pipe_subset("normal_merged_subset_neuronal_glial_subset_RG_OPC_NSC_lineage", "cleaned", "leiden_0_3",
#                 clusterIDs_removed=['11', '12'])


# 2020-11-09
# >>> sc.logging.print_header()
# scanpy==1.6.0 anndata==0.7.4 umap==0.3.10 numpy==1.19.2 scipy==1.5.2 pandas==1.1.3 scikit-learn==0.23.2 statsmodels==0.12.1 python-igraph==0.7.1 louvain==0.6.1 leidenalg==0.7.0

# run_pipe_subset("normal", "subset_neuronal_glial", "leiden_0_4",
#                 clusterIDs_removed=['0', '7', '17', '18', '19', '20', '23', '25', '27', '28'])

run_pipe_subset("normal_merged_subset_neuronal_glial", "cleaned", "leiden_0_2",
                clusterIDs_removed=['14'])


# without choroid plexus epithelia
run_pipe_subset("normal", "subset_neuronal_glial", "leiden_0_4",
                clusterIDs_removed=['0', '7', '17', '18', '19', '20', '23', '24', '25', '27', '28'])
run_pipe_subset("normal_subset_neuronal_glial", "cleaned", "leiden_0_2",
                clusterIDs_removed=['15'])

# 2020-11-13
run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_RG", "leiden_0_4",
                clusterIDs_selected=['2', '5'])

run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_RG_NSC_lineage", "leiden_0_4",
                clusterIDs_selected=['2', '4', '5', '8', '11'])
run_pipe_subset("normal_subset_neuronal_glial_cleaned_subset_RG_NSC_lineage", "cleaned", "leiden_0_3",
                clusterIDs_removed=['4', '11'])
                
run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_qNSCs", "leiden_0_4",
                clusterIDs_selected=['4', '11'])


# 2020-11-15
run_pipe_subset("normal", "subset_ependymal", "leiden_0_4",
                clusterIDs_selected=['16'])
run_pipe_subset("normal_subset_ependymal", "cleaned", "leiden_0_2",
                clusterIDs_removed=['5', '6'])


run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_ependymal", "leiden_0_4",
                clusterIDs_selected=['14'])
run_pipe_subset("normal_subset_neuronal_glial_cleaned_subset_ependymal", "cleaned", "leiden_0_1",
                clusterIDs_removed=['2', '3'])


run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_OPCs", "leiden_0_4",
                clusterIDs_selected=['6'])

run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_neuronal", "leiden_0_4",
                clusterIDs_selected=['0', '1', '9', '10', '12', '13', '15', '17'])
run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_RG_neuronal", "leiden_0_4",
                clusterIDs_selected=['0', '1', '2', '5', '9', '10', '12', '13', '15', '17'])
run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_RG_glial", "leiden_0_4",
                clusterIDs_selected=['2','3','4','5','6','7','8','18','16','11','14'])




run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_RG_NSC_OPC_lineage", "leiden_0_4",
                clusterIDs_selected=['2', '4', '5', '6', '7', '8', '11'])

# add_annotation("normal_subset_neuronal_glial_cleaned", "leiden_0_4", annot_subset_neuronal_glial)

run_pipe_subset("normal", "subset_embryonic_RG", "leiden_0_4",
                clusterIDs_selected=['6'])
run_pipe_subset("normal_subset_embryonic_RG", "cleaned", "leiden_0_7",
                clusterIDs_removed=['10', '11', '15', '17', '18', '19', '20' ,'21'])

run_pipe_subset("normal", "subset_juvenile_RG_TAPs", "leiden_0_4",
                clusterIDs_selected=['3'])
# run_pipe_subset("normal_subset_juvenile_RG_TAPs", "subset_juvenile_RG", "leiden_0_5",
#                 clusterIDs_selected=['0', '1', '6'])
# run_pipe_subset("normal_subset_juvenile_RG_TAPs", "subset_juvenile_RG", "leiden_1_0",
#                 clusterIDs_selected=['0', '2', '4', '12', '10', '3', '13', '5', '14'])
# 2020-12-01
run_pipe_subset("normal_subset_juvenile_RG_TAPs", "subset_juvenile_RG", "leiden_0_33",
                clusterIDs_selected=['1', '2', '4', '7'])



run_pipe_subset("normal", "subset_aNSCs", "leiden_0_4",
                clusterIDs_selected=['11'])
run_pipe_subset("normal", "subset_qNSCs", "leiden_0_4",
                clusterIDs_selected=['5', '14'])


run_pipe_subset("normal", "subset_RG", "leiden_0_4",
                clusterIDs_selected=['3', '6'])
run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_RG", "leiden_0_4",
                clusterIDs_selected=['2', '5'])


# 2020-11-23
run_pipe_subset("normal", "subset_OPCs", "leiden_0_4",
                clusterIDs_selected=['10'])


# 2020-12-08
run_pipe_subset("normal_subset_neuronal_glial_cleaned", "subset_NBs", "leiden_0_4",
                clusterIDs_selected=['0', '1', '9', '10', '12', '15'])

run_pipe_subset("normal_subset_neuronal_glial_cleaned_subset_NBs", "cleaned", "leiden_0_1",
                clusterIDs_removed=['6'])


run_pipe_subset("normal", "subset_ependymal", "leiden_0_4",
                clusterIDs_selected=['16'])

run_pipe_subset("normal_subset_ependymal", "cleaned", "leiden_0_2",
                clusterIDs_removed=['5', '6'])

# 2020-12-09
run_pipe_subset("normal_subset_embryonic_RG_cleaned", "wo_gliogenic", "leiden_0_33",
                clusterIDs_removed=['1', '5'])


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################




def run_pipe_merged_subset(subsetIDs, fileID):
        
        list_adata = []
        
        for subsetID in subsetIDs:
            # print(sampleID)
            file_path_main = "data/forebrain_"+subsetID+".h5ad"
            adata_temp = sc.read_h5ad(file_path_main)
            list_adata.append(adata_temp)
        
        adata_merged = list_adata[0].concatenate(list_adata[1:],
                                                 join='outer',
                                                 fill_value=0,
                                                 batch_key = 'subsetID',
                                                 batch_categories = subsetIDs,
                                                 index_unique=None)
        sampleIDs = adata_merged.obs.sampleID.unique()
        adata = get_merged_adata(sampleIDs)
        adata.obs.index = pd.Index(adata.obs['cellID'])
        
        adata = adata[adata_merged.obs_names.values, :].copy()
        filter_meta_cols = [not x.startswith(('annot_', 'leiden_')) for x in adata_merged.obs.columns]
        adata.obs = adata_merged.obs.loc[adata.obs_names.values, filter_meta_cols]
        
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



run_pipe_merged_subset(["normal_subset_embryonic_RG_cleaned",
                        "normal_subset_juvenile_RG_TAPs_subset_juvenile_RG",
                        "normal_subset_aNSCs",
                        "normal_subset_qNSCs"],
                        "normal_merged_RG_NSC_lineage")

run_pipe_merged_subset(["normal_subset_embryonic_RG_cleaned",
                        "normal_subset_juvenile_RG_TAPs_subset_juvenile_RG",
                        "normal_subset_aNSCs"],
                        "normal_merged_RG_aNSC_lineage")



run_pipe_merged_subset(["normal_E14_5_subset_neuronal_glial",
                        "normal_E16_5_subset_neuronal_glial",
                        "normal_E12_5_subset_neuronal_glial",
                        "normal_P0_subset_neuronal_glial",
                        "normal_P2_subset_neuronal_glial",
                        "normal_P7_subset_neuronal_glial",
                        "normal_P13_subset_neuronal_glial",
                        "normal_P19_subset_neuronal_glial",
                        "normal_P39_subset_neuronal_glial",
                        "normal_P111_subset_neuronal_glial",
                        "normal_P365_subset_neuronal_glial"],
                        "normal_merged_subset_neuronal_glial")

subsetIDs = ["normal_E14_5_subset_neuronal_glial",
                        "normal_E16_5_subset_neuronal_glial",
                        "normal_E12_5_subset_neuronal_glial",
                        "normal_P0_subset_neuronal_glial",
                        "normal_P2_subset_neuronal_glial",
                        "normal_P7_subset_neuronal_glial",
                        "normal_P13_subset_neuronal_glial",
                        "normal_P19_subset_neuronal_glial",
                        "normal_P39_subset_neuronal_glial",
                        "normal_P111_subset_neuronal_glial",
                        "normal_P365_subset_neuronal_glial"]

fileID = "normal_merged_subset_neuronal_glial"

# # 2020-07-10
# # ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
# run_pipe_subset("normal_E12_5", "subset_neuronal_glial", "leiden_0_1", ['0', '1', '2', '3', '5', '6', '7', '9'])
# run_pipe_subset("normal_E14_5", "subset_neuronal_glial", "leiden_0_05", ['0', '1', '2', '3'])
# run_pipe_subset("normal_E16_5", "subset_neuronal_glial", "leiden_0_05", ['0', '1', '4', '5', '6', '7', '8'])
# run_pipe_subset("normal_P0", "subset_neuronal_glial", "leiden_0_2", ['0', '3', '4', '5', '6', '8', '9', '11', '12'])
# run_pipe_subset("normal_P2", "subset_neuronal_glial", "leiden_0_1", ['0', '1', '2', '3', '4'])
# run_pipe_subset("normal_P7", "subset_neuronal_glial", "leiden_0_05", ['1', '2', '3', '4', '7'])
# run_pipe_subset("normal_P13", "subset_neuronal_glial", "leiden_0_2", ['1', '2', '3', '4', '6', '7'])
# run_pipe_subset("normal_P19", "subset_neuronal_glial", "leiden_0_2", ['1', '2', '3', '5', '6', '7', '8', '11', '14', '15'])
# run_pipe_subset("normal_P39", "subset_neuronal_glial", "leiden_0_2", ['0', '2', '4', '5', '6', '7', '8', '11'])
# run_pipe_subset("normal_P111", "subset_neuronal_glial", "leiden_0_2", ['1', '2', '3', '4', '6', '10', '13', '14'])
# run_pipe_subset("normal_P365", "subset_neuronal_glial", "leiden_0_4", ['1', '2', '3', '4', '5', '8', '9', '10', '12', '13', '17', '18'])





# run_pipe_subset("normal", "subset_ependymal", "leiden_0_45", ['17'])
# run_pipe_subset("normal_subset_ependymal", "cleaned", "leiden_0_2", ['0', '1', '2', '3', '4', '5'])
# run_pipe_subset("normal_subset_ependymal", "cleaned", "leiden_0_6", ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])



# run_pipe(sampleIDs=['HAM12604', 'HAM12605',
#                     'HAM13754', 'HAM13755',
#                     'HAM13072', 'HAM13073'], fileID='normal_prenatal')
# 
# 
# run_pipe(sampleIDs=['HAM11199', 'HAM11200',
#                     'HAM10572', 'HAM10573',
#                     'HAM10364', 'HAM10365'], fileID='normal_earlypostnatal')


# run_pipe(sampleIDs=['HAM11199', 'HAM11200',
#                     'HAM10572', 'HAM10573',
#                     'HAM10364', 'HAM10365',
#                     'HAM10288', 'HAM10289',
#                     'HAM11302', 'HAM11303'], fileID='normal_postnataldev')


# run_pipe(sampleIDs=['HAM10493', 'HAM10494',
#                     'HAM12214', 'HAM12215',
#                     'HAM13943', 'HAM13944'], fileID='normal_aging')



# ################################### not run 2020-07-20
# # ara-C
# run_pipe(sampleIDs=['HAM11657', 'HAM11658'], fileID='araC_T5')
# run_pipe(sampleIDs=['HAM11087', 'HAM11088'], fileID='araC_T5_2')
# run_pipe(sampleIDs=['HAM11197', 'HAM11198'], fileID='araC_T5_4')
# run_pipe(sampleIDs=['HAM11160', 'HAM11161'], fileID='araC_T5_7')
# run_pipe(sampleIDs=['HAM12271', 'HAM12272'], fileID='araC_T5_63')
# 
# run_pipe(sampleIDs=['HAM11657', 'HAM11658',
#                     'HAM11087', 'HAM11088',
#                     'HAM11197', 'HAM11198',
#                     'HAM11160', 'HAM11161',
#                     'HAM12271', 'HAM12272'], fileID='araC')
# 
# run_pipe(sampleIDs=['HAM11657', 'HAM11658',
#                     'HAM11087', 'HAM11088',
#                     'HAM11197', 'HAM11198',
#                     'HAM11160', 'HAM11161',
#                     'HAM12271', 'HAM12272',
#                     'HAM10493', 'HAM10494',
#                     'HAM12214', 'HAM12215',
#                     'HAM10811', 'HAM10812',
#                     'HAM10809', 'HAM10810'], fileID='araC_normal_saline')
# 
# 
# # saline
# run_pipe(sampleIDs=['HAM10811',
#                     'HAM10809',
#                     'HAM10493'], fileID='saline_T5_control_SOX2_pos')
# 
# # run_pipe(sampleIDs=['HAM10809', 'HAM10810',
# #                     'HAM10811', 'HAM10812',
# #                     'HAM10493', 'HAM10494'], fileID='saline_T5_control')
# 
# 
# # cancer
# 
# # P53-KO, PTEN-KO - majority of td tomato negative samples missing
# run_pipe(sampleIDs=['HAM11674',
#                     'HAM11675',
#                     'HAM11712',
#                     'HAM11713',
#                     'HAM13442',
#                     'HAM13797',
#                     'HAM13348'], fileID='P53-KO_PTEN-KO')
# 
# run_pipe(sampleIDs=['HAM11674',
#                     'HAM11712',
#                     'HAM13442',
#                     'HAM13797',
#                     'HAM13348'], fileID='P53-KO_PTEN-KO_tdTomato_positive')
# 
# 
# run_pipe(sampleIDs=['HAM13725',
#                     'HAM12214',
#                     'HAM11674',
#                     'HAM11712',
#                     'HAM13442',
#                     'HAM13797',
#                     'HAM13348'], fileID='P53-KO_PTEN-KO_tdTomato_positive_controls')
# 
# 
# run_pipe(sampleIDs=['HAM13725'], fileID='control_tdTomato_positive')
# 
# run_pipe(sampleIDs=['HAM11674', 'HAM11675'], fileID='P53-KO_PTEN-KO_preneoplastic-1')
# run_pipe(sampleIDs=['HAM11712', 'HAM11713'], fileID='P53-KO_PTEN-KO_preneoplastic-2')
# # run_pipe(sampleIDs=['HAM13442', 'HAM13443'], fileID='P53-KO_PTEN-KO_early-lesions-1')
# run_pipe(sampleIDs=['HAM13442'], fileID='P53-KO_PTEN-KO_early-lesions-1')
# run_pipe(sampleIDs=['HAM13797'], fileID='P53-KO_PTEN-KO_early-lesions-2')
# run_pipe(sampleIDs=['HAM13348'], fileID='P53-KO_PTEN-KO_endpoint-tumour-1')
# 
# 
# 
# 
# # run_pipe(sampleIDs=['HAM13725',
# #                     'HAM12279',
# #                     'HAM11674', 'HAM11712',
# #                     'HAM13442',
# #                     'HAM13797',
# #                     'HAM13348'], fileID='tumour_P53-KO_PTEN-KO_control_SOX2_lineage')
# # 
# # 
# # 
# # run_pipe(sampleIDs=['HAM13725',
# #                     'HAM12279',
# #                     'HAM12301'], fileID='tumour_P53-KO_injury_control_SOX2_lineage')
# 
# 
# # TO RUN TO RUN TO RUN TO RUN TO RUN TO RUN TO RUN TO RUN TO RUN TO RUN TO RUN TO RUN
# 
# # metadata_exp = pd.read_csv("data/metadata_experiment.tsv", sep='\t')
# # # filter_set1 = metadata_exp["set"] == "set1"
# # # filter_prenatal = metadata_exp['postnatal_day'] < 0
# # filter_set = metadata_exp.timepoint == 'E16.5'
# # sampleIDs = list(metadata_exp[filter_set].sampleID.sort_values())
# 
# ################################################################################
# ################################################################################
# # saline
# # run_pipe(sampleIDs=['HAM10811', 'HAM10812'], fileID='saline_T5')
# # run_pipe(sampleIDs=['HAM10809', 'HAM10810'], fileID='saline_control_T5')
# # 
# # run_pipe(sampleIDs=['HAM10878', 'HAM10879'], fileID='saline_T5_D7')
# # run_pipe(sampleIDs=['HAM10880', 'HAM10881'], fileID='saline_control_T5_D7')
# # 
# 
# # some samples are missing!
# # run_pipe(sampleIDs=['HAM12301', 'HAM12302'], fileID='lesions_saline_P53-KO')
# # run_pipe(sampleIDs=['HAM12279', 'HAM12280'], fileID='normal_P53-KO')
# # run_pipe(sampleIDs=['HAM11674', 'HAM11675'], fileID='normal_preneoplastic_P53-KO_PTEN-KO')
# # run_pipe(sampleIDs=['HAM13442', 'HAM13443'], fileID='lesions_P53-KO_PTEN-KO')
# # run_pipe(sampleIDs=['HAM13442'], fileID='lesions_P53-KO_PTEN-KO')
# # run_pipe(sampleIDs=['HAM13348', 'HAM13353'], fileID='tumour_P53-KO_PTEN-KO')
# # run_pipe(sampleIDs=['HAM13348'], fileID='tumour_P53-KO_PTEN-KO')
# # 
# # run_pipe(sampleIDs=['HAM11674', 'HAM13442', 'HAM13348'], fileID='tdTomato-pos_P53-KO_PTEN-KO')
# # 
# # run_pipe(sampleIDs=['HAM13684', 'HAM12279', 'HAM12301'], fileID='tdTomato-pos_P53-KO')
# # run_pipe(sampleIDs=['HAM13725', 'HAM13684', 'HAM12279', 'HAM12301'], fileID='tdTomato-WT-pos_P53-KO')
# # 
# # 
# # run_pipe(sampleIDs=['HAM12279', 'HAM12301', 'HAM11674', 'HAM13442', 'HAM13348'], fileID='tdTomato-pos_P53-KO_P53-KO_PTEN-KO')
# # 
# # # combinations
# # run_pipe(sampleIDs=['HAM12301', 'HAM10811'], fileID='tdTomato-pos_lesions_saline_P53-KO-SOX2-pos_saline_T5')
# 
# # 
# # 
# # run_pipe(sampleIDs=['HAM12604', 'HAM12605',
# #                     'HAM13754', 'HAM13755',
# #                     'HAM13072', 'HAM13073',
# #                     'HAM11199', 'HAM11200',
# #                     'HAM10572', 'HAM10573',
# #                     'HAM10364', 'HAM10365',
# #                     'HAM10288', 'HAM10289',
# #                     'HAM11302', 'HAM11303',
# #                     'HAM10493', 'HAM10494',
# #                     'HAM12214', 'HAM12215',
# #                     'HAM12301', 'HAM12302'], fileID='normal_cancer_D90')
# # 
# # 
# 
# 
# # 
# # 
# # run_pipe(sampleIDs=['HAM10811', 'HAM10812',
# #                     'HAM10878', 'HAM10879',
# #                     'HAM12301', 'HAM12302'], fileID='saline_T5_saline_T5_D7_cancer_D90')
# 
# 
# 
# 
# 
# 
# run_pipe_subset("normal_E12_5", "subset_RGCs", "leiden_0_1", ['0'])
# run_pipe_subset("normal_E14_5", "subset_RGCs", "leiden_0_1", ['2'])
# run_pipe_subset("normal_E16_5", "subset_RGCs", "leiden_0_1", ['3'])
# run_pipe_subset("normal_prenatal", "subset_RGCs", "leiden_0_2", ['0', '6'])
# run_pipe_subset("normal_P0", "subset_RGCs", "leiden_0_2", ['2', '3'])
# run_pipe_subset("normal_P2", "subset_RGCs", "leiden_0_2", ['2', '4'])
# run_pipe_subset("normal_P7", "subset_RGCs", "leiden_0_3", ['1', '7', '12'])
# run_pipe_subset("normal_P13", "subset_RGCs", "leiden_0_1", ['1', '2'])
# run_pipe_subset("normal_P19", "subset_RGCs", "leiden_0_2", ['1', '5', '6', '15'])
# run_pipe_subset("normal_P39", "subset_RGCs", "leiden_0_2", ['2', '3', '8'])
# run_pipe_subset("normal_P111", "subset_RGCs", "leiden_0_7", ['0', '5', '12'])
# run_pipe_subset("normal_P365", "subset_RGCs", "leiden_0_7", ['0', '5', '17'])
# 
# 
# # clean up subset RGCs and remove any debris
# run_pipe_subset("normal_E12_5_subset_RGCs", "wo_NCCs", "leiden_0_2", ['0', '1', '2', '3', '4', '5', '6', '7'])
# 
# run_pipe_subset("normal_E12_5_subset_RGCs", "cleaned", "leiden_0_2", ['0', '1', '2', '3', '5', '7'])
# run_pipe_subset("normal_E14_5_subset_RGCs", "cleaned", "leiden_1_0", ['0', '1', '3', '4', '5', '8', '9', '10', '14'])
# run_pipe_subset("normal_E16_5_subset_RGCs", "cleaned", "leiden_0_1", ['0'])
# run_pipe_subset("normal_P0_subset_RGCs", "cleaned", "leiden_0_3", ['0', '1', '3', '4'])
# run_pipe_subset("normal_P2_subset_RGCs", "cleaned", "leiden_1_0", ['0', '2', '3', '4', '6', '7', '9', '10', '11', '12', '13', '14'])
# run_pipe_subset("normal_P7_subset_RGCs", "cleaned", "leiden_1_0", ['0', '1', '2', '3', '4', '7', '10', '12', '14'])
# run_pipe_subset("normal_P13_subset_RGCs", "cleaned", "leiden_1_0", ['0', '1', '6', '7', '8', '11', '13'])
# run_pipe_subset("normal_P19_subset_RGCs", "cleaned", "leiden_0_2", ['0', '1', '3', '4'])
# run_pipe_subset("normal_P39_subset_RGCs", "cleaned", "leiden_0_2", ['0', '1'])
# run_pipe_subset("normal_P111_subset_RGCs", "cleaned", "leiden_0_2", ['0', '1'])
# run_pipe_subset("normal_P365_subset_RGCs", "cleaned", "leiden_1_0", ['0', '1', '2', '3', '4', '5', '6', '7', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27'])
# 
# 
# # 2020-07-10
# # ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
# run_pipe_subset("normal_E12_5", "subset_neuronal_glial", "leiden_0_1", ['0', '1', '2', '3', '5', '6', '7', '9'])
# run_pipe_subset("normal_E14_5", "subset_neuronal_glial", "leiden_0_05", ['0', '1', '2', '3'])
# run_pipe_subset("normal_E16_5", "subset_neuronal_glial", "leiden_0_05", ['0', '1', '4', '5', '6', '7', '8'])
# run_pipe_subset("normal_P0", "subset_neuronal_glial", "leiden_0_2", ['0', '3', '4', '5', '6', '8', '9', '11', '12'])
# run_pipe_subset("normal_P2", "subset_neuronal_glial", "leiden_0_1", ['0', '1', '2', '3', '4'])
# run_pipe_subset("normal_P7", "subset_neuronal_glial", "leiden_0_05", ['1', '2', '3', '4', '7'])
# run_pipe_subset("normal_P13", "subset_neuronal_glial", "leiden_0_2", ['1', '2', '3', '4', '6', '7'])
# run_pipe_subset("normal_P19", "subset_neuronal_glial", "leiden_0_2", ['1', '2', '3', '5', '6', '7', '8', '11', '14', '15'])
# run_pipe_subset("normal_P39", "subset_neuronal_glial", "leiden_0_2", ['0', '2', '4', '5', '6', '7', '8', '11'])
# run_pipe_subset("normal_P111", "subset_neuronal_glial", "leiden_0_2", ['1', '2', '3', '4', '6', '10', '13', '14'])
# run_pipe_subset("normal_P365", "subset_neuronal_glial", "leiden_0_4", ['1', '2', '3', '4', '5', '8', '9', '10', '12', '13', '17', '18'])
# 
# 
# # 2020-07-20
# run_pipe_subset("normal", "subset_neuroblasts", "leiden_0_52", ['1', '2', '8', '12', '14'])
# 
# # 2020-10-09
# run_pipe_subset("normal", "subset_neuronal_glial", "leiden_0_52", ['1', '2', '3', '4', '5', '7', '8', '9', '10', '11', '12', '13', '14', '15', '17', '21', '22', '25'])
# 
# 
# ################################################################################
# # remove debris cluster (no reclustering for now)
# fileID = "normal_subset_neuroblasts"
# subset_leidenID = "leiden_0_15"
# file_path_main = "data/brain-regeneration_"+fileID+".h5ad"
# adata = sc.read_h5ad(file_path_main)
# 
# filter_sub = [x in ['0', '1', '2', '3', '4', '5', '6'] for x in adata.obs[subset_leidenID]]
# 
# adata = adata[filter_sub,:].copy()
# 
# fileID = "normal_subset_neuroblasts_nodebris"
# # save results
# adata.write('data/brain-regeneration_'+fileID+'.h5ad')
# ################################################################################
# 
# 
# 
# run_pipe_subset("normal", "subset_ependymal", "leiden_0_7", ['18'])
# # run_pipe_subset("normal_subset_ependymal", "cleaned", "leiden_0_2", ['0', '1', '2', '3', '4', '5'])
# run_pipe_subset("normal_subset_ependymal", "cleaned", "leiden_0_6", ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10'])
# 
# 
# 
# run_pipe_subset("tumour_P53-KO_injury_control_SOX2_lineage", "subset_cancer", "leiden_0_7", ['0', '5', '8'])
# run_pipe_subset("tumour_P53-KO_PTEN-KO_control_SOX2_lineage", "subset_cancer", "leiden_0_3", ['4', '7', '10', '13'])
# 
# run_pipe_subset("saline_T5_control_SOX2_pos", "subset_injury_v1", "leiden_0_2", ['6'])
# run_pipe_subset("saline_T5_control_SOX2_pos", "subset_injury_v2", "leiden_0_7", ['14'])
# 
# 
# 
# # run_pipe_subset("normal_oligo", "subset_OPCs", "leiden_0_2", ['0', '5', '6'])
# # run_pipe_subset("normal_oligo_OLD", "subset_OPCs", "leiden2", ['0', '5', '6'])
# 
# # new
# # run_pipe_subset("normal", "subset_OPCs_Oligos", "leiden_0_7", ['3', '9', '17', '25', '29'])
# run_pipe_subset("normal", "subset_OPCs_oligos", "leiden_0_7", ['6', '30', '2', '27'])
# run_pipe_subset("normal_subset_OPCs_oligos", "cleaned", "leiden_0_25", ['0', '1', '2', '3', '4', '5', '6', '7'])
# 
# run_pipe_subset("normal_subset_OPCs_oligos_cleaned", "subset_OPCs", "leiden_0_2", ['0', '6', '7'])
# 
# # OLD
# run_pipe_subset("normal_subset_neuronal_glial", "subset_OPCs_Oligos", "leiden_0_3", ['4', '7', '16', '18'])
# 
# 
# # RGC subsets:
# run_pipe_subset("normal_RGCs_cleaned", "subset_early_embryonic", "leiden_0_13", ['1'])
# run_pipe_subset("normal_RGCs_cleaned", "subset_late_embryonic", "leiden_0_13", ['2'])
# run_pipe_subset("normal_RGCs_cleaned", "subset_juvenile", "leiden_0_13", ['4'])
# run_pipe_subset("normal_RGCs_cleaned", "subset_adult", "leiden_0_13", ['0', '3'])
# 
# run_pipe_subset("normal_RGCs_cleaned_subset_early_embryonic", "cleaned", "leiden_0_4", ['0', '1', '2', '3', '4', '5', '6', '7', '8'])
# 
# 
# 
# 
# 
# 
# 
# 
# # run_pipe_merged_subset(["normal_E12_5_subset_RGCs",
# #                         "normal_E14_5_subset_RGCs",
# #                         "normal_E16_5_subset_RGCs"],
# #                             "normal_prenatal_RGCs")
# 
# 
# # run_pipe_merged_subset(["normal_E12_5_subset_RGCs",
# #                         "normal_E14_5_subset_RGCs",
# #                         "normal_E16_5_subset_RGCs",
# #                         "normal_P0_subset_RGCs",
# #                         "normal_P2_subset_RGCs",
# #                         "normal_P7_subset_RGCs",
# #                         "normal_P13_subset_RGCs",
# #                         "normal_P19_subset_RGCs",
# #                         "normal_P39_subset_RGCs",
# #                         "normal_P111_subset_RGCs",
# #                         "normal_P365_subset_RGCs"],
# #                         "normal_RGCs")
# # run_pipe_subset("normal_RGCs", "cleaned", "leiden_0_3", ['0', '1', '2', '3', '4', '5', '6', '7', '8'])
# 
# run_pipe_merged_subset(["normal_E12_5_subset_RGCs_cleaned",
#                         "normal_E14_5_subset_RGCs_cleaned",
#                         "normal_E16_5_subset_RGCs_cleaned",
#                         "normal_P0_subset_RGCs_cleaned",
#                         "normal_P2_subset_RGCs_cleaned",
#                         "normal_P7_subset_RGCs_cleaned",
#                         "normal_P13_subset_RGCs_cleaned",
#                         "normal_P19_subset_RGCs_cleaned",
#                         "normal_P39_subset_RGCs_cleaned",
#                         "normal_P111_subset_RGCs_cleaned",
#                         "normal_P365_subset_RGCs_cleaned"],
#                         "normal_RGCs_NEW")
# run_pipe_subset("normal_RGCs_NEW", "cleaned", "leiden_0_3", ['0', '1', '2', '3', '4', '5', '6', '7', '8'])
# 
# 
# # to run:
# 
# run_pipe_merged_subset(["normal_subset_neuronal_glial",
#                         "tumour_P53-KO_injury_control_SOX2_lineage_subset_cancer",
#                         "tumour_P53-KO_PTEN-KO_control_SOX2_lineage_subset_cancer"],
#                         "integration-normal_subset_neuronal_glial-cancer_models")
# 
# # new cancer integration [2020-09-11]
# 
# run_pipe_merged_subset(["normal",
#                         "control_tdTomato_positive",
#                         "P53-KO_PTEN-KO"],
#                             "normal_development_control_cancer_P53-KO_PTEN-KO")
# 
# 
# 
# run_pipe_merged_subset(["normal_subset_neuronal_glial",
#                         "control_tdTomato_positive",
#                         "P53-KO_PTEN-KO"],
#                             "normal_subset_neuronal_glial_development_control_P53-KO_PTEN-KO")
# 
# run_pipe_merged_subset(["normal_subset_neuronal_glial",
#                         "control_tdTomato_positive",
#                         "P53-KO_PTEN-KO_endpoint-tumour-1"],
#                             "normal_subset_neuronal_glial_development_control_endpoint-tumour-1_P53-KO_PTEN-KO")
# 
# 
# run_pipe_merged_subset(["normal_E12_5_subset_neuronal_glial",
#                         "normal_E14_5_subset_neuronal_glial",
#                         "normal_E16_5_subset_neuronal_glial",
#                         "normal_P0_subset_neuronal_glial",
#                         "normal_P2_subset_neuronal_glial",
#                         "normal_P7_subset_neuronal_glial",
#                         "normal_P13_subset_neuronal_glial",
#                         "normal_P19_subset_neuronal_glial",
#                         "normal_P39_subset_neuronal_glial",
#                         "normal_P111_subset_neuronal_glial",
#                         "normal_P365_subset_neuronal_glial"],
#                         "normal_merged_subset_neuronal_glial")
# 
