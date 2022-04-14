import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import os
from sklearn.linear_model import LogisticRegression
from collections import Counter


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


def logreg_pred(adata_ext, adata_zeisel, grouping_zeisel):
    annot = list(adata_zeisel.obs[grouping_zeisel])
    
    logisticRegr = LogisticRegression(max_iter = 10000, n_jobs = -1, random_state = 0, C=0.2, solver='liblinear', multi_class='ovr')
    logisticRegr.fit(adata_zeisel.X.todense(), annot)
    
    # predict grouping_zeisel in adata_ext    
    predictions_class = logisticRegr.predict(adata_ext.X.todense())
    probabilities = logisticRegr.predict_proba(adata_ext.X.todense())
    probs_prediction = []
    for i,clus in enumerate(predictions_class):
        probs_prediction.append(probabilities[i,logisticRegr.classes_==clus][0])
    
    results = pd.DataFrame(data={'cellID': adata_ext.obs_names,
                                 'prediction_' + grouping_zeisel: predictions_class,
                                 'probability_' + grouping_zeisel: probs_prediction})
    
    return(results)



adata_zeisel = sc.read_h5ad('data/Zeisel_2018/Zeisel_2018-raw.h5ad')



#########
## Map ##
#########


import re
import fnmatch

adata_files = fnmatch.filter(os.listdir('data'), "brain-development_*.h5ad")

files_processed = fnmatch.filter(os.listdir('data/logreg'), "logreg_Zeisel_2018_*.csv.gz")
adata_processed = [re.sub('logreg_Zeisel_2018_', 'brain-development_', re.sub('\.csv.gz$', '.h5ad', sample)) for sample in files_processed]

adata_files = [x for x in adata_files if x not in adata_processed]
adata_files = ['forebrain_normal.h5ad']


for adata_file in adata_files:
    adata = sc.read_h5ad('data/'+adata_file)
    
    sampleIDs = np.unique(adata.obs.sampleID.values)
    cellIDs = np.concatenate([adata.obs_names.values])
    
    adata_raw = get_merged_adata(sampleIDs)
    adata_raw = adata_raw[cellIDs,:].copy()
    
    adata_zeisel = sc.read_h5ad('data/Zeisel_2018/Zeisel_2018-raw.h5ad')
    
    ## gene overlap
    gene_overlap = list(set(adata_zeisel.var_names) & set(adata_raw.var_names))
    adata_raw = adata_raw[:,gene_overlap]
    adata_zeisel = adata_zeisel[:,gene_overlap]
    
    
    ## normalise raw unprocessed counts 
    sc.pp.normalize_per_cell(adata_raw, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_raw)
    sc.pp.normalize_per_cell(adata_zeisel, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_zeisel)
    
    pred_Class = logreg_pred(adata_raw, adata_zeisel, "Class")
    pred_TaxonomyRank1 = logreg_pred(adata_raw, adata_zeisel, "TaxonomyRank1")
    
    dfs = [pred_Class,
           pred_TaxonomyRank1]
    
    
    dfs = [df.set_index('cellID') for df in dfs]
    dfs = dfs[0].join(dfs[1:])
    
    dfs.to_csv('data/logreg/logreg_Zeisel_2018_'+adata_file.split(".", 1)[0].split("_", 1)[1]+'.csv.gz', compression='gzip')


# run on 2020-11-01
# sc.logging.print_header()
## scanpy==1.6.0 anndata==0.7.4 umap==0.3.10 numpy==1.19.2 scipy==1.5.2 pandas==1.1.3 scikit-learn==0.23.2 statsmodels==0.12.1 python-igraph==0.7.1 louvain==0.6.1 leidenalg==0.7.0

