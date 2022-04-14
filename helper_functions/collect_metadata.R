library(readr)
library(dplyr)
library(tibble)


library(reticulate)
use_python("/usr/bin/python3")


fileID = "normal_raw"
leiden_res_annot = 'annot_leiden'

adata_path = paste0('data/forebrain_', fileID, '.h5ad')

py_run_string("import numpy as np")
py_run_string("import anndata as ad")
py_run_string("adata = ad.read_h5ad(r.adata_path)")

metadata_scanpy = as_tibble(py$adata$obs) %>%
                    select(cellID, sampleID)

metadata_experiment = read_tsv("data/metadata_experiment.tsv")
results_cyclone = read_csv("data/cyclone/cyclone_merged.csv.gz", col_types="ccddd") %>%
                    mutate(sampleID = sapply(cellID, function(cell) strsplit(cell, split="-")[[1]][3])) %>%
                    left_join(metadata_experiment, by="sampleID") %>%
                    mutate(cellID = gsub("-1", "", cellID)) %>%
                    select(-sampleID)
                    # select(cellID, tissue, gate, timepoint, postnatal_day, set, treatment, condition, cyclone, cyclone_G1, cyclone_S, cyclone_G2M)

results_logreg = read_csv(paste0("data/logreg/logreg_Zeisel_2018_normal.csv.gz"))


# filter QC
metadata = left_join(metadata_scanpy, results_cyclone, by="cellID") %>%
              mutate(cyclone = factor(cyclone, levels=c("G1", "S", "G2M"))) %>%
              mutate(timepoint = factor(timepoint, levels=unique(timepoint)[order(unique(postnatal_day))])) %>%
              left_join(results_logreg, by="cellID")


metadata = metadata %>%
            mutate(prediction_TaxonomyRank1 = gsub("Neurons", "Neuronal", prediction_TaxonomyRank1)) %>%
            mutate(prediction_TaxonomyRank1 = gsub("^Glia$", "Glial", prediction_TaxonomyRank1)) %>%
            mutate(gate_gfp = gsub("SOX2_-ve", "GFP-", gate)) %>%
            mutate(gate_gfp = gsub("SOX2_\\+ve", "GFP+ [SOX2]", gate_gfp)) %>%
            mutate(QC_pass = !is.na(prediction_TaxonomyRank1)) %>%
            relocate(QC_pass, gate_gfp, .after = sampleID) %>%
            relocate(cyclone, cyclone_G1, cyclone_S, cyclone_G2M, .after = last_col())


##############
## Overview ##
##############

fileID = "normal"
leiden_res_annot = 'annot_leiden'
adata_path = paste0('data/forebrain_', fileID, '.h5ad')


py_run_string("import numpy as np")
py_run_string("import anndata as ad")
py_run_string("adata = ad.read_h5ad(r.adata_path)")

metadata_scanpy = as_tibble(py$adata$obs, rownames="cellID") %>%
                    rename(mt_frac = percent_mito) %>%
                    mutate(UMAP_dim1 = py$adata$obsm[['X_umap']][,1],
                           UMAP_dim2 = py$adata$obsm[['X_umap']][,2]) %>%
                    mutate(annot_leiden = plyr::revalue(annot_leiden, c("[unclear/debris]" = NA))) %>%
                    select(-starts_with("leiden_"), -starts_with("annot_leiden_")) %>%
                    rename(UMAP_dim1_full = UMAP_dim1,
                           UMAP_dim2_full = UMAP_dim2,
                           annot_leiden_full = annot_leiden) %>%
                    relocate(annot_leiden_full, UMAP_dim1_full, UMAP_dim2_full, .after = last_col()) %>%
                    select(cellID, n_genes, n_counts, mt_frac, annot_leiden_full, UMAP_dim1_full, UMAP_dim2_full)

metadata = left_join(metadata, metadata_scanpy, by="cellID")


######################
## Neuronal & glial ##
######################

fileID = "normal_subset_neuronal_glial_cleaned"
label = "neuronal_glial"

adata_path = paste0('data/forebrain_', fileID, '.h5ad')

py_run_string("import numpy as np")
py_run_string("import anndata as ad")
py_run_string("adata = ad.read_h5ad(r.adata_path)")

metadata_scanpy = as_tibble(py$adata$obs, rownames="cellID") %>%
                    rename(mt_frac = percent_mito) %>%
                    mutate(UMAP_dim1 = py$adata$obsm[['X_umap']][,1],
                           UMAP_dim2 = py$adata$obsm[['X_umap']][,2]) %>%
                    select(cellID, annot_leiden, UMAP_dim1, UMAP_dim2) %>%
                    rename(!!paste0("UMAP_dim1_", label) := UMAP_dim1,
                           !!paste0("UMAP_dim2_", label) := UMAP_dim2,
                           !!paste0("annot_leiden_", label) := annot_leiden)


metadata = left_join(metadata, metadata_scanpy, by="cellID")


##########
## OPCs ##
##########

fileID = "normal_subset_OPCs"
label = "OPCs"

adata_path = paste0('data/forebrain_', fileID, '.h5ad')

py_run_string("import numpy as np")
py_run_string("import anndata as ad")
py_run_string("adata = ad.read_h5ad(r.adata_path)")

metadata_scanpy = as_tibble(py$adata$obs, rownames="cellID") %>%
                   rename(mt_frac = percent_mito) %>%
                   mutate(UMAP_dim1 = py$adata$obsm[['X_umap']][,1],
                          UMAP_dim2 = py$adata$obsm[['X_umap']][,2]) %>%
                   select(cellID, annot_leiden, UMAP_dim1, UMAP_dim2) %>%
                   rename(!!paste0("UMAP_dim1_", label) := UMAP_dim1,
                          !!paste0("UMAP_dim2_", label) := UMAP_dim2,
                          !!paste0("annot_leiden_", label) := annot_leiden)


metadata = left_join(metadata, metadata_scanpy, by="cellID")


##################
## Embryonic RG ##
##################

fileID = "normal_subset_embryonic_RG_cleaned"
label = "embryonic_RG"

adata_path = paste0('data/forebrain_', fileID, '.h5ad')

py_run_string("import numpy as np")
py_run_string("import anndata as ad")
py_run_string("adata = ad.read_h5ad(r.adata_path)")

metadata_scanpy = as_tibble(py$adata$obs, rownames="cellID") %>%
                   rename(mt_frac = percent_mito) %>%
                   mutate(UMAP_dim1 = py$adata$obsm[['X_umap']][,1],
                          UMAP_dim2 = py$adata$obsm[['X_umap']][,2]) %>%
                   select(cellID, annot_leiden, UMAP_dim1, UMAP_dim2) %>%
                   rename(!!paste0("UMAP_dim1_", label) := UMAP_dim1,
                          !!paste0("UMAP_dim2_", label) := UMAP_dim2,
                          !!paste0("annot_leiden_", label) := annot_leiden)


metadata = left_join(metadata, metadata_scanpy, by="cellID")


#################
## Neuroblasts ##
#################

fileID = "normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned"
label = "NBs"

adata_path = paste0('data/forebrain_', fileID, '.h5ad')

py_run_string("import numpy as np")
py_run_string("import anndata as ad")
py_run_string("adata = ad.read_h5ad(r.adata_path)")

metadata_scanpy = as_tibble(py$adata$obs, rownames="cellID") %>%
                   rename(mt_frac = percent_mito) %>%
                   mutate(UMAP_dim1 = py$adata$obsm[['X_umap']][,1],
                          UMAP_dim2 = py$adata$obsm[['X_umap']][,2]) %>%
                   select(cellID, annot_leiden, UMAP_dim1, UMAP_dim2) %>%
                   rename(!!paste0("UMAP_dim1_", label) := UMAP_dim1,
                          !!paste0("UMAP_dim2_", label) := UMAP_dim2,
                          !!paste0("annot_leiden_", label) := annot_leiden)


metadata = left_join(metadata, metadata_scanpy, by="cellID")



############
## Export ##
############

write_csv(metadata, "data/metadata_forebrain_atlas.csv.gz")

