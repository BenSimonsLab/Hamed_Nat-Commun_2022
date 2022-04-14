#################
## Cluster LFC ##
#################

library(readr)
library(dplyr)

source("code/core_settings.R")
source("code/core_functions.R")

library(reticulate)

add_cluster_LFC <- function(fileID, leiden_res){

  adata_path <<- paste0('data/forebrain_', fileID, '.h5ad')

  py_run_string("import numpy as np")
  py_run_string("import anndata as ad")

  py_run_string("adata = ad.read_h5ad(r.adata_path)")

  metadata_scanpy = as_tibble(py$adata$obs, rownames="cellID")

  # calculate cluster vs cluster LFC and add minimum LFC to DE file

  fpaths_DE = list.files(path=paste0("output/", fileID, "/DE-voom_limma"), pattern=paste0("DE-voom_limma_", fileID, "_DE_", leiden_res, "+"), full.names=T)



  sex_genes = c('Xist', 'Tsix', 'Eif2s3y', 'Ddx3y', 'Uty', 'Kdm5d')

  for (fpath in fpaths_DE){
    tbl_de = read_csv(fpath, col_types="cdddddd") %>%
                select(-matches("cluster_LFC")) # remove any previous runs
                
    fileID = basename(fpath)
    # leiden_res = stringr::str_extract(fileID, "leiden_[0-9]+_[0-9]+")
    clusterID = stringr::str_extract(stringr::str_extract(fileID, "cl_[0-9]+_vs_rest"), "[0-9]+")

    geneIDs_de = tbl_de %>%
                    filter(!(geneID %in% sex_genes)) %>%
                    pull(geneID)

    tbl_mean_expr = get_adata_raw_expr(geneIDs_de) %>%
                left_join(metadata_scanpy %>% select(cellID, !!leiden_res), by="cellID") %>%
                mutate(leiden_res = get(leiden_res)) %>%
                select(leiden_res, contains("expr_")) %>%
                group_by(leiden_res) %>%
                summarise_all(list(mean = mean))


    tbl_mean_expr_cl = tbl_mean_expr %>% filter(leiden_res == clusterID)
    tbl_mean_expr_rest = tbl_mean_expr %>% filter(leiden_res != clusterID)

    vec_mean_expr_cl = tbl_mean_expr_cl %>%
                        select(-leiden_res) %>%
                        unlist(., use.names=FALSE)
                        
    mat_expr_rest = tbl_mean_expr_rest %>%
                        select(-leiden_res) %>%
                        as.matrix() %>%
                        t()

    mat_cl_LFC = log2(mat_expr_rest) - log2(vec_mean_expr_cl)
    vec_cl_LFC = apply(mat_cl_LFC, 1, max)
    names(vec_cl_LFC) = gsub("^expr_", "", gsub("_mean$", "", names(vec_cl_LFC)))

    tbl_cl_LFC = tibble(geneID = names(vec_cl_LFC), cluster_LFC = vec_cl_LFC)

    de_results_merged = tbl_de %>%
                          left_join(tbl_cl_LFC, by="geneID")
                          
    write_csv(de_results_merged, fpath)
    print(paste0("Updated ",  basename(fpath)))
  }
  
return()
}



# add_cluster_LFC("normal_subset_embryonic_RG_cleaned", "leiden_0_33")

add_cluster_LFC("normal_merged_RG_NSC_lineage", "annot_merged_coarse_num")
add_cluster_LFC("normal_merged_RG_NSC_lineage", "annot_merged_fine_num")

# add_cluster_LFC("normal_merged_RG_aNSC_lineage", "annot_merged_num")
# add_cluster_LFC("normal_subset_qNSCs", "annot_superset_num")

# add_cluster_LFC("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", "leiden_0_13")

# add_cluster_LFC("normal_subset_OPCs", "leiden_0_1")
# add_cluster_LFC("normal_subset_ependymal_cleaned", "leiden_0_15")


