######################
## Run DE pipelines ##
######################

library(readr)
library(dplyr)
library(tibble)

Sys.setenv(RETICULATE_PYTHON = "/usr/bin/python3.6")
library(reticulate)
# use_python("/usr/bin/python3.6")

source("code/DE/voom_limma.R")
source_python("code/DE/get_raw_counts.py")


# filter counts
get_exprfrac_cl <- function(clusterID, leiden_res, counts, metadata){
  filter_cellIDs = metadata %>% dplyr::filter(get(leiden_res) == clusterID) %>% pull(cellID)
  expr_frac = rowSums(counts[, filter_cellIDs] >= 1)/dim(counts[, filter_cellIDs])[2]

  return(expr_frac)
}


get_filtered_counts <- function(leiden_res, counts, metadata, min_frac = 0.5){
  clusterIDs = metadata %>% pull(get(leiden_res)) %>% as.vector() %>% unique() %>% sort()
  expr_frac_cl = sapply(clusterIDs, function(clusterID) get_exprfrac_cl(clusterID, leiden_res, counts, metadata))

  # default: at least 50% in at least one cluster
  pass_cluster =  rowSums(expr_frac_cl > min_frac) > 0
  geneID_pass = names(pass_cluster)[pass_cluster]

  counts = counts[geneID_pass,]
  return(counts)
}


read_files <- function(fileID){
  adata_path = paste0('data/forebrain_', fileID, '.h5ad')
  assign("adata_path", adata_path, envir = .GlobalEnv)
  
  py_run_string("import numpy as np")
  py_run_string("import anndata as ad")
  py_run_string("adata = ad.read_h5ad(r.adata_path)")

  # use raw count matrix from scenic input
  # counts = t(read.csv(paste0("data/scenic/cmat_", fileID, ".csv.gz"), check.names=F, row.names=1))
  counts = t(get_raw_counts(fileID))

  metadata = as_tibble(py$adata$obs, rownames="cellID")

  # align counts and metadata
  counts = counts[, metadata$cellID]
  
  raw_input = list("count"=counts,
                   "metadata"=metadata)
  
  return(raw_input)
}


run_cl_rest <- function(raw_input, fileID, leiden_res, pipeline="voom_limma", clusterIDs=""){
  # run cluster vs rest DE
  if (clusterIDs == ""){
    clusterIDs = raw_input$metadata %>% pull(get(leiden_res)) %>% as.vector() %>% unique() %>% sort() %>% as.character()
  }
  
  if (pipeline == "voom_limma"){
    # prefilter geneIDs
    raw_input$count = get_filtered_counts(leiden_res, raw_input$count, raw_input$metadata)
  }
  
  for (clusterID in clusterIDs){
    print(paste0("Starting: ", leiden_res, ", cluster ", clusterID))
    raw_input$metadata = raw_input$metadata %>%
                mutate(condition = if_else(get(leiden_res) == clusterID, as.character(clusterID), "rest"))
    L = list("count"=raw_input$count,
             "condt"=raw_input$metadata %>% pull(condition, name="cellID"))
    if (pipeline == "voom_limma"){
      tbl_results = run_voomlimma(L)
    } else if (pipeline == "MAST"){
      tbl_results = run_MASTcpmDetRate(L)
      gc()
    }
    write_csv(tbl_results, paste0("output/", fileID, "/DE-", pipeline, "/DE-", pipeline, "_", fileID, "_DE_", leiden_res, "_cl_", gsub(" ", "_", clusterID), "_vs_rest.csv.gz"))
  }
  return()
}




run_DE_group1_vs_group2 <- function(raw_input, fileID, leiden_res, pipeline="voom_limma", clusterID_1="", clusterID_2=""){
  # run cluster vs rest DE
  
  raw_input$metadata = raw_input$metadata %>%
              filter(get(leiden_res) %in% c(clusterID_1, clusterID_2)) %>%
              mutate(condition = get(leiden_res))
  raw_input$count = raw_input$count[,raw_input$metadata$cellID]
  
  if (pipeline == "voom_limma"){
    # prefilter geneIDs
    raw_input$count = get_filtered_counts(leiden_res, raw_input$count, raw_input$metadata)
  }
  
  
  L = list("count"=raw_input$count,
           "condt"=raw_input$metadata %>% pull(condition, name="cellID") %>% as.character())
  if (pipeline == "voom_limma"){
    tbl_results = run_voomlimma(L)
  } else if (pipeline == "MAST"){
    tbl_results = run_MASTcpmDetRate(L)
    gc()
  }
  write_csv(tbl_results, paste0("output/", fileID, "/DE-", pipeline, "/DE-", pipeline, "_", fileID, "_DE_", leiden_res, "_cl_", clusterID_1, "_vs_", clusterID_2, ".csv.gz"))

  return()
}



##########
## OPCs ##
##########

# parameters
fileID = "normal_subset_OPCs"
leiden_res = 'leiden_0_1'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")




#####################
## Ependymal cells ##
#####################

# parameters
fileID = "normal_subset_ependymal_cleaned"
leiden_res = 'leiden_0_15'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")


#########
## NBs ##
#########

# parameters
fileID = "normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned"
leiden_res = 'leiden_0_13'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")



####################
## RG NSC lineage ##
####################

# parameters
fileID = "normal_merged_RG_NSC_lineage"
leiden_res = 'annot_merged_coarse_num'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")


# parameters
fileID = "normal_merged_RG_NSC_lineage"
leiden_res = 'annot_merged_fine_num'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")



####################
## RG aNSC lineage ##
####################

# parameters
fileID = "normal_merged_RG_aNSC_lineage"
leiden_res = 'annot_merged_num'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")




##################
## Embryonic RG ##
##################

# parameters
fileID = "normal_subset_embryonic_RG_cleaned"
leiden_res = 'annot_leiden'
leiden_res = 'leiden_0_33'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")




############################
## Compare qNSCs clusters ##
############################

# parameters
fileID = "normal_subset_qNSCs"
leiden_res = 'annot_superset_num'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")










#############################################################

# parameters
fileID = "normal_subset_neuronal_glial_cleaned_subset_qNSCs"
leiden_res = 'annot_superset'
leiden_res = 'annot_superset_num'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")


fileID = "normal_subset_neuronal_glial_cleaned_subset_qNSCs"
leiden_res = 'leiden_0_3'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")


# second comparison
# parameters
fileID = "normal"
leiden_res = 'leiden_0_4'

raw_input = read_files(fileID)
raw_input_store = raw_input

raw_input = raw_input_store


run_DE_group1_vs_group2(raw_input, fileID, leiden_res, pipeline="voom_limma", clusterID_1="5", clusterID_2="14")


############################
## Compare RG NSC lineage ##
############################

# parameters
fileID = "normal_subset_neuronal_glial_cleaned_subset_RG_NSC_lineage"
leiden_res = 'annot_superset'
leiden_res = 'annot_superset_num'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")





################################################################################
################################################################################
################################################################################

##################
## Radial Glial ##
##################

# parameters
fileID = "normal_RGCs_cleaned"
leiden_res = 'leiden_0_13'

raw_input = read_files(fileID)

# subset into developing and adult:
cellIDs_developing = raw_input$metadata %>% dplyr::filter(get(leiden_res) %in% c("1", "2", "4")) %>% pull(cellID)
raw_input_group_developing = list("count" = raw_input$count[,cellIDs_developing],
                                   "metadata" = raw_input$metadata %>% dplyr::filter(cellID %in% cellIDs_developing))

raw_input_group_developing$metadata %>% pull(get(leiden_res))

run_cl_rest(raw_input_group_developing, "normal_RGCs_cleaned_group_developing", leiden_res, pipeline="voom_limma")

cellIDs_adult = raw_input$metadata %>% dplyr::filter(get(leiden_res) %in% c("0", "3")) %>% pull(cellID)
raw_input_group_adult = list("count" = raw_input$count[,cellIDs_adult],
                                  "metadata" = raw_input$metadata %>% dplyr::filter(cellID %in% cellIDs_adult))

run_cl_rest(raw_input_group_adult, "normal_RGCs_cleaned_group_adult", leiden_res, pipeline="voom_limma")


# parameters
fileID = "normal_RGCs_cleaned"
leiden_res = 'leiden_0_13'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")


############################
## Embryonic Radial Glial ##
############################

# parameters
fileID = "normal_RGCs_cleaned_subset_early_embryonic_cleaned"
leiden_res = 'leiden_0_31'

raw_input = read_files(fileID)
run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")



###############
## Ependymal ##
###############

# parameters
fileID = "normal_subset_ependymal_cleaned"
leiden_res = 'leiden_0_15'

raw_input = read_files(fileID)

run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")



##########
## OPCs ##
##########

# parameters
fileID = "normal_subset_OPCs_oligos_cleaned_subset_OPCs"
leiden_res = 'leiden_0_1'

raw_input = read_files(fileID)

run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")



#################
## Neuroblasts ##
#################

# parameters
fileID = "normal_subset_neuroblasts"
leiden_res = 'leiden_0_15'

raw_input = read_files(fileID)

run_cl_rest(raw_input, fileID, leiden_res, pipeline="voom_limma")

