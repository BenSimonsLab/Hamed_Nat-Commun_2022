library(destiny)
library(Biobase)
library(dplyr)



calc_dm <- function(fileID=""){
  cmat_hvg = read.csv(file=paste0("data/URD/cmat_hvg_", fileID, ".csv.gz"), row.names=1, check.names=FALSE) %>%
                  as.matrix() %>%
                  t()

  # exprSet = as.ExpressionSet(cmat_hvg)
  exprSet = ExpressionSet(assayData=cmat_hvg)

  set.seed(23423)
  dm = DiffusionMap(exprSet, sigma="local")
  saveRDS(dm, file=paste0("data/URD/dm_sigma_local_", fileID, ".rds"))
  
  return()  
}



calc_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_10k")
calc_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_20k")
calc_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_30k")

calc_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_40k")
calc_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_50k")