
library(URD)
library(readr)
library(dplyr)


calc_URD_dm <- function(fileID){
  metadata = read_csv(paste0("data/URD/metadata_", fileID, ".csv.gz"))
  metadata = as.data.frame(metadata)
  rownames(metadata) = metadata$cellID

  cmat = read.csv(file=paste0("data/URD/cmat_hvg_", fileID, ".csv.gz"), row.names=1, check.names=FALSE) %>%
                  as.matrix() %>%
                  t()

  objectURD = createURD(count.data = cmat,
                        meta = metadata,
                        min.cells=1,
                        min.counts=1,
                        min.genes=1)

  set.seed(23423)

  objectURD = calcDM(objectURD, knn=round(sqrt(dim(objectURD@meta)[1])), sigma.use="local")
  dm = objectURD@dm
  saveRDS(dm, file=paste0("data/URD/dm_URD_sigma_local_", fileID, ".rds"))

return()
}



calc_URD_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_10k")
calc_URD_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_20k")
calc_URD_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_30k")
calc_URD_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_40k")
calc_URD_dm(fileID = "normal_subset_neuronal_glial_cleaned_subsample_50k")