##################################
## Cyclone Cell Cycle Inference ##
##################################

library(reticulate)
library(scran)
library(readr)
library(dplyr)


set.seed(73452)

sc = import("scanpy", convert = FALSE)


sampleIDs = list.dirs("data/samples/", full.names=FALSE, recursive=FALSE)

# get samplesIDs which have been run already
cyclone_sampleIDs = list.files("data/cyclone/", full.names=FALSE) %>%
                    {gsub("^cyclone_", "", .)} %>%
                    {gsub(".csv.gz$", "", .)}

sampleIDs = sampleIDs[!(sampleIDs %in% cyclone_sampleIDs)]



for (sampleID in sampleIDs){
  print(sampleID)
  adata = sc$read_10x_h5(paste0("data/samples/", sampleID, "/outs/filtered_feature_bc_matrix.h5"))
  adata$var_names_make_unique()

  cmat = t(as.matrix(adata$X$toarray()))
  colnames(cmat) = as.vector(adata$obs_names$values)
  rownames(cmat) = as.vector(adata$var['gene_ids']$values)

  sce = SingleCellExperiment(assays=list(counts=cmat))

  # mus musculus
  mm.pairs = readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

  assigned = cyclone(sce, pairs=mm.pairs)
  cyclone_tbl = tibble(cellID = paste0(colnames(sce), "-", sampleID),
                       cyclone = assigned$phases,
                       cyclone_G1 = assigned$scores$G1,
                       cyclone_S = assigned$scores$S,
                       cyclone_G2M = assigned$scores$G2M)

  write_csv(cyclone_tbl, paste0("data/cyclone/cyclone_", sampleID, ".csv.gz"))
}


###################
## Merge Results ##
###################

library(readr)
library(dplyr)


sampleIDs = list.files("data/cyclone/", full.names=FALSE, pattern = "^cyclone_HAM*.*\\.csv\\.gz$") %>%
                    {gsub("^cyclone_", "", .)} %>%
                    {gsub(".csv.gz$", "", .)}

cyclone_full = read_csv(paste0("data/cyclone/cyclone_", sampleIDs[1], ".csv.gz"))

for (i in 2:length(sampleIDs)){
  cyclone_sample = read_csv(paste0("data/cyclone/cyclone_", sampleIDs[i], ".csv.gz"))
  cyclone_full = bind_rows(cyclone_full, cyclone_sample)
}

write_csv(cyclone_full, "data/cyclone/cyclone_merged.csv.gz")
