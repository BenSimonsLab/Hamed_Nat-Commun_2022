#######################
## Monocle 2 DDRTree ##
#######################

library(monocle)
library(Matrix) # use sparse matrices
library(VGAM)

# library(ggplot2)
# library(viridis)
# library(patchwork)
# library(ggrepel)
# 
# library(readr)
# library(dplyr)
# library(tibble)
# 
# library(knitr)
# library(kableExtra)
# 
# library(reticulate)
# use_python("/usr/bin/python3.6")
# source("code/core_settings.R")



runMonocle_OLD <- function(counts, pd, fd, batch=FALSE){
  SPLT <- newCellDataSet(as(as.matrix(counts), "sparseMatrix"),
                  phenoData = pd,
                  featureData = fd,
                  expressionFamily=negbinomial.size())

  # processing for Monocle
  # adapted from Monocle vignette
  SPLT <- estimateSizeFactors(SPLT)
  SPLT <- estimateDispersions(SPLT)

  # select ordering genes based on variance
  # NB there are several alternatives to select genes
  disp_table <- dispersionTable(SPLT)
  ordering_genes <- subset(disp_table,
                    mean_expression >= 0.5 &
                    dispersion_empirical >= 1 * dispersion_fit)$gene_id
  # ordering_genes <- rownames(counts)
  
  SPLT <- setOrderingFilter(SPLT, ordering_genes)
  # reduce data dimensionality
  if (batch){
    SPLT <- reduceDimension(SPLT, max_components = 2, reduction_method = 'DDRTree', residualModelFormulaStr='~batch')
  } else {
    SPLT <- reduceDimension(SPLT, max_components = 2, method = 'DDRTree')
  }
  # order cells along the trajectory
  SPLT <- orderCells(SPLT)

  return(SPLT)
}


runMonocle <- function(counts, pd, fd, batch=FALSE){
  SPLT <- newCellDataSet(as(as.matrix(counts), "sparseMatrix"),
                  phenoData = pd,
                  featureData = fd,
                  expressionFamily=negbinomial.size())

  # processing for Monocle
  # adapted from Monocle vignette
  SPLT <- estimateSizeFactors(SPLT)
  SPLT <- estimateDispersions(SPLT)

  # select ordering genes based on variance
  # NB there are several alternatives to select genes
  disp_table <- dispersionTable(SPLT)
  ordering_genes <- subset(disp_table,
                    mean_expression >= 0.5 &
                    dispersion_empirical >= 1 * dispersion_fit)$gene_id
  # ordering_genes <- rownames(counts)
  
  SPLT <- setOrderingFilter(SPLT, ordering_genes)
  # reduce data dimensionality
  if (batch){
    SPLT <- reduceDimension(SPLT, max_components = 2, reduction_method = 'DDRTree', residualModelFormulaStr='~batch')
  } else {
    SPLT <- reduceDimension(SPLT)
  }
  # order cells along the trajectory
  SPLT <- orderCells(SPLT)

  return(SPLT)
}


fileID = "normal_subset_OPCs_oligos_cleaned"
fileID = "normal_subset_ependymal_cleaned"
fileID = "normal_RGCs_cleaned"
fileID = "normal_merged_subset_neuronal_glial_subset_RG_OPC_NSC_lineage_cleaned"

library(reticulate)



# use raw count matrix from scenic input
# counts = t(read.csv(paste0("data/scenic/cmat_", fileID, ".csv.gz"), check.names=F, row.names=1))

pd = new("AnnotatedDataFrame", data = data.frame(cellids=colnames(counts), row.names=colnames(counts)))
fd = new("AnnotatedDataFrame", data = data.frame(gene_short_name=rownames(counts), row.names=rownames(counts)))


SPLT = runMonocle(counts, pd, fd)

# save results
library(tidyverse)

features = attributes(SPLT)
# double check whether this is what is plotted by monocle 2 there exit other reduced dims as well "whitening"
dims_monocle = reducedDimS(SPLT)

results_monocle = tibble(cellID=SPLT$cellids,
                         monocle_state=SPLT$State,
                         monocle_state_2 = phenoData(SPLT)$State,
                         monocle_size_factor=SPLT$Size_Factor,
                         monocle_pseudotime=SPLT$Pseudotime,
                         monocle_dim1=dims_monocle[1,],
                         monocle_dim2=dims_monocle[2,])

write_csv(results_monocle, path=paste0("data/monocle2/pt_monocle2_", fileID, ".csv.gz"))
