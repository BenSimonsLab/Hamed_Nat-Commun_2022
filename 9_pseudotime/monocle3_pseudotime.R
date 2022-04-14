############################
## Monocle 3 - Pseudotime ##
############################

library(ggplot2)
library(viridis)
library(patchwork)
library(ggrepel)

library(tibble)
library(dplyr)

library(reticulate)
use_python("/usr/bin/python3")


theme_publication = theme_bw() +
                      theme(text=element_text(size=5), axis.text=element_text(size=5), axis.title=element_text(size=6), plot.title=element_text(size=6, hjust=0.5)) +
                      theme(strip.text= element_text(size=6)) +
                      theme(legend.title = element_text(size=5), legend.text = element_text(size=5)) +
                      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
                      theme(plot.tag=element_text(size=8, face="bold")) +
                      # remove facet box
                      theme(strip.background = element_blank())





fileID = "normal_subset_neuronal_glial_cleaned"

adata_path = paste0('data/forebrain_', fileID, '.h5ad')

 

py_run_string("import numpy as np")
py_run_string("import anndata as ad")

py_run_string("adata = ad.read_h5ad(r.adata_path)")


metadata_scanpy = as_tibble(py$adata$obs, rownames="cellID") %>%
                    rename(mt_frac = percent_mito) %>%
                    mutate(UMAP_dim1 = py$adata$obsm[['X_umap']][,1],
                           UMAP_dim2 = py$adata$obsm[['X_umap']][,2])



obs = as.matrix(metadata_scanpy)
rownames(obs) = metadata_scanpy$cellID
# umap = adata.obsm['X_umap']
umap = py$adata$obsm[['X_umap']]


library(monocle3)

X = matrix(0,nrow=1,ncol=dim(obs)[1])
colnames(X) = metadata_scanpy$cellID
mono = new_cell_data_set(X, cell_metadata=obs)


umap = as.matrix(umap)
# rownames(umap) = rownames(obs)
rownames(umap) = metadata_scanpy$cellID
# mono@reducedDims = SimpleList(UMAP=umap)
# mono@reduce_dim_aux = SimpleList(UMAP=umap)

# mono@int_colData@listData$reducedDims@listData[["UMAP"]] = umap
mono@int_colData@listData[["reducedDims"]][["UMAP"]] = umap




# cds_from_seurat@reduce_dim_aux@listData[["UMAP"]] 
# plot_cells(mono)



mono = cluster_cells(mono, reduction_method='UMAP')
# mono = learn_graph(mono)
mono = learn_graph(mono, use_partition = FALSE)



plot_cells(mono,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)



# mono = order_cells(mono)

# iroot = adata.obs_names[adata.uns['iroot']]
py_run_string("iroot = adata.obs_names[adata.uns['iroot']]")
iroot = py$iroot
iroot = names(which.max(umap[,1]))
mono = order_cells(mono, root_cells=iroot)



# plt_monocle3 = plot_cells(mono,
#                          color_cells_by = "pseudotime",
#                          label_cell_groups=FALSE,
#                          label_leaves=FALSE,
#                          label_branch_points=FALSE,
#                          graph_label_size=1.5)


saveRDS(mono, "data/monocle3/monocle3_neuronal_glial.rds")

plt_monocle3 = plot_cells(mono,
                         color_cells_by = "pseudotime",
                         label_cell_groups=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=FALSE,
                         rasterize = TRUE,
                         graph_label_size=1.5) +
                  theme_publication +
                  theme(axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank()) +
                  labs(colour="Monocle 3 /n pseudotime")


saveRDS(plt_monocle3, "data/monocle3/monocle3_neuronal_glial_ggplot.rds")


ggsave("test_monocle_3.pdf", plot = plt_monocle3, height=8, width=11, units="cm")

pseudotime = mono@principal_graph_aux[['UMAP']]$pseudotime
