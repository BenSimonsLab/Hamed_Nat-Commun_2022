
import scvelo as scv
import scanpy as sc



def calc_adata_velo(adata):
    # adata = sc.read_h5ad('data/'+adata_file)
    
    scv.settings.figdir = 'output/'+fileID+'/'
    
    ldata = scv.read('data/set1_normal_splicing.h5ad')
    ldata.var_names_make_unique()
    ldata = ldata[adata.obs_names.values, :].copy()
    
    adata_velo = scv.utils.merge(adata, ldata)
    
    scv.pp.filter_and_normalize(adata_velo)
    scv.pp.moments(adata_velo)
    
    scv.tl.velocity(adata_velo, mode='stochastic')
    scv.tl.velocity_graph(adata_velo)
        
    return adata_velo


############
## Normal ##
############

fileID = "normal"
adata = scv.read('data/forebrain_'+fileID+'.h5ad')

adata_velo = calc_adata_velo(adata)

# adata_velo.obs.annot_coarse_leiden_0_7
# adata_velo.uns["annot_coarse_leiden_0_7_colors"] = ['#fbada7', '#b0ce66', '#66d9dc', '#ddb0ff']
# leiden_res = "annot_coarse_leiden_0_7"

adata_velo.obs.logreg_Zeisel_2018
# adata_velo.uns["logreg_Zeisel_2018_colors"] = ['#fbada7', '#b0ce66', '#66d9dc', '#ddb0ff']
adata_velo.uns["logreg_Zeisel_2018_colors"] = ['#f8766d', '#7cae00', '#00bfc4', '#c77cff']
leiden_res = "logreg_Zeisel_2018"


scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, alpha=0.8, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, alpha=0.8, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)


scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, alpha=1, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_alpha_1_.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, alpha=1, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_alpha_1_dpi_600.png', dpi=600)




##################
## Embryonic RG ##
##################

fileID = "normal_subset_embryonic_RG_cleaned"
adata = scv.read('data/forebrain_'+fileID+'.h5ad')


adata_velo_embryonic_RG = calc_adata_velo(adata)
adata_velo_embryonic_RG.obs.timepoint


#adata.uns["{var}_colors"] = ['red', 'grey']
adata_velo_embryonic_RG.uns["timepoint_colors"] = ['#8e6698', '#7faebb', '#afe396', '#afe396', '#afe396']
adata_velo_embryonic_RG.uns["timepoint_colors"] = ['#fdcc8a', '#fc8d59', '#d7301f', '#bdc9e1', '#74a9cf', '#0570b0']



leiden_res = "timepoint"
scv.pl.velocity_embedding_stream(adata_velo_embryonic_RG, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo_embryonic_RG, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)



######################
## Neuronal & Glial ##
######################

fileID = "normal_subset_neuronal_glial_cleaned"
adata = scv.read('data/forebrain_'+fileID+'.h5ad')


adata_velo_neuronal_glial = calc_adata_velo(adata)
adata_velo_neuronal_glial.obs.timepoint


#adata.uns["{var}_colors"] = ['red', 'grey']
adata_velo_neuronal_glial.uns["timepoint_colors"] = ['#8e6698', '#7faebb', '#afe396', '#afe396', '#afe396']
adata_velo_neuronal_glial.uns["timepoint_colors"] = ["#fdcc8a", "#fc8d59", "#d7301f", "#bdc9e1", "#74a9cf", "#0570b0", "#d9d9d9", "#bdbdbd", "#969696", "#636363", "#252525"]


leiden_res = "timepoint"
scv.pl.velocity_embedding_stream(adata_velo_neuronal_glial, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo_neuronal_glial, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)


leiden_res = "annot_leiden"
adata_velo_neuronal_glial.uns["annot_leiden_colors"] = ['#8e6698', '#7faebb', '#afe396', '#afe396', '#afe396']
scv.pl.velocity_embedding_stream(adata_velo_neuronal_glial, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo_neuronal_glial, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)




#####################
## Ependymal cells ##
#####################

fileID = "normal_subset_ependymal_cleaned"
adata = scv.read('data/forebrain_'+fileID+'.h5ad')


adata_velo_ependymal = calc_adata_velo(adata)
adata_velo_ependymal.obs.timepoint


adata_velo_ependymal.uns["timepoint_colors"] = ["#d7301f", "#bdc9e1", "#74a9cf",
                                                "#0570b0", "#d9d9d9", "#bdbdbd",
                                                "#969696", "#636363", "#252525"]


leiden_res = "timepoint"
scv.pl.velocity_embedding_stream(adata_velo_ependymal, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo_ependymal, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)



#################
## Neuroblasts ##
#################

fileID = "normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned"
adata = scv.read('data/forebrain_'+fileID+'.h5ad')

adata_velo_ependymal = calc_adata_velo(adata)
adata_velo_ependymal.obs.timepoint


adata_velo_ependymal.uns["annot_leiden_colors"] = ["#C19F70", "#0F4A9C", "#139992",
"#8870AD", "#C594BF", "#9E6762", "#CC7818"]





leiden_res = "timepoint"
scv.pl.velocity_embedding_stream(adata_velo_ependymal, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo_ependymal, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)

leiden_res = "annot_leiden"
scv.pl.velocity_embedding_stream(adata_velo_ependymal, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo_ependymal, color=leiden_res, title="", legend_loc ="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)



####################
## RG NSC lineage ##
####################

fileID = "normal_merged_RG_NSC_lineage"
adata = scv.read('data/forebrain_'+fileID+'.h5ad')

adata_velo = calc_adata_velo(adata)


# "Embryonic RG"="#532C8A",
# "Juvenile RG"="#3F84AA",
# "Active NSCs"="#65A83E",
# "aNSCs"="#65A83E",
# "Quiescent NSCs [1]"="#DABE99",
# "Quiescent NSCs (dorsal)"="#DABE99",
# "Quiescent NSCs"="#9F8A70",
# "qNSCs"="#9F8A70",
# "Quiescent NSCs (ventral)"="#635547",
# "Quiescent NSCs [2]"="#635547"

adata_velo.uns["annot_merged_coarse_colors"] = ["#532c8a", "#3f84aa", "#65a83e", "#9f8a70"]
leiden_res = "annot_merged_coarse"

scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)


adata_velo.uns["annot_merged_fine_colors"] = ["#532c8a", "#3f84aa", "#dabe99", "#635547", "#65a83e"]
leiden_res = "annot_merged_fine"

scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)



# version 2
fileID = "normal_subset_neuronal_glial_cleaned_subset_RG_NSC_lineage_cleaned"
adata = scv.read('data/forebrain_'+fileID+'.h5ad')

adata_velo = calc_adata_velo(adata)

# adata_velo.obs.annot_coarse_leiden_0_7
# adata_velo.uns["annot_coarse_leiden_0_7_colors"] = ['#fbada7', '#b0ce66', '#66d9dc', '#ddb0ff']
# leiden_res = "annot_coarse_leiden_0_7"

# adata_velo.obs.logreg_Zeisel_2018
# # adata_velo.uns["logreg_Zeisel_2018_colors"] = ['#fbada7', '#b0ce66', '#66d9dc', '#ddb0ff']
# adata_velo.uns["logreg_Zeisel_2018_colors"] = ['#f8766d', '#7cae00', '#00bfc4', '#c77cff']
leiden_res = "annot_leiden_0_1"


scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, alpha=0.8, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, alpha=0.8, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_dpi_600.png', dpi=600)


scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, alpha=1, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_alpha_1_.png', dpi=300)
scv.pl.velocity_embedding_stream(adata_velo, color=leiden_res, alpha=1, title="", legend_loc="right margin", basis='umap', save='stochastic_'+fileID+'_'+leiden_res+'_alpha_1_dpi_600.png', dpi=600)
