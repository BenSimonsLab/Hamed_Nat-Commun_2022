####################
## Add annotation ##
####################

import numpy as np
import pandas as pd
import scanpy as sc


def add_annotation(fileID, leiden_res, annot_dict, annot_key="annot_"):
    # add annotation based on cluster numbers
    file_path = "data/forebrain_"+fileID+".h5ad"
    adata = sc.read_h5ad(file_path)
    
    adata.obs[annot_key+leiden_res] = adata.obs[leiden_res].values.rename_categories(annot_dict)
    adata.obs['annot_leiden'] = adata.obs[leiden_res].values.rename_categories(annot_dict)
    
    # save results
    adata.write('data/forebrain_'+fileID+'.h5ad')
    
    return


annot_complete = {'0': 'Adult microglia [1]',
                  '1': 'GE NBs [1]',
                  '2': 'Hippocampal NBs [2]',
                  '3': 'Juvenile RG & TAPs',
                  '4': 'Oligodendrocytes [2]',
                  '5': 'Quiescent NSCs (dorsal)',
                  '6': 'Embryonic RG',
                  '7': 'Juvenile microglia',
                  '8': 'Early EmNBs',
                  '9': 'Gliogenic precursors (APCs)',
                  '10': 'OPCs',
                  '11': 'Gliogenic precursors (aNSCs)',
                  '12': 'ImStNeurons',
                  '13': 'GE NBs [2]',
                  '14': 'Quiescent NSCs (ventral)',
                  '15': 'Hippocampal NBs [1]',
                  '16': 'Ependymal cells',
                  '17': 'Adult microglia [2]',
                  '18': 'VLMC',
                  '19': 'Macrophages',
                  '20': 'T-cells',
                  '21': 'GABAergic INs',
                  '22': 'Oligodendrocytes [1]',
                  '23': 'Erythrocytes',
                  '24': 'Choroid plexus epithelia',
                  '25': 'Endothelial cells',
                  '26': 'PreM-ODCs',
                  '27': 'Myeloid-DSCs',
                  '28': '[unclear/debris]'}

add_annotation("normal", "leiden_0_4", annot_complete)



annot_subset_neuronal_glial = {'0': 'GE NBs [1]',
                               '1': 'Hippocampal NBs [2]',
                               '2': 'Juvenile RG & TAPs',
                               '3': 'Oligodendrocytes [2]',
                               '4': 'Quiescent NSCs (dorsal)',
                               '5': 'Embryonic RG',
                               '6': 'OPCs',
                               '7': 'Gliogenic precursors (APCs)',
                               '8': 'Gliogenic precursors (aNSCs)',
                               '9': 'GE NBs [2]',
                               '10': 'Early EmNBs [1]',
                               '11': 'Quiescent NSCs (ventral)',
                               '12': 'Hippocampal NBs [1]',
                               '13': 'ImStNeurons',
                               '14': 'Ependymal cells',
                               '15': 'Early EmNBs [2]',
                               '16': 'Oligodendrocytes [1]',
                               '17': 'GABAergic INs',
                               '18': 'PreM-ODCs'}

add_annotation("normal_subset_neuronal_glial_cleaned", "leiden_0_4", annot_subset_neuronal_glial)




annot_embryonic_RG = {'0': 'Epithalamus',
                      '1': 'Gliogenic (cortical)',
                      '2': 'Ganglionic eminence',
                      '3': 'Cortical pallium',
                      '4': 'Thalamic eminence',
                      '5': 'Gliogenic (ganglionic)',
                      '6': 'Pretectum',
                      '7': 'Subthalamic nucleus',
                      '8': 'Cortical hem'}

add_annotation("normal_subset_embryonic_RG_cleaned", "leiden_0_33", annot_embryonic_RG)


annot_NBs = {'0': 'GE NBs [1]',
             '1': 'Hippocampal NBs [2]',
             '2': 'GE NBs [2]',
             '3': 'EmDienNBs',
             '4': 'Hippocampal NBs [1]',
             '5': 'EmCorNBs',
             '6': 'EmSthNBs'}

add_annotation("normal_subset_neuronal_glial_cleaned_subset_NBs_cleaned", "leiden_0_13", annot_NBs)


annot_ependymal = {'0': 'Adult ependymal cells',
                   '1': 'Juvenile ependymal cells',
                   '2': 'Juvenile ependymal cells [cycling]'}

add_annotation("normal_subset_ependymal_cleaned", "leiden_0_15", annot_ependymal)


annot_OPCs = {'0': 'Juvenile OPCs',
                   '1': 'Adult OPCs',
                   '2': 'Juvenile OPCs [cycling]'}

add_annotation("normal_subset_OPCs", "leiden_0_1", annot_OPCs)








annot_E12_5_subset_neuronal_glial = {'0': 'RG',
                                    '1': 'EEmCorNBs/CRCs',
                                    '2': 'EEmThNBs',
                                    '3': 'GE NBs',
                                    '4': 'GABAergic interneurons',
                                    '5': 'EEmSthNBs',
                                    '6': 'EEmDienNBs',
                                    '7': 'EEmPyrNBs',
                                    '8': '[TBD]'}
add_annotation("normal_E12_5_subset_neuronal_glial", "leiden_0_1", annot_E12_5_subset_neuronal_glial)


annot_E14_5_subset_neuronal_glial = {'0': '',
                                    '1': '',
                                    '2': '',
                                    '3': '',
                                    '4': '',
                                    '5': '',
                                    '6': '',
                                    '7': '',
                                    '8': ''}
add_annotation("normal_E14_5_subset_neuronal_glial", "leiden_0_", annot_E14_5_subset_neuronal_glial)


annot_E16_5_subset_neuronal_glial = {'0': '',
                                    '1': '',
                                    '2': '',
                                    '3': '',
                                    '4': '',
                                    '5': '',
                                    '6': '',
                                    '7': '',
                                    '8': ''}
add_annotation("normal_E16_5_subset_neuronal_glial", "leiden_0_", annot_E16_5_subset_neuronal_glial)


annot_P0_subset_neuronal_glial = {'0': 'Hippocampal NBs',
                                    '1': 'GE NBs [1]',
                                    '2': 'Gliogenic precursors (APCs)',
                                    '3': 'Gliogenic precursors (aNSCs)',
                                    '4': 'OPCs',
                                    '5': 'RG & TAPs',
                                    '6': 'ImStNeurons',
                                    '7': '[unclear/debris]',
                                    '8': 'GE NBs [2]'}
add_annotation("normal_P0_subset_neuronal_glial", "leiden_0_2", annot_P0_subset_neuronal_glial)


annot_P2_subset_neuronal_glial = {'0': 'OPCs',
                                    '1': 'Gliogenic precursors (APCs)',
                                    '2': 'Juvenile RG & TAPs',
                                    '3': 'GE NBs & ImStNeurons',
                                    '4': 'Gliogenic precursors (aNSCs)',
                                    '5': 'Hippocampal NBs',
                                    '6': 'Ependymal cells',
                                    '7': 'OPCs [cycling]',
                                    '8': 'ImPreMDs'}
add_annotation("normal_P2_subset_neuronal_glial", "leiden_0_2", annot_P2_subset_neuronal_glial)


annot_P7_subset_neuronal_glial = {'0': 'Gliogenic precursors (APCs)',
                                    '1': 'GE NBs',
                                    '2': 'OPCs',
                                    '3': 'TAPs',
                                    '4': 'RG',
                                    '5': 'Gliogenic precursors (aNSCs)',
                                    '6': 'Hippocampal NBs',
                                    '7': 'Ependymal cells',
                                    '8': 'ImPreMDs [+immune debris]',
                                    '9': '[unclear/debris]'}
add_annotation("normal_P7_subset_neuronal_glial", "leiden_0_3", annot_P7_subset_neuronal_glial)



annot_P13_subset_neuronal_glial = {'0': 'GE NBs',
                                    '1': 'Quiescent NSCs [2]',
                                    '2': 'Juvenile RG & TAPs',
                                    '3': 'OPCs',
                                    '4': 'Ependymal cells',
                                    '5': 'Quiescent NSCs [1]',
                                    '6': 'ImPreMDs',
                                    '7': 'Gliogenic precursors',
                                    '8': 'Hippocampal NBs',
                                    '9': 'TBD [Gliogenic precursors?]',
                                    '10': 'Debris [microglia]'}
add_annotation("normal_P13_subset_neuronal_glial", "leiden_0_3", annot_P13_subset_neuronal_glial)



annot_P19_subset_neuronal_glial = {'0': 'GE NBs',
                                    '1': 'Quiescent NSCs [1]',
                                    '2': 'Quiescent NSCs [2]',
                                    '3': 'Oligodendrocytes [1]',
                                    '4': 'Oligodendrocytes [2]',
                                    '5': 'TAPs',
                                    '6': 'Ependymal cells',
                                    '7': 'Gliogenic precursors',
                                    '8': 'ImPreMDs',
                                    '9': 'OPCs',
                                    '10': 'Hippocampal NBs'}
add_annotation("normal_P19_subset_neuronal_glial", "leiden_0_4", annot_P19_subset_neuronal_glial)


annot_P39_subset_neuronal_glial = {'0': 'Oligodendrocytes',
                                    '1': 'Quiescent NSCs [1]',
                                    '2': 'Quiescent NSCs [2]',
                                    '3': 'GE NBs',
                                    '4': 'Ependymal cells',
                                    '5': 'OPCs'}
add_annotation("normal_P39_subset_neuronal_glial", "leiden_0_1", annot_P39_subset_neuronal_glial)


annot_P111_subset_neuronal_glial = {'0': 'Oligodendrocytes',
                                    '1': 'Quiescent NSCs [1]',
                                    '2': 'Quiescent NSCs [2]',
                                    '3': 'GE NBs',
                                    '4': 'TAPs',
                                    '5': 'Ependymal cells',
                                    '6': 'OPCs'}
add_annotation("normal_P111_subset_neuronal_glial", "leiden_0_2", annot_P111_subset_neuronal_glial)


annot_P365_subset_neuronal_glial = {'0': 'Oligodendrocytes',
                                    '1': 'Quiescent NSCs [1]',
                                    '2': 'Quiescent NSCs [2]',
                                    '3': 'GE NBs',
                                    '4': 'OPCs'}
add_annotation("normal_P365_subset_neuronal_glial", "leiden_0_1", annot_P365_subset_neuronal_glial)


# 2020-11-04
annot_merged_subset_neuronal_glial = {'0': 'Quiescent NSCs [1]',
                                      '1': 'Juvenile RG & TAPs',
                                      '2': 'Embryonic RG',
                                      '3': 'OPCs',
                                      '4': 'Gliogenic precursors ()',
                                      '5': 'Gliogenic precursors (aNSCs)',
                                      '6': 'Quiescent NSCs [2]',
                                      '7': 'ImPreMDs'}
add_annotation("normal_merged_subset_neuronal_glial_subset_RG_OPC_NSC_lineage_cleaned", "leiden_0_2", annot_merged_subset_neuronal_glial)



# 2020-11-12
annot_map = {'0': 'Quiescent NSCs [1]',
             '1': 'Embryonic RG',
             '2': 'Juvenile RG',
             '3': 'Active NSCs',
             '4': 'Quiescent NSCs [2]'}
add_annotation("normal_subset_neuronal_glial_cleaned_subset_RG_NSC_lineage_v2_cleaned", "leiden_0_1", annot_map)




###########################################
## Add annotation based on other subsets ##
###########################################



fileID = "normal_merged_RG_NSC_lineage"
adata = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

# subsets this was merged from:
fileID = "normal_subset_embryonic_RG_cleaned"
adata_embryonic_RG = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_juvenile_RG_TAPs_subset_juvenile_RG"
adata_juvenile_RG = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_aNSCs"
adata_aNSCs = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_qNSCs"
adata_qNSCs = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")


annot_merged_coarse = [0]*len(adata)
annot_merged_fine = [0]*len(adata)

for i in range(len(adata)):
    cellID = adata.obs_names.values[i]
    if cellID in adata_embryonic_RG.obs_names.values:
        annot_merged_coarse[i] = "Embryonic RG"
        annot_merged_fine[i] = "Embryonic RG"
    elif cellID in adata_juvenile_RG.obs_names.values:
        annot_merged_coarse[i] = "Juvenile RG"
        annot_merged_fine[i] = "Juvenile RG"
    elif cellID in adata_aNSCs.obs_names.values:
        annot_merged_coarse[i] = "Active NSCs"
        annot_merged_fine[i] = "Active NSCs"
    elif cellID in adata_qNSCs.obs_names.values:
        annot_merged_coarse[i] = "Quiescent NSCs"
        annot_merged_fine[i] = adata_qNSCs.obs.annot_superset[cellID]
    else:
        annot_merged_coarse[i] = "NA"
        annot_merged_fine[i] = "NA"


adata.obs['annot_merged_coarse'] = annot_merged_coarse
adata.obs['annot_merged_fine'] = annot_merged_fine

adata.obs['annot_merged_coarse_num'] = pd.Series(annot_merged_coarse, index=adata.obs_names.values, dtype="category").cat.codes.astype('category')
adata.obs['annot_merged_fine_num'] = pd.Series(annot_merged_fine, index=adata.obs_names.values, dtype="category").cat.codes.astype('category')

# save results
fileID = "normal_merged_RG_NSC_lineage"
adata.write('data/forebrain_'+fileID+'.h5ad')



#############################


fileID = "normal_merged_RG_aNSC_lineage"
adata = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

# subsets this was merged from:
fileID = "normal_subset_embryonic_RG_cleaned"
adata_embryonic_RG = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_juvenile_RG_TAPs_subset_juvenile_RG"
adata_juvenile_RG = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_aNSCs"
adata_aNSCs = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")


annot_merged = [0]*len(adata)

for i in range(len(adata)):
    cellID = adata.obs_names.values[i]
    if cellID in adata_embryonic_RG.obs_names.values:
        annot_merged[i] = "Embryonic RG"
    elif cellID in adata_juvenile_RG.obs_names.values:
        annot_merged[i] = "Juvenile RG"
    elif cellID in adata_aNSCs.obs_names.values:
        annot_merged[i] = "Active NSCs"
    else:
        annot_merged[i] = "NA"


adata.obs['annot_merged'] = annot_merged

adata.obs['annot_merged_num'] = pd.Series(annot_merged, dtype="category").cat.codes.values

# save results
fileID = "normal_merged_RG_aNSC_lineage"
adata.write('data/forebrain_'+fileID+'.h5ad')




# # # # # # # #
# qNSC subset #
# # # # # # # #

fileID = "normal_subset_qNSCs"
adata_qNSCs = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal"
adata_normal = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

annot_superset = adata_normal.obs["annot_leiden"][adata_qNSCs.obs_names.values]

adata_qNSCs.obs['annot_superset'] = annot_superset
adata_qNSCs.obs['annot_superset_num'] = annot_superset.cat.codes.astype('category')

# save results
fileID = "normal_subset_qNSCs"
adata_qNSCs.write('data/forebrain_'+fileID+'.h5ad')



# # # # # # # # # # # # # # # #
# Embryonic RG w/o gliogenic  #
# # # # # # # # # # # # # # # #


fileID = "normal_subset_embryonic_RG_cleaned_wo_gliogenic"
adata_qNSCs = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_embryonic_RG_cleaned"
adata_normal = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

annot_superset = adata_normal.obs["annot_leiden"][adata_qNSCs.obs_names.values]

adata_qNSCs.obs['annot_superset'] = annot_superset
adata_qNSCs.obs['annot_superset_num'] = annot_superset.cat.codes.astype('category')

# save results
fileID = "normal_subset_embryonic_RG_cleaned_wo_gliogenic"
adata_qNSCs.write('data/forebrain_'+fileID+'.h5ad')



###############################################


fileID = "normal_subset_neuronal_glial_cleaned_subset_qNSCs"
adata_qNSCs = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_neuronal_glial_cleaned"
adata_neuronal_glial = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

annot_superset = adata_neuronal_glial.obs["annot_leiden"][adata_qNSCs.obs_names.values]

adata_qNSCs.obs['annot_superset'] = annot_superset
adata_qNSCs.obs['annot_superset_num'] = annot_superset.cat.codes.astype('category')

# save results
fileID = "normal_subset_neuronal_glial_cleaned_subset_qNSCs"
adata_qNSCs.write('data/forebrain_'+fileID+'.h5ad')


# # # # # # # # # # # # #
# RG NSC lineage subset #
# # # # # # # # # # # # #


fileID = "normal_subset_neuronal_glial_cleaned_subset_RG_NSC_lineage"
adata_RG_NSC = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_neuronal_glial_cleaned"
adata_neuronal_glial = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

annot_superset = adata_neuronal_glial.obs["annot_leiden"][adata_RG_NSC.obs_names.values]

adata_RG_NSC.obs['annot_superset'] = annot_superset
adata_RG_NSC.obs['annot_superset_num'] = annot_superset.cat.codes.astype('category')

# save results
fileID = "normal_subset_neuronal_glial_cleaned_subset_RG_NSC_lineage"
adata_RG_NSC.write('data/forebrain_'+fileID+'.h5ad')


# # # # # # # # # # # # # # # # # # # # # # # # # #
# RG lineage subset from neuronal & glial subset  #
# # # # # # # # # # # # # # # # # # # # # # # # # #


fileID = "normal_subset_neuronal_glial_cleaned_subset_RG"
adata_RG_NSC = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal_subset_neuronal_glial_cleaned"
adata_neuronal_glial = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

annot_superset = adata_neuronal_glial.obs["annot_leiden"][adata_RG_NSC.obs_names.values]

adata_RG_NSC.obs['annot_superset'] = annot_superset
adata_RG_NSC.obs['annot_superset_num'] = annot_superset.cat.codes.astype('category')

# save results
fileID = "normal_subset_neuronal_glial_cleaned_subset_RG"
adata_RG_NSC.write('data/forebrain_'+fileID+'.h5ad')


# # # # # # # # # # # 
# RG lineage subset #
# # # # # # # # # # #


fileID = "normal_subset_RG"
adata_RG_NSC = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

fileID = "normal"
adata_neuronal_glial = sc.read_h5ad("data/forebrain_"+fileID+".h5ad")

annot_superset = adata_neuronal_glial.obs["annot_leiden"][adata_RG_NSC.obs_names.values]

adata_RG_NSC.obs['annot_superset'] = annot_superset
adata_RG_NSC.obs['annot_superset_num'] = annot_superset.cat.codes.astype('category')

# save results
fileID = "normal_subset_RG"
adata_RG_NSC.write('data/forebrain_'+fileID+'.h5ad')


#######################
## Add logreg Zeisel ##
#######################

fileID = "normal"
file_path = "data/forebrain_"+fileID+".h5ad"
adata = sc.read_h5ad(file_path)

df_zeisel = pd.read_csv("data/logreg/logreg_Zeisel_2018_normal.csv.gz", index_col=0)

adata.obs["logreg_Zeisel_2018"] = df_zeisel["prediction_TaxonomyRank1"].loc[adata.obs_names.values].values

# save results
adata.write('data/forebrain_'+fileID+'.h5ad')