######################
## Merge DE results ##
######################

library(dplyr)
library(readr)
library(openxlsx)




###########
## qNSCs ##
###########

# fileID = "normal"
# file_path = paste0("output/", fileID, "/DE-voom_limma/")
# 
fpath = c("output/normal/DE-voom_limma/DE-voom_limma_normal_DE_leiden_0_4_cl_5_vs_14.csv.gz")
# 
# annot_mapping = c('4'='Thalamic eminence',
#                 )
# 
# 
list_tbl_de = list()
tbl_de = read_csv(fpath) %>% arrange(qval)
tableID = "qNSC_1 vs qNSC_2"
list_tbl_de[[tableID]] = tbl_de

write.xlsx(list_tbl_de, file = "data/supplementary_data/DE_voom_limma-qNSCs.xlsx")




################################################################################
################################################################################

####################
## Embryonic RGCs ##
####################

fileID = "normal_RGCs_cleaned_subset_early_embryonic_cleaned"
file_path = paste0("output/", fileID, "/DE-voom_limma/")

file_paths = list.files(file_path, full.names=T)

annot_mapping = c('0'='Epithalamus',
                  '1'='Gliogenic (1)',
                  '2'='Ganglionic eminence',
                  '3'='Cortical pallium',
                  '4'='Thalamic eminence',
                  '5'='Gliogenic (2)',
                  '6'='Cortical hem',
                  '7'='Pretectum',
                  '8'='Subthalamic nucleus')


list_tbl_de = list()

for (fpath in file_paths){
  tbl_de = read_csv(fpath) %>% arrange(qval)
  
  clusterID = strsplit(strsplit(fpath, "_cl_")[[1]][2], "_")[[1]][1]
  tableID = paste0(annot_mapping[clusterID], " vs rest")
  
  list_tbl_de[[tableID]] = tbl_de
}


write.xlsx(list_tbl_de, file = "data/supplementary_data/DE_voom_limma-embryonic_RG.xlsx")


##########
## RGCs ##
##########

# full
fileID = "normal_RGCs_cleaned"
file_path = paste0("output/", fileID, "/DE-voom_limma/")

file_paths = list.files(file_path, full.names=T)

annot_mapping = c('0'='qNSCs (1)',
                  '1'='Embryonic RG',
                  '2'='Juvenile RG',
                  '3'='qNSCs (2)',
                  '4'='aNSCs')


list_tbl_de = list()

for (fpath in file_paths){
  tbl_de = read_csv(fpath) %>% arrange(qval)
  
  clusterID = strsplit(strsplit(fpath, "_cl_")[[1]][2], "_")[[1]][1]
  tableID = paste0(annot_mapping[clusterID], " vs rest")
  
  list_tbl_de[[tableID]] = tbl_de
}


# developing
fileID = "normal_RGCs_cleaned_group_developing"
file_path = paste0("output/", fileID, "/DE-voom_limma/")

file_paths = list.files(file_path, full.names=T)

for (fpath in file_paths){
  tbl_de = read_csv(fpath) %>% arrange(qval)
  
  clusterID = strsplit(strsplit(fpath, "_cl_")[[1]][2], "_")[[1]][1]
  tableID = paste0(annot_mapping[clusterID], " vs dev.")
  
  list_tbl_de[[tableID]] = tbl_de
}


# adult
fileID = "normal_RGCs_cleaned_group_adult"
file_path = paste0("output/", fileID, "/DE-voom_limma/")

file_paths = list.files(file_path, full.names=T)

for (fpath in file_paths){
  tbl_de = read_csv(fpath) %>% arrange(qval)
  
  clusterID = strsplit(strsplit(fpath, "_cl_")[[1]][2], "_")[[1]][1]
  tableID = paste0(annot_mapping[clusterID], " vs adult")
  
  list_tbl_de[[tableID]] = tbl_de
}


write.xlsx(list_tbl_de, file = "data/supplementary_data/DE_voom_limma-RG.xlsx")


#################
## Neuroblasts ##
#################

fileID = "normal_subset_neuroblasts_nodebris"
file_path = paste0("output/", fileID, "/DE-voom_limma/")

file_paths = list.files(file_path, full.names=T)

annot_mapping = c('0'='GE NBs (1)',
                  '1'='Hippocampal NBs (1)',
                  '2'='GE NBs (2)',
                  '3'='Hippocampal NBs (2)',
                  '4'='Early EmDienNBs',
                  '5'='Early EmCorNBs',
                  '6'='Early EmSthNBs',
                  '7'='Early EmNBs')


list_tbl_de = list()

for (fpath in file_paths){
  tbl_de = read_csv(fpath) %>% arrange(qval)
  
  clusterID = strsplit(strsplit(fpath, "_cl_")[[1]][2], "_")[[1]][1]
  tableID = paste0(annot_mapping[clusterID], " vs rest")
  
  list_tbl_de[[tableID]] = tbl_de
}


write.xlsx(list_tbl_de, file = "data/supplementary_data/DE_voom_limma-neuroblasts.xlsx")


######################
## Ependymal Cells ##
######################

fileID = "normal_subset_ependymal_cleaned"
file_path = paste0("output/", fileID, "/DE-voom_limma/")

file_paths = list.files(file_path, full.names=T)

annot_mapping = c('0'='Adult ependymal',
                  '1'='Juvenile ependymal',
                  '2'='Jvnile ependymal (cycl)')


list_tbl_de = list()

for (fpath in file_paths){
  tbl_de = read_csv(fpath) %>% arrange(qval)
  
  clusterID = strsplit(strsplit(fpath, "_cl_")[[1]][2], "_")[[1]][1]
  tableID = paste0(annot_mapping[clusterID], " vs rest")
  
  list_tbl_de[[tableID]] = tbl_de
}


write.xlsx(list_tbl_de, file = "data/supplementary_data/DE_voom_limma-ependymal_cells.xlsx")


##########
## OPCs ##
##########

fileID = "normal_subset_OPCs_oligos_cleaned_subset_OPCs"
file_path = paste0("output/", fileID, "/DE-voom_limma/")

file_paths = list.files(file_path, full.names=T)

annot_mapping = c('0'='Juvenile OPCs',
                  '1'='Juvenile OPCs (cycl)',
                  '2'='Adult OPCs')


list_tbl_de = list()

for (fpath in file_paths){
  tbl_de = read_csv(fpath) %>% arrange(qval)
  
  clusterID = strsplit(strsplit(fpath, "_cl_")[[1]][2], "_")[[1]][1]
  tableID = paste0(annot_mapping[clusterID], " vs rest")
  
  list_tbl_de[[tableID]] = tbl_de
}


write.xlsx(list_tbl_de, file = "data/supplementary_data/DE_voom_limma-OPCs.xlsx")
