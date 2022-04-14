library(readr)
library(dplyr)

fpaths = list.files(path="data/URD/normal_subset_neuronal_glial_cleaned/", pattern="^tip_.*", full.names=T)

tbl_tips = tibble(cellID = character(),
                  clusterID = character())

for (fpath in fpaths){
    tbl_tip = read_csv(fpath, col_types="c") %>%
            mutate(clusterID = gsub(".csv", "", strsplit(basename(fpath), split="_")[[1]][2]))
    tbl_tips = bind_rows(tbl_tips, tbl_tip)
}


# table(table(tbl_tips$cellID))

# table_results = table(tbl_tips$cellID)
# duplicates = names(table_results)[table_results == 2]
# tbl_tips %>% filter(cellID %in% duplicates) %>% arrange(cellID)


write_csv(tbl_tips, "data/URD/normal_subset_neuronal_glial_cleaned/tips_merged.csv")