library(tidyverse)



## diffmeth sites files from RnBeads, named
tibble_list <- lapply(list.files(pattern = "^br_(.*)site(.*)"), FUN = readr::read_csv, col_names = T)

names_sites <- str_remove(list.files(pattern = "^br_(.*)site(.*)"), "_diff(.*)")



# Select only sensible CpG sites and add ID for Venn extraction
filtered_tibble_list <- tibble_list %>% 
  map(filter, (diffmeth.p.adj.fdr < 0.05) & abs(mean.diff) > 0.05) %>%
  map(mutate, ID = paste(Chromosome, Start, sep = "_"))

rm(tibble_list)



#### Add info on which comparison we are doing here ####
for (n in seq(from = 1, to = length(names_sites))){
  filtered_tibble_list[[n]]$comparison <- names_sites[n]
}

dir.create("adjpval0005sites")
for (n in seq(from = 1, to = length(names_sites))) {
  write_tsv(x = filtered_tibble_list[[n]], path = paste0("adjpval0005sites/adjpval0005sites_", names_sites[n], ".tsv"))
}
#### Add info on which comparison we are doing here ####



#### Here we create data for clustering via Venn ####
for_venn <- filtered_tibble_list %>%
  map(select, ID, comparison)

dir.create("for_venn")
for (n in seq(from = 1, to = length(names_sites))) {
  write_tsv(x = for_venn[[n]], path = paste0("for_venn/for_venn_", names_sites[n], ".tsv"))
}
#### Here we create data for clustering via Venn ####
