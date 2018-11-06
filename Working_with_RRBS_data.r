install.packages("Matrix", "nlme", "mgcv", "xml2", "tidyverse")
library(tidyverse)
setwd("/media/adrians/USB DISK1/Projekty/FE - Federica/Sequencing/Liver_data")

comparisons_ids <- c("b_vs_rest", "b_vs_c", "i_vs_rest", "i_vs_c", "e_vs_c")



############ GENE REGIONS #################################### 

### Here I upload comparisons
gene_tibble_list <- lapply(list.files(pattern = "*genes.csv"), FUN = read_csv, col_names = T, col_types = "ccnnccnnnnnnnnnnnnnn") 
### Here I checked is the comparisons are uploaded in order, and it seems that they indeed are
names <- list.files(pattern = "*genes.csv") 



### Initialize list, so I can work on it
annotated_tibble_list <- list() 

### Add comparison names to tables
for(n in 1:5){
  annotated_tibble_list[[n]] <- gene_tibble_list[[n]] %>%  ## As tibble
  mutate(
    Comparison = str_replace(string = Chromosome, pattern = ".*", replacement = comparisons_ids[n]) ### Add column of appropriate size and change all its values into name of the comparison
  ) %>%
  select(Comparison, everything()) ### Select all the columns, but put comparison id column at the begininng
}

### Select only p.val lower than 0.05.
gene_pval_list <- annotated_tibble_list  %>%  
  map(filter, comb.p.val <= 0.05) 



### Write the tables
for(n in 1:5){
  write.table(gene_pval_list[n], paste(gene_pval_list[[n]]$Comparison[1], "gene.txt", sep = "_"), sep="\t")
}
############ GENE REGIONS #################################### 
