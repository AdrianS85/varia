############ PREPARATION #####################################

library(tidyverse)
setwd("/media/adrians/USB DISK1/Projekty/FE - Federica/Sequencing/Liver_data")

### Names of comparisons
comparisons_ids <- c("b_vs_rest", "b_vs_c", "i_vs_rest", "i_vs_c", "e_vs_c")

############ PREPARATION #####################################



############ FUNCTIONS ####################################### 

###### Function for getting lists of regions with p_val lower than 0.05 ######
function_short_list <- function(tibble_list, analysis_name, comparison_number = 5, comparison_names = comparisons_ids){
### tibble_list = list of tibbles from RnBeads differential module results
### analysis_name = name of region, string
### comparison_names = string vector of the same number of elements as comparison_number
    
  ### Initialize list, so I can work on it
  annotated_tibble_list <- list()
  ### Establish appropriate amount of for loops
  loop_number <- seq(from = 1, to = comparison_number)
  
  ### Add comparison names to tables
  for(n in loop_number){
    annotated_tibble_list[[n]] <- tibble_list[[n]] %>%  ## As tibble
      mutate(
        Comparison = str_replace(string = Chromosome, pattern = ".*", replacement = comparison_names[n]) ### Add column of appropriate size and change all its values into name of the comparison
      ) %>%
      select(Comparison, everything()) ### Select all the columns, but put comparison id column at the begininng
  }
  
  ### Select only p.val lower than 0.05.
  pval_list <- annotated_tibble_list  %>%  
    map(filter, comb.p.val <= 0.05) 
  
  return(pval_list)  
}



###### Function for writing tables from lists of regions ######
### tibble_list = list of tibbles from RnBeads differential module results, analyzed as neccecary
### analysis_name = name of analysis, tissue, region, string
function_write_list <- function(annotated_tibble_list, analysis_name, comparison_number = 5){
  
  ### Establish appropriate amount of for loops
  loop_number <- seq(from = 1, to = comparison_number)
  
  for(n in loop_number){
    write.table(annotated_tibble_list[n], paste(annotated_tibble_list[[n]]$Comparison[1], "_", analysis_name, ".txt", sep = ""), sep="\t")
  } 
  
  return()  
}

############ FUNCTIONS ####################################### 



############ GENE REGIONS #################################### 

### Analysis name
this_analysis_name <- "liver_genes"

### Here I upload comparisons
gene_tibble_list <- lapply(list.files(pattern = "*genes.csv"), FUN = read_csv, col_names = T, col_types = "ccnnccnnnnnnnnnnnnnn") 

### Here I checked is the comparisons are uploaded in order, and it seems that they indeed are
names <- list.files(pattern = "*genes.csv") 

### Get p.val loerw than 5 tables
gene_pval_list <- function_short_list(tibble_list = gene_tibble_list, analysis_name = this_analysis_name)

### Write the tables
function_write_list(gene_pval_list, this_analysis_name)

############ GENE REGIONS ####################################
