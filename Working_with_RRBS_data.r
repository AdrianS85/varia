############ FUNCTIONS ####################################### 

###### Function for getting lists of regions with p_val lower than 0.05 ######
function_short_list <- function(tibble_list, comparison_number = 5, comparison_names = comparisons_ids){
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



############ PREPARATION #####################################

library(tidyverse)
setwd("/media/adrians/USB DISK1/Projekty/FE - Federica/Sequencing/Liver_data")

### Here I upload comparisons
gene_tibble_list <- lapply(list.files(pattern = "*genes.csv"), FUN = read_csv, col_names = T, col_types = "ccnnccnnnnnnnnnnnnnn") 

promoter_tibble_list <- lapply(list.files(pattern = "*promoters.csv"), FUN = read_csv, col_names = T, col_types = "ccnnccnnnnnnnnnnnnnn") 

cgi_tibble_list <- lapply(list.files(pattern = "*cpgislands.csv"), FUN = read_csv, col_names = T, col_types = "ccnnnnnnnnnnnnnnnn") 

tiling_tibble_list <- lapply(list.files(pattern = "*tiling.csv"), FUN = read_csv, col_names = T, col_types = "ccnnnnnnnnnnnnnnnn") 

tiling200_tibble_list <- lapply(list.files(pattern = "*tiling200bp.csv"), FUN = read_csv, col_names = T, col_types = "ccnnnnnnnnnnnnnnnn") 

#sites_tibble_list <- lapply(list.files(pattern = "*_site_*csv"), FUN = read_csv, col_names = T, col_types = "ccnnnnnnnnnnnnnnnn")

### Here I checked is the comparisons are uploaded in order, and it seems that they indeed are
list.files(pattern = "*genes.csv") 

### Names of comparisons
comparisons_ids <- c("b_vs_rest", "b_vs_c", "i_vs_rest", "i_vs_c", "e_vs_c")

############ PREPARATION #####################################



############ ANALYSIS #################################### 

### Analysis name
factor1 = "liver" #c("liver", "placenta", "brain")
factor2 = c("genes", "promoter", "cgi", "tiling", "tiling200", "sites")
this_analysis_name <- paste(factor1, factor2, sep = "_")
#xxx <- expand.grid(factor1, factor2)
#yyy2 <- paste(xxx[,2], xxx[,1], sep = "_")



### Get p.val less than 5 tables
gene_pval_list <- function_short_list(tibble_list = gene_tibble_list)
promoter_pval_list <- function_short_list(tibble_list = promoter_tibble_list)
cgi_pval_list <- function_short_list(tibble_list = cgi_tibble_list)
tiling_pval_list <- function_short_list(tibble_list = tiling_tibble_list)
tiling200_pval_list <- function_short_list(tibble_list = tiling200_tibble_list)
#sites_pval_list <- function_short_list(tibble_list = sites_tibble_list)

### Write the tables p.value
function_write_list(gene_pval_list, this_analysis_name[1])
function_write_list(promoter_pval_list, this_analysis_name[2])
function_write_list(cgi_pval_list, this_analysis_name[3])
function_write_list(tiling_pval_list, this_analysis_name[4])
function_write_list(tiling200_pval_list, this_analysis_name[5])
#function_write_list(sites_pval_list, this_analysis_name[6])



############ GET GENES FOR ENRICHMENT ###########

### Remove NA and select only relevant columns from gene and promoter analyses
gene_enrich <- gene_pval_list %>%
  map(filter, !is.na(symbol)) %>%
  map(select, Comparison, symbol)
  
promoter_enrich <- promoter_pval_list %>%
  map(filter, !is.na(symbol)) %>%
  map(select, Comparison, symbol)

### Merge gene and promoter data
gene_and_promoters <- list()
for(n in 1:5){
  gene_and_promoters[[n]] <- as.tibble(rbind2(gene_enrich[[n]], promoter_enrich[[n]]))
}

### Keep only unique names after merging
gene_and_promoters_unique <- lapply(X = gene_and_promoters, FUN = unique)

### Write tables
function_write_list(gene_and_promoters_unique, paste0("enrich_", this_analysis_name[1]))

############ GET GENES FOR ENRICHMENT ###########

############ CLUSTERING ###########
### Get all bedgraph files into single matrix
wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
tar jxvf bedtools-2.27.1.tar.gz
cp bedtools2/bin/* /usr/local/bin
gzip -dk *gz
ls *bedGraph > AllLongNames; sed 's/_R1.*//' AllLongNames > AllShortNames; paste AllShortNames AllLongNames > AllCompareNames
parallel --colsep '\t' "sort -k 1,1 -k2,2n {2} > {1}.clust.sort.bedGraph" :::: AllCompareNames
bedtools multiinter -header -i B*.clust.sort.bedGraph > BrainAllCpgs ##Requires that each interval file is sorted by chrom/start
bedtools multiinter -header -i L*.clust.sort.bedGraph > LiverAllCpgs
bedtools multiinter -header -i P*.clust.sort.bedGraph > PlacentaAllCpgs
############ CLUSTERING ###########


############ ANALYSIS ####################################
