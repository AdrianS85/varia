#########################################################
############ INSTALLS ###################################
#########################################################

# install.packages("BiocManager", "Matrix", "nlme", "mgcv", "xml2", "tidyverse", "matrixStats", "Rcpp", "bit", "rlist")
# install.packages()
# install.packages()
# BiocManager::install("GenomeInfoDb", version = "3.8")
# BiocManager::install("Biostrings", version = "3.8")
# BiocManager::install("Rsamtools", version = "3.8")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", version = "3.8")
# BiocManager::install("bumphunter", version = "3.8")

#########################################################
############ INSTALLS ###################################
#########################################################


##########################################################
############ FUNCTIONS ################################### 
##########################################################

###### FUNCTION FOR GETTING LISTS OF REGIONS WITH PVAL LOWER THAN O.O5 ######
### tibble_list = list of tibbles from RnBeads differential module results
### analysis_name = name of region, string
### comparison_names = string vector of the same number of elements as comparison_number
function_short_list <- function(tibble_list, comparison_number = 5, comparison_names = comparisons_ids){
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
      dplyr::select(Comparison, everything()) ### Select all the columns, but put comparison id column at the begininng
  }
  
  ### Select only p.val lower than 0.05.
  pval_list <- annotated_tibble_list  %>%  
    map(filter, comb.p.val <= 0.05) 
  
  return(pval_list)  
}
###### FUNCTION FOR GETTING LISTS OF REGIONS WITH PVAL LOWER THAN O.O5 ######



###### FUNCTIONS FOR ANNOTATIONG LISTS OF REGIONS WITH GENE NAMES ######
### pval_tibble_list = list of tibbles from RnBeads differential module results, analyzed as neccecary
function_annotate <- function(pval_tibble_list, comparison_number = 5){
  library(bumphunter)
  ### Change this if You want to use another organism
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  
  ### Initialize lists, so I can work on it
  base <- list()
  matched <- list()
  annotated_tibble_list_2 <- list()
  ### Establish appropriate amount of for loops
  loop_number <- seq(from = 1, to = comparison_number)
  ### Establish appropriate annotated genome
  annotations_for_analysis <- annotateTranscripts(txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annotationPackage = "org.Mm.eg.db", by = "gene", requireAnnotation = T)
  
  for (n in loop_number){
    ### Generate 
    base <- pval_tibble_list[[n]]
    matched <- matchGenes(x = pval_tibble_list[[n]], subject = annotations_for_analysis, type = "any")
    
    annotated_tibble_list_2[[n]] <- cbind(base, matched)
  }
  return(annotated_tibble_list_2)
}
###### FUNCTIONS FOR ANNOTATIONG LISTS OF REGIONS WITH GENE NAMES ######



###### FUNCTION FOR WRITING TABLES FROM LISTS OF REGIONS ######
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
###### FUNCTION FOR WRITING TABLES FROM LISTS OF REGIONS ######

##########################################################
############ FUNCTIONS ################################### 
##########################################################



##########################################################
############ PREPARATION #################################
##########################################################

library(tidyverse)
#setwd("/media/adrians/USB DISK1/Projekty/FE - Federica/Sequencing/Liver_data")
setwd("E:/Projekty/FE - Federica/Sequencing/Liver_data")

### Here I upload comparisons
gene_tibble_list <- lapply(list.files(pattern = "*genes.csv"), FUN = read_csv, col_names = T, col_types = "ccnnccnnnnnnnnnnnnnn") 

promoter_tibble_list <- lapply(list.files(pattern = "*promoters.csv"), FUN = read_csv, col_names = T, col_types = "ccnnccnnnnnnnnnnnnnn") 

cgi_tibble_list <- lapply(list.files(pattern = "*cpgislands.csv"), FUN = read_csv, col_names = T, col_types = "ccnnnnnnnnnnnnnnnn") 

tiling_tibble_list <- lapply(list.files(pattern = "*tiling.csv"), FUN = read_csv, col_names = T, col_types = "ccnnnnnnnnnnnnnnnn") 

tiling200_tibble_list <- lapply(list.files(pattern = "*tiling200bp.csv"), FUN = read_csv, col_names = T, col_types = "ccnnnnnnnnnnnnnnnn") 

#sites_tibble_list <- lapply(list.files(pattern = "*_site_*csv"), FUN = read_csv, col_names = T, col_types = "ccnnnnnnnnnnnnnnnn")

### Here I checked is the comparisons are uploaded in order, and it seems that they indeed are
#list.files(pattern = "*genes.csv") 

### Names of comparisons
comparisons_ids <- c("b_vs_rest", "b_vs_c", "i_vs_rest", "i_vs_c", "e_vs_c")
region_ids <- c("gene", "promoter", "cgi", "tiling", "tiling200") #sites

##########################################################
############ PREPARATION #################################
##########################################################



##########################################################
############ ANALYSIS ####################################
##########################################################

### Analysis name
factor1 = "liver" #c("liver", "placenta", "brain")
factor2 = c("genes", "promoter", "cgi", "tiling", "tiling200", "sites")
this_analysis_name <- paste(factor1, factor2, sep = "_")


############ (1) GET PVAL<0.05 TABLES ###########

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

############ (1) GET PVAL<0.05 TABLES ###########



############ (2) GET GENE NAMES FOR REGIONS ###########
### Use (1) as input here

gene_annotated_tibble_list <- function_annotate(gene_pval_list)
promoter_annotated_tibble_list <- function_annotate(promoter_pval_list)
cgi_annotated_tibble_list <- function_annotate(cgi_pval_list)
tiling_annotated_tibble_list <- function_annotate(tiling_pval_list)
tiling200_annotated_tibble_list <- function_annotate(tiling200_pval_list)
#sites_annotated_tibble_list <- function_annotate(sites_pval_list)

### Write the tables
function_write_list(gene_annotated_tibble_list, paste0("named_", this_analysis_name[1]))
function_write_list(promoter_annotated_tibble_list, paste0("named_", this_analysis_name[2]))
function_write_list(cgi_annotated_tibble_list, paste0("named_", this_analysis_name[3]))
function_write_list(tiling_annotated_tibble_list, paste0("named_", this_analysis_name[4]))
function_write_list(tiling200_annotated_tibble_list, paste0("named_", this_analysis_name[5]))
#function_write_list(sites_annotated_tibble_list, paste0("named_", this_analysis_name[6]))

############ (2) GET GENE NAMES FOR REGIONS ###########






  
  
#   for (n in loop_number){
#     gene_temp <- as.tibble(gene_annotated_tibble_list[[n]]$name) %>%
#       mutate(analysis_type = "gene", comparison = comparisons_ids[n])
#     promoter_temp <- as.tibble(promoter_annotated_tibble_list[[n]]$name) %>%
#       mutate(analysis_type = "promoter", comparison = comparisons_ids[n])
#     cgi_temp <- as.tibble(cgi_annotated_tibble_list[[n]]$name) %>%
#       mutate(analysis_type = "cgi", comparison = comparisons_ids[n])
#     tiling_temp <- as.tibble(tiling_annotated_tibble_list[[n]]$name) %>%
#       mutate(analysis_type = "tiling", comparison = comparisons_ids[n])
#     tiling200_temp <- as.tibble(tiling200_annotated_tibble_list[[n]]$name) %>%
#       mutate(analysis_type = "tiling200", comparison = comparisons_ids[n])
#     #sites_temp <- as.tibble(sites_annotated_tibble_list[[n]]$name) %>%
#     #  mutate(analysis_type = "sites")
#   
#     merged_table <- rbind(merged_tibble_list[n][[m]], merged_tibble_list[n][[m]], merged_tibble_list[n][[m]], merged_tibble_list[n][[m]], merged_tibble_list[n][[m]])
#     
#     merged_table_unique <- unique(merged_table)
#     merged_tibble_list[[n]] <- merged_table
#     rm(merged_table)
# 
# #annotation - tutaj są nowe nazwy genów. Dla pojedynczego porównania trzeba zebrać te annotations, zanotować skąd są i złożyć w jedną listę



############ (3) GET GENES FOR ENRICHMENT ###########
### Use (1) as input here

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

############ (3) GET GENES FOR ENRICHMENT ###########

##########################################################
############ ANALYSIS ####################################
##########################################################
