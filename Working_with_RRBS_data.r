#########################################################
############ INSTALLS ###################################
#########################################################

# install.packages("BiocManager", "Matrix", "nlme", "mgcv", "xml2", "tidyverse", "matrixStats", "Rcpp", "bit", "rlist", "enrichR")
# BiocManager::install("GenomeInfoDb", version = "3.8")
# BiocManager::install("Biostrings", version = "3.8")
# BiocManager::install("Rsamtools", version = "3.8")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene", version = "3.8")
# BiocManager::install("bumphunter", version = "3.8")
# BiocManager::install("org.Mm.eg.db", version = "3.8")

#########################################################
############ INSTALLS ###################################
#########################################################


##########################################################
############ FUNCTIONS ################################### 
##########################################################

###### FUNCTION FOR GETTING LISTS OF REGIONS WITH PVAL LOWER THAN O.O5 ######
### tibble_list = list of tibbles from RnBeads differential module results
function_short_list <- function(tibble_list){
  ### Initialize list, so I can work on it
  annotated_tibble_list <- list()

  
  ### Add comparison names to tables
  for(n in seq(from = 1, to = length(comparisons_ids))){
    annotated_tibble_list[[n]] <- tibble_list[[n]] %>%  ## As tibble
      dplyr::mutate(
        Comparison = stringr::str_replace(string = Chromosome, pattern = ".*", replacement = comparisons_ids[n]) ### Add column of appropriate size and change all its values into name of the comparison ##########CHANGED!!!!!!!!!!!!!!!
      ) %>%
      dplyr::select(Comparison, dplyr::everything()) ### Select all the columns, but put comparison id column at the begininng
  }
  
  ### Select only p.val lower than 0.05.
  pval_list <- annotated_tibble_list  %>%  
    purrr::map(filter, comb.p.val <= 0.05) 
  
  return(pval_list)  
}
###### FUNCTION FOR GETTING LISTS OF REGIONS WITH PVAL LOWER THAN O.O5 ######



###### FUNCTIONS FOR ANNOTATIONG LISTS OF REGIONS WITH GENE NAMES ######
### pval_tibble_list = list of tibbles from RnBeads differential module results, analyzed as neccecary
function_annotate <- function(pval_tibble_list){
  library(bumphunter)
  ### Change this if You want to use another organism
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  
  ### Initialize lists, so I can work on it
  base <- list()
  matched <- list()
  annotated_tibble_list_2 <- list()
  
  ### Establish appropriate annotated genome
  annotations_for_analysis <- annotateTranscripts(txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annotationPackage = "org.Mm.eg.db", by = "gene", requireAnnotation = T)
  
  for (n in seq(from = 1, to = length(comparisons_ids))){
    ### Generate 
    base <- pval_tibble_list[[n]]
    matched <- matchGenes(x = pval_tibble_list[[n]], subject = annotations_for_analysis, type = "any")
    
    annotated_tibble_list_2[[n]] <- cbind(base, matched)
  }
  return(annotated_tibble_list_2)
}
###### FUNCTIONS FOR ANNOTATIONG LISTS OF REGIONS WITH GENE NAMES ######



###### FUNCTION FOR GETTING ENRICHMENT OVERVIEW FROM ENTIRE DATA ######
function_enrich_overview <- function(annotated_tibble_list = list(gene_annotated_tibble_list, promoter_annotated_tibble_list, cgi_annotated_tibble_list, tiling_annotated_tibble_list, tiling200_annotated_tibble_list)){
  
  ### Initialize list, so I can work on it
  merged_tibble_list <- list()
  merged_number <- 1

  ### Double loop
  for (c in seq(from = 1, to = length(comparisons_ids))){
    for (r in seq(from = 1, to = length(region_ids))){
      merged_tibble_list[[merged_number]] <- as.tibble(unique(annotated_tibble_list[[r]][[c]]$name)) %>%
        mutate(Analysis_type = region_ids[r], Comparison = comparisons_ids[c]) ##########CHANGED!!!!!!!!!!!!!!! 
      merged_number <- merged_number + 1
    }
  }

  ### Merge list of tibbles into tibble
  merged_tibble <- rlist::list.rbind(merged_tibble_list)
  
  ### Create list of tibbles based on comparison names
  final_list <- list()
  for (c in seq(from = 1, to = length(comparisons_ids))){
    final_list[[c]] <- merged_tibble %>%
      filter(Comparison == comparisons_ids[c])
  }
  
  return (final_list)
}
###### FUNCTION FOR GETTING ENRICHMENT OVERVIEW FROM ENTIRE DATA ######



###### FUNCTION FOR GETTING ENRICHMENTS ######
### tibble_list - tibble list with tibbles for all 5 comparisons, 
### Currently I work with tibbles with 3 columns, with the first one being gene names and the third one being comparison_id. We should be able to use different ones
function_enrich_list <- function(tibble_list = enrich_overview_list, gene_names_column = 1){
  ### Create list of interesting databases
  dbs_Onto_Path <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "MGI_Mammalian_Phenotype_2017", "Human_Phenotype_Ontology", "KEGG_2016", "WikiPathways_2016", "Panther_2016", "Reactome_2016", "BioCarta_2016", "NCI-Nature_2016", "ARCHS4_Kinases_Coexp", "HumanCyc_2016", "BioPlex_2017", "SILAC_Phosphoproteomics")
  dbs_Regul <- c("Genome_Browser_PWMs", "TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs", "ChEA_2016",  "TF-LOF_Expression_from_GEO", "PPI_Hub_Proteins", "ENCODE_TF_ChIP-seq_2015", "ENCODE_Histone_Modifications_2015", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "CORUM", "Pfam_InterPro_Domains", "Phosphatase_Substrates_from_DEPOD", "TF_Perturbations_Followed_by_Expression", "ARCHS4_TFs_Coexp", "miRTarBase_2017", "TargetScan_microRNA_2017", "Enrichr_Submissions_TF-Gene_Coocurrence", "Epigenomics_Roadmap_HM_ChIP-seq")
  dbs_Drug_Tissue_Other <- c("ARCHS4_IDG_Coexp", "DrugMatrix", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "OMIM_Disease", "Jensen_DISEASES", "DSigDB",  "Jensen_COMPARTMENTS", "ARCHS4_Tissues", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB", "Mouse_Gene_Atlas", "ESCAPE", "Chromosome_Location", "MSigDB_Computational")
  dbs <- c(dbs_Onto_Path, dbs_Regul, dbs_Drug_Tissue_Other)
  
  ##  dbs_Removed <- c("Jensen_TISSUES", "dbGaP", "Genes_Associated_with_NIH_Grants", "GeneSigDB")
  ##  Error in scan(file = file, what = what, sep = sep, quote = quote, dec = dec,  : 
  ##  line 1 did not have 2 elements OR line 109 did not have 9 elements
  ##  https://www.biostars.org/p/349240/
  
  ### Initialize list, so I can work on it
  enriched <- list()

  for(n in seq(from = 1, to = length(comparisons_ids))){
    enriched[[n]] <- enrichR::enrichr(tibble_list[[n]][[gene_names_column]], dbs) %>%
      map(filter, P.value <= 0.05) %>%
      map(as.tibble)
  }
  
  ### Initialize list, so I can work on it
  enriched_names <- list()
  ### Initialize list counter for the loop
  merged_number_2 <- 1
  
  ### Add information on databaseses and comparisons used
  for(n in seq(from = 1, to = length(comparisons_ids))){
    for(d in seq(from = 1, to = length(dbs))){
      enriched_names[[merged_number_2]] <- enriched[[n]][[d]]%>%
        mutate(Database = dbs[d], Comparison = comparisons_ids[n])
      merged_number_2 <- merged_number_2 + 1  
    }
  }
  
  merged_enriched_names <- rlist::list.rbind(enriched_names)

  ### Create list of tibbles based on comparison names
  final_list_2 <- list()
  for (c in seq(from = 1, to = length(comparisons_ids))){
    final_list_2[[c]] <- merged_enriched_names %>%
      filter(Comparison == comparisons_ids[c])
  }
  
  return(final_list_2)
}
###### FUNCTION FOR GETTING ENRICHMENTS ######



###### FUNCTION FOR WRITING TABLES FROM LISTS OF REGIONS ######
### tibble_list = list of tibbles from RnBeads differential module results, analyzed as neccecary
### analysis_name = name of analysis, tissue, region, string
function_write_list <- function(annotated_tibble_list, analysis_name){

  for(n in seq(from = 1, to = length(comparisons_ids))){
    write.table(annotated_tibble_list[n], paste0(annotated_tibble_list[[n]]$Comparison[1], "_", analysis_name, ".txt"), sep="\t")
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

#setwd("/media/adrians/USB DISK1/Projekty/FE - Federica/Sequencing/Liver_data")
#setwd("E:/Projekty/FE - Federica/Sequencing/Liver_data")

library(tidyverse)

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
comparisons_ids <- c("b_vs_rest", "b_vs_c", "i_vs_rest", "i_vs_c", "e_vs_c") #
region_ids <- c("gene", "promoter", "cgi", "tiling", "tiling200") #sites

### Analysis name
factor1 = "liver" #c("liver", "placenta", "brain")
factor2 = c("genes", "promoter", "cgi", "tiling", "tiling200", "sites")
this_analysis_name <- paste(factor1, factor2, sep = "_")

##########################################################
############ PREPARATION #################################
##########################################################



##########################################################
############ ANALYSIS ####################################
##########################################################

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
#save(gene_annotated_tibble_list, promoter_annotated_tibble_list, cgi_annotated_tibble_list, tiling_annotated_tibble_list, tiling200_annotated_tibble_list, file = "annotated_tibble_list.Rdata")
############ (2) GET GENE NAMES FOR REGIONS ###########



############ (3) GET ENRICHMENT OVERVIEW ###########
### Use (2) as input here
### Get genes for enrichment overview from all regions
## For some reason liver, bvc needs special attension (doesnt work), hence:
## xxx <- enrich_overview_list
## xxx[2] <- NULL
## overview_enrichment <- function_enrich_list(tibble_list = xxx)
enrich_overview_list <- function_enrich_overview()

### Write genes for enrichment overview from all regions
function_write_list(enrich_overview_list, paste0("overview_", factor1, "_enrich_gene_list"))

### Get overview enrichment
overview_enrichment <- function_enrich_list()

### Write overview enrichment
function_write_list(overview_enrichment, paste0("overview_", factor1, "_enrichment"))

############ (3) GET GENES FOR ENRICHMENT OVERVIEW FROM ALL REGIONS ###########



############ (4) GET GENES FOR ENRICHMENT ###########
### Use (1) as input here

### Remove NA and select only relevant columns from gene and promoter analyses
gene_enrich <- gene_pval_list %>%
  map(filter, !is.na(symbol)) %>%
  map(dplyr::select, Comparison, symbol)
  
promoter_enrich <- promoter_pval_list %>%
  map(filter, !is.na(symbol)) %>%
  map(dplyr::select, Comparison, symbol)

### Merge gene and promoter data
gene_and_promoters <- list()
for(n in 1:5){
  gene_and_promoters[[n]] <- as.tibble(rbind2(gene_enrich[[n]], promoter_enrich[[n]]))
}

### Keep only unique names after merging
gene_and_promoters_unique <- lapply(X = gene_and_promoters, FUN = unique)

### Write tables
function_write_list(gene_and_promoters_unique, paste0("enrich_", factor1, "_"))

### Get enrichment?
enrichment <- function_enrich_list(tibble_list = gene_and_promoters_unique, gene_names_column = 2)

### Write overview enrichment
function_write_list(enrichment, paste0(factor1, "_enrichment"))

############ (4) GET GENES FOR ENRICHMENT ###########


############ CLUSTERING ###########
### Get all bedgraph files into single matrix
#In folder with all .bedgraph.gz files
### Get bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz
tar jxvf bedtools-2.27.1.tar.gz
cp bedtools2/bin/* /usr/local/bin

### Prepare sorted bedgraph files
gzip -dk *gz
ls *bedGraph > AllLongNames; sed 's/_R1.*//' AllLongNames > AllShortNames; paste AllShortNames AllLongNames > AllCompareNames
parallel --colsep '\t' "sort -k 1,1 -k2,2n {2} > {1}.clust.sort.bedGraph" :::: AllCompareNames ##sort a BED file by chromosome then by start position https://bedtools.readthedocs.io/en/latest/content/tools/sort.html


### PREPARE list of all captured CpG sites in all bedgraph files in given tissue
ls B*.clust.sort.bedGraph | sudo tee BrainMultiinterOrder
sudo bedtools multiinter -header -i B*.clust.sort.bedGraph | sudo tee BrainAllCpgs ##Requires that each interval file is sorted by chrom/start
ls L*.clust.sort.bedGraph | sudo tee LiverMultiinterOrder
bedtools multiinter -header -i L*.clust.sort.bedGraph > LiverAllCpgs
ls L*.clust.sort.bedGraph | sudo tee PlacentaMultiinterOrder
bedtools multiinter -header -i P*.clust.sort.bedGraph > PlacentaAllCpgs

### INTO R:
x <- readr::read_tsv(file = "BrainAllCpgs", col_types = "cnnnc") ### Read only relevant columns
x <- x %>% dplyr::()
x2 <- x[, 4:6]

### PREPARE .bedgraph files
https://stackoverflow.com/questions/14096814/merging-a-lot-of-data-frames

bedgraph_tibble_list <- lapply(list.files(pattern = "*sort.bedGraph"), FUN = readr::read_tsv, col_names = F, col_types = "cnnn") 
a_bedgraph_tibble_list <- bedgraph_tibble_list %>% map(mutate, ID = paste(X1, X2, sep = "_"))
b_bedgraph_tibble_list <- a_bedgraph_tibble_list %>% map(select, X4, ID)
## Name bedgraph columns with filenames
bedgraph_names <- str_remove(list.files(pattern = "*sort.bedGraph"), "_S(.*)") 
for (n in seq(from = 1, to = length(bedgraph_names))) { colnames(b_bedgraph_tibble_list[[n]]) <- c(bedgraph_names[n], "ID") }

c_bedgraph_tibble_list <- Reduce(function(x, y) merge(x, y, by = "ID", all=TRUE), b_bedgraph_tibble_list)

write.table(c_bedgraph_tibble_list, file = "all_bedgraphs.tsv", sep = "\t", dec = ".")

### INTO R

### PREPARE DIFFERENTIALLY METHYLATED CPG FILES
### Move diffmeth cpg files into the working folder. Rename them so that they are recognizable
rename 's/^diffMethTable/b1_diffMethTable/' diffMethTable*


############ CLUSTERING ###########

##########################################################
############ ANALYSIS ####################################
##########################################################



##########################################################
############ bvc, liver ##################################
##########################################################

# removed unknown genes, pseudogenes, Gm genes

one <- read_csv("xxx.txt")

two <- enrichR::enrichr(one[[1]], dbs)

three <- rlist::list.rbind(two)

write.table(three, "bvc_liver_overvew_enrichment.txt")

##########################################################
############ bvc, liver ##################################
##########################################################
