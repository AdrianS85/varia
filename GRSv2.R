# library(tidyverse)
# 
# setwd("E:/Projekty/GRS - GJt Review Stress")
# 
# GRS <- read_tsv("GRS.txt", col_types = "ccccnnnnnccccc")
# 
# one <- read.table("GPL6887.txt", header = T)
# two <- read.table("GPL6887_MA.txt", header = T)
# 
# three <- merge(one, two, all.x = T, by = "ProbeID")
# write.table(three, "xxx6.txt", sep = "\t")





library(tidyverse)

COMPARISONS <- readr::read_tsv("comparisons.txt", col_types = "ncccccccccccccccc")

#PRE_ORG_DATA <- readr::read_tsv("test_dataset.txt", col_types = "ccccnnnnncccc", locale = locale(decimal_mark = ","))
#ORG_DATA <- PRE_ORG_DATA %>%
#  dplyr::mutate(ID_NUM = row.names(PRE_ORG_DATA))

DATA <- readr::read_tsv("G1.txt", col_types = "nccccccnnnnncccc", locale = locale(decimal_mark = ","))

LIST_DATA <- split(DATA, f = DATA$Paper)





###### EXTRACT ANNOTATIONS FROM BIOMART ######

library(biomaRt)
#listMarts()

usedMartRAT <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl")
usedMartMUS <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")

#attributesRAT <- listAttributes(usedMartRAT)
#attributesMUS <- listAttributes(usedMartMUS)

###### Here we are trying to choose proper microarrays ###### 
filtersRAT <- listFilters(usedMartRAT)
filtersMUS <- listFilters(usedMartMUS)

PLATFORMS <- read.csv("platforms.txt", header = F)

# Name ensembl databases
PLATFORMS$Biomart <- NA
PLATFORMS$Biomart[1] <- subset(filtersMUS, filtersMUS$description == "ILLUMINA MouseWG 6 V2 probe ID(s) [e.g. ILMN_1240829]")[[1]]
PLATFORMS$Biomart[2] <- subset(filtersRAT, filtersRAT$description == "AFFY RaGene 1 0 st v1 probe ID(s) [e.g. 10930560]")[[1]]
PLATFORMS$Biomart[3] <- subset(filtersMUS, filtersMUS$description == "AGILENT SurePrint G3 GE 8x60k probe ID(s) [e.g. A_65_P05358]")[[1]]
PLATFORMS$Biomart[4] <- subset(filtersMUS, filtersMUS$description == "AFFY MoGene 2 1 st v1 probe ID(s) [e.g. 17532593]")[[1]]
PLATFORMS$Biomart[5] <- subset(filtersMUS, filtersMUS$description == "AGILENT SurePrint G3 GE 8x60k probe ID(s) [e.g. A_65_P05358]")[[1]]
PLATFORMS$Biomart[6] <- subset(filtersRAT, filtersRAT$description == "AFFY Rat230 2 probe ID(s) [e.g. 1375651_at]")[[1]]
PLATFORMS$Biomart[7] <- subset(filtersMUS, filtersMUS$description == "AGILENT WholeGenome 4x44k v1 probe ID(s) [e.g. A_51_P323880]")[[1]]
PLATFORMS$Biomart[8] <- subset(filtersMUS, filtersMUS$description == "AFFY MG U74Av2 probe ID(s) [e.g. 96290_f_at]")[[1]] 
PLATFORMS$Biomart[9] <- subset(filtersMUS, filtersMUS$description == "AFFY MoGene 1 0 st v1 probe ID(s) [e.g. 10598025]")[[1]]
PLATFORMS$Biomart[10] <- NA
PLATFORMS$Biomart[11] <- subset(filtersMUS, filtersMUS$description == "AFFY MoGene 1 0 st v1 probe ID(s) [e.g. 10598025]")[[1]]
PLATFORMS$Biomart[12] <- NA

# Add species information
PLATFORMS$Species <- c("m", "r", "m", "m", "m", "r", "m", "m", "m", "m", "m", "m")

# Make platform information
PLATFORMS <- PLATFORMS %>%
  mutate(PL_ID = str_remove(V1, "[/ ].*"), )
  

###### Here we are trying to choose proper microarrays ###### 



###### Here we are annotating experiments with biomart platform descriptions ###### 

for_annot_COMPARISONS <- COMPARISONS %>%
  dplyr::select(Paper, PL_ID) %>%
  unique()

FOR_ANNOT_PLATFORMS <- merge(for_annot_COMPARISONS, PLATFORMS, by = "PL_ID") %>%
  arrange(Paper)

###### Here we are annotating experiments with biomart platform descriptions ###### 



###### Here we need to remove experiments with microarrays not captured in ensembl ######
WHICH_EXP_TO_ANAL <- seq(1,nrow(FOR_ANNOT_PLATFORMS))[c(-13, -16)]

## Here we are doing probe ensembl annotation
ANNOT_LIST_DATA <- list()

for (n in WHICH_EXP_TO_ANAL){
  ANNOT_LIST_DATA[[n]] <- getBM(attributes = c(FOR_ANNOT_PLATFORMS$Biomart[n], "external_gene_name"), 
                      filters = FOR_ANNOT_PLATFORMS$Biomart[n], 
                      values = LIST_DATA[[n]]$Probe_ID,
                      uniqueRows = F,
                      mart = 
                        if(FOR_ANNOT_PLATFORMS$Species[n] == "m"){ usedMartMUS } else if(FOR_ANNOT_PLATFORMS$Species[n] == "r"){ usedMartRAT }
  )
  }

## HERE WE WILL COLLAPSE duplicated probe-genename pairs WITHIN ANNOTATED GENE LIST (and also change column names to standarized ones)
UNIQ_ANNOT_LIST_DATA <- list()

for (n in WHICH_EXP_TO_ANAL){
  UNIQ_ANNOT_LIST_DATA[[n]] <- unique(ANNOT_LIST_DATA[[n]])
  colnames(UNIQ_ANNOT_LIST_DATA[[n]]) <- c("Probe_ID", "ensembl_gene_name")
}

## Here we will collapse all gene names for given probe
for (n in WHICH_EXP_TO_ANAL){
  UNIQ_ANNOT_LIST_DATA[[n]] <- aggregate(ensembl_gene_name~Probe_ID, data = UNIQ_ANNOT_LIST_DATA[[n]], FUN = str_c) ## Here we aggregate the genenames into single row
  UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- lapply(X = UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name, FUN = paste, collapse = "; ") ##
  UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name <- as.character(UNIQ_ANNOT_LIST_DATA[[n]]$ensembl_gene_name) ## Here we make sure that collapsed genename column is not a list (need for writing function)
}

# Je?eli w ensemble nie ma przypisanego genu (oznaczenie NA w przys?anym pliku "100_genow") to program pokazuje we wszystkich rubrykach oznaczenie NA chocia? te kom?rki zawiera?y wcze?niej dane (plik AAAreview11). W pe?nych danych kt?re przys?a?e? mi przed ?wi?tami (plik 1 znajduj?cy si? w skopresowanych danych "TEST_ANNOTATION") tych sond w og?le nie ma. Wygl?da na to, ?e to co przys?a?e? wczoraj odnosi si? chyba do wcze?niejszego pliku w kt?rym dane dla sond nie maj?cych przypisanego genu w ensemble zosta?y ca?kowicie usuni?te a to nie jest dobre :( - AMS W LIST_DATA S? TE DANE, ALE W TEST_ANNOTATION JUZ NIE

## Here we marge probes annotated with ensembl with probes annotated with other methods
TEST_ANNOTATION <- list()
for (n in WHICH_EXP_TO_ANAL){
  TEST_ANNOTATION[[n]] <- merge(LIST_DATA[[n]], UNIQ_ANNOT_LIST_DATA[[n]], by = "Probe_ID", all.x = T)
}  

# Here we write output test annotation files
dir.create("TEST_ANNOTATION")
for (n in WHICH_EXP_TO_ANAL){
  write.table(TEST_ANNOTATION[[n]], file = paste0("TEST_ANNOTATION/", n, ".txt"), sep = "\t", dec = ",")
} 


### Here write probe extraction function, but for now we will have this ###
xxx <- lapply(list.files(pattern = "*.txt"), FUN = read.delim, header = T, sep = "\t", dec = ",")
### Here we merge list of tables into single table
EXTRACTION_TABLE <- rlist::list.rbind(TEST_ANNOTATION) #####!!!!!!!!!!!!

TO_BE_EXTRACTED <- read.delim(file = "test adnotacji wg ensemble.txt", header = T) ### Here we load GJts probes to be extracted

POST_EXTRACTION <- merge(x = TO_BE_EXTRACTED, y = EXTRACTION_TABLE, by = "Probe_ID", all.x = T) ### Here we extract the probes...

write_tsv(x = POST_EXTRACTION, path = "100_genow_V2.txt", na = "NA", append = FALSE)






### ADDITIONAL
# MISSING PROBES
LIST_DATA[[11]]$Probe_ID[!(LIST_DATA[[11]]$Probe_ID %in% TEST_ANNOTATION[[11]]$Probe_ID)]



####### To trzeba da? gdzie? do g??wnego flowa
### Here we get all the written output tables

### Here we get all the written output tables



### Here we get all the written output tables
zzz <- read_tsv(file = "100_genow.txt") 
xxx <- aggregate(ensembl_gene_name~Probe_ID, data = zzz, FUN = str_c) ## Here we aggregate
xxx$upgraded_ens_names <- map(.x = xxx$ensembl_gene_name, .f = unique, collapse = ' ') # Here we remove duplicate gene names
xxx$upgraded_ens_names <- lapply(X = xxx$upgraded_ens_names, FUN = paste, collapse = "; ") ##!!!!!!!! Here we collapse character vector into single string



paste((d), )


yyy <- as.data.frame(xxx$ensembl_gene_name)


xxx$ensembl_gene_name <- str_remove(string = xxx$ensembl_gene_name, pattern = '[c][(][\\]')

write.table(x = xxx, file = "test.txt", sep = "\t", )
### HERE WE EXTRACT SPECIFIC PROBES FROM DATASET FOR GRZEGORZ ###

head(zzz)
###### Here we need to remove experiments with microarrays not captured in ensembl ######






                                
                        
                                
                                
                                

#### TESTING MULTIPLE-PROBE ANNOATION ####
# Add: make the analysis on 100 randomly selected IDs, make a table annotating which filters were actually used


exp_species <- c("m", "r", "m")
#exp_species <- c("m", "r", "m", "m", "m", "r", "m", "m", "m", "m", "m", "m", "m", "m", "m", "m")
#SHORT_LIST_DATA

SHORT_LIST_DATA <- SHORT_LIST_DATA[1:3]


final_list <- list()
ANNOT_SHORT_LIST_DATA <- list(list(), list(), list())
##### This needs to be looped for all data files ##### 
for(n in seq_along(SHORT_LIST_DATA)){
  
  # Here we establish which mart(species) we are using in this given dataset
  usedMart <- switch(exp_species[n], 
                     "m" = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl"), 
                     "r" = useMart("ENSEMBL_MART_ENSEMBL", dataset = "rnorvegicus_gene_ensembl"))
  filters <- listFilters(usedMart)
  
  # Here we extract all the potential gene identifiers
  potental_identifiers <- c(filters[grep(pattern = "^ensembl(.*)", filters[[1]]) , 1], 
                             filters[grep(pattern = "^refseq(.*)", filters[[1]]) , 1], 
                             filters[grep(pattern = "^affy(.*)", filters[[1]]) , 1], 
                             filters[grep(pattern = "^agilent(.*)", filters[[1]]) , 1], 
                             filters[grep(pattern = "^illumina(.*)", filters[[1]]) , 1])
  
  # Here we are annotating given datasets with data from all the relevant databases
  for (m in seq_along(potental_identifiers)){
    ANNOT_SHORT_LIST_DATA[[n]][[m]] <- getBM(attributes = c(potental_identifiers[[m]], "external_gene_name"), 
                                  filters = potental_identifiers[[m]], 
                                  values = SHORT_LIST_DATA[[n]]$Probe_ID,
                                  uniqueRows = F,
                                  mart = usedMart
    )
  }
}

ANNOT_SHORT_LIST_DATA <- list()

### Here we are getting all relevant possible databases
potentalNames <- c(filtersMUS[grep(pattern = "^ensembl(.*)", filtersMUS[[1]]) , 1], filtersMUS[grep(pattern = "^refseq(.*)", filtersMUS[[1]]) , 1], filtersMUS[grep(pattern = "^affy(.*)", filtersMUS[[1]]) , 1], filtersMUS[grep(pattern = "^agilent(.*)", filtersMUS[[1]]) , 1], filtersMUS[grep(pattern = "^illumina(.*)", filtersMUS[[1]]) , 1])



### Here we are extracting annotations from all the relevant databases
test_list <- list()
for(n in seq_along(potentalNames)) {
  test_list[[n]] <- getBM(attributes = c(potentalNames[[n]], "external_gene_name"), 
                          filters = potentalNames[[n]], 
                          values = test[[1]],
                          uniqueRows = T,
                          mart = usedMartMUS) 
}

# Here I get vector with length of each relevant database
for(n in seq_along(test_list)) {
  tables[n] <- length(test_list[[n]][[1]])
}
# Here I find out at which list the highest number of hits and put it into final data extraction
final_list[[1]] <- test_list[[which(tables == max(tables))]] ##
##### This needs to be looped for all data files ##### 


#### TESTING MULTIPLE-PROBE ANNOATION ####
                                
                                
                                





                                
                                



#Tutaj trzeba zrobić tak: zrobić .txt files z każdego z eksperymentów. Przefiltrować odpowiednią bazę danych przez nazwy sond z tego eksperymentu.
#Potem jeszcze trzeba je sprowadzić do jednolitego nazewnictwa (mus or homo?)

featureDat <- getBM(attributes = c("agilent_wholegenome_4x44k_v1", "affy_mg_u74av2", "affy_mogene_1_0_st_v1"), 
                    mart = usedMartMUS)

featureDat <- getBM(attributes = c("affy_ragene_2_1_st_v1", "external_gene_name"), 
                    filters = "affy_ragene_2_1_st_v1", 
                    values = rownames(dataMatrixMW),
                    uniqueRows = TRUE,
                    mart = usedMart) #This is from Biomart. Here we download the names
xxx <- aggregate(external_gene_name~affy_ragene_2_1_st_v1, data = featureDat, FUN = c) #Here we collapse names, cause row names have to be unique and many probes correspond to more than 1 gene name in biomart

yyy <- as.data.frame(rownames(dataMatrixMW))#Here we produce full list of features, cause some of them dont have annotations in biomart and they are not returned when annotations

colnames(yyy) <- "affy_ragene_2_1_st_v1" #This is for the merge to work
xxxx <- merge(yyy, xxx, by = "affy_ragene_2_1_st_v1", all.x = TRUE) #Here we produce annotated full list of features
rownames(xxxx) <- rownames(dataMatrixMW) #Here we name the features
xxxx$affy_ragene_2_1_st_v1 <- as.character(xxxx$affy_ragene_2_1_st_v1) #Not sure if thats nececcary
xxxx$external_gene_name <- as.character(xxxx$external_gene_name) #Not sure if thats nececcary

featureMW <- new("AnnotatedDataFrame", data = xxxx) #Feature data is in Annotated Dataframe format
dim(featureDat)
rm(usedMart, featureDat, xxx, xxxx, yyy)

###### EXTRACT ANNOTATIONS FROM BIOMART ######



### HERE WE CHECK IF ALL VALUES IN INPUT TABLE ARE EITHER NONES, NAS OR CORRECT VALUES
FUNCTION__1__check_rm_of_unid_val <- function(){
# Table containing all NONEs
ORG_DATA__1__NONEs <- ORG_DATA %>%
  dplyr::filter(Gene_symbol == "NONE") 

# Table containing all NAs
ORG_DATA__2__NAs <- ORG_DATA %>%
  dplyr::filter(is.na(Gene_symbol)) 

# Number of removed NONEs and NAs
ORG_DATA__3__sum_of_NONEs_and_NAs <- nrow(ORG_DATA__1__NONEs) + nrow(ORG_DATA__2__NAs)

# Here is number of unidentified rows sliped during filtering out bad values
ORG_DATA__x__missing_processed_rows <- nrow(ORG_DATA) - (ORG_DATA__3__sum_of_NONEs_and_NAs + nrow(NO_UNIDS_ORG_DATA))

if (ORG_DATA__x__missing_processed_rows != 0){ 
  stop("Hey, buddy! You have some wierd values in Your raw data, guy! Better check whats happening, or Your results will smell of farts!")
}

return(ORG_DATA__x__missing_processed_rows)
}
### HERE WE CHECK IF ALL VALUES IN INPUT TABLE ARE EITHER NONES, NAS OR CORRECT VALUES



## Subset of table without "NONE" gene ids.
NO_UNIDS_ORG_DATA <- ORG_DATA %>%
  dplyr::filter(Gene_symbol != "NONE" & !is.na(Gene_symbol))

FUNCTION__1__check_rm_of_unid_val() # Check if it went well


###### WHOLE DATASET ANALYSIS ######

# Tutaj liczymy ile razy geny wyst?puj? w oryginalnym dataset, nie patrz?c czy s? up czy down
WHOLE_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(Gene_symbol) %>%
  group_by(Gene_symbol) %>%
  summarise(number = n())



# Tutaj liczymy ile razy geny wyst?puj? w oryginalnym dataset, patrz?c czy s? up czy down
UorDWHOLE_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(Gene_symbol, logFC) %>%
  mutate(Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) %>%
  mutate(Symbol_direction = paste(Gene_symbol, Symbol_direction, sep = "_")) %>%
  group_by(Symbol_direction) %>%
  summarise(number = n()) %>% 
  mutate(Gene_symbol2 = str_remove(Symbol_direction, "_.*"))

###### WHOLE DATASET ANALYSIS ######



###### COMPARISONS-CENTERED ANALYSIS ######



### Here we set whether we want to analyze papers or comparisons
P_or_C = quo(Paper) #" GroupID OR Paper "



# Tutaj liczymy ile razy geny wyst?puj? W KA?DYM Z POR?WNA?, nie patrz?c czy s? up czy down
COMP_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(!!P_or_C, Gene_symbol) %>%
  group_by(!!P_or_C, Gene_symbol) %>%
  summarise(number = n())



# Tutaj liczymy ile razy geny wyst?puj? W KA?DYM Z POR?WNA?, patrz?c czy s? up czy down    
UorDCOMP_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(!!P_or_C, Gene_symbol, logFC) %>%
  mutate(Symbol_direction = ifelse(logFC > 0, "UP", "DOWN")) %>%
  mutate(Symbol_direction = paste(Gene_symbol, Symbol_direction, sep = "_")) %>%
  group_by(!!P_or_C, Symbol_direction) %>%
  summarise(Sym_dir_number = n()) %>%
  mutate(Gene_symbol2 = str_remove(Symbol_direction, "_.*"))



#Divide data into genes expressed in single direction in given comparison, vs genes expressed in different direction (bad genes)
nonUNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UorDCOMP_NO_UNIDS_ORG_DATA %>%
  group_by(!!P_or_C) %>%
  filter(duplicated(Gene_symbol2, fromLast = T) | duplicated(Gene_symbol2))

UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UorDCOMP_NO_UNIDS_ORG_DATA %>%
  group_by(!!P_or_C) %>%
  filter(!duplicated(Gene_symbol2, fromLast = T) & !duplicated(Gene_symbol2))



# Check if unique/duplicated division went well           
if (nrow(UorDCOMP_NO_UNIDS_ORG_DATA) - (nrow(nonUNIQ_UorDCOMP_NO_UNIDS_ORG_DATA) + nrow(UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA)) != 0){ 
  stop("Hey, fwend! You have some wierd values in Your counted data, buddy! Better check whats happening, or Your results will smell of moose scrotum!")
}



# Here we make a table only with genes that were replicated in few comparisons
REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA %>%
  filter(Sym_dir_number >= 3)


#Annotate base on Paper OR GroupID
ANNO_REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA <- merge(REPL_UNIQ_UorDCOMP_NO_UNIDS_ORG_DATA, COMPARISONS, by = "Paper")


###### COMPARISONS-CENTERED ANALYSIS ######



             
      
    
    
    

  
  
  
  
  
  
KAJA:

setwd("E:/Projekty/Kaja Review LDH")
  
library(enrichR)

dbs_Onto_Path <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018", "MGI_Mammalian_Phenotype_2017", "Human_Phenotype_Ontology", "KEGG_2016", "WikiPathways_2016", "Panther_2016", "Reactome_2016", "BioCarta_2016", "NCI-Nature_2016", "ARCHS4_Kinases_Coexp", "HumanCyc_2016", "BioPlex_2017", "SILAC_Phosphoproteomics")
dbs_Regul <- c("Genome_Browser_PWMs", "TRANSFAC_and_JASPAR_PWMs", "Transcription_Factor_PPIs", "ChEA_2016",  "TF-LOF_Expression_from_GEO", "PPI_Hub_Proteins", "ENCODE_TF_ChIP-seq_2015", "ENCODE_Histone_Modifications_2015", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X", "CORUM", "Pfam_InterPro_Domains", "Phosphatase_Substrates_from_DEPOD", "TF_Perturbations_Followed_by_Expression", "ARCHS4_TFs_Coexp", "miRTarBase_2017", "TargetScan_microRNA_2017", "Enrichr_Submissions_TF-Gene_Coocurrence", "Epigenomics_Roadmap_HM_ChIP-seq")
dbs_Drug_Tissue_Other <- c("Jensen_TISSUES", "ARCHS4_IDG_Coexp", "DrugMatrix", "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO", "OMIM_Disease", "Jensen_DISEASES", "DSigDB",  "Jensen_COMPARTMENTS", "ARCHS4_Tissues", "Tissue_Protein_Expression_from_Human_Proteome_Map", "Tissue_Protein_Expression_from_ProteomicsDB", "Mouse_Gene_Atlas", "ESCAPE", "Chromosome_Location", "MSigDB_Computational", "dbGaP", "Genes_Associated_with_NIH_Grants", "GeneSigDB")
dbs <- c(dbs_Onto_Path, dbs_Regul, dbs_Drug_Tissue_Other)

genes <- c("ldha", "ldhb")

ALL_ENRICHR <- enrichR::enrichr(genes, dbs)
ALL_ENRICHR_DATA <- rlist::list.rbind(ALL_ENRICHR)
write.table(ALL_ENRICHR_DATA, "ALL_ENRICHR_DATA.txt", sep="\t")

getwd()

tp <- read.table("TF_PPI_FROM_ENRICHR.txt", header = T, stringsAsFactors = F)

TP_ENRICHR <- enrichR::enrichr(tp$Enrichr_TFs_PPIs_unique, dbs_Onto_Path)
ALL_TP_ENRICHR_DATA <- rlist::list.rbind(TP_ENRICHR)
write.table(ALL_TP_ENRICHR_DATA, "ALL_TP_ENRICHR_DATA.txt", sep="\t")
