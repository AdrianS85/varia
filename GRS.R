library(tidyverse)

setwd("E:/Projekty/GRS - GJt Review Stress")

GRS <- read_tsv("GRS.txt", col_types = "ccccnnnnnccccc")

one <- read.table("GPL6887.txt", header = T)
two <- read.table("GPL6887_MA.txt", header = T)

three <- merge(one, two, all.x = T, by = "ProbeID")
write.table(three, "xxx6.txt", sep = "\t")





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
PLATFORMS$Biomart[8] <- subset(filtersMUS, filtersMUS$description == "AFFY MG U74Av2 probe ID(s) [e.g. 102126_at]")[[1]]
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

# Error in getBM(attributes = c(FOR_ANNOT_PLATFORMS$Biomart[n], "external_gene_name"),  :     
#                  Invalid attribute(s): NA 
#                Please use the function 'listAttributes' to get valid attribute names
#                
#                REMOVE NA PLATFORMS


for (n in seq(1,nrow(FOR_ANNOT_PLATFORMS))){
  featureDat <- getBM(attributes = c(FOR_ANNOT_PLATFORMS$Biomart[n], "external_gene_name"), 
                      filters = FOR_ANNOT_PLATFORMS$Biomart[n], 
                      values = LIST_DATA[[n]]$Probe_ID,
                      uniqueRows = F,
                      mart = 
                        if(FOR_ANNOT_PLATFORMS$Species[n] == "m"){ usedMartMUS } else if(FOR_ANNOT_PLATFORMS$Species[n] == "r"){ usedMartRAT }
  )
  }

for (n in seq(1,nrow(FOR_ANNOT_PLATFORMS))){

}

ifelse(FOR_ANNOT_PLATFORMS$Species[n]), yes, no)
FOR_ANNOT_PLATFORMS$Species[1]


#Tutaj trzeba zrobić tak: zrobić .txt files z każdego z eksperymentów. Przefiltrować odpowiednią bazę danych przez nazwy sond z tego eksperymentu.
#Potem jeszcze trzeba je sprowadzić do jednolitego nazewnictwa (mus or homo?)

featureDat <- getBM(attributes = c("agilent_wholegenome_4x44k_v1", "affy_mg_u74av2", "affy_mogene_1_0_st_v1"), 
                    mart = usedMartMUS)

featureDat <- getBM(attributes = c("affy_ragene_2_1_st_v1", "external_gene_name"), 
                    filters = "affy_ragene_2_1_st_v1", 
                    values = rownames(dataMatrixMW),
                    uniqueRows = TRUE,
                    mart = usedMart) #This is from Biomart. Here we download the names
xxx <- aggregate(external_gene_name~affy_ragene_2_1_st_v1, data = featureDat, c) #Here we collapse names, cause row names have to be unique and many probes correspond to more than 1 gene name in biomart

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
