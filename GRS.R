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

# Jeżeli w ensemble nie ma przypisanego genu (oznaczenie NA w przysłanym pliku "100_genow") to program pokazuje we wszystkich rubrykach oznaczenie NA chociaż te komórki zawierały wcześniej dane (plik AAAreview11). W pełnych danych które przysłałeś mi przed świętami (plik 1 znajdujący się w skopresowanych danych "TEST_ANNOTATION") tych sond w ogóle nie ma. Wygląda na to, że to co przysłałeś wczoraj odnosi się chyba do wcześniejszego pliku w którym dane dla sond nie mających przypisanego genu w ensemble zostały całkowicie usunięte a to nie jest dobre :( - AMS W LIST_DATA SĄ TE DANE, ALE W TEST_ANNOTATION JUZ NIE

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


