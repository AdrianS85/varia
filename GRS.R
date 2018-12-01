library(tidyverse)

COMPARISONS <- readr::read_tsv("comparisons.txt", col_types = "cccccccccccccccc")

PRE_ORG_DATA <- readr::read_tsv("test_dataset.txt", col_types = "ccccnnnnncccc", locale = locale(decimal_mark = ","))

ORG_DATA <- PRE_ORG_DATA %>%
  dplyr::mutate(ID_NUM = row.names(PRE_ORG_DATA))



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

# Tutaj liczymy ile razy geny występują w oryginalnym dataset, nie patrząc czy są up czy down
WHOLE_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(Gene_symbol) %>%
  group_by(Gene_symbol) %>%
  summarise(number = n())



# Tutaj liczymy ile razy geny występują w oryginalnym dataset, patrząc czy są up czy down
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



# Tutaj liczymy ile razy geny występują W KAŻDYM Z PORÓWNAŃ, nie patrząc czy są up czy down
COMP_NO_UNIDS_ORG_DATA <- NO_UNIDS_ORG_DATA %>%
  select(!!P_or_C, Gene_symbol) %>%
  group_by(!!P_or_C, Gene_symbol) %>%
  summarise(number = n())



# Tutaj liczymy ile razy geny występują W KAŻDYM Z PORÓWNAŃ, patrząc czy są up czy down    
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
