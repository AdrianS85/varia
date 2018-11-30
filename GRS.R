library(tidyverse)

PRE_ORG_DATA <- readr::read_tsv("test_dataset.txt", col_types = "ccccnnnnncccc", locale = locale(decimal_mark = ","))

ORG_DATA <- test %>%
  mutate(ID_NUM = row.names(test))



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
ORG_DATA__x__missing_processed_rows <- nrow(ORG_DATA) - (ORG_DATA__3__sum_of_NONEs_and_NAs + nrow(ORG_DATA_NO_UNIDS))

if (ORG_DATA__x__missing_processed_rows != 0){ 
  stop("Hey, buddy! You have some wierd values in Your raw data, guy! Better check whats happening, or Your results will smell of farts!")
}

return(ORG_DATA__x__missing_processed_rows)
}
### HERE WE CHECK IF ALL VALUES IN INPUT TABLE ARE EITHER NONES, NAS OR CORRECT VALUES



## Subset of table without "NONE" gene ids.
ORG_DATA_NO_UNIDS <- ORG_DATA %>%
  dplyr::filter(Gene_symbol != "NONE" & !is.na(Gene_symbol))

FUNCTION__1__check_rm_of_unid_val() # Check if it went well


after_filter <- as.tibble(c(test_NONE$ID_NUM, test_NONErm$ID_NUM))
colnames(after_filter) <- "ID_NUM"

xxx <- anti_join(testmutate, after_filter, by = "ID_NUM")

test_N <- replace_na(data = test$Gene_symbol, replace = "NONE")
length(unique(test$Probe_ID))  
xx <- is.na(test$Gene_symbol)

  dplyr::filter(Gene_symbol == is.na(Gene_symbol)) %>%
