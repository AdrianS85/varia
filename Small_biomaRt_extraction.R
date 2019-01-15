#library(tidyverse)
#library(biomaRt)

###### Load input (single column with header) ######
PROBES <- readr::read_tsv("test adnotacji wg ensemble.txt")


###### Extract annotations from ensembl ######
ANNOT_PROBES <- biomaRt::getBM(attributes = c("illumina_mousewg_6_v2", "external_gene_name"), 
                              filters = "illumina_mousewg_6_v2", 
                              values = PROBES,
                              uniqueRows = F,
                              mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
)

###### Remove duplicates and aggregate data ######
UNIQ_ANNOT_PROBES <- aggregate(external_gene_name~illumina_mousewg_6_v2, data = unique(ANNOT_PROBES), FUN = c)
UNIQ_ANNOT_PROBES$external_gene_name <- as.character(lapply(X = UNIQ_ANNOT_PROBES$external_gene_name, FUN = paste, collapse = "; ")) 

###### Write output ######
readr::write_tsv(x = UNIQ_ANNOT_PROBES, path = "output.txt")
