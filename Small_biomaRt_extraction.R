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
colnames(UNIQ_ANNOT_PROBES) <- c("Probe_ID", "ens_biomart_annot")


###### Write output ######
readr::write_tsv(x = UNIQ_ANNOT_PROBES, path = "output.txt")


###### Load annotation files from gemma and illumina ######
GEMMA_ANNOT <- readr::read_tsv(file = "GPL6887_noParents.an.txt", col_names = c("Probe_ID", "gemma_annot"), col_types = "cc", skip = 8)
ILLUM_ANNOT <- readr::read_tsv(file = "MouseWG-6_V2_0_R3_11278593_A.txt", col_names = c("illum_annot", "Probe_ID"), col_types = "-----------c-c", skip = 9)


###### Merge all annotation datasets ######
LIST_ANNOT <- list(PROBES, UNIQ_ANNOT_PROBES, GEMMA_ANNOT, ILLUM_ANNOT)
test <- Reduce(function(x, y) merge(x, y, by = "Probe_ID", all.x = TRUE), LIST_ANNOT) #https://stackoverflow.com/questions/14096814/merging-a-lot-of-data-frames
