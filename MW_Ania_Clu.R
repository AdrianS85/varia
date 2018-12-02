https://bioconductor.org/packages/release/bioc/html/multiClust.html   
clustComp
ClusterJudge
clusterSeq
cogena
ConsensusClusterPlus
debrowser
geneRecommender
iterClust
M3C
MantelCorr
MCbiclust
QUBIC
rqubic
runibic
TTMap


FGNet
clusterProfiler
meshes

https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-0002-introduction-to-computational-thinking-and-data-science-fall-2016/lecture-videos/index.htm

BiocManager::install("multiClust", version = "3.8")    
BiocManager::install("runibic", version = "3.8")    

library(multiClust)
library(tidyverse)

setwd("E:/Projekty/W - MWieczorek/Results for MW/For Anna new clustering vagotomy/Nowe od Anny 1118/Clustering")



COMPARISONS <- readr::read_tsv("OnlyVagxClusteringF.txt", col_types = "cnnnnnnnnnnnnnnnnnnnnn")
## Tutaj nadal są zduplikowane geny
COMP <- COMPARISONS %>%
  mutate(Gene_Symbol = str_remove(Genes, ".* ")) %>%
  mutate(ID = str_remove(Genes, " .*")) #%>%
  #filter(Gene_Symbol != "NA")


## Tutaj mamy plik z sondami z nowymi anotacjami
NEW_ANNOT_PROBE <- readr::read_tsv("From OnlyVagotomy v ania k.txt", col_types = "cc")


## Niestety trzeba jeszcze zanotować COMPARISONS with newer probe names
TEST_CLU <- readr::read_tsv("CIV.txt", col_types = "cc")
# Tutaj dodamy sondy do genów
T_CLU <- merge(TEST_CLU, NEW_ANNOT_PROBE, by = "Gene_Symbol", all.x = T)


COMP_check_TEST_CLU <- merge(T_CLU, COMP, by = "ID", all.x = T) %>%
  mutate(ID_Symbol = paste(Gene_Symbol.x, ID, sep = "_"))

FINAL_COMP_check_TEST_CLU <- COMP_check_TEST_CLU %>%
  mutate(
    ID = NULL,
    Gene_Symbol.x = NULL,
    Gene_Symbol.y = NULL,
    Cluster = NULL,
    Genes = NULL
  ) %>%
  select(ID_Symbol, everything())
