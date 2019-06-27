library(tidyverse)

setwd("/media/adrians/USB DISK1/Projekty/W - MWieczorek/qPCR/Analysis")
setwd("E:/Projekty/W - MWieczorek/qPCR/Analysis")

test <- readr::read_tsv("gene_tests.tsv")

melt_test <- tidyr::gather(test, key = "Repeat", value = "CT", rep_1, rep_2, rep_3)

# Just a quick look at the density plot of CTs between repeats
melt_test %>% 
  ggplot(aes(x = CT, color = Repeat))+
  geom_density()+
  facet_wrap(~Gene)




CT_minus_CTmean_melt_test <- melt_test %>% 
  group_by(Gene, Repeat) %>% 
  dplyr::mutate(x = mean(CT, na.rm = T)) %>%
  ungroup() %>%
  dplyr::mutate(CT_minus_CTmean = CT - x)


## "Normalized" density plot of CTs between repeats
CT_minus_CTmean_melt_test %>% 
  ggplot(aes(x = CT_minus_CTmean, color = Repeat))+
  geom_density()+
  facet_wrap(~Gene)
