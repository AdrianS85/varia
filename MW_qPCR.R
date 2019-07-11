setwd("/media/adrians/USB DISK1/Projekty/W - MWieczorek/qPCR/Analysis")

library(tidyverse)

qpcrResults <- readr::read_tsv("qpcr_results.tsv")

tidy_qpcrResults <- qpcrResults %>%
  tidyr::gather(key = "Gene", value = "Cq", -c(Group, Name) )

tidy_qpcrResults %>%
  ggplot(aes(x = Group, y = Cq, fill = Group))+
  geom_boxplot()+
  facet_wrap(~Gene)+
  theme_light()





### Lets look at means and SDs
meanSD_tidy_qpcrResults <- tidy_qpcrResults %>%
  select(-Name) %>%
  group_by(Group, Gene) %>%
  nest() %>%
  mutate(
    Srednia = data %>% map_dbl(map_dbl, mean),
    Od_std = data %>% map_dbl(~.[[1]] %>% sd)) %>% ## THIS WILL BE MY MAPPING METHOD OF CHOICE FOR NOW
  unnest()



### Here we make SW test for normality
SW_tidy_qpcrResults <- tidy_qpcrResults %>%
  select(-Name) %>%
  dplyr::group_by(Gene, Group) %>%
  nest() 

library(magrittr)
SW_results <- SW_tidy_qpcrResults %$%
  map(data, ~ shapiro.test(.x[[1]]))

SW_tidy_qpcrResults <- SW_tidy_qpcrResults %>%
  dplyr::mutate(Shapiro_Wilk_p_val = 
                  map_chr(SW_results, ~ .x$p.value ) 
                )



### Here we find out whether any SW analysis gave significant value, meaning that distribution of Cqs in given group in given gene is not normal. 
is_gene_normal <- SW_tidy_qpcrResults %>%
  select( -c(data, Group) ) %>%
  group_by(Gene) %>%
  nest() %>%
  mutate( gene_is_normal = data %>% map(~ case_when(
    .x < 0.05  ~ FALSE,
    TRUE ~ TRUE
  ) ) ) %>%
  mutate( gene_is_normal_short = 
            gene_is_normal %>% map( ~ !any(.x) ) 
        )


### Statistics:
###### Rank Rs within given gene
ranked_tidy_qpcrResults <- tidy_qpcrResults %>%
  select(-Name) %>%
  group_by(Gene) %>%
  mutate( ranked_Cq = rank(Cq) )

###### ANOVA this shit
nested_ranked_tidy_qpcrResults <- ranked_tidy_qpcrResults %>%
  group_by(Gene) %>%
  nest()

ANOVA_ <- map(.x = nested_ranked_tidy_qpcrResults$data, 
              ~ broom::tidy(aov(ranked_Cq ~ Group, data = .x)) )

TUKEY_ <- map(.x = nested_ranked_tidy_qpcrResults$data, 
              ~ broom::tidy(TukeyHSD(aov(ranked_Cq ~ Group, data = .x))) )

ANOVAd_nested_ranked_tidy_qpcrResults <- nested_ranked_tidy_qpcrResults %>% 
  mutate(ANOVA = ANOVA_)

###### Should I perform post-hoc test?
#Sprawdza czy pvalue w ANOVIE jest większe od 0.05 albo nie istnieje i daje FALSE to tych elementów
ShouldIDunnet <- ANOVA_ %>% map_lgl(~ if (is.null(.x)) { F } else if(.x[[1,6]] < 0.05) { T } else { F })
ShouldIDunnet_ANOVAd_nested_ranked_tidy_qpcrResults <- ANOVAd_nested_ranked_tidy_qpcrResults %>% 
  mutate(Dunnet = ShouldIDunnet)



ANOVA_ <-  map(.x = Dawid4_filtered, ~ tidy(TukeyHSD(aov(Stosunek ~ Stadium, data = .x))) )
Dawid4_filtered_A_D$xxx <- Dawid4_filtered_A_D$ANOVA %>% map(tidy)





kruskal_tidy_qpcrResults <- tidy_qpcrResults %>%
  select(-Name) %>%
  group_by(Gene) %>%
  nest() %>%
  dplyr::mutate(Shapiro_Wilk_p_val = data %>%
                  map( ~ kruskal.test( .x$Cq ~ as.factor(.x$Group) ) ) 
  )




rm(xxx)



### Here we need to do correlation analysis, perhaps something good comes out of it
maResults <- readr::read_tsv("ma_results.tsv")


meanedRab1_maResults <- maResults %>%
  mutate( Rab1b = 
            as.numeric( purrr::map2( .x = Rab1b, .y = Rab1b_1, .f = ~ mean(c(.x, .y)) ) ),
          Rab1b_1 = NULL
  )


tidy_meanedRab1_maResults <- meanedRab1_maResults %>%
  select(-Name_2) %>%
  tidyr::gather(key = "Gene", value = "Exp_value", -c(Group, Name) )


mergin_tidy_meanedRab1_maResults <- tidy_meanedRab1_maResults %>%
  mutate( ID =  paste0(tidy_meanedRab1_maResults$Name, "_", tidy_meanedRab1_maResults$Gene)) %>%
  select( -c(Group, Name, Gene) )


merged_qpcr_ma <- merge(mergin_tidy_qpcrResults, mergin_tidy_meanedRab1_maResults, by = "ID") %>%
  select(-ID)


forCorrelations_merged_qpcr_ma <- merged_qpcr_ma %>%
  select(Gene, Cq, Exp_value) %>%
  group_by(Gene) %>%
  nest()
  

Correlations <- map(forCorrelations_merged_qpcr_ma$data, ~ corrr::correlate( .x[,1:2]) )  


bleble_forCorrelations_merged_qpcr_ma <- forCorrelations_merged_qpcr_ma %>%
  mutate(Cor_value = 
    map_dbl(Correlations, ~ .x[[2,2]])
  )



merged_qpcr_ma %>%
  ggplot(aes(x = Exp_value, y = Cq, color = Gene))+
  geom_point()+
  theme_light()+
  facet_wrap(~Gene)
