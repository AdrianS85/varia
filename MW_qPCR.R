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



### Lets make Mann-Whitney-Wilcoxon
#wilcox.test(formula, data, subset) - needs 2 levels of factor for comparison

MWW_test <- wilcox.test(Actn2 ~ Group, data = qpcrResults, subset = Group == c("C", "S_1"))


MWW_test <- qpcrResults


MWW_test$C_vs_S_1 <-   lapply(X = qpcrResults[3:9], FUN = function(X) { wilcox.test(X ~ qpcrResults$Group, subset = qpcrResults$Group == c("C", "S_1")) })

MWW_test$C_vs_S_3 <-   lapply(X = qpcrResults[3:9], FUN = function(X) { wilcox.test(X ~ qpcrResults$Group, subset = qpcrResults$Group == c("C", "S_3")) })

MWW_test$C_vs_V_1 <-   lapply(X = qpcrResults[3:9], FUN = function(X) { wilcox.test(X ~ qpcrResults$Group, subset = qpcrResults$Group == c("C", "V_1")) })

MWW_test$C_vs_V_3 <-   lapply(X = qpcrResults[3:9], FUN = function(X) { wilcox.test(X ~ qpcrResults$Group, subset = qpcrResults$Group == c("C", "V_3")) })

MWW_test$S_1_vs_S_3 <-   lapply(X = qpcrResults[3:9], FUN = function(X) { wilcox.test(X ~ qpcrResults$Group, subset = qpcrResults$Group == c("S_1", "S_3")) })

MWW_test$S_1_vs_V_1 <-   lapply(X = qpcrResults[3:9], FUN = function(X) { wilcox.test(X ~ qpcrResults$Group, subset = qpcrResults$Group == c("S_1", "V_1")) })

MWW_test$S_3_vs_V_3 <-   lapply(X = qpcrResults[3:9], FUN = function(X) { wilcox.test(X ~ qpcrResults$Group, subset = qpcrResults$Group == c("S_3", "V_3")) })

MWW_test$V_1_vs_V_3 <-   lapply(X = qpcrResults[3:9], FUN = function(X) { wilcox.test(X ~ qpcrResults$Group, subset = qpcrResults$Group == c("V_1", "V_3")) })

# Write-out resulting analysis
write_lines(x = rjson::toJSON(MWW_test), path = "MWW_test.json")
### Lets make Mann-Whitney-Wilcoxon



### Lets make a nice list o pvalues from Mann-Whitney-Wilcoxon
MWW_pval <- list()
nb <- 1
for(n_comp in seq_along(MWW_test))
  {

  for(n_gene in seq_along(MWW_test[[n_comp]]))
    {
    MWW_pval$Comparison[[nb]] <- names(MWW_test[n_comp])
    MWW_pval$Gene[[nb]] <- names(MWW_test[[n_comp]][n_gene])
    MWW_pval$pval[[nb]] <- MWW_test[[n_comp]][[n_gene]]$p.value
    nb = nb + 1
    }
  
  }

df_MWW_pval <- as_tibble(MWW_pval)
### Lets make a nice list o pvalues from Mann-Whitney-Wilcoxon



### Lets statistically adjust Mann-Whitney-Wilcoxon tests
geneList_df_MWW_pval <- df_MWW_pval %>%
  group_by(Gene) %>%
  nest()

geneList_df_MWW_pval$data <- map(.x = geneList_df_MWW_pval$data, .f = 
    ~ mutate(.data = .x, MWW_adjusted = p.adjust(as.numeric(.x$pval), method = "BH")
    )
)

geneList_df_MWW_pval <- geneList_df_MWW_pval %>%
  unnest()
### Lets statistically adjust Mann-Whitney-Wilcoxon tests









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


ANOVA_ <- map(.x = nested_ranked_tidy_qpcrResults$data, 
              ~ broom::tidy(aov(ranked_Cq ~ Group, data = .x)) )




### Try contrasts - https://stats.stackexchange.com/questions/89021/how-to-get-only-desirable-comparisons-from-post-hoc
factored_ranked_tidy_qpcrResults <- ranked_tidy_qpcrResults
factored_ranked_tidy_qpcrResults$Group <- as.factor(factored_ranked_tidy_qpcrResults$Group)

levels(factored_ranked_tidy_qpcrResults$Group)

factors_ <- as.matrix(read.delim(file = "factors.txt", header = F, sep = "\t"))

contrasts(factored_ranked_tidy_qpcrResults$Group) <- factors_

factored_ranked_tidy_qpcrResults <- factored_ranked_tidy_qpcrResults %>%
  group_by(Gene) %>%
  nest()




for(n in seq_along(factored_nested_ranked_tidy_qpcrResults$data)){
  factored_nested_ranked_tidy_qpcrResults$data[[n]][[1]] <- as.factor(factored_nested_ranked_tidy_qpcrResults$data[[n]][[1]])
}


for(n in seq_along(factored_nested_ranked_tidy_qpcrResults$data)){
  contrasts(factored_nested_ranked_tidy_qpcrResults$data[[n]]$Group) <- factors_
}

p.adjust


contrasts(factored_nested_ranked_tidy_qpcrResults$data[[1]]$Group)
class(factored_nested_ranked_tidy_qpcrResults$data[[1]]$Group)
levels(factored_nested_ranked_tidy_qpcrResults$data[[1]]$Group)

LM_ <- map(.x = factored_nested_ranked_tidy_qpcrResults$data, 
           ~ summary( lm(ranked_Cq ~ Group, data = .x) ) )

ANOVA_ <- map(.x = factored_ranked_tidy_qpcrResults$data, 
              ~ broom::tidy(aov(ranked_Cq ~ Group, data = .x)) )

TUKEY_ <- map(.x = factored_nested_ranked_tidy_qpcrResults$data, 
              ~ broom::tidy(TukeyHSD(aov(ranked_Cq ~ Group, data = .x))) )

### Try contrasts - https://stats.stackexchange.com/questions/89021/how-to-get-only-desirable-comparisons-from-post-hoc




TUKEY_ <- map(.x = nested_ranked_tidy_qpcrResults$data, 
              ~ broom::tidy(TukeyHSD(aov(ranked_Cq ~ Group, data = .x))) )

ANOVAd_nested_ranked_tidy_qpcrResults <- nested_ranked_tidy_qpcrResults %>% 
  mutate(ANOVA = ANOVA_)

###### Should I perform post-hoc test?
#Sprawdza czy pvalue w ANOVIE jest wiÄ™ksze od 0.05 albo nie istnieje i daje FALSE to tych elementĂłw
ShouldIDunnet <- ANOVA_ %>% map_lgl(~ if (is.null(.x)) { F } else if(.x[[1,6]] < 0.05) { T } else { F })
ShouldIDunnet_ANOVAd_nested_ranked_tidy_qpcrResults <- ANOVAd_nested_ranked_tidy_qpcrResults %>% 
  mutate(Dunnet = ShouldIDunnet)



ANOVA_ <-  map(.x = ShouldIDunnet_ANOVAd_nested_ranked_tidy_qpcrResults$data, ~ broom::tidy(TukeyHSD(aov(ranked_Cq ~ Group, data = .x))) )
Dawid4_filtered_A_D$xxx <- Dawid4_filtered_A_D$ANOVA %>% map(tidy)





kruskal_tidy_qpcrResults <- tidy_qpcrResults %>%
  select(-Name) %>%
  group_by(Gene) %>%
  nest() %>%
  dplyr::mutate(Shapiro_Wilk_p_val = data %>%
                  map( ~ kruskal.test( .x$Cq ~ as.factor(.x$Group) ) ) 
  )







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

#############################################################################
testingActn2 <- merged_qpcr_ma %>%
  filter(Group %in% c("C", "S_1", "V_1", "V_3") ) %>%
  filter(Gene == "Actn2")

cor_testingActn2 <- corrr::correlate( testingActn2[,4:5], method = "spearman")


testingCxcl14 <- merged_qpcr_ma %>%
  filter(Group %in% c("S_3", "V_1", "V_3") ) %>%
  filter(Gene == "Cxcl14")

cor_testingCxcl14 <- corrr::correlate( testingCxcl14[,4:5], method = "spearman")


testingGabrg1 <- merged_qpcr_ma %>%
  filter(Group %in% c("C", "S_1", "V_1") ) %>%
  filter(Gene == "Gabrg1")

cor_testingGabrg1 <- corrr::correlate( testingGabrg1[,4:5], method = "spearman")


testingGria1 <- merged_qpcr_ma %>%
  filter(Group %in% c("C", "V_3", "S_1", "V_1") ) %>%
  filter(Gene == "Gria1")

cor_testingGria1 <- corrr::correlate( testingGria1[,4:5], method = "spearman")


testingRab1b <- merged_qpcr_ma %>%
  filter(Group %in% c("V_3", "S_3", "V_3", "C") ) %>%
  filter(Gene == "Rab1b")

cor_testingRab1b <- corrr::correlate( testingRab1b[,4:5], method = "spearman")


testingOxt <- merged_qpcr_ma %>%
  filter(Group %in% c("C", "V_3", "S_1", "V_1") ) %>%
  filter(Gene == "Oxt")

cor_testingOxt <- corrr::correlate( testingOxt[,4:5], method = "spearman")


testingVim <- merged_qpcr_ma %>%
  filter(Group %in% c("V_3", "S_3", "V_3", "C") ) %>%
  filter(Gene == "Vim")

cor_testingVim <- corrr::correlate( testingVim[,4:5], method = "spearman")
#############################################################################





forCorrelations_merged_qpcr_ma <- merged_qpcr_ma %>%
  select(Gene, Cq, Exp_value) %>%
  group_by(Gene) %>%
  nest()
  

Correlations <- map(forCorrelations_merged_qpcr_ma$data, ~ corrr::correlate( .x[,1:2], method = "spearman") )


bleble_forCorrelations_merged_qpcr_ma <- forCorrelations_merged_qpcr_ma %>%
  mutate(Cor_value = 
    map_dbl(Correlations, ~ .x[[2,2]])
  )



merged_qpcr_ma %>%
  ggplot(aes(x = Exp_value, y = Cq, color = Gene))+
  geom_point()+
  theme_light()+
  facet_wrap(~Gene)


