library(magrittr) ## For %>% operator
library(ggplot2)

global <- readr::read_tsv("methylation percent fede experiment2.txt", 
                   col_names = T, 
                   col_types = "cccn", 
                   locale = readr::locale(decimal_mark = ".")) ## GOOD FUNCTION locale

global %>% 
  dplyr::mutate(Group = factor(Group, levels = c("Control-NM", "EmbrioTransfer", "IVculture", "Biopsy"))) %>% 
  ggplot(aes(x = Group, y = percent_mCpG, color = Group))+
  geom_boxplot()+
  facet_wrap(~Tissue)+
  theme(axis.text.x = element_text(angle = 90))
scale_x_discrete()

library(tidyverse)



# Get mean and SD
mean.SD_global <- global %>%
  select(-Sample_Name) %>%
  group_by(Tissue, Group) %>%
  nest() %>%
  mutate(
    Srednia = data %>% map(map, mean),
    Od_std = as.numeric(data %>% map(~.[[1]] %>% sd)),
    Shap = data %>% map(map, shapiro.test) 
  ) %>%
  mutate(
    Mean = as.numeric(Srednia %>% map(~.[[1]], unlist)),
    p_Shapiro = Shap %>% map(map, ~.[[2]])
  ) %>% 
  mutate(
    p_Shapiro = as.numeric(p_Shapiro %>% map(~.[[1]], unlist))
  ) %>% 
  select(-Srednia, -Shap) 


# Get mean and SD
xxx_mean.SD_global <- mean.SD_global %>%
  unnest() %>%
  group_by(Tissue) %>%
  nest()

# Do ANOVA
ANOVA_ <- map(.x = xxx_mean.SD_global$data, ~ broom::tidy(aov(percent_mCpG ~ Group, data = .x)))
Tukey_ <- map(.x = xxx_mean.SD_global$data, ~ broom::tidy(TukeyHSD(aov(percent_mCpG ~ Group, data = .x))) )

xxx_mean.SD_global <- xxx_mean.SD_global %>%
  mutate(
    ANOVA = ANOVA_,
    Tukey = Tukey_
  )

rlist::list.save(x = xxx_mean.SD_global, file = "pvalues_global_meth.json")
