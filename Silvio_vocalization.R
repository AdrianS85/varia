#setwd("E:/Projekty/SV - Silvio Vocalization")
setwd("/media/adrians/USB DISK1/Projekty/SV - Silvio Vocalization")

library(Hmisc)
library(readxl)
library(tidyverse)
library(rlist)


#path = "USVs BTBR-C57 - SONATA.xlsx"
path = "rep2_USV-DATA_male_intruder.xls"

xl_list <- lapply(readxl::excel_sheets(path), readxl::read_excel, path = path)

named_xl_list <- purrr::map2(.x = xl_list, .y = readxl::excel_sheets(path), .f = ~ dplyr::mutate(.x, mouse = .y ) )

named_xl_df <- rlist::list.rbind(.data = named_xl_list)

grouped_named_xl_df <- named_xl_df %>%
  dplyr::mutate( group = stringr::str_replace(mouse, "[0-9]+", "") )



# We should eliminate all calls having „freq centre” equal to zero. These are noises.
fcNotZero_grouped_named_xl_df <- grouped_named_xl_df %>%
  dplyr::filter(`freq centre` %nin% 0)



# Codify each call basing on its freq start/centre/end, intended as less (<), more (>) or equal (=) frequencies compared to the previous.
# For example: freq start=73700, centre=73200, end=71200, this means that this call is descending, so called „downward”. Other calls can be ascending (upward), or „chevron” when first going up and then down. The specific parameters is frequency modulation. We should find a numer which express the ability of the pup to modulate the frequencies of its call. Such as a modulation score? I dont know if I explain properly. But prlease, let’s deep on this frequency modulation.
coded_fcNotZero_grouped_named_xl_df <- fcNotZero_grouped_named_xl_df %>%
  mutate( code = case_when(
    `freq start` > `freq centre` & `freq centre` > `freq end` ~ "downward",
    `freq start` < `freq centre` & `freq centre` < `freq end` ~ "upward",
    `freq start` < `freq centre` & `freq centre` > `freq end` ~ "chevron",
    `freq start` > `freq centre` & `freq centre` < `freq end` ~ "reversed_chevron",
    `freq start` == `freq centre` & `freq centre` == `freq end` ~ "flat",
    TRUE ~ "other"
  ) )

counting_coded_for_sample <- coded_fcNotZero_grouped_named_xl_df %>%
  group_by(group, mouse, code) %>%
  dplyr::count()

coded_fcNotZero_grouped_named_xl_df %>%
  mutate(code = factor(code, levels = c("downward", "upward", "chevron", "reversed_chevron", "flat", "other"))) %>% 
  group_by(group, mouse, code) %>%
  ggplot(aes(x = mouse, fill = group))+
  geom_bar()+
  facet_grid(code ~ .)+
  theme_light()



# Duration of the calls: short long?
call_duration_characteristics <- coded_fcNotZero_grouped_named_xl_df %>%
  group_by(group, mouse) %>%
  dplyr::summarize(mean_duration = mean(duration), duration_sd = sd(duration))

coded_fcNotZero_grouped_named_xl_df %>%
  group_by(group, mouse) %>%
  ggplot(aes(x = mouse, y = duration, fill = group))+
  geom_boxplot()+
  theme_light()


readr::write_tsv(grouped_named_xl_df, path = paste0("grouped_named_xl_df_", path, ".tsv") )





## Here, we compare our call names with manually established call names
compare_coded_fcNotZero_grouped_named_xl_df <- coded_fcNotZero_grouped_named_xl_df

#to merge the same names inputed differently
compare_coded_fcNotZero_grouped_named_xl_df$`class of call` <- tolower(compare_coded_fcNotZero_grouped_named_xl_df$`class of call`)


counting_coded_compare_coded_fcNotZero_grouped_named_xl_df <- compare_coded_fcNotZero_grouped_named_xl_df %>%
  group_by(`class of call`, code) %>%
  dplyr::count()


counting_coded_compare_coded_fcNotZero_grouped_named_xl_df %>%
  #mutate(code = factor(code, levels = c("downward", "upward", "chevron", "reversed_chevron", "flat", "other"))) %>% 
  #group_by(group, mouse, code) %>%
  ggplot(aes(x = code, y = n, fill = `class of call`))+
  geom_col()+
  facet_wrap(~`class of call`)+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
## Here, we compare our call names with manually established call names
