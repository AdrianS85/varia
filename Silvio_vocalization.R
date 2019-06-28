setwd("E:/Projekty/SV - Silvio Vocalization")

library(readxl)
library(tidyverse)
library(rlist)

path = "USVs BTBR-C57 - SONATA.xlsx"

xl_list <- lapply(excel_sheets(path), read_excel, path = path)

named_xl_list <- map2(.x = xl_list, .y = excel_sheets(path), .f = ~ mutate(.x, mouse = .y ) )

named_xl_df <- rlist::list.rbind(.data = named_xl_list)

grouped_named_xl_df <- named_xl_df %>%
  mutate( group = str_replace(mouse, "[0-9]+", "") )

readr::write_tsv(grouped_named_xl_df, path = "grouped_named_xl_df_USVs_BTB_C57_SONATA.tsv")

grouped_named_xl_df %>%
  ggplot(aes(x = `Call #`, y = mouse, color = group))+
  geom_point()

