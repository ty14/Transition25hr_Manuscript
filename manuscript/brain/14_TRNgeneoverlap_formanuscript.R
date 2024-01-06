library(tidyverse)

#TRN vs DOM 
td <- read_csv("manuscript/brain/results_tables/TRNgenesvsDOM_bothtimepoints.csv")
td

td70 <- td %>% filter(time == 70)

td70_up <- td70 %>% filter(dd_logFC > 0.2)
td70_up$symbol %>% unique(.) # 62 


td70_down <- td70 %>% filter(dd_logFC < 0.2)
td70_down$symbol %>% unique(.) # 66


td25 <- td %>% filter(time == 25)

td25_up <- td25 %>% filter(dd_logFC > 0.2)
td25_up$symbol %>% unique(.) # 28

td25_down <- td25 %>% filter(dd_logFC < 0.2)
td25_down$symbol %>% unique(.) # 32


td70x <- td70 %>% select(symbol, logFC = dd_logFC)
td25x <- td25 %>% select(symbol, logFC = dd_logFC)
source("functions/gettop10GO.R")

gettop10GO(td70x, my_showCategory) %>% 
  mutate(comparison = "TRN-DOM 70min") -> top10go1

gettop10GO(td25x, my_showCategory ) %>% 
  mutate(comparison = "TRN-DOM 25hr") -> top10go2

rbind(top10go1,top10go2) -> top10_GOterms

write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_TRNvsDOM_both_time_points.csv", row.names = F)


# TRN vs SUB 

ts <- read_csv("manuscript/brain/results_tables/TRNgenesvsSUB_bothtimepoints.csv")
ts

ts70 <- ts %>% filter(time == 70)

ts70_up <- ts70 %>% filter(ds_logFC > 0.2)
ts70_up$symbol %>% unique(.) # 150



ts70_down <- ts70 %>% filter(ds_logFC < 0.2)
ts70_down$symbol %>% unique(.) # 140


ts25 <- ts %>% filter(time == 25)

ts25_up <- ts25 %>% filter(ds_logFC > 0.2)
ts25_up$symbol %>% unique(.) # 105

ts25_down <- ts25 %>% filter(ds_logFC < 0.2)
ts25_down$symbol %>% unique(.) # 73


ts70x <- ts70 %>% select(symbol, logFC = ds_logFC)
ts25x <- ts25 %>% select(symbol, logFC = ds_logFC)
source("functions/gettop10GO.R")

gettop10GO(ts70x, my_showCategory) %>% 
  mutate(comparison = "TRN-SUB 70min") -> top10go1

gettop10GO(ts25x, my_showCategory ) %>% 
  mutate(comparison = "TRN-SUB 25hr") -> top10go2

rbind(top10go1,top10go2) -> top10_GOterms

write.csv(top10_GOterms,"manuscript/brain/results_tables/topBP_GOterms_TRNvsSUBth_time_points.csv", row.names = F)

