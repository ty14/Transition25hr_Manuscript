library(tidyverse)


t_mea <- read_csv("manuscript/brain/results_tables/TransitionGenes_MEA.csv")
  
t_d <- read_csv("manuscript/brain/results_tables/TRNgenesvsDOM_bothtimepoints.csv")

t_s <- read_csv("manuscript/brain/results_tables/TRNgenesvsSUB_bothtimepoints.csv")


t_mea$symbol[t_mea$symbol %in% t_d$symbol]
# [1] "Dennd1a"
t_mea$symbol[t_mea$symbol %in% t_s$symbol]
# [1] "Lamb1"  "Alkbh1" "Klhl4" 


t_mea %>% filter(symbol %in% c("Dennd1a", "Lamb1", "Alkbh1", "Klhl4")) # up in TRN, expect "Klhl4" 
t_d %>% filter(symbol %in% c("Dennd1a")) # down in TRN vs DOM 70
t_s %>% filter(symbol %in% c("Lamb1", "Alkbh1", "Klhl4")) #Alkbh1, Klh14 Up TRN , lamb1 down in TRN
                 