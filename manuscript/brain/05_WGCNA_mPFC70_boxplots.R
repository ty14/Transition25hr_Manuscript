#boxplots for 70 min WGCNA.

library(WGCNA)
library(tidyverse)

coldata <- read.csv("manuscript/brain/results_tables/coldata.csv", row.names = 1)
# get rid of controls
coldata <- coldata %>% filter(condition != "control")


MEs <- readRDS("manuscript/brain/results/WGCNA_MEs_mPFC70_Power5.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:26, names_to = "Module") %>% 
  full_join(coldata) %>% filter(Module != "MEgrey") %>% 
  mutate(status = condition1)

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)




#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
source("functions/newggtheme.R")