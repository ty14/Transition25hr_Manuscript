#boxplots for 70 min WGCNA.

library(WGCNA)
library(tidyverse)

coldata <- read.csv("manuscript/brain/results_tables/coldata.csv")
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



ME_df %>% 
  ggplot(aes(status, value, fill = status, color = status))+
  geom_boxjitter(outlier.color = NA, jitter.shape = 21, jitter.color = NA,
                 alpha = 0.3,
                 jitter.height = 0.02, jitter.width = 0.07, errorbar.draw = TRUE,
                 position = position_dodge(0.85)) +
  # facet_wrap(~Module, ncol =5)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Social condition",
       y = "Module eigengene")+
  facet_wrap(~Module) +newggtheme_with_legends+ 
  theme(legend.position = "none", text =element_text(size =10))




ME_df %>% 
  ggplot(aes(post_Ncort, value, fill = status, color = status))+
 geom_point()+   stat_smooth(method = "lm", se =F)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Normalized CORT",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol =4) +newggtheme_with_legends+ 
  theme(text =element_text(size =10))

ME_df %>% 
  ggplot(aes(received1, value, fill = status, color = status))+
  geom_point()+   stat_smooth(method = "lm", se =F)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "rate of agg. given",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol =4) +newggtheme_with_legends+ 
  theme(text =element_text(size =10))


ME_df <- ME_df %>% mutate(tot_agg = given1+received1)

ME_df %>% 
  ggplot(aes(tot_agg, value, fill = status, color = status))+
  geom_point()+   stat_smooth(method = "lm", se =F)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "Rate of Aggression",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol =4) +newggtheme_with_legends+ 
  theme(text =element_text(size =10))

