#boxplots for 70 min WGCNA.

library(WGCNA)
library(tidyverse)

coldata <- read.csv("manuscript/brain/results_tables/coldata_ALLXAGG.csv")
# get rid of controls
coldata <- coldata  %>% filter(time != 70)

fix <- coldata %>% filter(condition1 == "SUB")

fix$CORT <- ifelse(is.na(fix$CORT), mean(fix$CORT, na.rm = T), fix$CORT)
fix$wt_d8 <- ifelse(is.na(fix$wt_d8), mean(fix$wt_d8, na.rm = T), fix$wt_d8)
fix$post.given1 <- ifelse(is.na(fix$post.given1), mean(fix$post.given1, na.rm = T), fix$post.given1)
fix$post.received1 <- ifelse(is.na(fix$post.received1), mean(fix$post.received1, na.rm = T), fix$post.received1)
amy1 <- coldata %>% full_join(fix) 
coldatax <- amy1[c(1:15,17:23),]


# b2.5.1


MEs <- readRDS("manuscript/brain/results/WGCNA_MEs_mPFC25_Power4.RDS")

ME_df <-MEs%>% data.frame() %>% 
  tibble::rownames_to_column(var = "SampleID") %>%
  pivot_longer(cols = 2:17, names_to = "Module") %>% 
  full_join(coldatax) %>% filter(Module != "MEgrey") %>% 
  mutate(status = condition1) %>% filter(SampleID != "b2.5.1")

head(ME_df)

ME_df$Module <- gsub("ME", "", ME_df$Module)


table(ME_df$condition1)

ME_df$condition2 <- ifelse(ME_df$condition1 == "SUB", "SD", ME_df$condition1)
ME_df$condition2 <- ifelse(ME_df$condition2 == "DOM", "SD", ME_df$condition2)

#boxplot of eigengenes and status 
source("functions/geom_boxjitter.R")
source("functions/newggtheme.R")



ME_df %>% filter(Module %in% c("purple", "salmon")) %>% 
  ggplot(aes(condition2, value, fill = condition2, color = condition2))+
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


# TRN: black, blue, pink, purple, red, salmon,turquoise
lm.df <- ME_df %>% pivot_wider(names_from = Module, values_from = value)
lm.df$condition3 <- factor(lm.df$condition2, c("DES", "SD", "ASC"))
colnames(lm.df)
#DES: grey60 
g6 <- lmer(turquoise ~condition2 +(1|batch), data = lm.df)
summary(g6)

g62 <- lmer(turquoise~condition3 +(1|batch), data = lm.df)
summary(g62)

#purple asc, 
#salmon des
ME_df %>% 
  ggplot(aes(CORT, value, fill = status, color = status))+
  geom_point()+   stat_smooth(method = "lm", se =F)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "CORT",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol =4) +newggtheme_with_legends+ 
  theme(text =element_text(size =10))

colnames(ME_df)
ME_df %>% filter(condition1 %in% c("DES", "SUB")) %>% 
  ggplot(aes(post.received1, value, fill = status, color = status))+
  geom_point()+   stat_smooth(method = "lm", se =F)+
  scale_color_manual(values = viridis::viridis(4))+
  scale_fill_manual(values = viridis::viridis(4))+         
  labs(x = "rate of agg. received",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol =4) +newggtheme_with_legends+ 
  theme(text =element_text(size =10))




ME_df %>% filter(condition1 %in% c("ASC", "DOM")) %>% 
  ggplot(aes(post.given1, value, fill = status, color = status))+
  geom_point()+   stat_smooth(method = "lm", se =F)+
  scale_color_manual(values = viridis::viridis(2))+
  scale_fill_manual(values = viridis::viridis(2))+         
  labs(x = "rate of agg. given",
       y = "Module eigengene")+
  facet_wrap(~Module, ncol =4) +newggtheme_with_legends+ 
  theme(text =element_text(size =10))





## linear models and dot plot
library(WGCNA)
library(tidyverse)

set.seed(312)
datExpr <- readRDS("manuscript/brain/results/WGCNA_datExpr25.RDS") 
net <- readRDS("manuscript/brain/results/WGCNA_net_mPFC25_Power4_minmod75.RDS")
MEs <- readRDS("manuscript/brain/results/WGCNA_MEs_mPFC25_Power4.RDS")

MEs 

MEsx <- MEs[,c(1:15)]

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
moduleNumber = length(unique(moduleColors))
rownames(datExpr)
colnames(MEsx)<- gsub("ME", "", colnames(MEsx))

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes

# get linear model data 

coldatax$SampleID

#agg recieved 
cd <- coldatax %>% filter(condition1 != "ASC") %>% filter(condition1 != "DOM")


cd$condition1 <- factor(cd$condition1, levels = c("DES", "SUB"))
# cd <- cd %>% rownames_to_column(., var = "SampleID")
cd$pre_id <- str_sub(cd$pre_idbatchcage,1,1)
cd$post_id <- str_sub(cd$post_idbatch, 1,1)

MEs0
orderMEs(MEs0) %>% 
  rownames_to_column("SampleID") %>%
  full_join(cd) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  relocate(condition1, batch, pre_id, post_id, CORT, Postds,post.given1, post.received1 ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -pre_idbatchcage,
                -time, -Prerank, -wt_d8, -Postrank) %>% na.omit(.)-> ME_df

lm_result_list <- list()

ME_df$condition1
library(lme4)
library(lmerTest)

for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  lmer(module ~ condition1*post.received1 +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod1
  summary(mod1)
  df
  rbind(summary(lm(df[,9] ~ df[,1]*df[,8]))$coefficients[2,],
        summary(lm(df[,9] ~ df[,1]*df[,8]))$coefficients[3,],
        summary(lm(df[,9] ~ df[,1]*df[,8]))$coefficients[4,]) %>%
    as.data.frame() %>%
    cbind(key = c("DES-SUB","post.recieved1", "interaction")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}


lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey")  %>% format(., scientific = F) -> d

d %>% filter(`Pr(>|t|)`< 0.08)

#agg
#yellow condition and agg
# darkgrey royalblue interaction with agg


#agg give
cs <- coldatax %>% filter(condition1 != "DES") %>% filter(condition1 != "SUB")


cs$condition1 <- factor(cs$condition1, levels = c("ASC", "DOM"))
# cd <- cd %>% rownames_to_column(., var = "SampleID")
cs$pre_id <- str_sub(cs$pre_idbatchcage,1,1)
cs$post_id <- str_sub(cs$post_idbatch, 1,1)

MEs0
orderMEs(MEs0) %>% 
  rownames_to_column("SampleID") %>%
  full_join(cs) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  relocate(condition1, batch, pre_id, post_id, CORT, Postds,post.given1, post.received1 ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -pre_idbatchcage,
                -time, -Prerank, -wt_d8, -Postrank) %>% na.omit(.)-> ME_df

lm_result_list <- list()

ME_df$condition1
library(lme4)
library(lmerTest)

for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  lmer(module ~ condition1*post.given1 +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod1
  summary(mod1)
  df
  rbind(summary(lm(df[,9] ~ df[,1]*df[,7]))$coefficients[2,],
        summary(lm(df[,9] ~ df[,1]*df[,7]))$coefficients[3,],
        summary(lm(df[,9] ~ df[,1]*df[,7]))$coefficients[4,]) %>%
    as.data.frame() %>%
    cbind(key = c("ASC-DOM","post.given1", "interaction")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}


lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey")  %>% format(., scientific = F) -> d
d %>% filter(`Pr(>|t|)`< 0.08)

# a


#CORT linear models 
data <- coldata 


data$condition1 <- factor(data$condition1, levels = c("ASC", "DOM", "DES", "SUB"))
data$condition2 <- factor(data$condition1, levels = c("DOM", "ASC", "DES", "SUB"))
data$condition3 <- factor(data$condition1, levels = c("DES", "ASC", "DOM", "SUB"))
# cd <- cd %>% rownames_to_column(., var = "SampleID")
data$pre_id <- str_sub(data$pre_idbatchcage,1,1)
data$post_id <- str_sub(data$post_idbatch, 1,1)

MEs0
orderMEs(MEs0) %>% 
  rownames_to_column("SampleID") %>%
  full_join(data) %>% 
  mutate(batch = as.factor(batch)) %>% 
  mutate_if(is.numeric,scale) %>% 
  relocate(condition1, batch, pre_id, post_id, CORT,condition2, condition3, Preds ) %>% 
  dplyr::select(-SampleID, -post_idbatch, -pre_idbatchcage,
                -time, -Prerank, -wt_d8, -Postrank) %>% na.omit(.)-> ME_df

lm_result_list <- list()

ME_df$condition1
library(lme4)
library(lmerTest)

for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  lmer(module ~ condition1*CORT +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod1
  lmer(module ~ condition2*CORT +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod2
  lmer(module ~ condition3*CORT +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod3
  summary(mod1)
  summary(mod2)
  summary(mod3)
  df
  rbind(summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[5,],
        summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[6,],
        summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[7,],
        summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[8,],
        summary(lm(df[,9] ~ df[,6]*df[,5]))$coefficients[7,],
        summary(lm(df[,9] ~ df[,6]*df[,5]))$coefficients[8,],
        summary(lm(df[,9] ~ df[,7]*df[,5]))$coefficients[8,]) %>%
    as.data.frame() %>%
    cbind(key = c("CORT","ASC-DOM:CORT", "ASC-DES:CORT", 'ASC-SUB:CORT', "DOM-DES:CORT", "DOM-SUB:CORT", "DES-SUB:CORT")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}


lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey")  %>% format(., scientific = F) -> d
d %>% filter(`Pr(>|t|)` < 0.05)



#condition 

data
for(x in 1:length(MEs0)){
  k = x + 8
  ME_df[,c(1:8,k)] -> df
  md <- gsub("ME","",colnames(df)[9])
  colnames(df)[9] <- "module"
  lmer(module ~ condition1 +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod1
  lmer(module ~ condition2 +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod2
  lmer(module ~ condition3 +(1|batch)+(1|pre_id)+(1|post_id) , data = df) -> mod3
  summary(mod1)
  summary(mod2)
  summary(mod3)
  df
  rbind(summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[2,],
        summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[3,],
        summary(lm(df[,9] ~ df[,1]*df[,5]))$coefficients[4,],
        summary(lm(df[,9] ~ df[,6]*df[,5]))$coefficients[3,],
        summary(lm(df[,9] ~ df[,6]*df[,5]))$coefficients[4,],
        summary(lm(df[,9] ~ df[,7]*df[,5]))$coefficients[4,]) %>%
    as.data.frame() %>%
    cbind(key = c("ASC-DOM", "ASC-DES", 'ASC-SUB', "DOM-DES", "DOM-SUB", "DES-SUB")) %>%
    mutate(module = md) -> lm_result_list[[x]]
}


lm_result_list %>% 
  do.call(rbind,.)%>% 
  filter(.,module != "grey")  %>% format(., scientific = F) -> d
d %>% filter(`Pr(>|t|)` < 0.06)


