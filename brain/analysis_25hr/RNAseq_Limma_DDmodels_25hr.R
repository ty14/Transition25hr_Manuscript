# libraries 
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(tidyverse)

# my data 25hrs Descenders and Dominant

mPFC <- read_csv("manuscript/brain/PFC_25hcounts.csv")

my_tissue = "mPFC"
mPFC  -> dlNorm
colnames(dlNorm)[1] <- "engene"

#Getting metadata ready 
coldata <- read_csv("manuscript/brain/sample25hr_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(coldata$post_Ncort)

##getting condition1
# Dom-Dom to Descenders (DOM to SUB)(4->1)
# Sub-Sub to Ascenders (Sub to DOM)  (1->4)
coldata <- coldata %>% filter(groupEX != "control")
table(coldata$condition)

coldata$condition <- ifelse(coldata$Postrank == 4 & coldata$Prerank == 4, "same", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "Dominant", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "Subordinate", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "Descenders", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "Ascenders", coldata$condition1)
table(coldata$condition1)


coldata %>% 
  filter(Postrank != 3) %>% 
  filter(Prerank != 3) %>% 
  filter(condition1 != "Ascenders") %>%
  filter(condition1 != "Subordinate") %>%
  filter(region != "A") %>% 
  dplyr::select(SampleNames, condition1,AggGiven70min, AggRec70min, post_Ncort) -> var_info

ifelse(var_info$condition1 == "Dominant", 1,-1) -> var_info$condition1
table(var_info$condition1)

#just get samples I want
row.names <- var_info$SampleNames

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

# Create DGEList object
do <- DGEList(dlNorm)
do
d0 <- calcNormFactors(do)
d0

d <- d0
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d) # number of genes left 14770    12

ExprMatrix <- d$counts
d$samples

samplenames<- substring(colnames(d),1, nchar(colnames(d)))
samplenames
colnames(d) <- samplenames


d$samples <-cbind(SampleNames=rownames(d$samples), d$samples)

d$samples <- d$samples %>% full_join(var_info)
rownames(d$samples) <- d$samples$SampleNames

d$samples <- d$samples %>% dplyr::select(-SampleNames)


##Model Comparison

# What selection model is the best. 
d1 <- model.matrix(~0+condition1,data=d$samples) #just condition bic3, aic 1
d2 <- model.matrix(~0+post_Ncort,data=d$samples)
d3 <- model.matrix(~0+AggGiven70min,data=d$samples)
d4 <- model.matrix(~0+AggRec70min,data=d$samples)
d5 <- model.matrix(~0+condition1+post_Ncort,data=d$samples) # bic 3, aic 12
d6 <- model.matrix(~0+condition1+post_Ncort+condition1*post_Ncort,data=d$samples) # best model 
d7 <- model.matrix(~0+condition1+post_Ncort+AggGiven70min,data=d$samples) 
d8 <- model.matrix(~0+condition1+post_Ncort+AggGiven70min+ condition1*post_Ncort,data=d$samples)
dlist <- list(d1,d2,d3,d4,d5,d6,d7,d8)

sm1 <- selectModel(ExprMatrix,dlist,criterion="aic")
sm2 <- selectModel(ExprMatrix,dlist,criterion="bic")
sm1
sm2

barplot(table(sm1$pref), main = "AIC")
barplot(table(sm2$pref), main = "BIC")
table(sm1$pref)
table(sm2$pref)


# my data 25hrs Descenders and Dominant

MeA <- read_csv("manuscript/brain/AMY_25hcounts.csv")

my_tissue = "MeA"
MeA  -> dlNorm
colnames(dlNorm)[1] <- "engene"

#Getting metadata ready 
coldata <- read_csv("manuscript/brain/sample25hr_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(coldata$post_Ncort)

##getting condition1
# Dom-Dom to Descenders (DOM to SUB)(4->1)
# Sub-Sub to Ascenders (Sub to DOM)  (1->4)
coldata <- coldata %>% filter(groupEX != "control")
table(coldata$condition)

coldata$condition <- ifelse(coldata$Postrank == 4 & coldata$Prerank == 4, "same", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "Dominant", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "Subordinate", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "Descenders", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "Ascenders", coldata$condition1)
table(coldata$condition1)


coldata %>% 
  filter(Postrank != 3) %>% 
  filter(Prerank != 3) %>% 
  filter(condition1 != "Ascenders") %>%
  filter(condition1 != "Subordinate") %>%
  filter(region == "A") %>% 
  dplyr::select(SampleNames, condition1,AggGiven70min, AggRec70min, post_Ncort) -> var_info

ifelse(var_info$condition1 == "Dominant", 1,-1) -> var_info$condition1
table(var_info$condition1)

#just get samples I want
row.names <- var_info$SampleNames

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

# Create DGEList object
do <- DGEList(dlNorm)
do
d0 <- calcNormFactors(do)
d0

d <- d0
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
dim(d) # number of genes left 14770    12

ExprMatrix <- d$counts
d$samples

samplenames<- substring(colnames(d),1, nchar(colnames(d)))
samplenames
colnames(d) <- samplenames


d$samples <-cbind(SampleNames=rownames(d$samples), d$samples)

d$samples <- d$samples %>% full_join(var_info)
rownames(d$samples) <- d$samples$SampleNames

d$samples <- d$samples %>% dplyr::select(-SampleNames)


##Model Comparison

# What selection model is the best. 
d1 <- model.matrix(~0+condition1,data=d$samples) #bic 0, aic 0 
d2 <- model.matrix(~0+post_Ncort,data=d$samples)
d3 <- model.matrix(~0+AggGiven70min,data=d$samples)
d4 <- model.matrix(~0+AggRec70min,data=d$samples)
d5 <- model.matrix(~0+condition1+post_Ncort,data=d$samples) #bic 8, aic 12
d6 <- model.matrix(~0+condition1+post_Ncort+condition1*post_Ncort,data=d$samples) # best model 
d7 <- model.matrix(~0+condition1+post_Ncort+AggGiven70min,data=d$samples) 
d8 <- model.matrix(~0+condition1+post_Ncort+AggGiven70min+ condition1*post_Ncort,data=d$samples)
dlist <- list(d1,d2,d3,d4,d5,d6,d7,d8)

sm1 <- selectModel(ExprMatrix,dlist,criterion="aic")
sm2 <- selectModel(ExprMatrix,dlist,criterion="bic")
sm1
sm2

barplot(table(sm1$pref), main = "AIC")
barplot(table(sm2$pref), main = "BIC")
table(sm1$pref)
table(sm2$pref)

