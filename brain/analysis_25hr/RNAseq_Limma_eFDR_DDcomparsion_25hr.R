# libraries 
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(pheatmap)
library(annotables)
grcm38 # mouse genes


library(tidyverse)

# my data 70 mins Descenders and Dominant


# Bring in count data for mPFC
a_countdata <- read_csv("manuscript/brain/PFC_25hcounts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
row.names(a_count) <- nrows
dlNorm <-as.data.frame(a_count)
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
  dplyr::select(SampleNames, condition1, post_Ncort) -> var_info  

#just get samples I want
row.names <- var_info$SampleNames

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code

ifelse(var_info$condition1 == "Dominant", 1,-1) -> group.dl
var_info$post_Ncort %>% 
  scale %>% 
  as.numeric -> cort.dl 



#limma modeling
d0= DGEList(dlNorm, group = group.dl)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
#  14770    12 # how many we have left

design.dl <- model.matrix(~ group.dl + cort.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]



saveRDS(v.dl, "manuscript/brain/results/limma_vdl_mPFC_25hr_DD.RDS")
# How many random sampling
R = 5000

p.dl.rand = vector('list',length = R)

for(i in 1 : R){
  print(paste("Starting on Permutation", i))

  # Randomize the traits
  cort.dl.rand = sample(cort.dl)
  group.dl.rand = sample(group.dl)

  # Model
  design.dl.rand = model.matrix(~group.dl.rand + cort.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  efit.dl.rand = eBayes(vfit.dl.rand)
  p.dl.rand[[i]] = efit.dl.rand[["p.value"]]
}

q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))

for(i in 1 : R){
  print(paste("Calculating Permutation", i))

  temp = p.dl.rand[[i]]

  for(c in 1 : 3){
    for(r in 1 : nrow(p.dl.limma)){
      if(temp[r, c] <= p.dl.limma[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}

q.dl = q.dl / R
colnames(q.dl) <- mycolnames
q.dl = as.data.frame(q.dl)
row.names(q.dl) <- rownames(dge.dl)

saveRDS(q.dl,("manuscript/brain/results/limma_vdl_cutoff5_2000_tworand_mPFC_25hr_DD.RDS"))
q.dl <- readRDS("manuscript/brain/results/limma_vdl_cutoff5_2000_tworand_mPFC_25hr_DD.RDS")


head(round(p.dl.limma,5))
head(q.dl)
hist(q.dl-p.dl.limma)

p.dl.limma %>% 
  as.data.frame() %>% 
  filter(group.dl < 0.05) %>% 
  nrow()

q.dl %>% 
  as.data.frame() %>% 
  filter(group.dl < 0.05) %>% 
  nrow()

p.dl.limma %>% 
  as.data.frame() %>% 
  filter(cort.dl < 0.05) %>% 
  nrow()

q.dl %>% 
  as.data.frame() %>% 
  distinct() %>% 
  filter(cort.dl < 0.05) %>% 
  nrow()


# ========================================================================
# replace p value with q value 
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))

colnames(efit.dl)
tmp_status <- contrasts.fit(efit.dl, coef = 2) # group
tmp_status <- eBayes(tmp_status)

tmp_cort <- contrasts.fit(efit.dl, coef = 3) # cort
tmp_cort <- eBayes(tmp_cort)

limma_list <- list()

topTable(tmp_status, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$status


topTable(tmp_cort, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cort






saveRDS(limma_list,"manuscript/brain/results/limma_mPFC_25hrDD.RDS")


my_logFC_threshold = 0.2 ## what should this be?

limma_list<- readRDS("manuscript/brain/results/limma_mPFC_25hrDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) 

#colnames(limma_list$status)
limma_list %>% 
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_status_DEG <- limma_list$status 
limma_cort_DEG <- limma_list$cort 


limma_list %>% map(~hist(.$logFC))

limma_cort_DEG %>% 
  arrange(logFC) %>% 
  head(10) 


limma_status_DEG %>% 
  arrange(logFC) %>% 
  head(10) 


#################################################

# Bring in count data for MEA
a_countdata <- read_csv("manuscript/brain/AMY_25hcounts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
row.names(a_count) <- nrows
dlNorm <-as.data.frame(a_count)
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
  dplyr::select(SampleNames, condition1, post_Ncort) -> var_info  

#just get samples I want
row.names <- var_info$SampleNames

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code

ifelse(var_info$condition1 == "Dominant", 1,-1) -> group.dl
var_info$post_Ncort %>% 
  scale %>% 
  as.numeric -> cort.dl 



#limma modeling
d0= DGEList(dlNorm, group = group.dl)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
#  14883    12 # how many we have left

design.dl <- model.matrix(~ group.dl + cort.dl)
colnames(design.dl) -> mycolnames


v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]

saveRDS(v.dl, "manuscript/brain/results/limma_vdl_MeA_25hr_DD.RDS")
# How many random sampling
R = 5000

p.dl.rand = vector('list',length = R)

for(i in 1 : R){
  print(paste("Starting on Permutation", i))
  
  # Randomize the traits
  cort.dl.rand = sample(cort.dl)
  group.dl.rand = sample(group.dl)
  
  # Model
  design.dl.rand = model.matrix(~group.dl.rand + cort.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  efit.dl.rand = eBayes(vfit.dl.rand)
  p.dl.rand[[i]] = efit.dl.rand[["p.value"]]
}

q.dl = matrix(0, nrow = nrow(p.dl.limma), ncol = ncol(p.dl.limma))

for(i in 1 : R){
  print(paste("Calculating Permutation", i))
  
  temp = p.dl.rand[[i]]
  
  for(c in 1 : 3){
    for(r in 1 : nrow(p.dl.limma)){
      if(temp[r, c] <= p.dl.limma[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}

q.dl = q.dl / R
colnames(q.dl) <- mycolnames
q.dl = as.data.frame(q.dl)
row.names(q.dl) <- rownames(dge.dl)

saveRDS(q.dl,("manuscript/brain/results/limma_vdl_cutoff5_2000_tworand_MeA_25hr_DD.RDS"))
q.dl <- readRDS("manuscript/brain/results/limma_vdl_cutoff5_2000_tworand_MeA_25hr_DD.RDS")


head(round(p.dl.limma,5))
head(q.dl)
hist(q.dl-p.dl.limma)

p.dl.limma %>% 
  as.data.frame() %>% 
  filter(group.dl < 0.05) %>% 
  nrow()

q.dl %>% 
  as.data.frame() %>% 
  filter(group.dl < 0.05) %>% 
  nrow()

p.dl.limma %>% 
  as.data.frame() %>% 
  filter(cort.dl < 0.05) %>% 
  nrow()

q.dl %>% 
  as.data.frame() %>% 
  distinct() %>% 
  filter(cort.dl < 0.05) %>% 
  nrow()


# ========================================================================
# replace p value with q value 
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))


colnames(efit.dl)
tmp_status <- contrasts.fit(efit.dl, coef = 2) # group
tmp_status <- eBayes(tmp_status)

tmp_cort <- contrasts.fit(efit.dl, coef = 3) # cort
tmp_cort <- eBayes(tmp_cort)

limma_list <- list()

topTable(tmp_status, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$status


topTable(tmp_cort, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cort






saveRDS(limma_list,"manuscript/brain/results/limma_MeA_25hrDD.RDS")


my_logFC_threshold = 0.2 ## what should this be?

limma_list<- readRDS("manuscript/brain/results/limma_MeA_25hrDD.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) 

#colnames(limma_list$status)
limma_list %>% 
  map(~summarise(.,Up = sum(logFC>0),
                 Down = sum(logFC<0))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_status_DEG <- limma_list$status 
limma_cort_DEG <- limma_list$cort 


limma_list %>% map(~hist(.$logFC))

