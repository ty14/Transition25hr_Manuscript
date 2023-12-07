library(limma)
library(edgeR)
library(Mus.musculus)
organism = 'org.Mm.eg.db'
library(organism, character.only = TRUE)
library(biomaRt)
library(AnnotationDbi)
library(annotables)
grcm38 # mouse genes
library(tidyverse)



# Expression values
dlNorm <-  read.csv("brain/PFC_counts.csv", row.names = 1)
#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)

#Group traits
#Getting metadata ready 
coldata <- read_csv("manuscript/brain/results_tables/coldata_ALLxAGG.csv")
head(coldata)
str(coldata)

coldata <- coldata %>% filter(time == 70)
coldata <- coldata[c(1:25,27:29),]
coldata <- coldata %>% column_to_rownames(., var ='SampleID')

#check before normalizing 
dlNorm<- dlNorm[, rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))

#normalize and filter with all groups 
dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = coldata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)


d0 <- calcNormFactors(d0)
d0$samples
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# 14734    28

# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("ASC", "DOM","DES", "SUB")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata  %>% 
  dplyr::select(condition1,CORT, post.received1) -> var_info  


dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##get groups I want


var_info$CORT%>% 
  scale %>% 
  as.numeric -> cort.dl 


var_info$post.received1 %>% 
  scale %>% 
  as.numeric -> agg.dl 


design.dl <- model.matrix(~cort.dl*agg.dl)
colnames(design.dl) -> mycolnames

dim(dge.dl)

v.dl = voom(dge.dl, design.dl, plot = F)

# Run limma model over each gene: lmFit
vfit.dl = lmFit(v.dl, design.dl)


v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]

saveRDS(v.dl, "manuscript/brain/results/limma_vdl_mPFC_CORT_AGGREC.RDS")



#### Permutation analysis
# How many random sampling
R = 5000
set.seed(312)

#to store pvalues in
p.dl.rand = vector('list',length = R)

# to store "t" values (coefficients)
p.dl.rand.t = vector('list',length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  
  agg.dl.rand = sample(agg.dl)
  cort.dl.rand = sample(cort.dl)
  
  # Model
  design.dl.rand = model.matrix(~cort.dl.rand*agg.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  efit.dl.rand = eBayes(vfit.dl.rand)
  p.dl.rand[[g]] = efit.dl.rand[["p.value"]]
  
}


#Counting how many observations are above/below observed values

q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma)
}))

q.dl

q.dl = q.dl / R
colnames(q.dl) = mycolnames
q.dl = as.data.frame(q.dl)


## Replacing pvalues from limma with permutation p-values from pvalues?
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))



#save eFDR values 
saveRDS(q.dl,("manuscript/brain/results/limma_vdl_cutoff5_2000_tworand_mPFC_CORT_AGGREC.RDS"))


#########  READ IN RESULTS
q.dl <- readRDS("manuscript/brain/results/limma_vdl_cutoff5_2000_tworand_mPFC_CORT_AGGREC.RDS")


##### Analysis pulling genes out for each contrast 
tmp1 <- contrasts.fit(efit.dl, coef = 2) # cort
tmp2 <- contrasts.fit(efit.dl, coef = 3) # agg
tmp3 <- contrasts.fit(efit.dl, coef = 4) # cort*agg

limma_list <- list()

topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cort

topTable(tmp2, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$aggrec


topTable(tmp3, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cort_aggrec


saveRDS(limma_list,"manuscript/brain/results/limma_mPFC_CORT_AGGREC.RDS")


#quick look at number of genes
limma_list <- readRDS("manuscript/brain/results/limma_mPFC_CORT_AGGREC.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

limma_list %>% map(~filter(., P.Value<0.05)) %>% 
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_list %>% map(~hist(.$logFC))



###### now agg given 

# Expression values
dlNorm <-  read.csv("brain/PFC_counts.csv", row.names = 1)
#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)

#Group traits
#Getting metadata ready 
coldata <- read_csv("manuscript/brain/results_tables/coldata_ALLxAGG.csv")
head(coldata)
str(coldata)

coldata <- coldata %>% filter(time == 70)
coldata <- coldata[c(1:25,27:29),]
coldata <- coldata %>% column_to_rownames(., var ='SampleID')

#check before normalizing 
dlNorm<- dlNorm[, rownames(coldata)]
all(rownames(coldata) == colnames(dlNorm))

#normalize and filter with all groups 
dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = coldata$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)


d0 <- calcNormFactors(d0)
d0$samples
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# 14734    28

# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("ASC", "DOM","DES", "SUB")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata  %>% 
  dplyr::select(condition1,CORT, post.given1) -> var_info  


dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##get groups I want


var_info$CORT%>% 
  scale %>% 
  as.numeric -> cort.dl 


var_info$post.given1 %>% 
  scale %>% 
  as.numeric -> agg.dl 


design.dl <- model.matrix(~cort.dl*agg.dl)
colnames(design.dl) -> mycolnames

dim(dge.dl)

v.dl = voom(dge.dl, design.dl, plot = F)

# Run limma model over each gene: lmFit
vfit.dl = lmFit(v.dl, design.dl)


v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)
efit.dl = eBayes(vfit.dl)
p.dl.limma = efit.dl[["p.value"]]

saveRDS(v.dl, "manuscript/brain/results/limma_vdl_mPFC_CORT_AGGiven.RDS")



#### Permutation analysis
# How many random sampling
R = 5000
set.seed(312)

#to store pvalues in
p.dl.rand = vector('list',length = R)

# to store "t" values (coefficients)
p.dl.rand.t = vector('list',length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  
  agg.dl.rand = sample(agg.dl)
  cort.dl.rand = sample(cort.dl)
  
  # Model
  design.dl.rand = model.matrix(~cort.dl.rand*agg.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  efit.dl.rand = eBayes(vfit.dl.rand)
  p.dl.rand[[g]] = efit.dl.rand[["p.value"]]
  
}


#Counting how many observations are above/below observed values

q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma)
}))

q.dl

q.dl = q.dl / R
colnames(q.dl) = mycolnames
q.dl = as.data.frame(q.dl)


## Replacing pvalues from limma with permutation p-values from pvalues?
efit.dl[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl$coefficients)))



#save eFDR values 
saveRDS(q.dl,("manuscript/brain/results/limma_vdl_cutoff5_2000_tworand_mPFC_CORT_AGGiven.RDS"))


#########  READ IN RESULTS
q.dl <- readRDS("manuscript/brain/results/limma_vdl_cutoff5_2000_tworand_mPFC_CORT_AGGiven.RDS")


##### Analysis pulling genes out for each contrast 
tmp1 <- contrasts.fit(efit.dl, coef = 2) # cort
tmp2 <- contrasts.fit(efit.dl, coef = 3) # agg
tmp3 <- contrasts.fit(efit.dl, coef = 4) # cort*agg

limma_list <- list()

topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cort

topTable(tmp2, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$agiven


topTable(tmp3, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cort_agiven


saveRDS(limma_list,"manuscript/brain/results/limma_mPFC_CORT_AGGgiven.RDS")


#quick look at number of genes
limma_list <- readRDS("manuscript/brain/results/limma_mPFC_CORT_AGGgiven.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= 0.2)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

limma_list %>% map(~filter(., P.Value<0.05)) %>% 
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_list %>% map(~hist(.$logFC))



x <- limma_list$cort_agiven
x %>% arrange(logFC) %>%  head(.,10)
