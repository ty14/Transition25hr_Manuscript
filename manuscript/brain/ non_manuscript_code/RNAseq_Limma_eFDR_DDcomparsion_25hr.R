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
a_countdata <- read_csv("brain/AMY_25hcounts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]

dlNorm <-as.data.frame(a_count)
row.names(dlNorm) <- nrows

colnames(dlNorm) <- substr(colnames(dlNorm),6,11)

#Group traits
#Getting metadata ready 
coldata <- read_csv("brain/sample25hr_table.csv")
head(coldata)
str(coldata)

coldata <- coldata %>% filter(region == 'P')

coldata$condition <- ifelse(is.na(coldata$condition), "same", coldata$condition)
# fix sub values
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)

coldata$SampleNames <- substr(coldata$SampleNames, 6,11)
coldata1 <- coldata %>% column_to_rownames(., var = 'SampleNames')


#check before normalizing 
dlNorm<- dlNorm[,rownames(coldata1)]
all(rownames(coldata1) == colnames(dlNorm))

#normalize and filter with all groups 
dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = coldata1$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)


d0 <- calcNormFactors(d0)
d0$samples
cutoff <- 5
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 15127    24

# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("DOM","DES")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata1  %>% filter(condition1 != "ASC")%>% filter(condition1 != "SUB") %>% 
  dplyr::select(condition1,mean_con_ng_ul) -> var_info  


dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##get groups I want


var_info$mean_con_ng_ul%>% 
  scale %>% 
  as.numeric -> cort.dl 


agg.dl <- ifelse(var_info$condition1 == "DOM", 1,-1)


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

saveRDS(v.dl, "limma_vdl_AMY_DD25.RDS")



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
saveRDS(q.dl,("limma_vdl_Q.dl_AMY_DD25.RDS"))


#########  READ IN RESULTS
q.dl <- readRDS("limma_vdl_Q.dl_AMY_DD25.RDS")


##### Analysis pulling genes out for each contrast 
tmp1 <- contrasts.fit(efit.dl, coef = 2) # cort
tmp2 <- contrasts.fit(efit.dl, coef = 3) # status
tmp3 <- contrasts.fit(efit.dl, coef = 4) # cort*status

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
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$status


topTable(tmp3, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cort_status


saveRDS(limma_list,"limma_AMY_DD25.RDS")


#quick look at number of genes
limma_list  %>% 
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


