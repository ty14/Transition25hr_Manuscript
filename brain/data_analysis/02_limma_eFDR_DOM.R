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
dlNorm <-  read.csv("brain/AMY_counts.csv", row.names = 1)
#remove zeros
dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)

#Group traits
#Getting metadata ready 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#fixing things
coldata$region <- gsub("mPF", "PFC", coldata$region)
coldata$groupEX <- coldata$group


# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)


#just get samples I want
coldata$SampleID <- substr(coldata$SampleName, 7, 13)

coldata <- coldata %>%  dplyr::select(-SampleName) %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != 'ascenders') %>% 
  pivot_wider(
    names_from = region, 
    values_from = region)

row.names <- coldata$SampleID
row.names(coldata) <- row.names #Assigning row names from as sample names  
head(coldata)



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
#15184    40

# Now take out groups that you want
#DOMs first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CDOM", "DOM", "DES")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

coldata %>% 
  filter(condition1 != "ASC") %>%
  filter(condition1 != "SUB") %>%
  filter(condition1 != "CSUB") %>%
  dplyr::select(SampleID, condition1) -> var_info  

row.names <- var_info$SampleID

row.names(var_info) <- row.names #Assigning row names from as sample names  
head(var_info )

dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>%
  factor(.,levels = c("CDOM","DOM","DES")) -> group.dl


design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)

# Run limma model over each gene: lmFit
vfit.dl = lmFit(v.dl, design.dl)

contr.matrix <- makeContrasts(group.dlCDOM-group.dlDOM,
                                 group.dlCDOM-group.dlDES, 
                                 group.dlDOM-group.dlDES,
                                 levels=design.dl)

vfit.dl2 <- contrasts.fit(vfit.dl, contr.matrix)

efit.dl2 = eBayes(vfit.dl2)

p.dl.limma2 = efit.dl2[["p.value"]]
p.dl.limma2.t = efit.dl2[["t"]]
head(p.dl.limma2)

# saveRDS(v.dl, "brain/results/RDS/limma_vdl_MEA_CDOM.RDS")



#### Permutation analysis

## Observed Data from Above
p.dl.limma2   #pvalues
p.dl.limma2.t  #t-values


# How many random sampling
R = 5000
# R = 5

set.seed(312)

#to store pvalues in
p.dl.rand = vector('list',length = R)

# to store "t" values (coefficients)
p.dl.rand.t = vector('list',length = R)

for( g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  
  group.dl.rand = sample(group.dl)
  
  # Model
  design.dl.rand = model.matrix(~0 + group.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contr.matrix)
  
  efit.dl.rand2 = eBayes(vfit.dl.rand2)
  
  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  p.dl.rand.t[[g]] = efit.dl.rand2[["t"]]

}



####
#Counting how many observations are above/below observed values

# # create empty matrix  
# q.dl = matrix(0, nrow = nrow(p.dl.limma2), ncol = ncol(p.dl.limma2))
# 
# 
# for(h in 1 : R){
#   print(paste("Calculating Permutation", h))
#   
#   temp = p.dl.rand[[h]]
#   
#   for(c in 1 : 3){
#     for(r in 1 : nrow(p.dl.limma2)){
#       if(temp[r, c] <= p.dl.limma2[r, c]){
#         q.dl[r, c] = q.dl[r, c] + 1
#       }
#     }
#   }
# }

## Faster Way than above
q.dl <- Reduce(`+`, lapply(p.dl.rand, \(x) {
  (x < p.dl.limma2)
}))


q.dl

q.dl = q.dl / R
colnames(q.dl) = mycolnames
q.dl = as.data.frame(q.dl)



#### Counting permutations more positive than t or more negative than t
xx<-Reduce(`+`, lapply(p.dl.rand.t, \(x) {
  (abs(x) > abs(p.dl.limma2.t) & sign(x) == sign(p.dl.limma2.t))
}))

xx = xx / R
colnames(xx) = mycolnames
xx = as.data.frame(xx)



## Duplicate object just for below.
efit.dl2x <- efit.dl2

## Replacing pvalues from limma with permutation p-values from t-coefficients?
efit.dl2x[["p.value"]] <- xx
row.names(xx) <- NULL
sum(duplicated(row.names(efit.dl2x$coefficients)))
saveRDS(xx,("brain/results/RDS/jamesxx.RDS"))


## Replacing pvalues from limma with permutation p-values from pvalues?
efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))
saveRDS(q.dl,("brain/results/RDS/jamesqdl.RDS"))

# saveRDS(q.dl,("brain/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_CDOM.RDS"))

xx

#########  READ IN RESULTS

q.dl <- readRDS("brain/results/RDS/limma_vdl_cutoff5_2000_tworand_MEA_CDOM.RDS")
efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))



##### Analysis

tmp1 <- contrasts.fit(efit.dl2, coef = 1) # 
tmp2 <- contrasts.fit(efit.dl2, coef = 2) # 
tmp3 <- contrasts.fit(efit.dl2, coef = 3) # 

limma_list <- list()

topTable(tmp1, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cdom

topTable(tmp2, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cdes


topTable(tmp3, sort.by = "P", n = Inf) %>% 
  rownames_to_column('ensgene') %>% 
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$domdes


saveRDS(limma_list,"brain/results/RDS/limma_MEA_CDOM.RDS")

        limma_list <- readRDS("brain/results/RDS/limma_MEA_CDOM.RDS") %>% 
          map(~distinct(.)) %>% 
          map(~filter(.,abs(logFC) >= 0.2)) %>%
          map(~filter(.,P.Value <0.05)) %>%
          map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
          map(~filter(.,!is.na(entrez))) 
        
        
        
limma_list <- 
        limma_list %>% 
          map(~distinct(.)) %>% 
          map(~filter(.,abs(logFC) >= 0.2)) %>%
          map(~filter(.,P.Value <0.05)) %>%
          map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
          map(~filter(.,!is.na(entrez))) 
        
        
        #quick look at number of genes
        limma_list %>% map(~filter(., P.Value<0.05)) %>% 
          map(~summarise(.,Up = sum(logFC>0.2),
                         Down = sum(logFC<0.2))) %>% 
          map(~mutate(.,Total = Up + Down))
        
        limma_list %>% map(~hist(.$logFC))
        
        
        cdom <- limma_list$cdom
        cdes <- limma_list$cdes 
        domdes <- limma_list$domdes 
        
        domdes %>% arrange(logFC) %>% head(.,10)
        
        

#### redoing with coefficients pvalues
        
        tmp1 <- contrasts.fit(efit.dl2x, coef = 1) # 
        tmp2 <- contrasts.fit(efit.dl2x, coef = 2) # 
        tmp3 <- contrasts.fit(efit.dl2x, coef = 3) # 
        
        limma_list <- list()
        
        topTable(tmp1, sort.by = "P", n = Inf) %>% 
          rownames_to_column('ensgene') %>% 
          left_join(grcm38) %>%
          filter(!is.na(symbol)) %>% 
          dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cdom
        
        topTable(tmp2, sort.by = "P", n = Inf) %>% 
          rownames_to_column('ensgene') %>% 
          left_join(grcm38) %>%
          filter(!is.na(symbol)) %>% 
          dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cdes
        
        
        topTable(tmp3, sort.by = "P", n = Inf) %>% 
          rownames_to_column('ensgene') %>% 
          left_join(grcm38) %>%
          filter(!is.na(symbol)) %>% 
          dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$domdes
        
        
        #saveRDS(limma_list,"brain/results/RDS/limma_MEA_CDOM.RDS")
        
        # limma_list <- readRDS("brain/results/RDS/limma_MEA_CDOM.RDS") %>% 
        #   map(~distinct(.)) %>% 
        #   map(~filter(.,abs(logFC) >= 0.2)) %>%
        #   map(~filter(.,P.Value <0.05)) %>%
        #   map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
        #   map(~filter(.,!is.na(entrez))) 
        # 
        
        
        limma_list <- 
          limma_list %>% 
          map(~distinct(.)) %>% 
          map(~filter(.,abs(logFC) >= 0.2)) %>%
          map(~filter(.,P.Value <0.05)) %>%
          map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
          map(~filter(.,!is.na(entrez))) 
        
        
        #quick look at number of genes
        limma_list %>% map(~filter(., P.Value<0.05)) %>% 
          map(~summarise(.,Up = sum(logFC>0.2),
                         Down = sum(logFC<0.2))) %>% 
          map(~mutate(.,Total = Up + Down))
        
        limma_list %>% map(~hist(.$logFC))
        
        
        cdom <- limma_list$cdom
        cdes <- limma_list$cdes 
        domdes <- limma_list$domdes 
        
        domdes %>% arrange(logFC) %>% head(.,10)
        
        
        
        
        
        
        
        
         ##########################################################################################       
        
       # Expression values
        dlNorm <-  read.csv("brain/PFC_counts.csv", row.names = 1)
        #remove zeros
        dlNorm <- dlNorm[apply(dlNorm[], 1, function(x) !all(x==0)),]
        #trim sample ids
        colnames(dlNorm)[c(1:67)] <- substr(colnames(dlNorm)[c(1:67)], 7, 13)
        
        #Group traits
        #Getting metadata ready 
        coldata <- read_csv("brain/sample70min_table.csv")
        head(coldata)
        str(coldata)
        
        coldata <- coldata %>% 
          filter(post_idbatch != "3-4Batch12") 
        
        #fixing things
        coldata$region <- gsub("mPF", "PFC", coldata$region)
        coldata$groupEX <- coldata$group
        
        
        # Normalizing cort data
        # df <- transform(df, N = (N - min(N)) / (max(N) - min(N))
        
        coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
        coldata <- coldata %>%  dplyr::select(-group, -period)
        
        #getting condition1
        table(coldata$condition)
        # CDOM, RDOM to Descenders (DOM to SUB)(4->1)
        # CSUB, SUB to Ascenders (Sub to DOM)  (1->4)
        
        coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
        coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
        coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
        coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)
        coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 4, "CSUB", coldata$condition1)
        coldata$condition1 <- ifelse(coldata$condition == "control" & coldata$Postrank == 1, "CDOM", coldata$condition1)
        
        
        #just get samples I want
        coldata$SampleID <- substr(coldata$SampleName, 7, 13)
        
        coldata <- coldata %>%  dplyr::select(-SampleName) %>% 
          filter(Postrank != 3) %>% 
          filter(condition1 != 'ascenders') %>% 
          pivot_wider(
            names_from = region, 
            values_from = region)
        
        row.names <- coldata$SampleID
        row.names(coldata) <- row.names #Assigning row names from as sample names  
        head(coldata)
        
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
        
        cutoff <- 5
        drop <- which(apply(cpm(d0), 1, max) < cutoff)
        dge.dl <- d0[-drop,]
        dim(dge.dl)
        # 14991    39
        
        # Now take out groups that you want
        #DOMs first 
        dge.dl$samples$group
        
        dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("CDOM", "DOM", "DES")]
        dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
        dge.dl_dom$samples$group
        dge.dl<- dge.dl_dom
        dge.dl$samples$group
        
        coldata %>% 
          filter(condition1 != "ASC") %>%
          filter(condition1 != "SUB") %>%
          filter(condition1 != "CSUB") %>%
          dplyr::select(SampleID, condition1) -> var_info  
        
        row.names <- var_info$SampleID
        
        row.names(var_info) <- row.names #Assigning row names from as sample names  
        head(var_info )
        
        dlNorm<- dlNorm[, rownames(var_info)]
        all(rownames(var_info) == colnames(dlNorm)) #check
        
        ##following Won's code
        var_info$condition1 %>%
          factor(.,levels = c("CDOM","DOM","DES")) -> group.dl
        
        
        design.dl <- model.matrix(~ 0 + group.dl)
        colnames(design.dl) -> mycolnames
        
        v.dl = voom(dge.dl, design.dl, plot = F)
        vfit.dl = lmFit(v.dl, design.dl)
        
        contr.matrix <- makeContrasts(group.dlCDOM-group.dlDOM,
                                      group.dlCDOM-group.dlDES, 
                                      group.dlDOM-group.dlDES,
                                      levels=design.dl)
        
        vfit.dl2 <- contrasts.fit(vfit.dl, contr.matrix)
        
        efit.dl2 = eBayes(vfit.dl2)
        
        p.dl.limma2 = efit.dl2[["p.value"]]
        head(p.dl.limma2)
        
        saveRDS(v.dl, "brain/results/RDS/limma_vdl_mPFC_CDOM.RDS")
        
      
        
        # How many random sampling
        R = 5000
        
        set.seed(312)
        
        p.dl.rand = vector('list',length = R)
        
        for( g in 1 : R){
          print(paste("Starting on Permutation", g))
          
          # Randomize the traits
          
          group.dl.rand = sample(group.dl)
          
          # Model
          design.dl.rand = model.matrix(~0 + group.dl.rand)
          colnames(design.dl.rand) <- mycolnames
          
          # Calculate p-values based on randomized traits
          v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
          vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
          
          vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contr.matrix)
          
          efit.dl.rand2 = eBayes(vfit.dl.rand2)
          
          p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
          head(p.dl.rand[[g]])
        }
        q.dl = matrix(0, nrow = nrow(p.dl.limma2), ncol = ncol(p.dl.limma2))
        
        for(h in 1 : R){
          print(paste("Calculating Permutation", h))
          
          temp = p.dl.rand[[h]]
          
          for(c in 1 : 3){
            for(r in 1 : nrow(p.dl.limma2)){
              if(temp[r, c] <= p.dl.limma2[r, c]){
                q.dl[r, c] = q.dl[r, c] + 1
              }
            }
          }
        }
        
        
        
        q.dl = q.dl / R
        colnames(q.dl) = mycolnames
        q.dl = as.data.frame(q.dl)
        
        efit.dl2[["p.value"]] <- q.dl
        row.names(q.dl) <- NULL
        sum(duplicated(row.names(efit.dl2$coefficients)))
        
        saveRDS(q.dl,("brain/results/RDS/limma_vdl_cutoff5_2000_tworand_mPFC_CDOM.RDS"))
        
        efit.dl2[["p.value"]] <- q.dl
        row.names(q.dl) <- NULL
        sum(duplicated(row.names(efit.dl2$coefficients)))
        
        tmp1 <- contrasts.fit(efit.dl2, coef = 1) # 
        tmp2 <- contrasts.fit(efit.dl2, coef = 2) # 
        tmp3 <- contrasts.fit(efit.dl2, coef = 3) # 
        
        limma_list <- list()
        
        topTable(tmp1, sort.by = "P", n = Inf) %>% 
          rownames_to_column('ensgene') %>% 
          left_join(grcm38) %>%
          filter(!is.na(symbol)) %>% 
          dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$cdom
        
        topTable(tmp2, sort.by = "P", n = Inf) %>% 
          rownames_to_column('ensgene') %>% 
          left_join(grcm38) %>%
          filter(!is.na(symbol)) %>% 
          dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$cdes
        
        
        topTable(tmp3, sort.by = "P", n = Inf) %>% 
          rownames_to_column('ensgene') %>% 
          left_join(grcm38) %>%
          filter(!is.na(symbol)) %>% 
          dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$domdes
        
        
        saveRDS(limma_list,"brain/results/RDS/limma_mPFC_CDOM.RDS")
        
        limma_list <- readRDS("brain/results/RDS/limma_mPFC_CDOM.RDS") %>% 
          map(~distinct(.)) %>% 
          map(~filter(.,abs(logFC) >= 0.2)) %>%
          map(~filter(.,P.Value <0.05)) %>%
          map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
          map(~filter(.,!is.na(entrez))) 
        
        #quick look at number of genes
        limma_list %>% map(~filter(., P.Value<0.05)) %>% 
          map(~summarise(.,Up = sum(logFC>0.2),
                         Down = sum(logFC<0.2))) %>% 
          map(~mutate(.,Total = Up + Down))
        
        limma_list %>% map(~hist(.$logFC))
        
        
        cdom <- limma_list$cdom
        cdes <- limma_list$cdes 
        domdes <- limma_list$domdes
        
        
        cdom %>% arrange(-logFC) %>% head(.,10)
        