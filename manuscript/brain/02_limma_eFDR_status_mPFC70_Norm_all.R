# libraries 
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


#Getting metadata ready 
coldata <- read_csv("brain/sample25hr_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
# coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(post_Ncort)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)

coldata$condition1  <- coldata$condition1  %>% replace_na('SUB')
coldata$condition1  <- paste(coldata$condition1, coldata$time)

coldata <- coldata %>% 
  filter(region == "P")

coldata$Sampleid <- substr(coldata$SampleNames,6,12)

coldata25 <- coldata %>% select(condition1, Sampleid) 

#remove outlier
#b1.2.1 outlier DES
pfc_data25 <- coldata25 %>% filter(Sampleid != 'b1.2.1.')

#now 70 min data 
coldata <- read_csv("brain/sample70min_table.csv")
head(coldata)
str(coldata)

#group thing 
coldata$groupEX <- coldata$group
# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
# coldata$post_Ncort[is.na(coldata$post_Ncort)] = mean(post_Ncort)

#getting condition1
table(coldata$condition)
# CDOM, RDOM to Descenders (DOM to SUB)(4->1)
# CSUB, SUB to Ascenders (Sub to DOM)  (1->4)

coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "DOM", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "SUB", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "DES", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "ASC", coldata$condition1)

coldata$condition1  <- coldata$condition1  %>% replace_na('SUB')
coldata$condition1 <- paste(coldata$condition1, coldata$time)
pfc_data1<- coldata %>% 
  filter(region != "AMY") %>% 
  filter(Postrank != 3) %>% 
  filter(condition1 != "control 1 hr") %>% 
  filter(condition1 != "ascenders 1 hr")

table(pfc_data1$condition1)
pfc_data1$Sampleid <- substr(pfc_data1$SampleName, 7,13)
pfc_data70 <- pfc_data1 %>% select(condition1, Sampleid)


pfc_data <- pfc_data70 %>% rbind(pfc_data25) 
pfc_data <- pfc_data %>% column_to_rownames(., var = "Sampleid")

### 70 min data 
# Expression values
dlNorm70 <-  read.csv("brain/PFC_counts.csv", row.names = 1)
#remove zeros
dlNorm70 <- dlNorm70[apply(dlNorm70[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm70)[c(1:67)] <- substr(colnames(dlNorm70)[c(1:67)], 7, 13)


# Bring in count data for mPFC 25 hr 
dlNorm25 <-  read.csv("brain/PFC_25hcounts.csv", row.names = 1)
#remove zeros
dlNorm25 <- dlNorm25[apply(dlNorm25[], 1, function(x) !all(x==0)),]
#trim sample ids
colnames(dlNorm25)[c(1:24)] <- substr(colnames(dlNorm25)[c(1:24)], 6, 12)
#remove outlier B1.1.2.
dlNorm25 <- dlNorm25[,c(1:19,21:24)]


#combinf counts to normalize together 
dl70 <- dlNorm70 %>% rownames_to_column(., var = "gene")
dl25 <- dlNorm25 %>% rownames_to_column(., var = "gene")

dlNorm <- dl70 %>% full_join(dl25)
dlNorm <- dlNorm %>% column_to_rownames(., var= "gene")


colnames(dlNorm)
rownames(pfc_data)

#check before normalizing 
dlNorm<- dlNorm[ ,rownames(pfc_data)]
all(rownames(pfc_data) == colnames(dlNorm))

#normalize and filter with all groups 
# row.names(dlNorm) <- nrows
dlNorm <- dlNorm[!is.na(rowSums(dlNorm)),]

d = apply(dlNorm, 2, as.numeric)
dim(d)

d0= DGEList(d, group = pfc_data$condition1)
dim(d0)
rownames(d0) <- rownames(dlNorm)
d0 <- calcNormFactors(d0)

cutoff <- 10
drop <- which(apply(cpm(d0), 1, max) < cutoff)
dge.dl <- d0[-drop,]
dim(dge.dl)
# [1] 13367    51
# Now take out groups that you want
#70min first 
dge.dl$samples$group

dge.dl_dom <- dge.dl[, dge.dl$samples$group %in% c("ASC 1 hr", "DES 1 hr", "DOM 1 hr", "SUB 1 hr")]
dge.dl_dom$samples$group <- droplevels(dge.dl_dom$samples$group)
dge.dl_dom$samples$group
dge.dl<- dge.dl_dom
dge.dl$samples$group

pfc_data %>% 
  filter( condition1 %in% c("ASC 1 hr", "DES 1 hr", "DOM 1 hr", "SUB 1 hr") ) -> var_info  



dlNorm<- dlNorm[, rownames(var_info)]
all(rownames(var_info) == colnames(dlNorm)) #check

##following Won's code
var_info$condition1 %>% 
  factor(.,levels = c("DOM 1 hr","ASC 1 hr","DES 1 hr", "SUB 1 hr")) -> group.dl

group.dl<- gsub(" 1 hr", "", group.dl)

design.dl <- model.matrix(~ 0 + group.dl)
colnames(design.dl) -> mycolnames

v.dl = voom(dge.dl, design.dl, plot = F)
vfit.dl = lmFit(v.dl, design.dl)

contrast.matrix <- makeContrasts(group.dlASC-group.dlDOM,
                                 group.dlDOM-group.dlSUB,
                                 group.dlDES-group.dlDOM, 
                                 group.dlDES-group.dlASC,
                                 group.dlASC-group.dlSUB,
                                 group.dlDES-group.dlSUB,
                                 levels=design.dl)

vfit.dl2 <- contrasts.fit(vfit.dl, contrast.matrix)

efit.dl2 = eBayes(vfit.dl2)

p.dl.limma2 = efit.dl2[["p.value"]]
head(p.dl.limma2)

saveRDS(v.dl, "manuscript/brain/results_use/limma_vdl_PFC70min_NormRG.RDS")

# How many random sampling
R = 5000

p.dl.rand = vector('list',length = R)

for(g in 1 : R){
  print(paste("Starting on Permutation", g))
  
  # Randomize the traits
  
  group.dl.rand = sample(group.dl)
  
  # Model
  design.dl.rand = model.matrix(~0 + group.dl.rand)
  colnames(design.dl.rand) <- mycolnames
  
  # Calculate p-values based on randomized traits
  v.dl.rand = voom(dge.dl, design.dl.rand, plot = F)
  vfit.dl.rand = lmFit(v.dl.rand, design.dl.rand)
  
  vfit.dl.rand2 <- contrasts.fit(vfit.dl.rand, contrast.matrix)
  
  efit.dl.rand2 = eBayes(vfit.dl.rand2)
  
  p.dl.rand[[g]] = efit.dl.rand2[["p.value"]]
  head(p.dl.rand[[g]])
}

q.dl = matrix(0, nrow = nrow(p.dl.limma2), ncol = ncol(p.dl.limma2))

for(h in 1 : R){
  print(paste("Calculating Permutation", h))
  
  temp = p.dl.rand[[h]]
  
  for(c in 1 : 6){
    for(r in 1 : nrow(p.dl.limma2)){
      if(temp[r, c] <= p.dl.limma2[r, c]){
        q.dl[r, c] = q.dl[r, c] + 1
      }
    }
  }
}

q.dl = q.dl / R
q.dl = as.data.frame(q.dl)
row.names(q.dl) <- rownames(dge.dl)
colnames(q.dl) <- mycolnames

saveRDS(q.dl,("manuscript/brain/results_use/limma_vdl_cutoff5_2000_tworand_PFC70min_ReorganizedGroups.RDS"))

q.dl <- readRDS("manuscript/brain/results_use/limma_vdl_cutoff5_2000_tworand_PFC70min_ReorganizedGroups.RDS")


efit.dl2[["p.value"]] <- q.dl
row.names(q.dl) <- NULL
sum(duplicated(row.names(efit.dl2$coefficients)))

tmp1 <- contrasts.fit(efit.dl2, coef = 1) # DOM_ASC

tmp2 <- contrasts.fit(efit.dl2, coef = 2) #DOM-SUB

tmp3 <- contrasts.fit(efit.dl2, coef = 3) #DES-DOM

tmp4 <- contrasts.fit(efit.dl2, coef = 4) #DES-ASC

tmp5 <- contrasts.fit(efit.dl2, coef = 5) #ASC-SUB

tmp6 <- contrasts.fit(efit.dl2, coef = 6) #DES-SUB

limma_list <- list()


topTable(tmp1, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>% 
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) ->limma_list$ascdom



topTable(tmp2, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val)  -> limma_list$domsub



topTable(tmp3, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$desdom


topTable(tmp4, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$desasc


topTable(tmp5, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$ascsub


topTable(tmp6, sort.by = "P", n = Inf) %>%
  rownames_to_column('ensgene') %>%
  left_join(grcm38) %>%
  filter(!is.na(symbol)) %>%
  dplyr::select(symbol,logFC,P.Value,adj.P.Val) -> limma_list$dessub


saveRDS(limma_list,"manuscript/brain/results_use/limma_PFC70min_Norm_RG.RDS")


#quick look at number of genes
limma_list %>% map(~filter(., P.Value<0.05)) %>% 
  map(~summarise(.,Up = sum(logFC>0.2),
                 Down = sum(logFC<0.2))) %>% 
  map(~mutate(.,Total = Up + Down))

limma_list %>% map(~hist(.$logFC))

