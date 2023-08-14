library(limma)
library(Mus.musculus)
library(tidyverse)
library(scales)


library(tidyverse)
# library(edgeR)
library(DESeq2)
library(pheatmap)
library(annotables)
grcm38 # mouse genes



rawcount_list <- list()
temp = list.files(path = "manuscript/brain",pattern="*_25hcounts.csv")
temp0 <- paste0("manuscript/brain/",temp)
rawcount_list = lapply(temp0, read_csv)
regions <- gsub("_25hcounts.csv","",temp)

# behavior data cleaned from rna_seq github repo 
names(rawcount_list) <- regions
lapply(rawcount_list, head)

rawcount_list<- rawcount_list %>% map(~mutate(., ensgene = X)) %>% map(~dplyr::select(., -X))

counts <- rawcount_list %>%
  reduce(full_join, by = 'ensgene')

# total gene counts per sample ========================================

counts %>% 
  summarize_if(is.numeric,sum,na.rm = T) %>% 
  t() %>% 
  as.data.frame %>% 
  rename(genecounts = V1) %>% 
  rownames_to_column(var = 'sampleID') -> genecounts


genecounts %>% 
  ggplot(aes(genecounts)) +
  geom_histogram(bins = 40,alpha =0.5,color = 'grey') +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_minimal() -> p_genecounts

print(p_genecounts)
 # ggsave(p_genecounts,filename = 'results_figures/25hp_genecounts.png',width = 8, height = 5)
 
#Getting metadata ready 
coldata <- read_csv("manuscript/brain/sample25hr_table.csv")
head(coldata)
str(coldata)

##getting condition1
# Dom-Dom to Descenders (DOM to SUB)(4->1)
# Sub-Sub to Ascenders (Sub to DOM)  (1->4)
coldata$condition <- ifelse(coldata$Postrank == 4 & coldata$Prerank == 4, "same", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "Dominant", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "Subordinate", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "Descenders", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "Ascenders", coldata$condition1)
##getting condition2
coldata$condition2 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "HeldRank", coldata$condition)
coldata$condition2 <- ifelse(coldata$condition == "same" & coldata$Prerank == 4, "HeldRank", coldata$condition2)
coldata$condition2 <- ifelse(coldata$condition == "descenders"| coldata$condition == "ascenders", "ChangedRank", coldata$condition2)


# factor
coldata$condition1 <- factor(coldata$condition1, c("Ascenders", "Descenders", "Subordinate", "Dominant"))
coldata$condition2 <- factor(coldata$condition2, c("ChangedRank", "HeldRank"))
# First AMY data 
amy_data<- coldata %>% filter(region == "A")%>% arrange(condition) %>% arrange(group) 

row.names <- amy_data$SampleNames

amy_data1 <- amy_data%>% dplyr::select(-SampleNames)


row.names(amy_data1) <- row.names #Assigning row names from as sample names  
head(amy_data1)


# Bring in count data for mAMY
a_countdata <- read_csv("manuscript/brain/AMY_25hcounts.csv")
a_countdata$X -> nrows
a_count <- a_countdata[,-1]
row.names(a_count) <- nrows
a_count <-as.data.frame(a_count)
a_count <- a_count[rowSums(a_count >= 10) > round((length(amy_data1))*0.9), ]
#checks
all(row.names(amy_data1) %in% colnames(a_count)) #check 

a_count <- a_count[,rownames(amy_data1)]
all(rownames(amy_data1) == colnames(a_count)) #check



# Second mpfc data data 
pfc_data<- coldata %>% filter(region == "P")%>% arrange(condition) %>% arrange(group) 

row.names <- pfc_data$SampleNames

pfc_data1 <- pfc_data%>% dplyr::select(-SampleNames)

row.names(pfc_data1) <- row.names #Assigning row names from as sample names  
head(pfc_data1)


# Bring in count data for mpfc
p_countdata <- read_csv("manuscript/brain/PFC_25hcounts.csv")
p_countdata$X -> nrows
p_count <- p_countdata[,-1]
row.names(p_count) <- nrows
p_count <-as.data.frame(p_count)

p_count <- p_count[rowSums(p_count >= 10) > round((length(pfc_data1))*0.9), ]
all (row.names(pfc_data1) %in% colnames(p_count)) #check 

p_count <- p_count[, rownames(pfc_data1)]
all(rownames(pfc_data1) == colnames(p_count)) #check


#MeA
dag <- DESeqDataSetFromMatrix(countData = a_count,
                              colData = amy_data1,
                              design = ~condition1)



nrow(dag) # 13097 lost about half

#mPFC
dpg <- DESeqDataSetFromMatrix(countData = p_count,
                              colData = pfc_data1,
                              design = ~condition1)

nrow(dpg) # 12870

#dds
da_dom <- DESeq(dag)
dp_dom <- DESeq(dpg)

DESeq_domlist<- list(MeA=da_dom, mPFC=dp_dom)

# adjusted pval and fold change threshold to 'define' DEG
adj_pval_cutoff = 0.1
threshold = 0.9  # 90% of the subject
LFC_threshold = 0.2
pvalue_threshold = 0.05 # didn't use it for now



get_LFC_results <- function(dds)
{ 
  res1 <- results(dds, contrast=c("condition1","Dominant","Descenders")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "Dominant - Descenders") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(x=., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  res2 <- results(dds, contrast=c("condition1","Dominant","Descenders")) %>% 
    as.data.frame() %>% 
    rownames_to_column('ensgene') %>% 
    as_tibble() %>% 
    dplyr::select(ensgene, log2FoldChange, pvalue, padj) %>% 
    mutate(contrast = "Dominant - Descenders") %>% 
    filter(abs(log2FoldChange)>=LFC_threshold) %>%
    filter(padj <=adj_pval_cutoff ) %>% filter(pvalue <=pvalue_threshold) %>%
    arrange(desc(abs(log2FoldChange))) %>% 
    left_join(., y =grcm38[,c("ensgene", "symbol", "description")],by = "ensgene") %>% 
    unique()
  
  
  rbind(res1, res2) %>% 
    arrange(desc(abs(log2FoldChange))) -> result
  
  return(result)
}

results_domlist <- DESeq_domlist %>% map(get_LFC_results)
results_dom <- results_domlist %>% map2_df(.,names(.), ~mutate(.x,region = .y))%>%  arrange(desc(log2FoldChange)) %>% unique(.)
head(results_dom)
tail(results_dom)
# write.csv(results_dom, 'manuscript/brain/results/DD25_DEG_table.csv' ,row.names = F)
