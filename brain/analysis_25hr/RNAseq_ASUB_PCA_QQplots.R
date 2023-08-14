library(limma)
library(Mus.musculus)
library(DESeq2)
library(tidyverse)
global_size = 10


coldata <- read_csv("manuscript/brain/sample25hr_table.csv")
head(coldata)
str(coldata)

coldata <- coldata %>% filter(condition != "control")
table(coldata$condition)

# Normalizing cort data
# df <- transform(df, N = (N - min(N)) / (max(N) - min(N))

coldata <-coldata %>% mutate(post_Ncort =  (mean_con_ng_ul - min(mean_con_ng_ul,na.rm=T))/(max(mean_con_ng_ul,na.rm=T)-min(mean_con_ng_ul,na.rm=T)))
coldata <- coldata %>%  dplyr::select(-group, -period)
coldata$post_Ncort[is.na(coldata$post_Ncort)] = 0


coldata$condition <- ifelse(coldata$Postrank == 4 & coldata$Prerank == 4, "same", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition == "same" & coldata$Prerank == 1, "Dominant", coldata$condition)
coldata$condition1 <- ifelse(coldata$condition1 == "same" & coldata$Prerank == 4, "Subordinate", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "descenders" & coldata$Postrank == 4 & coldata$Prerank == 1, "Descender", coldata$condition1)
coldata$condition1 <- ifelse(coldata$condition == "ascenders" & coldata$Prerank == 4 & coldata$Postrank == 1, "Ascender", coldata$condition1)
table(coldata$condition1)


coldata <- coldata %>% 
  filter(condition1 != "Descender") %>%
  filter(condition1 != "Dominant") %>%
  dplyr::select(SampleNames, region, condition1,AggGiven70min, AggRec70min, post_Ncort)


# First AMY data 
amy_data<- coldata %>% filter(region == "A")%>% arrange(condition1)

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

#filter
a_count <- a_count[rowSums(a_count >= 20) > round((length(amy_data1))*0.9), ]

#checks
all (row.names(amy_data1) %in% colnames(a_count)) #check 

a_count <- a_count[,rownames(amy_data1)]
all(rownames(amy_data1) == colnames(a_count)) #check



# Second mpfc data data 
pfc_data<- coldata %>% filter(region != "A")%>% arrange(condition1) 

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

#filter
p_count <- p_count[rowSums(p_count >= 20) > round((length(pfc_data1))*0.9), ]
# pkeep <- rowSums(p_count >= 10)
# p_count <- p_count[pkeep,]


p_count <- p_count[, rownames(pfc_data1)]
all(rownames(pfc_data1) == colnames(p_count)) #check


#MeA
dag <- DESeqDataSetFromMatrix(countData = a_count,
                              colData = amy_data1,
                              design = ~condition1)



nrow(dag) # 12776 lost about half

dag_sub <- dag[, dag$condition1 %in% c("Subordinate", "Ascender")]
dag_sub$condition1 <- droplevels(dag_sub$condition1)
dag_sub$condition1
dag<- dag_sub
dag$condition1



#mPFC
dpg <- DESeqDataSetFromMatrix(countData = p_count,
                              colData = pfc_data1,
                              design = ~condition1)

nrow(dpg) # 14875

dpg_sub <- dpg[, dpg$condition1 %in% c("Subordinate", "Ascender")]
dpg_sub$condition1 <- droplevels(dpg_sub$condition1)
dpg_sub$condition1
dpg<- dpg_sub
dpg$condition1


#dds
da_dom <- DESeq(dag)
dp_dom <- DESeq(dpg)


vsd <- vst(da_dom, blind=FALSE) 
plotPCA(vsd, intgroup = c("condition1"), returnData = T) %>%
  left_join(coldata %>% 
              rename(name = SampleNames))-> d
ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_manual(values = c("#29AF7FFF", "#DCE319FF"))+
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       color = "Status",
       fill = "Status") +theme_bw()+
  theme(legend.position = c(0.83,0.15),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        text=element_text(size=15))+
  ggtitle('25hr MeA') -> p
ggsave("manuscript/brain/results/results_figures/ASUB25_MeA_PCA.png",p, dpi=300)
dev.off()
limma_list<- readRDS("manuscript/brain/results/limma_MeA_25hrASUB.RDS") 

limma_list %>%
  map(~arrange(., P.Value)) %>%
  map(~head(.,20))
p <- limma_list$status$P.Value
nn = length(p)
xx =  -log10((1:nn)/(nn+1))
p_cort <- limma_list$cort$P.Value
length(p) == length(p_cort) # TRUE
sum(p == p_cort)
df <- cbind(exp = xx,
            Status = -sort(log10(p)),
            CORT = -sort(log10(p_cort))) %>%
  as.data.frame() %>%
  gather(key,value,2:3) %>%
  mutate(key = factor(key, levels = c("Status","CORT")))
lim = log10(5000)


# png(filename = glue("results_figures/QQplot_{my_tissue}.png"),
#     width = 8, height = 9, units = "cm", res = 600)


df %>%
  filter(value < lim) %>%
  ggplot(aes(exp, value, color = key))+
  geom_point(shape  = 21,size =1)+
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype ="dashed")+
  theme_bw()+
  labs(x = expression(-log[10]~p-values~(expected)),
       y = expression(-log[10]~p-values~(observed)),
       color = "",
       fill = "")+
  theme(legend.position = c(0.8,0.2),
        legend.key.height = unit(0,"mm"),
        legend.title = element_blank(),
        text=element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 3.3)))+
  scale_color_manual(values = c("#E7B800","#00AFBB" ))+
  scale_fill_manual(values = c("#E7B800","#00AFBB" ))+
  scale_x_continuous(limits = c(0, 4), expand=expansion(mult=c(0,0.0)))+
  scale_y_continuous(limits = c(0, 4), expand=expansion(mult=c(0,0.0))) -> p2


ggsave("manuscript/brain/results/results_figures/ASUB25MeA_QQ.png",p2, dpi=300)
# "#29AF7FFF", "#FDE725FF"

dev.off()

vsd <- vst(dp_dom, blind=FALSE)
plotPCA(vsd, intgroup = c("condition1"), returnData = T) %>%
  left_join(coldata %>% 
              rename(name = SampleNames))-> d

ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) +
  geom_point(size = 3, alpha = 0.5) + 
  # geom_text(aes(label = subjectID), vjust = -0.5)+
  scale_color_manual(values = c("#29AF7FFF", "#DCE319FF")) + 
  labs(x = paste0("PC1: ", 
                  round(attributes(d)$percentVar[1] * 100), 
                  "% variance"),
       y = paste0("PC2: ", 
                  round(attributes(d)$percentVar[2] * 100), 
                  "% variance"),
       color = "Status",
       fill = "Status") + theme_bw()+
  theme(legend.position = c(0.20,0.9),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        text=element_text(size=15))+ 
  ggtitle('25hr mPFC') -> p

ggsave("manuscript/brain/results/results_figures/ASUB25_mPFCPCA.png",p, dpi=300)
dev.off()
limma_list<- readRDS("manuscript/brain/results/limma_mPFC_25hrASUB.RDS") 

limma_list %>%
  map(~arrange(., P.Value)) %>%
  map(~head(.,20))
p <- limma_list$status$P.Value
nn = length(p)
xx =  -log10((1:nn)/(nn+1))
p_cort <- limma_list$cort$P.Value
length(p) == length(p_cort) # TRUE
sum(p == p_cort)
df <- cbind(exp = xx,
            Status = -sort(log10(p)),
            CORT = -sort(log10(p_cort))) %>%
  as.data.frame() %>%
  gather(key,value,2:3) %>%
  mutate(key = factor(key, levels = c("Status","CORT")))
lim = log10(5000)


# png(filename = glue("results_figures/QQplot_{my_tissue}.png"),
#     width = 8, height = 9, units = "cm", res = 600)


df %>%
  filter(value < lim) %>%
  ggplot(aes(exp, value, color = key))+
  geom_point(shape  = 21,size =1)+
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype ="dashed")+
  theme_bw()+
  labs(x = expression(-log[10]~p-values~(expected)),
       y = expression(-log[10]~p-values~(observed)),
       color = "",
       fill = "")+
  theme(legend.position = c(0.8,0.2),
        legend.key.height = unit(0,"mm"),
        legend.title = element_blank(),
        text=element_text(size=15))+
  guides(color = guide_legend(override.aes = list(size = 3.3)))+
  scale_color_manual(values = c("#E7B800","#00AFBB" ))+
  scale_fill_manual(values = c("#E7B800","#00AFBB" ))+
  scale_x_continuous(limits = c(0, 4), expand=expansion(mult=c(0,0.0)))+
  scale_y_continuous(limits = c(0, 4), expand=expansion(mult=c(0,0.0))) -> p2


ggsave("manuscript/brain/results/results_figures/ASUB25mPFC_QQ.png",p2, dpi=300)
