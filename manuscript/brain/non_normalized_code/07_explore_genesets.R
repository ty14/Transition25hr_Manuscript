library(annotables)
library(tidyverse)
grcm38 # mouse genes



chol <- c('Chrm2', 'Chrm4','Gnai1','Gnai2','Chrm1','Chrm3','Chrm5','Gna11', 'Gna14','Gnaq','Ache',
          'Chat','Grk2', 'Grk5', 'Rgs2', 'Rgs4', 'Rgs6', 'Slc18a3', 'Slc5a7', 'Nat1', 'Lhx8', 'Slc10a4', 
          'Gbx1','Chrna2', 'Chrna3', 'Chrna6', 'Chrna7', 'Chrnb4', 'Chrnb3', 'Agrn', 'Chrna1', 'Chrna10', 
          'Chrna4', 'Chrna5', 'Chrna9', 'Chrnb1', 'Chrnb2', 'Chrnd', 'Chrne', 'Chrng', 'Dok7', 'Lrp4', 'Musk',
          'Rapsn')


my <- c('Cnp', 'Mobp', 'Mbp', 'Myrf', 'Mal', 'Bcas1', 'Mog', 'Mag', 'Lpar1',  'Plp1', 'Tspan2', 'Cntn2', 'Opalin', 'Arhgef10')

ar <- c("Maoa", "Erbb4", "Gria3", 'Mecp2', 'Prnp', 'Avpra1', 'Chrmp2b', "En2", "Fgf14", "Hdac4", "kcnj18", "Lrrc7", "Ache", "Cadm1", 
        "Crhr1", "Dnajb5","Ecm1","Eef1a2", "Ehmt1","Gad2", "Gdl1",
        "Grid1", "Grn", "Gsk3a", "Hsf1", "Lama2", "Mapk15", "Mme", 
        "Nfkb1", "Npy1r", "Osmr", "Pnoc", "Rbfox1", "Spast","Syn1", "Wdr62")



wnt <- c("Fosl1", "Wnt1", "T", "Wnt16", "Wisp1", "Fgf4", "Wnt6",
         "Fzd8", "Pitx2", "Wnt11", "Wnt7b", "Pygo1", "Wnt9a", "Ctbp2",
         "Myc", "Wnt10a", "Wif1", "Wnt3", "Apc", "Fbxw4", 
         "Wnt5a", "Ccd2", "Gsk3b", "Frat1", "Tcf3", "Sox17", "Ctbp1", 
         "Csnk1d", "Wnt2b", "Daam1", "Bcl9", "Sfrp1", "Tcf7", "Frzb", 
         "Jun", "Nkd1", "Wnt7a", "Fbxw2", "Dvl2", "Fzd3", "Aes", "Ctnnbip1",
         "Wnt4", "Ep300", "Tle1", "Wnt5b", "Sfrp2", "Fzd1", "Dvl1", "Lrp5", 
         "Sfrp4", "Kremen1", "Csnk1a1", "Ctnnb1", "Porcn", "Axin1", 
         "Btrc", "Senp2", "Csn2a1", "Ppp2ca", "Nlk", "Ccnd1", "Ppp2r1a", 
         "Ppp2r5d", "Lrp6","Rhou", "Fbxw11", "Slc9a3r1", "Wnt8a", "Foxn1", 
         "Ccnd3", "Lef1", "Fzd6", "Dixdc1", "Wnt2", "Dkk1", "Fzd4", "Fzd7",
         "Fshb", "Wnt3a", "Tle2", "Fzd5")

ser <- read_csv("manuscript/brain/gene_sets/Serotonergic_genes.csv")
 head(ser)
ser %>%
as_tibble() %>%.$symbol -> serx
serx <- str_to_title(serx)  


ser <- read_csv("manuscript/brain/gene_sets/Dopaminergic_genes.csv")
head(ser)
ser %>%
  as_tibble() %>%.$symbol -> serx
serx <- str_to_title(serx)  


mt <- read_csv("manuscript/brain/gene_sets/Mouse.MitoCarta3.csv")
head(mt)
mt %>%
  as_tibble() %>%.$Symbol -> mtx

#70 min 
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_PFC70min_ReorganizedGroup.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>%
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>%
  map(~filter(.,!is.na(entrez))) 

dd <- limma_list$desdom

ds <- limma_list$dessub

ad <- limma_list$ascdom

as <- limma_list$ascsub



dd %>% filter(symbol %in% chol)#0
ds %>% filter(symbol %in% chol)#2
ad %>% filter(symbol %in% chol)#1
as %>% filter(symbol %in% chol)#2



dd %>% filter(symbol %in% my)#0
ds %>% filter(symbol %in% my)#10
ad %>% filter(symbol %in% my)#0
as %>% filter(symbol %in% my)#5


dd %>% filter(symbol %in% ar)#0
ds %>% filter(symbol %in% ar)#6
ad %>% filter(symbol %in% ar)#0
as %>% filter(symbol %in% ar)#1
#Nfkb1 - transition gene 


dd %>% filter(symbol %in% wnt)#1
ds %>% filter(symbol %in% wnt)#9
ad %>% filter(symbol %in% wnt)#1 
as %>% filter(symbol %in% wnt)#3



dd %>% filter(symbol %in% serx)#0
ds %>% filter(symbol %in% serx)#2
ad %>% filter(symbol %in% serx)#1
as %>% filter(symbol %in% serx)#2

dd %>% filter(symbol %in% mtx)#57
ds %>% filter(symbol %in% mtx)#101
ad %>% filter(symbol %in% mtx)#63
as %>% filter(symbol %in% mtx)#66


#25 hour
my_logFC_threshold = 0.2

limma_list<- readRDS("manuscript/brain/results/limma_PFC_ReorganizedGroups_outlierRemoved.RDS") %>% 
  map(~distinct(.)) %>% 
  map(~filter(.,abs(logFC) >= my_logFC_threshold)) %>%
  map(~filter(.,P.Value <0.05)) %>% 
  map(~ left_join(., grcm38 %>% dplyr::select(symbol, entrez))) %>% 
  map(~filter(.,!is.na(entrez))) 

dd25 <- limma_list$desdom

ds25 <- limma_list$dessub

ad25 <- limma_list$ascdom

as25 <- limma_list$ascsub


dd25 %>% filter(symbol %in% mtx)#22
ds25 %>% filter(symbol %in% mtx)#53
ad25 %>% filter(symbol %in% mtx)#37
as25 %>% filter(symbol %in% mtx)#53

dd25 %>% filter(symbol %in% chol)#0
ds25 %>% filter(symbol %in% chol)#4
ad25 %>% filter(symbol %in% chol)#1
as25 %>% filter(symbol %in% chol)#1



dd25 %>% filter(symbol %in% my)#0
ds25 %>% filter(symbol %in% my)#2
ad25 %>% filter(symbol %in% my)#0
as25 %>% filter(symbol %in% my)#4


dd25 %>% filter(symbol %in% ar)#1
ds25 %>% filter(symbol %in% ar)#1
ad25 %>% filter(symbol %in% ar)#1
as25 %>% filter(symbol %in% ar)#1
#Nfkb1 - transition gene 


dd25 %>% filter(symbol %in% wnt)#3
ds25 %>% filter(symbol %in% wnt)#4
ad25 %>% filter(symbol %in% wnt)#1 
as25 %>% filter(symbol %in% wnt)#2

dd25 %>% filter(symbol %in% serx)#0
ds25 %>% filter(symbol %in% serx)#2
ad25 %>% filter(symbol %in% serx)#1
as25 %>% filter(symbol %in% serx)#2




# look at endrocrine pathway ways, kegg pathways
# how do I know C7 is immune hallmark gene set? 
# Go to; http://www.gsea-msigdb.org/gsea/msigdb/collection_details.jsp#C7
library(msigdbr)
msigdbr_collections() %>% as.data.frame()
#immune =  C7     IMMUNESIGDB
# biological pathways = C2 CP:WIKIPATHWAYS  
#hallmark = H
#reactome = immune, metabolic 


msigdbr(species = "Mus musculus", subcategory ="CP:REACTOME") -> reactome_sets
unique(reactome_sets$gs_name)

#IMMUNE 
reactome_sets %>% 
  filter(grepl("immune",gs_name,ignore.case = T)) -> my_gs_sets

unique(my_gs_sets$gs_name)

#adaptive immune = 1263 genes 
reactome_sets %>% 
  filter(grepl("REACTOME_ADAPTIVE_IMMUNE_SYSTEM",gs_name,ignore.case = T)) %>% select(gene_symbol) -> a_immune

dd %>% filter(symbol %in% a_immune$gene_symbol)#27
ds %>% filter(symbol %in% a_immune$gene_symbol)#57
ad %>% filter(symbol %in% a_immune$gene_symbol)#11 
as %>% filter(symbol %in% a_immune$gene_symbol)#15


dd25 %>% filter(symbol %in% a_immune$gene_symbol)#14
ds25 %>% filter(symbol %in% a_immune$gene_symbol)#25
ad25 %>% filter(symbol %in% a_immune$gene_symbol)#16 
as25 %>% filter(symbol %in% a_immune$gene_symbol)#24


# "REACTOME_INNATE_IMMUNE_SYSTEM" #827 genes 
 reactome_sets %>% 
  filter(grepl("REACTOME_INNATE_IMMUNE_SYSTEM",gs_name,ignore.case = T)) %>% select(gene_symbol) -> in_immune

dd %>% filter(symbol %in% in_immune$gene_symbol)#25
ds %>% filter(symbol %in% in_immune$gene_symbol)#57 
ad %>% filter(symbol %in% in_immune$gene_symbol)#11 
as %>% filter(symbol %in% in_immune$gene_symbol)#15

dd25 %>% filter(symbol %in% in_immune$gene_symbol)#11
ds25 %>% filter(symbol %in% in_immune$gene_symbol)#27
ad25 %>% filter(symbol %in% in_immune$gene_symbol)#9
as25 %>% filter(symbol %in% in_immune$gene_symbol)#27

# "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM" 
#827 genes 
reactome_sets %>% 
  filter(grepl("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",gs_name,ignore.case = T)) -> cyt_immune

dd %>% filter(symbol %in% cyt_immune$gene_symbol)#21
ds %>% filter(symbol %in% cyt_immune$gene_symbol)#41
ad %>% filter(symbol %in% cyt_immune$gene_symbol)#20
as %>% filter(symbol %in% cyt_immune$gene_symbol)#20

dd25 %>% filter(symbol %in% cyt_immune$gene_symbol)#11
ds25 %>% filter(symbol %in% cyt_immune$gene_symbol)#27 
ad25 %>% filter(symbol %in% cyt_immune$gene_symbol)#9 
as25 %>% filter(symbol %in% cyt_immune$gene_symbol)#17


# [5] "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE"  
# # "REACTOME_PRC2_METHYLATES_HISTONES_AND_DNA"
# 11] "REACTOME_DNA_METHYLATION"                                                                       
# [12] "REACTOME_DNA_REPAIR"
# # [57] "REACTOME_TRANSCRIPTIONAL_REGULATION_BY_SMALL_RNAS"
# [27] "REACTOME_NEGATIVE_EPIGENETIC_REGULATION_OF_RRNA_EXPRESSION"                                                    
# [28] "REACTOME_PIWI_INTERACTING_RNA_PIRNA_BIOGENESIS"                                                                
# [29] "REACTOME_POSITIVE_EPIGENETIC_REGULATION_OF_RRNA_EXPRESSION" 
# [1] "REACTOME_EPIGENETIC_REGULATION_OF_GENE_EXPRESSION"  
#87 genes 
reactome_sets %>% 
  filter(grepl("REACTOME_DNA_REPAIR" ,gs_name,ignore.case = T)) -> cs

dd %>% filter(symbol %in% cs$gene_symbol)#2
ds %>% filter(symbol %in% cs$gene_symbol)#8
ad %>% filter(symbol %in% cs$gene_symbol)#1
as %>% filter(symbol %in% cs$gene_symbol)#2

dd25 %>% filter(symbol %in% cs$gene_symbol)#2
ds25 %>% filter(symbol %in% cs$gene_symbol)#6 
ad25 %>% filter(symbol %in% cs$gene_symbol)#2
as25 %>% filter(symbol %in% cs$gene_symbol)#2


# [6] "REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES"
#30 genes 

reactome_sets %>% 
  filter(grepl("REACTOME_FOXO_MEDIATED_TRANSCRIPTION_OF_OXIDATIVE_STRESS_METABOLIC_AND_NEURONAL_GENES",gs_name,ignore.case = T)) -> fox

dd %>% filter(symbol %in% fox$gene_symbol)#1
ds %>% filter(symbol %in% fox$gene_symbol)#5
ad %>% filter(symbol %in% fox$gene_symbol)#3
as %>% filter(symbol %in% fox$gene_symbol)#2

dd25 %>% filter(symbol %in% fox$gene_symbol)#1
ds25 %>% filter(symbol %in% fox$gene_symbol)#4 
ad25 %>% filter(symbol %in% fox$gene_symbol)#1
as25 %>% filter(symbol %in% fox$gene_symbol)#0

# [7] "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE"

reactome_sets %>% 
  filter(grepl("REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE",gs_name,ignore.case = T)) -> oxi #132

dd %>% filter(symbol %in% oxi$gene_symbol)#4
ds %>% filter(symbol %in% oxi$gene_symbol)#9
ad %>% filter(symbol %in% oxi$gene_symbol)#5
as %>% filter(symbol %in% oxi$gene_symbol)#4

dd25 %>% filter(symbol %in% oxi$gene_symbol)#2
ds25 %>% filter(symbol %in% oxi$gene_symbol)#6 
ad25 %>% filter(symbol %in% oxi$gene_symbol)#2
as25 %>% filter(symbol %in% oxi$gene_symbol)#2

 
reactome_sets %>% 
  filter(grepl("oxida",gs_name,ignore.case = T)) -> my_gs_sets
unique(my_gs_sets$gs_name)
#to many 
# [47] "REACTOME_METABOLISM_OF_RNA"      
reactome_sets %>% 
  filter(grepl("REACTOME_METABOLISM_OF_RNA",gs_name,ignore.case = T)) -> mr #662

dd %>% filter(symbol %in% mr$gene_symbol)#29
ds %>% filter(symbol %in% mr$gene_symbol)#26
ad %>% filter(symbol %in% mr$gene_symbol)#36
as %>% filter(symbol %in% mr$gene_symbol)#34

dd25 %>% filter(symbol %in% mr$gene_symbol)#13
ds25 %>% filter(symbol %in% mr$gene_symbol)#27
ad25 %>% filter(symbol %in% mr$gene_symbol)#27
as25 %>% filter(symbol %in% mr$gene_symbol)#47



# [49] "REACTOME_METABOLISM_OF_STEROIDS"
reactome_sets %>% 
  filter(grepl("REACTOME_METABOLISM_OF_STEROIDS",gs_name,ignore.case = T)) -> ms #160

dd %>% filter(symbol %in% ms$gene_symbol)#6
ds %>% filter(symbol %in% ms$gene_symbol)#13
ad %>% filter(symbol %in% ms$gene_symbol)#6
as %>% filter(symbol %in% ms$gene_symbol)#5

dd25 %>% filter(symbol %in% ms$gene_symbol)#1
ds25 %>% filter(symbol %in% ms$gene_symbol)#4
ad25 %>% filter(symbol %in% ms$gene_symbol)#6
as25 %>% filter(symbol %in% ms$gene_symbol)#2




# [30] "REACTOME_INTEGRATION_OF_ENERGY_METABOLISM"  
reactome_sets %>% 
  filter(grepl("REACTOME_INTEGRATION_OF_ENERGY_METABOLISM",gs_name,ignore.case = T)) -> e#108

dd %>% filter(symbol %in% e$gene_symbol)#8
ds %>% filter(symbol %in% e$gene_symbol)#17
ad %>% filter(symbol %in% e$gene_symbol)#7
as %>% filter(symbol %in% e$gene_symbol)#3

dd25 %>% filter(symbol %in% e$gene_symbol)#1
ds25 %>% filter(symbol %in% e$gene_symbol)#4
ad25 %>% filter(symbol %in% e$gene_symbol)#1
as25 %>% filter(symbol %in% e$gene_symbol)#0



# [19] "REACTOME_FATTY_ACID_METABOLISM"
fat$gs_name
reactome_sets %>% 
  filter(grepl("REACTOME_GLYCOLYSIS",gs_name,ignore.case = T)) -> fat#186

dd %>% filter(symbol %in% fat$gene_symbol)#10
ds %>% filter(symbol %in% fat$gene_symbol)#23
ad %>% filter(symbol %in% fat$gene_symbol)#2
as %>% filter(symbol %in% fat$gene_symbol)#5

dd25 %>% filter(symbol %in% fat$gene_symbol)#3
ds25 %>% filter(symbol %in% fat$gene_symbol)#6
ad25 %>% filter(symbol %in% fat$gene_symbol)#3
as25 %>% filter(symbol %in% fat$gene_symbol)#8



# [1] "REACTOME_BIOLOGICAL_OXIDATIONS"     
reactome_sets %>% 
  filter(grepl("REACTOME_BIOLOGICAL_OXIDATIONS",gs_name,ignore.case = T)) -> bo#238

dd %>% filter(symbol %in% bo$gene_symbol)#7
ds %>% filter(symbol %in% bo$gene_symbol)#14
ad %>% filter(symbol %in% bo$gene_symbol)#8
as %>% filter(symbol %in% bo$gene_symbol)#8

dd25 %>% filter(symbol %in% bo$gene_symbol)#4
ds25 %>% filter(symbol %in% bo$gene_symbol)#9
ad25 %>% filter(symbol %in% bo$gene_symbol)#4
as25 %>% filter(symbol %in% bo$gene_symbol)#5





reactome_sets %>% 
  filter(grepl("methylation",gs_name,ignore.case = T)) -> my_gs_sets
# [1] "REACTOME_METHYLATION"   
reactome_sets %>% 
  filter(grepl("REACTOME_DNA_METHYLATION",gs_name,ignore.case = T)) -> meth#73

dd %>% filter(symbol %in% meth$gene_symbol)#0
ds %>% filter(symbol %in% meth$gene_symbol)#5
ad %>% filter(symbol %in% meth$gene_symbol)#0
as %>% filter(symbol %in% meth$gene_symbol)#1

dd25 %>% filter(symbol %in% meth$gene_symbol)#4
ds25 %>% filter(symbol %in% meth$gene_symbol)#9
ad25 %>% filter(symbol %in% meth$gene_symbol)#4
as25 %>% filter(symbol %in% meth$gene_symbol)#5

                        

#halmark sets 
msigdbr(species = "Mus musculus", category ="H") -> hallmark_sets
unique(hallmark_sets$gs_name)

# "HALLMARK_ANDROGEN_RESPONSE" 
    
hallmark_sets %>% 
  filter(grepl("HALLMARK_OXIDATIVE_PHOSPHORYLATION",gs_name,ignore.case = T)) -> ando # 126


dd %>% filter(symbol %in% ando$gene_symbol)#8
ds %>% filter(symbol %in% ando$gene_symbol)#21
ad %>% filter(symbol %in% ando$gene_symbol)#3
as %>% filter(symbol %in% ando$gene_symbol)#6

dd25 %>% filter(symbol %in% ando$gene_symbol)#2
ds25 %>% filter(symbol %in% ando$gene_symbol)#2
ad25 %>% filter(symbol %in% ando$gene_symbol)#3
as25 %>% filter(symbol %in% ando$gene_symbol)#2


# "HALLMARK_APOPTOSIS"   
hallmark_sets %>% 
  filter(grepl("HALLMARK_APOPTOSIS",gs_name,ignore.case = T)) -> apop # 126


dd %>% filter(symbol %in% apop$gene_symbol)#6
ds %>% filter(symbol %in% apop$gene_symbol)#20
ad %>% filter(symbol %in% apop$gene_symbol)#3
as %>% filter(symbol %in% apop$gene_symbol)#5

dd25 %>% filter(symbol %in% apop$gene_symbol)#3
ds25 %>% filter(symbol %in% apop$gene_symbol)#9
ad25 %>% filter(symbol %in% apop$gene_symbol)#6
as25 %>% filter(symbol %in% apop$gene_symbol)#8

# "HALLMARK_CHOLESTEROL_HOMEOSTASIS"  
hallmark_sets %>% 
  filter(grepl("HALLMARK_CHOLESTEROL_HOMEOSTASIS",gs_name,ignore.case = T)) -> ch # 75


dd %>% filter(symbol %in% ch$gene_symbol)#6
ds %>% filter(symbol %in% ch$gene_symbol)#9
ad %>% filter(symbol %in% ch$gene_symbol)#4
as %>% filter(symbol %in% ch$gene_symbol)#3

dd25 %>% filter(symbol %in% ch$gene_symbol)#4
ds25 %>% filter(symbol %in% ch$gene_symbol)#2
ad25 %>% filter(symbol %in% ch$gene_symbol)#3
as25 %>% filter(symbol %in% ch$gene_symbol)#3

# "HALLMARK_COAGULATION"
hallmark_sets %>% 
  filter(grepl("HALLMARK_COAGULATION",gs_name,ignore.case = T)) -> coa # 75


dd %>% filter(symbol %in% coa$gene_symbol)#3
ds %>% filter(symbol %in% coa$gene_symbol)#14
ad %>% filter(symbol %in% coa$gene_symbol)#4
as %>% filter(symbol %in% coa$gene_symbol)#4

dd25 %>% filter(symbol %in% coa$gene_symbol)#2
ds25 %>% filter(symbol %in% coa$gene_symbol)#5
ad25 %>% filter(symbol %in% coa$gene_symbol)#0
as25 %>% filter(symbol %in% coa$gene_symbol)#4

# "HALLMARK_CHOLESTEROL_HOMEOSTASIS"  
hallmark_sets %>% 
  filter(grepl("HALLMARK_CHOLESTEROL_HOMEOSTASIS",gs_name,ignore.case = T)) -> ch # 75


dd %>% filter(symbol %in% ch$gene_symbol)#6
ds %>% filter(symbol %in% ch$gene_symbol)#9
ad %>% filter(symbol %in% ch$gene_symbol)#4
as %>% filter(symbol %in% ch$gene_symbol)#3

dd25 %>% filter(symbol %in% ch$gene_symbol)#4
ds25 %>% filter(symbol %in% ch$gene_symbol)#2
ad25 %>% filter(symbol %in% ch$gene_symbol)#3
as25 %>% filter(symbol %in% ch$gene_symbol)#3

# "HALLMARK_DNA_REPAIR"   
hallmark_sets %>% 
  filter(grepl("HALLMARK_DNA_REPAIR",gs_name,ignore.case = T)) -> dre # 150


dd %>% filter(symbol %in% dre$gene_symbol)#7
ds %>% filter(symbol %in% dre$gene_symbol)#5
ad %>% filter(symbol %in% dre$gene_symbol)#6
as %>% filter(symbol %in% dre$gene_symbol)#5

dd25 %>% filter(symbol %in% dre$gene_symbol)#5
ds25 %>% filter(symbol %in% dre$gene_symbol)#4
ad25 %>% filter(symbol %in% dre$gene_symbol)#1
as25 %>% filter(symbol %in% dre$gene_symbol)#4


# "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
hallmark_sets %>% 
  filter(grepl("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",gs_name,ignore.case = T)) -> emt # 206


dd %>% filter(symbol %in% emt$gene_symbol)#13
ds %>% filter(symbol %in% emt$gene_symbol)#34
ad %>% filter(symbol %in% emt$gene_symbol)#5
as %>% filter(symbol %in% emt$gene_symbol)#6

dd25 %>% filter(symbol %in% emt$gene_symbol)#2
ds25 %>% filter(symbol %in% emt$gene_symbol)#2
ad25 %>% filter(symbol %in% emt$gene_symbol)#2
as25 %>% filter(symbol %in% emt$gene_symbol)#3

# "HALLMARK_FATTY_ACID_METABOLISM"  
hallmark_sets %>% 
  filter(grepl("HALLMARK_FATTY_ACID_METABOLISM",gs_name,ignore.case = T)) -> emt # 158


dd %>% filter(symbol %in% emt$gene_symbol)#6
ds %>% filter(symbol %in% emt$gene_symbol)#22
ad %>% filter(symbol %in% emt$gene_symbol)#6
as %>% filter(symbol %in% emt$gene_symbol)#9

dd25 %>% filter(symbol %in% emt$gene_symbol)#2
ds25 %>% filter(symbol %in% emt$gene_symbol)#6
ad25 %>% filter(symbol %in% emt$gene_symbol)#3
as25 %>% filter(symbol %in% emt$gene_symbol)#5

# "HALLMARK_GLYCOLYSIS" 
hallmark_sets %>% 
  filter(grepl("HALLMARK_GLYCOLYSIS",gs_name,ignore.case = T)) -> emt # 200


dd %>% filter(symbol %in% emt$gene_symbol)#10
ds %>% filter(symbol %in% emt$gene_symbol)#23
ad %>% filter(symbol %in% emt$gene_symbol)#17
as %>% filter(symbol %in% emt$gene_symbol)#12

dd25 %>% filter(symbol %in% emt$gene_symbol)#2
ds25 %>% filter(symbol %in% emt$gene_symbol)#15
ad25 %>% filter(symbol %in% emt$gene_symbol)#3
as25 %>% filter(symbol %in% emt$gene_symbol)#8

# "HALLMARK_INFLAMMATORY_RESPONSE" 
hallmark_sets %>% 
  filter(grepl("HALLMARK_INFLAMMATORY_RESPONSE",gs_name,ignore.case = T)) -> emt # 202


dd %>% filter(symbol %in% emt$gene_symbol)#9
ds %>% filter(symbol %in% emt$gene_symbol)#17
ad %>% filter(symbol %in% emt$gene_symbol)#1
as %>% filter(symbol %in% emt$gene_symbol)#6

dd25 %>% filter(symbol %in% emt$gene_symbol)#2
ds25 %>% filter(symbol %in% emt$gene_symbol)#7
ad25 %>% filter(symbol %in% emt$gene_symbol)#4
as25 %>% filter(symbol %in% emt$gene_symbol)#5

# "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
hallmark_sets %>% 
  filter(grepl("HALLMARK_OXIDATIVE_PHOSPHORYLATION",gs_name,ignore.case = T)) -> emt # 199


dd %>% filter(symbol %in% emt$gene_symbol)#3
ds %>% filter(symbol %in% emt$gene_symbol)#16
ad %>% filter(symbol %in% emt$gene_symbol)#16
as %>% filter(symbol %in% emt$gene_symbol)#10

dd25 %>% filter(symbol %in% emt$gene_symbol)#3
ds25 %>% filter(symbol %in% emt$gene_symbol)#2
ad25 %>% filter(symbol %in% emt$gene_symbol)#8
as25 %>% filter(symbol %in% emt$gene_symbol)#1
# "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"  

hallmark_sets %>% 
  filter(grepl("HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" ,gs_name,ignore.case = T)) -> emt # 49


dd %>% filter(symbol %in% emt$gene_symbol)#3
ds %>% filter(symbol %in% emt$gene_symbol)#6
ad %>% filter(symbol %in% emt$gene_symbol)#6
as %>% filter(symbol %in% emt$gene_symbol)#1

dd25 %>% filter(symbol %in% emt$gene_symbol)#0
ds25 %>% filter(symbol %in% emt$gene_symbol)#0
ad25 %>% filter(symbol %in% emt$gene_symbol)#5
as25 %>% filter(symbol %in% emt$gene_symbol)#2
 
# "HALLMARK_SPERMATOGENESIS"    
hallmark_sets %>% 
  filter(grepl("HALLMARK_SPERMATOGENESIS" ,gs_name,ignore.case = T)) -> emt # 135


dd %>% filter(symbol %in% emt$gene_symbol)#3
ds %>% filter(symbol %in% emt$gene_symbol)#13
ad %>% filter(symbol %in% emt$gene_symbol)#0
as %>% filter(symbol %in% emt$gene_symbol)#3

dd25 %>% filter(symbol %in% emt$gene_symbol)#2
ds25 %>% filter(symbol %in% emt$gene_symbol)#2
ad25 %>% filter(symbol %in% emt$gene_symbol)#3
as25 %>% filter(symbol %in% emt$gene_symbol)#6


# "HALLMARK_P53_PATHWAY" 
hallmark_sets %>% 
  filter(grepl("HALLMARK_P53_PATHWAY" ,gs_name,ignore.case = T)) -> emt # 54


dd %>% filter(symbol %in% emt$gene_symbol)#8
ds %>% filter(symbol %in% emt$gene_symbol)#18
ad %>% filter(symbol %in% emt$gene_symbol)#7
as %>% filter(symbol %in% emt$gene_symbol)#7

dd25 %>% filter(symbol %in% emt$gene_symbol)#7
ds25 %>% filter(symbol %in% emt$gene_symbol)#5
ad25 %>% filter(symbol %in% emt$gene_symbol)#11
as25 %>% filter(symbol %in% emt$gene_symbol)#5

msigdbr_collections() %>% as.data.frame()
msigdbr(species = "Mus musculus", subcategory ='CP:WIKIPATHWAYS' ) -> wp_sets
unique(wp_sets$gs_name)

# "WP_ANDROGEN_RECEPTOR_SIGNALING_PATHWAY"
# "WP_APOPTOSIS_MODULATION_AND_SIGNALING" 
# "WP_CORTICOTROPINRELEASING_HORMONE_SIGNALING_PATHWAY" 
                                                                         
# "WP_DOPAMINERGIC_NEUROGENESIS"   
# [190] "WP_FATTY_ACID_BIOSYNTHESIS"   
# [229] "WP_GLUCOCORTICOID_RECEPTOR_PATHWAY"                                                                
# [265] "WP_HISTONE_MODIFICATIONS" 
# [398] "WP_MONOAMINE_GPCRS"                                                                                
# [399] "WP_MONOAMINE_TRANSPORT"  
# [426] "WP_NEUROINFLAMMATION" 


# # WP_TNFALPHA_SIGNALING_PATHWAY
# 88] "WP_MITOCHONDRIAL_COMPLEX_I_ASSEMBLY_MODEL_OXPHOS_SYSTEM"                                           
# [389] "WP_MITOCHONDRIAL_COMPLEX_II_ASSEMBLY"                                                              
# [390] "WP_MITOCHONDRIAL_COMPLEX_III_ASSEMBLY"                                                             
# [391] "WP_MITOCHONDRIAL_COMPLEX_IV_ASSEMBLY"                      

# "WP_MITOCHONDRIAL_COMPLEX_" 
# WP_GPCRS_CLASS_C_METABOTROPIC_GLUTAMATE_PHEROMONE"  



wp_sets %>% 
  filter(grepl("WP_TNFALPHA_SIGNALING_PATHWAY" ,gs_name,ignore.case = T))-> glu
dd %>% filter(symbol %in% glu$gene_symbol)#8
ds %>% filter(symbol %in% glu$gene_symbol)#18
ad %>% filter(symbol %in% glu$gene_symbol)#7
as %>% filter(symbol %in% glu$gene_symbol)#7

dd25 %>% filter(symbol %in% glu$gene_symbol)#7
ds25 %>% filter(symbol %in% glu$gene_symbol)#5
ad25 %>% filter(symbol %in% glu$gene_symbol)#11
as25 %>% filter(symbol %in% glu$gene_symbol)#5
# ] "WP_BRAINDERIVED_NEUROTROPHIC_FACTOR_BDNF_SIGNALING_PATHWAY"   
# [427] "WP_NEUROINFLAMMATION_AND_GLUTAMATERGIC_SIGNALING" 
wp_sets %>% 
  filter(grepl( "WP_BRAINDERIVED_NEUROTROPHIC_FACTOR_BDNF_SIGNALING_PATHWAY" ,gs_name,ignore.case = T))-> glu
dd %>% filter(symbol %in% glu$gene_symbol)#8
ds %>% filter(symbol %in% glu$gene_symbol)#18
ad %>% filter(symbol %in% glu$gene_symbol)#7
as %>% filter(symbol %in% glu$gene_symbol)#7

dd25 %>% filter(symbol %in% glu$gene_symbol)#7
ds25 %>% filter(symbol %in% glu$gene_symbol)#5
ad25 %>% filter(symbol %in% glu$gene_symbol)#11
as25 %>% filter(symbol %in% glu$gene_symbol)#5

# [473] "WP_OXIDATIVE_DAMAGE_RESPONSE"                                                                      
# [474] "WP_OXIDATIVE_PHOSPHORYLATION"                                                                      
# [475] "WP_OXIDATIVE_STRESS_RESPONSE"                                                               
                                                               
# [478] "WP_P53_TRANSCRIPTIONAL_GENE_NETWORK"  
# 62] "WP_SEROTONIN_AND_ANXIETY"                                                                          
# [563] "WP_SEROTONIN_AND_ANXIETYRELATED_EVENTS"   
# [564] "WP_SEROTONIN_HTR1_GROUP_AND_FOS_PATHWAY"  
wp_sets %>% 
  filter(grepl("WP_p38" ,gs_name,ignore.case = T))-> ser
dd %>% filter(symbol %in% ser$gene_symbol)#8
ds %>% filter(symbol %in% ser$gene_symbol)#18
ad %>% filter(symbol %in% ser$gene_symbol)#7
as %>% filter(symbol %in% ser$gene_symbol)#7

dd25 %>% filter(symbol %in% ser$gene_symbol)#7
ds25 %>% filter(symbol %in% ser$gene_symbol)#5
ad25 %>% filter(symbol %in% ser$gene_symbol)#11
as25 %>% filter(symbol %in% ser$gene_symbol)#5


# [565] "WP_SEROTONIN_RECEPTOR_2_AND_ELKSRFGATA4_SIGNALING"                                                 
# [566] "WP_SEROTONIN_RECEPTOR_467_AND_NR3C_SIGNALING"                                                      
# [567] "WP_SEROTONIN_TRANSPORTER_ACTIVITY" 
wp_sets %>% 
  filter(grepl("WP_STEROID_BIOSYNTHESIS" ,gs_name,ignore.case = T))-> ser
dd %>% filter(symbol %in% ser$gene_symbol)#8
ds %>% filter(symbol %in% ser$gene_symbol)#18
ad %>% filter(symbol %in% ser$gene_symbol)#7
as %>% filter(symbol %in% ser$gene_symbol)#7

dd25 %>% filter(symbol %in% ser$gene_symbol)#7
ds25 %>% filter(symbol %in% ser$gene_symbol)#5
ad25 %>% filter(symbol %in% ser$gene_symbol)#11
as25 %>% filter(symbol %in% ser$gene_symbol)#5
# [569] "WP_SLEEP_REGULATION"    

# [659] "WP_WNT_SIGNALING"                                                                                  
# [661] "WP_WNT_SIGNALING_PATHWAY" 

wp_sets %>% 
  filter(grepl("WP_WNT_SIGNALING_PATHWAY" ,gs_name,ignore.case = T))-> ser
dd %>% filter(symbol %in% ser$gene_symbol)#8
ds %>% filter(symbol %in% ser$gene_symbol)#18
ad %>% filter(symbol %in% ser$gene_symbol)#7
as %>% filter(symbol %in% ser$gene_symbol)#7

dd25 %>% filter(symbol %in% ser$gene_symbol)#7
ds25 %>% filter(symbol %in% ser$gene_symbol)#5
ad25 %>% filter(symbol %in% ser$gene_symbol)#11
as25 %>% filter(symbol %in% ser$gene_symbol)#5

msigdbr(species = "Mus musculus", subcategory ='GO:BP' ) -> bp_sets
unique(x$gs_name)

x <- bp_sets %>% filter(grepl("MITO", gs_name))
#GLUTaMATE
# "GOBP_SYNAPTIC_SIGNALING"
# # [25] "GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC" 
# [3] "GOBP_REGULATION_OF_LONG_TERM_NEURONAL_SYNAPTIC_PLASTICITY"         
# [4] "GOBP_REGULATION_OF_NEURONAL_SYNAPTIC_PLASTICITY"                   
# [5] "GOBP_REGULATION_OF_SHORT_TERM_NEURONAL_SYNAPTIC_PLASTICITY"        
# [6] "GOBP_REGULATION_OF_SYNAPSE_STRUCTURAL_PLASTICITY"                  
# [7] "GOBP_REGULATION_OF_SYNAPTIC_PLASTICITY"                            
#



# "GOBP_CELLULAR_RESPONSE_TO_STRESS"    
# [58] "GOBP_MITOCHONDRION_ORGANIZATION"  
                                                                            
bp_sets %>% 
  filter(grepl("GOBP_CELLULAR_RESPONSE_TO_STRESS" ,gs_name,ignore.case = T))-> glu
dd %>% filter(symbol %in% glu$gene_symbol)#8
ds %>% filter(symbol %in% glu$gene_symbol)#18
ad %>% filter(symbol %in% glu$gene_symbol)#7
as %>% filter(symbol %in% glu$gene_symbol)#7

dd25 %>% filter(symbol %in% glu$gene_symbol)#7
ds25 %>% filter(symbol %in% glu$gene_symbol)#5
ad25 %>% filter(symbol %in% glu$gene_symbol)#11
as25 %>% filter(symbol %in% glu$gene_symbol)#5


# [7] "GOBP_GLUTAMATE_METABOLIC_PROCESS"  
bp_sets %>% 
  filter(grepl("GOBP_GLUTAMATE_METABOLIC_PROCESS" ,gs_name,ignore.case = T))-> glu
dd %>% filter(symbol %in% glu$gene_symbol)#8
ds %>% filter(symbol %in% glu$gene_symbol)#18
ad %>% filter(symbol %in% glu$gene_symbol)#7
as %>% filter(symbol %in% glu$gene_symbol)#7

dd25 %>% filter(symbol %in% glu$gene_symbol)#7
ds25 %>% filter(symbol %in% glu$gene_symbol)#5
ad25 %>% filter(symbol %in% glu$gene_symbol)#11
as25 %>% filter(symbol %in% glu$gene_symbol)#5

# [8] "GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY"  
bp_sets %>% 
  filter(grepl("GOBP_GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY" ,gs_name,ignore.case = T))-> glu
dd %>% filter(symbol %in% glu$gene_symbol)#8
ds %>% filter(symbol %in% glu$gene_symbol)#18
ad %>% filter(symbol %in% glu$gene_symbol)#7
as %>% filter(symbol %in% glu$gene_symbol)#7

dd25 %>% filter(symbol %in% glu$gene_symbol)#7
ds25 %>% filter(symbol %in% glu$gene_symbol)#5
ad25 %>% filter(symbol %in% glu$gene_symbol)#11
as25 %>% filter(symbol %in% glu$gene_symbol)#5


#SERTONIN
# [3] "GOBP_REGULATION_OF_SEROTONIN_SECRETION"  
bp_sets %>% 
  filter(grepl("GOBP_REGULATION_OF_SEROTONIN_SECRETION" ,gs_name,ignore.case = T))-> ser
dd %>% filter(symbol %in% ser$gene_symbol)#8
ds %>% filter(symbol %in% ser$gene_symbol)#18
ad %>% filter(symbol %in% ser$gene_symbol)#7
as %>% filter(symbol %in% ser$gene_symbol)#7

dd25 %>% filter(symbol %in% ser$gene_symbol)#7
ds25 %>% filter(symbol %in% ser$gene_symbol)#5
ad25 %>% filter(symbol %in% ser$gene_symbol)#11
as25 %>% filter(symbol %in% ser$gene_symbol)#5

# [6] "GOBP_SEROTONIN_SECRETION"   
bp_sets %>% 
  filter(grepl("GOBP_SEROTONIN_SECRETION" ,gs_name,ignore.case = T))-> ser
dd %>% filter(symbol %in% ser$gene_symbol)#8
ds %>% filter(symbol %in% ser$gene_symbol)#18
ad %>% filter(symbol %in% ser$gene_symbol)#7
as %>% filter(symbol %in% ser$gene_symbol)#7

dd25 %>% filter(symbol %in% ser$gene_symbol)#7
ds25 %>% filter(symbol %in% ser$gene_symbol)#5
ad25 %>% filter(symbol %in% ser$gene_symbol)#11
as25 %>% filter(symbol %in% ser$gene_symbol)#5


# [7] "GOBP_SEROTONIN_TRANSPORT"  
bp_sets %>% 
  filter(grepl("GOBP_SEROTONIN_TRANSPORT" ,gs_name,ignore.case = T))-> ser
dd %>% filter(symbol %in% ser$gene_symbol)#8
ds %>% filter(symbol %in% ser$gene_symbol)#18
ad %>% filter(symbol %in% ser$gene_symbol)#7
as %>% filter(symbol %in% ser$gene_symbol)#7

dd25 %>% filter(symbol %in% ser$gene_symbol)#7
ds25 %>% filter(symbol %in% ser$gene_symbol)#5
ad25 %>% filter(symbol %in% ser$gene_symbol)#11
as25 %>% filter(symbol %in% ser$gene_symbol)#5

# [8] "GOBP_SEROTONIN_UPTAKE" 
bp_sets %>% 
  filter(grepl("GOBP_SEROTONIN_UPTAKE" ,gs_name,ignore.case = T))-> ser
dd %>% filter(symbol %in% ser$gene_symbol)#8
ds %>% filter(symbol %in% ser$gene_symbol)#18
ad %>% filter(symbol %in% ser$gene_symbol)#7
as %>% filter(symbol %in% ser$gene_symbol)#7

dd25 %>% filter(symbol %in% ser$gene_symbol)#7
ds25 %>% filter(symbol %in% ser$gene_symbol)#5
ad25 %>% filter(symbol %in% ser$gene_symbol)#11
as25 %>% filter(symbol %in% ser$gene_symbol)#5


#dopamine

# [4] "GOBP_DOPAMINE_METABOLIC_PROCESS"                   
bp_sets %>% 
  filter(grepl("GOBP_DOPAMINE_METABOLIC_PROCESS" ,gs_name,ignore.case = T))-> dop
dd %>% filter(symbol %in% dop$gene_symbol)#8
ds %>% filter(symbol %in% dop$gene_symbol)#18
ad %>% filter(symbol %in% dop$gene_symbol)#7
as %>% filter(symbol %in% dop$gene_symbol)#7

dd25 %>% filter(symbol %in% dop$gene_symbol)#7
ds25 %>% filter(symbol %in% dop$gene_symbol)#5
ad25 %>% filter(symbol %in% dop$gene_symbol)#11
as25 %>% filter(symbol %in% dop$gene_symbol)#5

# [5] "GOBP_DOPAMINE_RECEPTOR_SIGNALING_PATHWAY"           
bp_sets %>% 
  filter(grepl("GOBP_DOPAMINE_RECEPTOR_SIGNALING_PATHWAY" ,gs_name,ignore.case = T))-> dop
dd %>% filter(symbol %in% dop$gene_symbol)#8
ds %>% filter(symbol %in% dop$gene_symbol)#18
ad %>% filter(symbol %in% dop$gene_symbol)#7
as %>% filter(symbol %in% dop$gene_symbol)#7

dd25 %>% filter(symbol %in% dop$gene_symbol)#7
ds25 %>% filter(symbol %in% dop$gene_symbol)#5
ad25 %>% filter(symbol %in% dop$gene_symbol)#11
as25 %>% filter(symbol %in% dop$gene_symbol)#5

# [6] "GOBP_DOPAMINE_SECRETION"                            
bp_sets %>% 
  filter(grepl("GOBP_DOPAMINE_SECRETION" ,gs_name,ignore.case = T))-> dop
dd %>% filter(symbol %in% dop$gene_symbol)#8
ds %>% filter(symbol %in% dop$gene_symbol)#18
ad %>% filter(symbol %in% dop$gene_symbol)#7
as %>% filter(symbol %in% dop$gene_symbol)#7

dd25 %>% filter(symbol %in% dop$gene_symbol)#7
ds25 %>% filter(symbol %in% dop$gene_symbol)#5
ad25 %>% filter(symbol %in% dop$gene_symbol)#11
as25 %>% filter(symbol %in% dop$gene_symbol)#5

# [7] "GOBP_DOPAMINE_TRANSPORT" 
bp_sets %>% 
  filter(grepl("GOBP_DOPAMINE_TRANSPORT" ,gs_name,ignore.case = T))-> dop
dd %>% filter(symbol %in% dop$gene_symbol)#8
ds %>% filter(symbol %in% dop$gene_symbol)#18
ad %>% filter(symbol %in% dop$gene_symbol)#7
as %>% filter(symbol %in% dop$gene_symbol)#7

dd25 %>% filter(symbol %in% dop$gene_symbol)#7
ds25 %>% filter(symbol %in% dop$gene_symbol)#5
ad25 %>% filter(symbol %in% dop$gene_symbol)#11
as25 %>% filter(symbol %in% dop$gene_symbol)#5

# [19] "GOBP_RESPONSE_TO_DOPAMINE"  
bp_sets %>% 
  filter(grepl("GOBP_RESPONSE_TO_DOPAMINE" ,gs_name,ignore.case = T))-> dop

dd %>% filter(symbol %in% dop$gene_symbol)#8
ds %>% filter(symbol %in% dop$gene_symbol)#18
ad %>% filter(symbol %in% dop$gene_symbol)#7
as %>% filter(symbol %in% dop$gene_symbol)#7

dd25 %>% filter(symbol %in% dop$gene_symbol)#7
ds25 %>% filter(symbol %in% dop$gene_symbol)#5
ad25 %>% filter(symbol %in% dop$gene_symbol)#11
as25 %>% filter(symbol %in% dop$gene_symbol)#5
# [8] "GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION"     

bp_sets %>% 
  filter(grepl("GOBP_DOPAMINERGIC_NEURON_DIFFERENTIATION" ,gs_name,ignore.case = T))-> dop

dd %>% filter(symbol %in% dop$gene_symbol)#8
ds %>% filter(symbol %in% dop$gene_symbol)#18
ad %>% filter(symbol %in% dop$gene_symbol)#7
as %>% filter(symbol %in% dop$gene_symbol)#7

dd25 %>% filter(symbol %in% dop$gene_symbol)#7
ds25 %>% filter(symbol %in% dop$gene_symbol)#5
ad25 %>% filter(symbol %in% dop$gene_symbol)#11
as25 %>% filter(symbol %in% dop$gene_symbol)#5
# [20] "GOBP_SYNAPTIC_TRANSMISSION_DOPAMINERGIC"                                            

bp_sets %>% 
  filter(grepl("GOBP_SYNAPTIC_TRANSMISSION_DOPAMINERGIC" ,gs_name,ignore.case = T))-> dop

dd %>% filter(symbol %in% dop$gene_symbol)#8
ds %>% filter(symbol %in% dop$gene_symbol)#18
ad %>% filter(symbol %in% dop$gene_symbol)#7
as %>% filter(symbol %in% dop$gene_symbol)#7

dd25 %>% filter(symbol %in% dop$gene_symbol)#7
ds25 %>% filter(symbol %in% dop$gene_symbol)#5
ad25 %>% filter(symbol %in% dop$gene_symbol)#11
as25 %>% filter(symbol %in% dop$gene_symbol)#5



  #gaba

# [5] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC" 
bp_sets %>% 
  filter(grepl("GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC" ,gs_name,ignore.case = T))-> ga

dd %>% filter(symbol %in% ga$gene_symbol)#8
ds %>% filter(symbol %in% ga$gene_symbol)#18
ad %>% filter(symbol %in% ga$gene_symbol)#7
as %>% filter(symbol %in% ga$gene_symbol)#7

dd25 %>% filter(symbol %in% ga$gene_symbol)#7
ds25 %>% filter(symbol %in% ga$gene_symbol)#5
ad25 %>% filter(symbol %in% ga$gene_symbol)#11
as25 %>% filter(symbol %in% ga$gene_symbol)#5

# [7] "GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC" 
  
bp_sets %>% 
  filter(grepl("GOBP_SYNAPTIC_TRANSMISSION_GABAERGIC" ,gs_name,ignore.case = T))-> ga

dd %>% filter(symbol %in% ga$gene_symbol)#8
ds %>% filter(symbol %in% ga$gene_symbol)#18
ad %>% filter(symbol %in% ga$gene_symbol)#7
as %>% filter(symbol %in% ga$gene_symbol)#7

dd25 %>% filter(symbol %in% ga$gene_symbol)#7
ds25 %>% filter(symbol %in% ga$gene_symbol)#5
ad25 %>% filter(symbol %in% ga$gene_symbol)#11
as25 %>% filter(symbol %in% ga$gene_symbol)#5

#chol

# [4] "GOBP_CHOLINE_METABOLIC_PROCESS"  
bp_sets %>% 
  filter(grepl("GOBP_CHOLINE_METABOLIC_PROCESS" ,gs_name,ignore.case = T))-> ache

dd %>% filter(symbol %in% ache$gene_symbol)#8
ds %>% filter(symbol %in% ache$gene_symbol)#18
ad %>% filter(symbol %in% ache$gene_symbol)#7
as %>% filter(symbol %in% ache$gene_symbol)#7

dd25 %>% filter(symbol %in% ache$gene_symbol)#7
ds25 %>% filter(symbol %in% ache$gene_symbol)#5
ad25 %>% filter(symbol %in% ache$gene_symbol)#11
as25 %>% filter(symbol %in% ache$gene_symbol)#5
#Slc44a1, Enpp6


# [13] "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_CHOLINERGIC"   
bp_sets %>% 
  filter(grepl("GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_CHOLINERGIC" ,gs_name,ignore.case = T))-> ache

dd %>% filter(symbol %in% ache$gene_symbol)#8
ds %>% filter(symbol %in% ache$gene_symbol)#18
ad %>% filter(symbol %in% ache$gene_symbol)#7
as %>% filter(symbol %in% ache$gene_symbol)#7

dd25 %>% filter(symbol %in% ache$gene_symbol)#7
ds25 %>% filter(symbol %in% ache$gene_symbol)#5
ad25 %>% filter(symbol %in% ache$gene_symbol)#11
as25 %>% filter(symbol %in% ache$gene_symbol)#5
# [14] "GOBP_RESPONSE_TO_ACETYLCHOLINE" 
bp_sets %>% 
  filter(grepl("GOBP_RESPONSE_TO_ACETYLCHOLINE"  ,gs_name,ignore.case = T))-> ache

dd %>% filter(symbol %in% ache$gene_symbol)#8
ds %>% filter(symbol %in% ache$gene_symbol)#18
ad %>% filter(symbol %in% ache$gene_symbol)#7
as %>% filter(symbol %in% ache$gene_symbol)#7

dd25 %>% filter(symbol %in% ache$gene_symbol)#7
ds25 %>% filter(symbol %in% ache$gene_symbol)#5
ad25 %>% filter(symbol %in% ache$gene_symbol)#11
as25 %>% filter(symbol %in% ache$gene_symbol)#5

#norepinephrine


# [4] "GOBP_NOREPINEPHRINE_SECRETION"       
bp_sets %>% 
  filter(grepl("GOBP_NOREPINEPHRINE_SECRETION" ,gs_name,ignore.case = T))-> n

dd %>% filter(symbol %in% n$gene_symbol)#8
ds %>% filter(symbol %in% n$gene_symbol)#18
ad %>% filter(symbol %in% n$gene_symbol)#7
as %>% filter(symbol %in% n$gene_symbol)#7

dd25 %>% filter(symbol %in% n$gene_symbol)#7
ds25 %>% filter(symbol %in% n$gene_symbol)#5
ad25 %>% filter(symbol %in% n$gene_symbol)#11
as25 %>% filter(symbol %in% n$gene_symbol)#5
as25 %>% filter(symbol %in% n$gene_symbol)#5


#CATECHOLAMINE
# [1] "GOBP_CATECHOLAMINE_SECRETION"     
bp_sets %>% 
  filter(grepl("GOBP_CATECHOLAMINE_SECRETION" ,gs_name,ignore.case = T))-> cat

dd %>% filter(symbol %in% cat$gene_symbol)#8
ds %>% filter(symbol %in% cat$gene_symbol)#18
ad %>% filter(symbol %in% cat$gene_symbol)#7
as %>% filter(symbol %in% cat$gene_symbol)#7

dd25 %>% filter(symbol %in% cat$gene_symbol)#7
ds25 %>% filter(symbol %in% cat$gene_symbol)#5
ad25 %>% filter(symbol %in% cat$gene_symbol)#11
as25 %>% filter(symbol %in% cat$gene_symbol)#5
# [7] "GOBP_REGULATION_OF_CATECHOLAMINE_METABOLIC_PROCESS"
bp_sets %>% 
  filter(grepl("GOBP_REGULATION_OF_CATECHOLAMINE_METABOLIC_PROCESS" ,gs_name,ignore.case = T))-> cat

dd %>% filter(symbol %in% cat$gene_symbol)#8
ds %>% filter(symbol %in% cat$gene_symbol)#18
ad %>% filter(symbol %in% cat$gene_symbol)#7
as %>% filter(symbol %in% cat$gene_symbol)#7

dd25 %>% filter(symbol %in% cat$gene_symbol)#7
ds25 %>% filter(symbol %in% cat$gene_symbol)#5
ad25 %>% filter(symbol %in% cat$gene_symbol)#11
as25 %>% filter(symbol %in% cat$gene_symbol)#5

# "GOBP_RESPONSE_TO_CATECHOLAMINE"
bp_sets %>% 
  filter(grepl("GOBP_RESPONSE_TO_CATECHOLAMINE" ,gs_name,ignore.case = T))-> cat

dd %>% filter(symbol %in% cat$gene_symbol)#8
ds %>% filter(symbol %in% cat$gene_symbol)#18
ad %>% filter(symbol %in% cat$gene_symbol)#7
as %>% filter(symbol %in% cat$gene_symbol)#7

dd25 %>% filter(symbol %in% cat$gene_symbol)#7
ds25 %>% filter(symbol %in% cat$gene_symbol)#5
ad25 %>% filter(symbol %in% cat$gene_symbol)#11
as25 %>% filter(symbol %in% cat$gene_symbol)#5

