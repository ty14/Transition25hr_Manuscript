# make a central nervous system gene set 
library(msigdbr)
msigdbr_collections() %>% as.data.frame()
library(tidyverse)


#genes from the MeA paper 
my <- c('Cnp', 'Mobp', 'Mbp', 'Myrf', 'Mal', 'Bcas1', 'Mog', 'Mag', 
        'Lpar1',  'Plp1', 'Tspan2', 'Cntn2', 'Opalin', 'Arhgef10')

myx <- as.data.frame(my)
colnames(myx)[1] <- "gene_symbol"

#genes from bpsets for central nervous system myelin 
msigdbr(species = "Mus musculus", subcategory ='GO:BP' ) -> bp_sets
x <- bp_sets %>% filter(grepl("myelin", gs_name,ignore.case = T))
unique(x$gs_name)
# [1] "GOBP_CENTRAL_NERVOUS_SYSTEM_MYELIN_FORMATION"      "GOBP_CENTRAL_NERVOUS_SYSTEM_MYELIN_MAINTENANCE"   
# [3] "GOBP_MYELIN_ASSEMBLY"                              "GOBP_MYELIN_MAINTENANCE"                          
# [5] "GOBP_NEGATIVE_REGULATION_OF_MYELINATION"           "GOBP_PERIPHERAL_NERVOUS_SYSTEM_MYELIN_MAINTENANCE"
# [7] "GOBP_POSITIVE_REGULATION_OF_MYELINATION"           "GOBP_REGULATION_OF_MYELINATION"                   
# [9] "GOBP_SPHINGOMYELIN_BIOSYNTHETIC_PROCESS"           "GOBP_SPHINGOMYELIN_CATABOLIC_PROCESS"             
# [11] "GOBP_SPHINGOMYELIN_METABOLIC_PROCESS"             


#1 "GOBP_CENTRAL_NERVOUS_SYSTEM_MYELIN_FORMATION" 
#2 "GOBP_CENTRAL_NERVOUS_SYSTEM_MYELIN_MAINTENANCE" 
my1 <- x %>% 
  filter(grepl("CENTRAL_NERVOUS_SYSTEM_MYELIN",gs_name,ignore.case = T)) %>% select(gene_symbol) %>%  unique(.)
#3 "GOBP_REGULATION_OF_MYELINATION"            
my2 <- x %>% 
  filter(grepl("GOBP_REGULATION_OF_MYELINATION",gs_name,ignore.case = T)) %>% select(gene_symbol) %>%  unique(.)

#4 "GOBP_MYELIN_ASSEMBLY" 
my3 <- x %>% 
 filter(grepl("GOBP_MYELIN_ASSEMBLY",gs_name,ignore.case = T)) %>% select(gene_symbol) %>%  unique(.)


myelin_genes <- my1 %>% rbind(my2, my3, myx) %>% unique(.)
 

write.csv(myelin_genes, "manuscript/brain/gene_sets/myelin.csv", row.names =F)
