#Making venndiagram using Differentially expressed genes within QTL, expressed in PGCs, expressed in adult ovaries, genes with In-phase SNPs (synonymous and non-synonymous)
library(VennDiagram)
library(dplyr)

setwd("~/Documents/Paper/fig5Analysis")

#uploading list of differentially expressed genes within the QTLs (QTL-21d and QTL-3d) 
QTLdiffExp <- read.csv("QTL_diffGenes.csv", sep = ",", header = T)
#Selecting QTL-21d (QTL from 21 day old datasets)
QTL21_diff <- QTLdiffExp %>% filter(arm == "QTL-21d")
QTL21_diff <- QTL21_diff[,c(1,2)]
colnames(QTL21_diff)<- c("FBgn", "gene_symbol")
#Selecting QTL-3d (QTL from 3 day old datasets)
QTL3_diff <- QTLdiffExp %>% filter(arm == "QTL-3d")
QTL3_diff <- QTL3_diff[,c(1,2)]
colnames(QTL3_diff)<- c("FBgn", "gene_symbol")

#Extracting genes that are expressed in primordial germ cells, PGCs (publicly available from single cell RNA-seq data)
PGC_exp <- read.csv("pgc-expressed.txt", sep = ",", header = T)

setwd("~/Documents/RNAseq/Count_file/TE-count/unaligned-gene-counts/kallisto/abundance.tsv") 
#This is dm6 annotation.
B <- read.csv("coordi_Gene2.csv", header = T)
colnames(PGC_exp) <- "gene_symbol"

PGC_exp_fbgn <- merge(PGC_exp, B, by="gene_symbol")
extra <- anti_join(PGC_exp, PGC_exp_fbgn)

setwd("~/Documents/Paper/fig5Analysis")
#write.csv(extra, file="extra_PGC_nonmatched.csv")
extra_fbgn <- read.csv("extra_PGC_flybaseFbgn.csv", header = T)
extra_fbgn <- extra_fbgn[,c(1,2)]
PGC_exp_fbgn <- PGC_exp_fbgn[,c(1,5)]
PGC_exp_fbgn <- rbind(extra_fbgn,PGC_exp_fbgn)

#Extracting genes that are expressed in adult ovaries (source: flybase)
ovary_exp <- read.csv("ovary-germline-expressed.txt", sep = ",", header = T)
colnames(ovary_exp) <- "gene_symbol"
ovary_exp_fbgn <- merge(ovary_exp, B, by="gene_symbol")
extra <- anti_join(ovary_exp, ovary_exp_fbgn)

setwd("~/Documents/Paper/fig5Analysis")
#write.csv(extra, file="extra_ovary_nonmatched.csv")
extra_fbgn <- read.csv("extra_ovary_flybaseFbgn.csv", header = T)
ovary_exp_fbgn <- ovary_exp_fbgn[,c(1,5)]
ovary_exp_fbgn <- rbind(extra_fbgn,ovary_exp_fbgn)
PGC_ovary_exp <- rbind(PGC_exp_fbgn,ovary_exp_fbgn) #11536
PGC_ovary_exp <- PGC_ovary_exp[,2] 
PGC_ovary_exp <- unique(PGC_ovary_exp) %>% data.frame() #11529
#merge <-  merge(ovary_exp_fbgn, PGC_exp_fbgn, by = "FBgn") #9796 (9796+1129+597= 11522)

#Extracting all the genes with non-synonymouse Inphase SNPs
non_synonym <- read.csv("non-synonymouseSnps.csv", sep = ",", header = T)
non_synonym_fbgn <- merge(non_synonym, B, by="gene_symbol")
non_synonym_fbgn <- non_synonym_fbgn[,c(1,2,9)]
extra <- anti_join(non_synonym, non_synonym_fbgn) #null
QTL21_non_synom <- non_synonym_fbgn %>% filter(QTL == "QTL-21d") #4 genes
QTL3_non_synom <- non_synonym_fbgn %>% filter(QTL == "QTL-3d") #38 genes

#Now moving on to the in-phase SNPs that includes non-synonymous ones for QTL21-d and QTL3-d
setwd("~/Documents/Paper/Updated_2L2RQTL/")
QTL21_synonym <- read.csv("QTL1_allpredictedUniqueGenes.csv", sep = ",", header = T)
QTL21_synonym <- QTL21_synonym[,2]  %>% data.frame()
colnames(QTL21_synonym) <- "gene_symbol" #42 in-phase SNPs
QTL21_synonym_fbgn <- merge(QTL21_synonym, B, by="gene_symbol")
extra <- anti_join(QTL21_synonym, QTL21_synonym_fbgn)

setwd("~/Documents/Paper/fig5Analysis")
extra_fbgn <- read.csv("Extra_QTL21d_in-phaseSNP_gene.csv", sep = ",", header = T)
colnames(extra_fbgn) <- c("gene_symbol", "FBgn")
QTL21_synonym_fbgn <- QTL21_synonym_fbgn[, c(1,5)]
QTL21_synonym_fbgn <- rbind(QTL21_synonym_fbgn,extra_fbgn)

setwd("~/Documents/Paper/Updated_2L2RQTL/")
QTL3_2R_synonym <- read.csv("QTL2_2R_predictedUniqueGenes.csv", sep = ",", header = T)
QTL3_2L_synonym <- read.csv("QTL2_allpredictedUniqueGenes.csv", sep = ",", header = T)
QTL3_2L_synonym <- QTL3_2L_synonym[,2] %>% data.frame()
QTL3_2R_synonym <-  QTL3_2R_synonym[,2] %>% data.frame()
QTL3_synonym <- rbind(QTL3_2R_synonym,QTL3_2L_synonym) #258 in-phase SNPs
colnames(QTL3_synonym) <- "gene_symbol" 
setwd("~/Documents/Paper/fig5Analysis")
#write.csv(QTL3_synonym, "QTL3_in-phaseSNPs")

#I manually added Fbgn Id to the gene_symbols from flybase.
QTL3_synonym <- read.csv("QTL3_InphaseSnp_fbgn.csv", sep = ",", header = T)


#Lets start with QTL-21d
A <- as.character(QTL21_diff[,1])  # 30
B <- as.character(PGC_ovary_exp[,1])  #11529
C <- as.character(QTL21_non_synom[,3])   # 38 genes
D <- as.character(QTL21_synonym_fbgn[,2])    # 258 genes

#making a list
x <- list(A,B,D,C)
#making the venn diagram
library(ggVennDiagram)
ggVennDiagram(x,c("Diff_exp","PGC_ovary","in-phaseSNP","non-synonym"), label = "count", lwd=0.8, color= "black", lty =1) +  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")


#Lets look at most important interactions 
A <- (QTL21_diff[,1]) %>% data.frame() # 30
B <- (PGC_ovary_exp[,1])  %>% data.frame() #11529
C <- (QTL21_non_synom[,3])  %>% data.frame()  # 38 genes
D <- (QTL21_synonym_fbgn[,2]) %>% data.frame() # 258 genes

colnames(A) <- "FBgn" #Differentially expressed in QTL-21d
colnames(B) <- "FBgn" #Genes expressed in PGC and ovary germline (from Kelleher)
colnames(C) <- "FBgn" #Non-synonymous in-phase Snp
colnames(D) <- "FBgn" #In-phase SNP

a_b <- merge(A,B, by= "FBgn") #diff_expressed + PGC_ovary expressed = 13 GENES
a_d <- merge(A,D, by= "FBgn") #diff expressed + in-phase SNP = 5 GENES
a_b_d <- merge(a_b,D, by= "FBgn") #diff_expressed + PGC_ovary expressed + in-phase SNP = 5 Genes
b_c <- merge(B,C, by= "FBgn") #non-synonymouse SNPs + PGC_ovary expressed = 3 genes

setwd("~/Documents/RNAseq/Count_file/TE-count/unaligned-gene-counts/kallisto/abundance.tsv") 
#This is dm6 annotation.
Cordi <- read.csv("coordi_Gene2.csv", header = T)
Cordi <- Cordi[,c(4,5)]

a_b<- merge(Cordi,a_b, by="FBgn")
a_d<- merge(Cordi,a_d, by="FBgn")
a_b_d<- merge(Cordi,a_b_d, by="FBgn")
b_c<- merge(Cordi,b_c, by="FBgn")

a_b <- a_b[,2] %>% data.frame()
a_d <- a_d[,2] %>% data.frame()
a_b_d <- a_b_d[,2] %>% data.frame()
b_c <-b_c[,2] %>% data.frame()


colnames(a_b) <- "DiffExp+PGC/ovaryExp" #diff_expressed + PGC_ovary expressed = 29 GENES
colnames(a_d) <- "DiffExp+in-phaseSNP" #diff expressed + in-phase SNP = 12 GENES
colnames(a_b_d) <- "DiffExp+PGC/ovaryExp+in-phaseSNP" #diff_expressed + PGC_ovary expressed + in-phase SNP = 12 Genes
colnames(b_c) <- "Non-synonymouseSNP+PGC"

setwd("~/Documents/Paper/fig5Analysis")

#write.csv(a_b, "QTL-21_DEG_PGC")
#write.csv(a_d, "QTL-21_DEG_in-phase")
#write.csv(a_b_d, "QTL-21_DEG_PGC_in_phase")
#write.csv(b_c, "QTL-21_Non-synonym_PGC")



##################################################################
#Lets start with QTL-3d
A <- as.character(QTL3_diff[,1])  # 30
B <- as.character(PGC_ovary_exp[,1])  #11529
C <- as.character(QTL3_non_synom[,3])   # 38 genes
D <- as.character(QTL3_synonym[,2])    # 258 genes

#making a list
x <- list(A,B,D,C)

library(ggVennDiagram)
ggVennDiagram(x,c("Diff_exp","PGC_ovary","in-phaseSNP","non-synonym"), label = "count", lwd=0.8, color= "black", lty =1) +  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")

# find overlap between all possible combinations (11 plus 4 unique to each list = 15 total)
#overlap <- calculate.overlap(x)
# get the lengths of each list - these will be the numbers used for the Venn diagram
#overlap_values <- lapply(overlap, function(x) length(x))

#venn.diagram(x, category.names = c("Diff_exp","PGC_ovary","non-synonym","in-phaseSNP"), filename = 'QTL3d.png', output=TRUE, imagetype = "png", scaled =FALSE )
A <- QTL3_diff[,1]  %>%  data.frame()# 30
B <- PGC_ovary_exp %>%  data.frame() #11529
C <- QTL3_non_synom [,3] %>%  data.frame() # 38 genes
D <- QTL3_synonym  [,2] %>%  data.frame() # 258 genes

colnames(A) <- "FBgn" #Differentially expressed in QTL-21d
colnames(B) <- "FBgn" #Genes expressed in PGC and ovary germline (from Kelleher)
colnames(C) <- "FBgn" #non-synonymous in-phase SNPs
colnames(D) <- "FBgn" # all in-phase SNPs

#Lets look at most important interactions.
A_B <- merge(A,B, by= "FBgn") #diff_expressed + PGC_ovary expressed = 29 GENES
A_D <- merge(A,D, by= "FBgn") #diff expressed + in-phase SNP = 12 GENES
A_B_D <- merge(A_B,D, by= "FBgn") #diff_expressed + PGC_ovary expressed + in-phase SNP = 12 Genes
B_C <- merge(B,C, by= "FBgn") #non-synonymouse SNPs + PGC_ovary expressed = 33 genes

setwd("~/Documents/RNAseq/Count_file/TE-count/unaligned-gene-counts/kallisto/abundance.tsv") 
#This is dm6 annotation.
Cordi <- read.csv("coordi_Gene2.csv", header = T)
Cordi <- Cordi[,c(4,5)]
A_B<- merge(Cordi,A_B, by="FBgn")
A_D<- merge(Cordi,A_D, by="FBgn")
A_B_D<- merge(Cordi,A_B_D, by="FBgn")
B_C<- merge(Cordi,B_C, by="FBgn")

A_B <- A_B[,2] %>% data.frame()
A_D <- A_D[,2] %>% data.frame()
A_B_D <- A_B_D[,2] %>% data.frame()
B_C <-B_C[,2] %>% data.frame()

colnames(A_B) <- "DiffExp+PGC/ovaryExp" #diff_expressed + PGC_ovary expressed = 29 GENES
colnames(A_D) <- "DiffExp+in-phaseSNP" #diff expressed + in-phase SNP = 12 GENES
colnames(A_B_D) <- "DiffExp+PGC/ovaryExp+in-phaseSNP" #diff_expressed + PGC_ovary expressed + in-phase SNP = 12 Genes
colnames(B_C) <- "Non-synonymouseSNP+PGC"

setwd("~/Documents/Paper/fig5Analysis")

#write.csv(a_b, "QTL-3_DEG_PGC")
#write.csv(a_d, "QTL-3_DEG_in-phase")
#write.csv(a_b_d, "QTL-3_DEG_PGC_in_phase")
#write.csv(b_c, "QTL-3_Non-synonym_PGC")

