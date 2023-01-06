#Larval and Ovary Expressing Inphase SNPs

library(reshape2)
library(dplyr)
library(data.table)

#opening the file where the SNP position effecting each genes are seperated in a column.
genes <- read.csv("QTL2_genes.csv", header=T, na.strings = c("", "NA"))
genes_copy <- cbind(genes) #copying dataframe 

#replacing duplicated genes within each row by NA
genes[which(t(apply(genes,1,function(genes) duplicated(genes))), arr.ind = T)] <- NA

#To compare the genes in every position, first melt the genes in different rows and sort it on the basis of position in ascending order
A<- melt(genes, id = c("CHROM", "POS"))
sorted1 <- A[order(A$POS), ]
sorted1$duplicate <- duplicated(sorted1$POS) #checking the duplicated postions to check whether sorting has occured
summary(sorted1)
A1<- sorted1[!is.na(sorted1$value), ] #removing the NA values.
A2 <- setDT(A1)[,.(value=paste(value,collapse=",")),by="POS"] #all the genes in every postitions with NA values removed.

#Comparing with the larval expressed genes (L3) by merging them by the column with gene names
larval_genes <- read.delim("larval_Genes.txt", header=F, col.names = c("value"))
merged_genes <- merge(A, larval_genes, by = "value") 
sorted <- merged_genes[order(merged_genes$POS), ] 
sorted$duplicate <- duplicated(sorted$POS) #checking if sorting has occured by looking at duplicated rows
sorted2 <- setDT(sorted)[,.(value=paste(value,collapse=",")),by="POS"] #merging the row with same Position into one with the genes seperated by comma
sorted2$duplicate <- duplicated(sorted2$POS)
summary(sorted2) #checking number of FALSE and TRUE 
Larval_sortGene <- select(sorted2, POS, value) # selecting columns without the duplication column


genes$larval_expressed <- genes$POS %in% sorted2$POS #gives a TRUE or FALSE value for matched/unmatched POS
summary(genes) #1414 larval expressed matched genes, 1414 TRUE values in larval_expressed column
Larval_SNPs <- select(genes, CHROM, POS, larval_expressed)

#Comparing with the adult ovary expressed genes by merging them by the column with gene names
Ovary_genes <- read.delim("Ovary_Genes.txt", header=F, col.names = c("value"))
merged_genes2 <- merge(A, Ovary_genes, by = "value") 
Ovary_sorted <- merged_genes2[order(merged_genes2$POS), ] 
Ovary_sorted$duplicate <- duplicated(Ovary_sorted$POS) #checking if sorting has occured looking for duplicated rows

Ovary_sorted2 <- setDT(Ovary_sorted)[,.(value=paste(value,collapse=",")),by="POS"]
Ovary_sorted2$duplicate <- duplicated(Ovary_sorted2$POS)
summary(Ovary_sorted2) #checking number of FALSE and TRUE 
Ovary_sort <- select(Ovary_sorted2, POS, value)


genes$Ovary_expressed <- genes$POS %in% Ovary_sorted2$POS #gives a TRUE or FALSE value for matched/unmatched POS
summary(genes) #4102 ovary expressed matched genes, 4102 TRUE values in larval_expressed column
Ovary_SNPs <- select(genes, CHROM, POS, Ovary_expressed)

#merging all SNP Positions with only the inphase SNP postions for QTL2 2L criteria a
Inphase_QTL2a <- read.delim("Inphase_SNP_QTL2.a.txt", header=T)
Inphase_QTL2a_larval <- merge(Larval_SNPs, Inphase_QTL2a, by= 'POS')
Inphase_QTL2a_ovaryLarval <- merge(Ovary_SNPs, Inphase_QTL2a_larval, by= 'POS')

Inphase_QTL2a_ovaryLarval1 <- merge(Ovary_sort, Inphase_QTL2a_ovaryLarval, by ='POS', all.y = TRUE) #also merging the ovary expressed gene name without deleting the extra SNP positions not expressing it.
Inphase_QTL2a_ovaryLarval2 <- merge(Larval_sortGene, Inphase_QTL2a_ovaryLarval1, by ='POS', all.y = TRUE)
write.table(Inphase_QTL2a_ovaryLarval2, file = "Inphase_QTL2a_ovaryLarval.txt", sep = "\t")

#merging all SNP Positions with only the inphase SNP postions for QTL2 2L criteria b
Inphase_QTL2b <- read.delim("Inphase_SNP_QTL2.b.txt", header=T)
Inphase_QTL2b_larval <- merge(Larval_SNPs, Inphase_QTL2b, by= 'POS')
Inphase_QTL2b_ovaryLarval <- merge(Ovary_SNPs, Inphase_QTL2b_larval, by= 'POS')
Inphase_QTL2b_ovaryLarval1 <- merge(Ovary_sort, Inphase_QTL2b_ovaryLarval, by ='POS', all.y = TRUE) #also merging the ovary expressed gene name without deleting the extra SNP positions not expressing it.
Inphase_QTL2b_ovaryLarval2 <- merge(Larval_sortGene, Inphase_QTL2b_ovaryLarval1, by ='POS', all.y = TRUE)
write.table(Inphase_QTL2b_ovaryLarval2, file = "Inphase_QTL2b_ovaryLarval.txt", sep = "\t")

#counting the total number of Genes affected by inphase SNPs 
library(splitstackshape)
OvaryGene_number <- select(Inphase_QTL2a_ovaryLarval2, POS, value.y)
OvaryGene_number1 <- cSplit(OvaryGene_number, "value.y", ",", "long")
OvaryGene_number2 <- select(OvaryGene_number1, value.y)
OvaryGene_number3 <- OvaryGene_number2[complete.cases(OvaryGene_number2), ]
OvaryGene_number4 <- OvaryGene_number3[!duplicated(OvaryGene_number3),]

LarvalGene_number <- select(Inphase_QTL2a_ovaryLarval2, POS, value.x)
LarvalGene_number1 <- cSplit(LarvalGene_number, "value.x", ",", "long")
LarvalGene_number2 <- select(LarvalGene_number1, value.x)
LarvalGene_number3 <- LarvalGene_number2[complete.cases(LarvalGene_number2), ]
LarvalGene_number4 <- LarvalGene_number3[!duplicated(LarvalGene_number3),]

#counting the duplicated Genes in all In-phase SNP positions.
A3 <- merge(A2, Inphase_QTL2a, by= 'POS')
Gene_number <- select(A3, POS, value)
Gene_number1 <- cSplit(Gene_number, "value", ",", "long")
#Gene_number1$duplicate <- duplicated(Gene_number1$POS)
Gene_number2 <- select(Gene_number1, value)
Gene_number3 <- Gene_number2[complete.cases(Gene_number2), ]
Gene_number4 <- Gene_number3[!duplicated(Gene_number3),]
#161 total genes

#counting the duplicated Genes in all In-phase SNP positions.
B3 <- merge(A2, Inphase_QTL2b, by= 'POS')
B3_Gene_number <- select(B3, POS, value)
B3_Gene_number1 <- cSplit(B3_Gene_number, "value", ",", "long")
B3_Gene_number2 <- select(B3_Gene_number1, value)
B3_Gene_number3 <- B3_Gene_number2[complete.cases(B3_Gene_number2), ]
B3_Gene_number4 <- B3_Gene_number3[!duplicated(B3_Gene_number3),]

OvaryGene_number <- select(Inphase_QTL2b_ovaryLarval2, POS, value.y)
OvaryGene_number1 <- cSplit(OvaryGene_number, "value.y", ",", "long")
OvaryGene_number2 <- select(OvaryGene_number1, value.y)
OvaryGene_number3 <- OvaryGene_number2[complete.cases(OvaryGene_number2), ]
OvaryGene_number4 <- OvaryGene_number3[!duplicated(OvaryGene_number3),]'

#According to Gene'
#Larvae 
#Gene_number1-- all positions with SNPs common with Inphase SNPs
#sorted -- all positions with genes common with larval genes
Gene_number1$larvae <- Gene_number1$value  %in% sorted$value #common gene values assigned TRUE, rest given FALSE
#Ovary
Gene_number1$ovary <- Gene_number1$value  %in% Ovary_sorted$value 

Filter_gene <- Gene_number1[order(Gene_number1$value), ] #sorting based on gene name 
Filter_gene1 <- data.frame( table(Filter_gene$value)) #calculating frequency of repeated genes for No. of SNP positions for each gene
colnames(Filter_gene1) <- c('value', '#ofSNPs') #changing the colnames.
Filter_gene2 <- Filter_gene[!duplicated(Filter_gene$value), ]
Filter_gene3 <- merge(Filter_gene2, Filter_gene1, by='value')
CHROM <- rep("Chr2L",121)
CHROM <- data.frame(x)
QTL_peak <- rep("QTL2", 121)
QTL_peak <- data.frame(QTL_peak) 
Filter_Gene4 <- cbind(Filter_gene3, CHROM, QTL_peak)
write.table(Filter_Gene4, file = "QTL2_2L_Genefrequency.txt", sep= '\t')

#According to Gene'CRITERIA --b
#Larvae
B3_Gene_number1$larvae <- B3_Gene_number1$value  %in% sorted$value #common gene values assigned TRUE, rest given FALSE
#Ovary
B3_Gene_number1$ovary <- B3_Gene_number1$value  %in% Ovary_sorted$value 

B3_Filter_gene <- B3_Gene_number1[order(B3_Gene_number1$value), ] #sorting based on gene name 
B3_Filter_gene1 <- data.frame( table(B3_Filter_gene$value)) #calculating frequency of repeated genes for No. of SNP positions for each gene
colnames(Filter_gene1) <- c('value', '#ofSNPs') #changing the colnames.
B3_Filter_gene2 <- B3_Filter_gene[!duplicated(B3_Filter_gene$value), ]
B3_Filter_gene3 <- merge(B3_Filter_gene2, B3_Filter_gene1, by='value')
CHROM <- rep("Chr2L",121)
CHROM <- data.frame(x)
QTL_peak <- rep("QTL2", 121)
QTL_peak <- data.frame(QTL_peak) 
B3_Filter_Gene4 <- cbind(B3_Filter_gene3, CHROM, QTL_peak)
write.table(B3_Filter_Gene4, file = "QTL2_2L_Genefrequency.txt", sep= '\t')
