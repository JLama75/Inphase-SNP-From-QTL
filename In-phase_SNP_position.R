#Release 5 drosophgila genome dm3
# Table of SNP positions for each QTL peaks 
#install.packages('vcfR') 

#Analysis for QTL1
library(vcfR)
library(dplyr)
QTL1 <- read.vcfR("DSPRsnps_(chr2L_19,010,000-20,000,000)].vcf")
head(QTL1)

polymorphic <- (is.polymorphic(QTL1, na.omit = TRUE))
summary (polymorphic)
biallelic <- (is.biallelic(QTL1))
summary(biallelic)

chrom <- create.chromR(name='RAD_data', vcf = QTL1)
plot(chrom) # plot the data 

#quick check genotype distribution per strain
dp1 <- extract.gt(QTL1, element='GT', as.numeric=TRUE) %>% data.frame (dp1) #extracting the Genotype of each SNP of all strains 

library(data.table)
dp2 <- setDT(dp1, keep.rownames = TRUE)[] #converting the rowname to the 1st column

#selecting the B founder columns. excluding B5
dp3 <- dp1[ ,c(1, 9:12, 14, 18:19) ]

summary(dp3) #check for NA values
dp4 <- dp3[complete.cases(dp3),] #removing NA from the data
summary (dp4)

#converting the data table dp4 into data frame.
dp4 <- as.data.frame(dp4)

#Screening for the inphase SNPs using ifelese criteria 2. B6 and B2  not equal to any of the B founders.
dp4[ ,"inphase2"] <- NA 
#putting the inphase conditions B6=B2 and B7
dp4$inphase2 <- ifelse ((dp4$B6 != dp4$B1 & dp4$B6 == dp4$B2 & dp4$B6 != dp4$B3 & dp4$B6 != dp4$B4 & dp4$B6 != dp4$B7 & dp4$B6 != dp4$AB8), "Pass", "Fail")

#selecting only thoes rows that pass the set criteria into a different file. 
Inphase_QTL1b <- subset(dp4, grepl("Pass", dp4$inphase2))
Inphase_QTL1b[, "POS"] <- NA
#sepearting the position from Chr. 
Inphase_QTL1b$POS <-  substr(Inphase_QTL1b$rn, 7, 14)  
#write.table(Inphase_QTL1b, file = "Inphase_SNP_QTL1.2.txt", sep = "\t") 

#extracting the information for each SNP position 
Test <- vcfR2tidy(QTL1)
Test1 <- Test$fix
View(Test1)
Info <- Test1[, c("POS","REF","ALT","EFF")]
lapply(Info, class)
Info$POS <- as.character(Info$POS)
lapply(Inphase_QTL1b, class)

POS <- Inphase_QTL1b$rn
library(dplyr)
merge <- full_join(Info, Inphase_QTL1b, by ="POS")
summary(merge)
QTL1_SNP <- merge[complete.cases(merge), ]
write.table(QTL2_SNP, file = "Inphase_SNP_QTL2.a.txt", sep = "\t")

QTL1_SNP$Chrom <- "2L"
QTL1_SNP <- QTL1_SNP[,c(1,14,2,3,4,7:12,6)]
write.table(QTL1_SNP, file = "Inphase_SNP_QTL1.txt", sep="\t")
######################################################################################################
######################################################################################################

#QTL2

#install.packages('vcfR')
library("vcfR")
setwd("~/Documents/Paper/QTL_merge SNPs/")
QTL2 <- read.vcfR("QTL2(chr2L_20,550,000-23,011,544)].vcf")
head(QTL2)

setwd("~/Documents/Paper/Updated_2L2RQTL/")
#QTL2.1 <- QTL2@fix[ , ]
#QTL2@fix[ ,1:7]
polymorphic <- (is.polymorphic(QTL2, na.omit = TRUE))
summary (polymorphic)
biallelic <- (is.biallelic(QTL2))
summary(biallelic)
QTL2 <- QTL2[is.biallelic(QTL2), ] #only selecting the biallelic ones.
Biallelic <- is.biallelic(QTL2)
summary(Biallelic)
chrom <- create.chromR(name='RAD_data', vcf = QTL2)
plot(chrom) # plot the data 

#quick check genotype distribution per strain
dp1 <- extract.gt(QTL2, element='GT', as.numeric=TRUE) #etracting the Genotype of each SNP of all strains
dp1 <- data.frame (dp1)
ncol(dp1)
colnames(dp1)
rownames(dp1)
library(data.table)
dp2 <- setDT(dp1, keep.rownames = TRUE)[] #converting the rowname to the 1st column
ncol(dp2)
dp3 <- dp1[ ,c(1, 9:12, 14, 18:19) ]
summary(dp3) #check last row for NA values
dp4 <- dp3[complete.cases(dp3),] #removing NA from the data
summary (dp4)
#converting the data table dp4 into data frame.
dp4 <- as.data.frame(dp4)
#adding a new column for inphase SNPs
dp4[ ,"inphase"] <- NA 

#Screening for the inphase SNPs using ifelese Criteria 1. B6 not equal to any of the B founders.
dp4$inphase <- ifelse ((dp4$B6 != dp4$B1 & dp4$B6 != dp4$B2 & dp4$B6 != dp4$B3 & dp4$B6 != dp4$B4 & dp4$B6 != dp4$B7 & dp4$B6 != dp4$AB8), "Pass", "Fail")
#selecting only thoes rows that pass the set criteria into a different file. 
Inphase_QTL2a <- subset(dp4, grepl("Pass", dp4$inphase))
Inphase_QTL2a[, "POS"] <- NA
Inphase_QTL2a$POS <-  substr(Inphase_QTL2a$rn, 7, 14)  
#grep <- subset(Inphase_QTL2a, grep(2,Inphase_QTL2a))
#write.table(Inphase_QTL2a, file = "Inphase_SNP_QTL2.txt", sep = "\t") 


#extracting the information for each SNP position 
Test <- vcfR2tidy(QTL2)
Test1 <- Test$fix
View(Test1)
Info <- Test1[, c("POS","REF","ALT","EFF")]
lapply(Info, class)
Info$POS <- as.character(Info$POS)
lapply(Inphase_QTL2a, class)

POS <- Inphase_QTL2a$rn
library(dplyr)
merge <- full_join(Info, Inphase_QTL2a, by ="POS")
summary(merge)
QTL2_SNP <- merge[complete.cases(merge), ]

QTL2_SNP <- subset(QTL2_SNP, POS >= 20820000)
setwd("~/Documents/Paper/Supplementals")
QTL2_SNP$Chrom <- "2L"
QTL2_SNP <- QTL2_SNP[,c(1,14,2,3,4,7:12,6)]

write.table(QTL2_SNP, file = "Inphase_SNP_QTL2.txt", sep = "\t")

######################################################################################################
######################################################################################################

#QTL2-2R
setwd("~/Documents/Paper/Updated_2L2RQTL/")
QTL2_3_2R <- read.vcfR("Galaxy23-[UCSC_Main_on_D._melanogaster__hub_55499_DSPRsnps_(chr2R_1-6,942,495)].vcf") # note that these are dm3 positions of QTL
head(QTL2_3_2R)

#but the positions are in dm6. converting to dm3-- the QTL position becomes 1-2830000. So we have to subset thoes positions
#QTL2_3_2R.1 <- QTL2_3_2R@fix[ , ]
#QTL2_3_2R@fix[ ,1:7]
polymorphic <- (is.polymorphic(QTL2_3_2R, na.omit = TRUE))
summary (polymorphic)
biallelic <- (is.biallelic(QTL2_3_2R))
summary(biallelic)
QTL2_3_2R <- QTL2_3_2R[is.biallelic(QTL2_3_2R), ] #only selecting the biallelic ones.
Biallelic <- is.biallelic(QTL2_3_2R)
summary(Biallelic)

chrom <- create.chromR(name='RAD_data', vcf = QTL2_3_2R)
plot(chrom) # plot the data 

#quick check read depth distribution per individual
dp <- extract.gt(QTL2_3_2R, element='DP', as.numeric=TRUE) 
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)",
        las=2, cex=0.5, cex.axis=1)
abline(h=8, col="red")
dev.off() 

#quick check genotype distribution per strain
dp1 <- extract.gt(QTL2_3_2R, element='GT', as.numeric=TRUE) #etracting the Genotype of each SNP of all strains
dp1 <- data.frame (dp1)


ncol(dp1)
colnames(dp1)
rownames(dp1)
library(data.table)
dp2 <- setDT(dp1, keep.rownames = TRUE)[] #converting the rowname to the 1st column
ncol(dp2)

dp3 <- dp1[ ,c(1, 9:12, 14, 18:19) ]
summary(dp3) #check last row for NA values
dp4 <- dp3[complete.cases(dp3),] #removing NA from the data
summary (dp4)
#converting the data table dp4 into data frame.
dp4 <- as.data.frame(dp4)
#adding a new column for inphase SNPs
dp4[ ,"inphase"] <- NA 

dp3 <- dp4 #making a copy
#Screening for the inphase SNPs using ifelese Criteria 1. B6 not equal to any of the B founders.
dp4$inphase <- ifelse ((dp4$B6 != dp4$B1 & dp4$B6 != dp4$B2 & dp4$B6 != dp4$B3 & dp4$B6 != dp4$B4 & dp4$B6 != dp4$B7 & dp4$B6 != dp4$AB8), "Pass", "Fail")
#selecting only thoes rows that pass the set criteria into a different file. 
Inphase_QTL2_3_2Ra <- subset(dp4, grepl("Pass", dp4$inphase))
Inphase_QTL2_3_2Ra[, "POS"] <- NA

Inphase_QTL2_3_2Ra$POS <-  substr(Inphase_QTL2_3_2Ra$rn, 7, 13)  #605

#extracting the information for each SNP position 
Test <- vcfR2tidy(QTL2_3_2R)
Test1 <- Test$fix
View(Test1)
Info <- Test1[, c("POS","REF","ALT","EFF")]
Info$POS <- as.character(Info$POS)
lapply(Inphase_QTL2_3_2Ra, class)

POS <- Inphase_QTL2_3_2Ra$rn
library(dplyr)
merge <- full_join(Info, Inphase_QTL2_3_2Ra, by ="POS")
summary(merge)
QTL2_3_2R_SNP <- merge[complete.cases(merge), ]

#subseting the QTL dm3 positions as the inphase Data is according to R5/dm3
QTL2_3_2R_SNP <- subset(QTL2_3_2R_SNP, POS <= 2830000)

setwd("~/Documents/Paper/Supplementals")
QTL2_3_2R_SNP$Chrom <- "2R"
QTL2_3_2R_SNP <- QTL2_3_2R_SNP[,c(1,14,2,3,4,7:12,6)]

write.table(QTL2_3_2R_SNP, file = "Inphase_SNP_QTL2_2R.txt", sep = "\t")

#write.table(QTL2_3_2R_SNP, file = "Inphase_SNP_QTL2_2R.txt", sep = "\t")

cluster42AB_distal <- subset(QTL2_3_2R_SNP, POS >= 2336719 & POS <=2386719) #50 kb
cluster42AB_distal <- subset(QTL2_3_2R_SNP, POS >= 2380000 & POS <=2383000) #near the promoter region

cluster42AB <- subset(QTL2_3_2R_SNP, POS >= 2144200 & POS <=2386719) #11
write.csv(cluster42AB, "cluster42AB_In-phaseSNP.csv")
##########################################################

#QTL3 both narrowed peak

QTL3_both <- subset(QTL2_3_2R_SNP, POS >= 2030000 & POS <=2310000)
setwd("~/Documents/Paper/Supplementals")

write.table(QTL3_both, file = "Inphase_SNP_QTL3_both.txt", sep = "\t")

######################################################################################################################
######################################################################################################################

