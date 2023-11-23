setwd("~/Documents/HMM_R4_QTL/QTL2")

#Aim: To look at expression of pericentromeric genes among the 8 founder alleles using microarray gene expression data from head tissue of Drosophila melanogaster

library("dplyr")
library("rlist")
library("ggplot2")
library("tidyverse")

#Following is the csv file with 864 unique RILs (recombinant inbred line) and information on the founder allele in 347 genomic positions
#HMM release 4 (used dm3 co-ordinates of QTL2) : greater than 20550000(2L) and less than 1550000 (2R)
HMM1<- read.csv("p95_HMM_2L_founderAssigned_sorted_dm3.csv")
HMM2 <- read.csv("p95_HMM_2R_founderAssigned_sorted_dm3.csv")
HMM1 <- subset(HMM1, select= -X)
HMM2 <- subset(HMM2, select= -X)

#merge 2l and 2R hmm data
HMM <- rbind(HMM1,HMM2)
HMM<- HMM[order(HMM$RIL),]
rownames(HMM) <- seq(1:299808) 
#347 --number of positions

#creating a list of RILs
list1 <- unique(HMM$RIL) 
length(list1) #864 RILs
HMM$continguent <- NA

#Only selecting thoes RILs that have single founder allele in all 347 positions

for(i in seq(1,(dim(HMM)[1])-346, 347 ))
{ #print(i)
  n = i
  loop <- HMM[i:(i+346),] 
  if (length(unique(loop$RIL)) == 1)
  {
    if (length(unique(loop$founder)) == 1)
    {
      HMM[i:(i+346),4] <- "True"
    }
    else {HMM[i:(i+346),4] <- "False"}
  }
  else message("Error in ", n)
  
}

#How many RILs not continguent across QTL2?
(dim(subset(HMM, continguent == "True")) [1])/347 #736 RILs out of 864 RILs
False <- subset(HMM, continguent == "False")
True <- subset(HMM, continguent == "True")

#subsetting all the continguent genotypes and removing NA values
HMM3 <- subset(HMM, HMM$continguent == "True")
list2 <- unique(HMM3$RIL)
length(list1) #864 RILs
length(list2) #736 RILs

HMM4 <- HMM3[!is.na(HMM3$founder),]
list3 <- unique(HMM4$RIL)
length(list3) #719 RILs

#now looking at the founder alleles at QTL2 locus
final_list <- unique(subset(HMM4, select = c(RIL, founder)))
final_list <- final_list[order(final_list$founder),]

#subsetting individual founders
BB1 <- subset(final_list, final_list$founder == "BB1", select = RIL)
colnames(BB1) <- "BB1"
BB2 <- subset(final_list, final_list$founder == "BB2", select = RIL)
colnames(BB2) <- "BB2"
BB3 <- subset(final_list, final_list$founder == "BB3", select = RIL)
colnames(BB3) <- "BB3"
BB4 <- subset(final_list, final_list$founder == "BB4", select = RIL)
colnames(BB4) <- "BB4"
BB5 <- subset(final_list, final_list$founder == "BB5", select = RIL)
colnames(BB5) <- "BB5"
BB6 <- subset(final_list, final_list$founder == "BB6", select = RIL)
colnames(BB6) <- "BB6"
BB7 <- subset(final_list, final_list$founder == "BB7", select = RIL)
colnames(BB7) <- "BB7"
BB8 <- subset(final_list, final_list$founder == "BB8", select = RIL)
colnames(BB8) <- "BB8"

#extracting the gene expression of head microarray data of the founders in the list 
setwd("~/Documents/RNAseq/DSPR_HMM/HMMregB_R2")
RIL<- read.csv("FemaleHeadExpression.csv", sep = ",", header = T)
setwd("~/Documents/HMM_R4_QTL/QTL2")

#Now merging RIL and final list by paternal RIL to find the common RILs among them.
colnames(final_list) <- c("patRIL", "founder")
RIL2 <- merge(final_list, RIL, by = "patRIL") #out of 719 we have 501 RILs
RIL2 <- RIL2[order(RIL2$founder),]

#looking at the genes in head expression data.
list_genes <- colnames(RIL2) %>% data.frame
colnames(list_genes) <- "genes"
list_genes <- list_genes[order(list_genes$genes),] %>% data.frame()
list_genes <- data.frame(list_genes)
#removing patRIL, matRIL and founders from the list of genes
list_genes <- list_genes[-c(11065:11067),] 
list_genes <- data.frame(list_genes)

#the list has splice variants too. therefore, replacing the name of splice variant by common gene name
list_alleles <- list_genes
list_alleles$gene <- list_alleles$list_genes
list_alleles$list_genes <- gsub("(\\w+\\d+)\\.(\\w+)", "\\1",list_genes$list_genes)
more <- c("patRIL","founder","matRIL") %>% data.frame
colnames(more) <- "list_genes"
more$gene <- c("patRIL","founder","matRIL")
list_alleles <- rbind(list_alleles,more)
list_alleles <- list_alleles[,c(2,1)] 
rownames(list_alleles)
rownames(list_alleles) <- list_alleles[,1]

list_alleles <- subset(list_alleles, select = -gene)
#write.csv(list_alleles, "listalleles.csv")
list_all <- data.frame(t(list_alleles))

#adding the new gene names (replaced) into a new row in the expression data RIL2
RIL3 <- rbind(RIL2[,which(colnames(RIL2) %in% colnames(list_all))],list_all[,which(colnames(list_all) %in% colnames(RIL2))])
#moving the new row up
RIL3 <- RIL3[c(502,1:501),]

#-----------------------------------------------------------------------------------------------------------------------
#Subsetting Pericentromeric genes
setwd("~/Documents/HMM_R4_QTL/QTL2")
perigenes <- read.csv("perigenes.csv", header=T, sep=",") #664 genes
pat_mat_founder <- more$gene %>% data.frame()
colnames(pat_mat_founder) <- "gene"
B_peri <- rbind(pat_mat_founder,B_peri)

#now subseting the pericentromeric genes in the RIL3 dataset.
RIL4 <- RIL3[,which(RIL3[1,] %in% B_peri$gene)] 
#is the number equal?
check <- list_alleles[which(list_alleles$list_genes %in% B_peri$gene),] %>% data.frame 
RIL5 <- RIL4[-1,]
write.csv(RIL5, "RIL5.csv")

#----------------------------------------------------------------------------------------------------------------
#Plotting expression of pericentromeric genes for all 8 founder alleles
RIL5 <- read.csv("RIL5.csv", sep=",", header = T)
#RIL5 has the RILs we want and the expression of pericentromeric genes in these genes
BB1_new <- subset(RIL5, RIL5$founder == "BB1")
BB2_new <- subset(RIL5, RIL5$founder == "BB2")
BB3_new <- subset(RIL5, RIL5$founder == "BB3")
BB4_new <- subset(RIL5, RIL5$founder == "BB4")
BB5_new <- subset(RIL5, RIL5$founder == "BB5")
BB6_new <- subset(RIL5, RIL5$founder == "BB6")
BB7_new <- subset(RIL5, RIL5$founder == "BB7")
BB8_new <- subset(RIL5, RIL5$founder == "BB8")

#obtaining the mean of all pericentromeric gene expression for each RILs with BB1 QTL2 allele.
BB1_rowCountPeri <- BB1_new %>% mutate(m=rowMeans(BB1_new[,5:579]))
BB1_rowCountPeri <- BB1_rowCountPeri[,c(2:4,580)]
#now taking the mean of all RIL expression and calculating SE.
se <- function(x) {sd(x)/sqrt(length(x))}
BB1Peri_plot <- data.frame("founder" = c("BB1"), "mean.expression" = mean(BB1_rowCountPeri$m),"se.expression" = se(BB1_rowCountPeri$m))

BB2_rowCountPeri <- BB2_new %>% mutate(m=rowMeans(BB2_new[,5:579]))
BB2_rowCountPeri <- BB2_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB2Peri_plot <- data.frame("founder" = c("BB2"), "mean.expression" = mean(BB2_rowCountPeri$m),"se.expression" = se(BB2_rowCountPeri$m))

BB3_rowCountPeri <- BB3_new %>% mutate(m=rowMeans(BB3_new[,5:579]))
BB3_rowCountPeri <- BB3_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB3Peri_plot <- data.frame("founder" = c("BB3"), "mean.expression" = mean(BB3_rowCountPeri$m),"se.expression" = se(BB3_rowCountPeri$m))

BB4_rowCountPeri <- BB4_new %>% mutate(m=rowMeans(BB4_new[,5:579]))
BB4_rowCountPeri <- BB4_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB4Peri_plot <- data.frame("founder" = c("BB4"), "mean.expression" = mean(BB4_rowCountPeri$m),"se.expression" = se(BB4_rowCountPeri$m))

BB6_rowCountPeri <- BB6_new %>% mutate(m=rowMeans(BB6_new[,5:579]))
BB6_rowCountPeri <- BB6_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB6Peri_plot <- data.frame("founder" = c("BB6"), "mean.expression" = mean(BB6_rowCountPeri$m),"se.expression" = se(BB6_rowCountPeri$m))

BB7_rowCountPeri <- BB7_new %>% mutate(m=rowMeans(BB7_new[,5:579]))
BB7_rowCountPeri <- BB7_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB7Peri_plot <- data.frame("founder" = c("BB7"), "mean.expression" = mean(BB7_rowCountPeri$m),"se.expression" = se(BB7_rowCountPeri$m))

BB8_rowCountPeri <- BB8_new %>% mutate(m=rowMeans(BB8_new[,5:579]))
BB8_rowCountPeri <- BB8_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB8Peri_plot <- data.frame("founder" = c("BB8"), "mean.expression" = mean(BB8_rowCountPeri$m),"se.expression" = se(BB8_rowCountPeri$m))

Peri_plot <- rbind(BB1Peri_plot,BB2Peri_plot,BB3Peri_plot,BB4Peri_plot,BB6Peri_plot,BB7Peri_plot,BB8Peri_plot)
Peri_plot$bais <- c("A","A","A","A","C","A","A")
Peri_plot$founder <- c("B1", "B2","B3","B4","B6","B7","B8")

#Combining all the founders to perform anova test
rowCount_Peri <- rbind(BB1_rowCountPeri,BB2_rowCountPeri,BB3_rowCountPeri,BB4_rowCountPeri,BB6_rowCountPeri,BB7_rowCountPeri,BB8_rowCountPeri) 
resAnova_new_rowCountPeri <- aov(m~founder, rowCount_Peri)
summary(resAnova_new_rowCountPeri)
#F = 5.66 P<1.06e-05 ***
# Compare with alternative method by calculating mean RILs instead of mean genes: F=23.82 P<2e-16 ***
R2 <- 6/(6+494) #Df(founder/Residuals)
R2 #0.012
#1.2% variation in expression explained by founders

#Performing Tukey test
TukeyHSD(resAnova_new_rowCountPeri,conf.level=0.95)

#Linear Regression
linearMod_eu_new_rowCountPeri<- lm(m~founder, rowCount_Peri)
summary(linearMod_eu_new_rowCountPeri)
anova(linearMod_eu_new_rowCountPeri)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (20), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_text( colour="black",family = "Helvetica", size = (18)), 
                      axis.title = element_text(family = "Helvetica", size = (24), colour = "black"),
                      axis.text = element_text(family = "Helvetica", colour = "black", size = (20)))


p <- ggplot(Peri_plot, aes(as.factor(founder),mean.expression, color=bais)) +
  geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2) + xlab("") +
  ylab("Adjusted Expression\n") + theme_bw() +  mynamestheme + theme(legend.position="none") + 
  scale_colour_manual(values=c("red3","pink1")) + geom_hline(aes(yintercept=0), alpha=.5, size=1) + ylim(-0.25, 0.2)

p+annotate("text", x = 1, y= 0.04, label="C", size = 7, color= "red3", fontface="bold") + annotate("text", x = 2, y= 0.06, label="AC", size = 7, color= "red3", fontface="bold") +
  annotate("text", x = 3, y= -0.07, label="C", size = 7, color= "red3", fontface="bold") + annotate("text", x = 4, y= -0.13, label="ABC", size = 7, color= "red3", fontface="bold")  +
  annotate("text", x = 5, y= 0.14, label="A", size = 7, color= "pink1", fontface="bold") + annotate("text", x = 6, y= 0.09, label="ABC", size = 7, color= "red3", fontface="bold") +
  annotate("text", x = 7, y= -0.13, label="B", size = 7, color= "red3", fontface="bold")

