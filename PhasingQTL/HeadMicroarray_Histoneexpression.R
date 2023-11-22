setwd("~/Documents/RNAseq/DSPR_HMM/HMMregB_R2")
library("dplyr")
library("rlist")
library(ggplot2)
#Aim: To look at expression of histone genes among the 8 founder alleles using microarray gene expression data from head tissue of Drosophila melanogaster

#csv file with 864 unique RILs (recombinant inbred line) and information on the founder allele in 19 genomic positions
HMM<- read.csv("HMM_sorted_release4.csv") 
HMM <- subset(HMM, select= -X)
colnames(HMM) <- c("pos", "RIL", "founder")
list1 <- unique(HMM$RIL) %>% data.frame
HMM$continguent <- NA

#Only selecting thoes RILs that have single founder allele in all 19 positions

for(i in seq(1,(dim(HMM)[1])-18, 19 )){ 
  #print(c(i, HMM[i,1], HMM[i,2]))
  n = i
  loop <- HMM[i:(i+18),] 
  if (length(unique(loop$RIL)) == 1)
    {
      if (length(unique(loop$founder)) == 1)
        {
        HMM[i:(i+18),4] <- "True"
      }
       else {HMM[i:(i+18),4] <- "False"}
    }
  else message("Error in ", n)
}

#subsetting all the continguent genotypes and removing NA values
HMM3 <- subset(HMM, HMM$continguent == "True")
HMM3 <- HMM3[!is.na(HMM3$founder),]
final_list <- unique(subset(HMM3, select = c(RIL, founder)))
final_list <- final_list[order(final_list$founder),]

#subsetting individual founders.
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
RIL<- read.csv("FemaleHeadExpression.csv", sep = ",", header = T)
#Now merging RIL and final list by paternal RIL
colnames(final_list) <- c("patRIL", "founder")
RIL2 <- merge(final_list, RIL, by = "patRIL") 
RIL2 <- RIL2[order(RIL2$founder),]

#also adding the founder genotype info in the data frame, for that first subsetting the founder and patRIL in a sepearate table.
founders <- subset(RIL2, select = c(founder, patRIL))
rownames(RIL2)<- RIL2$patRIL
# Total gene expression of (563 from 597 microarray data)

#extracting only the 111 histone genes 
gene_names <- read.csv("Histones.csv", header = F)
colnames(gene_names) = "X"
rownames(gene_names) <- gene_names$X
gene_names1 <- data.frame(t(gene_names)) 

RIL3 <- rbind(RIL2[,which(colnames(RIL2) %in% colnames(gene_names1))],gene_names1[,which(colnames(gene_names1) %in% colnames(RIL2))])
RIL3$patRIL <- rownames(RIL3)
RIL3 <- merge(RIL3,founders, by = "patRIL") #Adding the founder genotype column
RIL3 <- RIL3[,c(113,1:112)]
RIL3 <- RIL3[order(RIL3$founder),]

#converting CG to geneSymbol and merging all the copies from each of the five histone genes
CGG <- read.table("Histone_CG_FBgn.txt", header = TRUE) #contains gene names in CG and FBgn id
GeneSymbol <- read.table("Histone_FBgn_geneSymbol.txt", header = TRUE) #contains FBgn id and gene symbol
merged_gene <- merge(CGG,GeneSymbol, by = "FBgn")
merged_gene$histoneGene <- merged_gene$Gene_name

#substituting the His1:CG number to simply His1 and so on.
merged_gene$histoneGene <- gsub("His1:CG\\d\\d\\d\\d\\d", "His1", merged_gene$Gene_name, perl= TRUE)
merged_gene$histoneGene <- gsub("His2A:CG\\d\\d\\d\\d\\d", "His2A", merged_gene$histoneGene, perl= TRUE)
merged_gene$histoneGene <- gsub("His2B:CG\\d\\d\\d\\d\\d", "His2B", merged_gene$histoneGene, perl= TRUE)
merged_gene$histoneGene <- gsub("His3:CG\\d\\d\\d\\d\\d", "His3", merged_gene$histoneGene, perl= TRUE)
merged_gene$histoneGene <- gsub("His4:CG\\d\\d\\d\\d\\d", "His4", merged_gene$histoneGene, perl= TRUE)
merged_gene <- merged_gene[order(merged_gene$histoneGene),]

#only selecting gene name (CG...) and histone names 
Histone_Genes <- subset(merged_gene, select = c(CG, histoneGene))
rownames(Histone_Genes) <- Histone_Genes$CG
Histone_Genes <- subset(Histone_Genes, select = histoneGene)

#Converting row into column to facilitate merging with the RIL3 data set.
Histone_Genes <-data.frame(t(Histone_Genes))
Histone_Genes$founder <- NA
Histone_Genes$patRIL <- NA
Histone_Genes <- Histone_Genes[,c(112:113, 1:111)] #rearranging the column
#merging based on colnames.
RIL4 <- rbind(RIL3[,which(colnames(RIL3) %in% colnames(Histone_Genes))],Histone_Genes[,which(colnames(Histone_Genes) %in% colnames(RIL3))])
RIL4 <- RIL4[c(564, 1:563),]

#Now taking the mean of all the copies of individual classes of histone genes for each founders (each row)

#Listing the column number containing copies of histone1.
His1.data <- list(0)
for(i in 3:113) {
  
  if (RIL4[1,i]== "His1") {
    col_num <- i
    His1.data <- list.append(His1.data, col_num)
  }
}

#Calculating mean expression of histone 1 gene copies for all rows.
RIL4$His1 <- NA
His1_list <- list(0) #another list for keeping the values of  histone 1 
for(k in 2:dim(RIL4)[1])
{
  for(i in 2:length(His1.data)){
  #creating a list of gene expression values of histone 1 gene 
  His1_list <- list.append(His1_list,RIL4[k, His1.data[[i]]])
  }

  RIL4[k,114] <- mean(as.numeric(His1_list[-1]))
  His1_list <- list(0) #resetting the list
}

#For checking the above code works or not:
mean(as.numeric(RIL4[2,c(5,13,24,26,31,35,37,40,41,47,48,58,59,62,70,73,79,84,88,90,100,108,113)])) 
RIL4[2,114]

#For His2A, His2B, His3, His4:
His2A.data <- list(0)
His2B.data <- list(0)
His3.data <- list(0)
His4.data <- list(0)

for(i in 3:113) {
  
  if (RIL4[1,i]== "His2A") {
    col_num <- i
    His2A.data <- list.append(His2A.data, col_num)
    print(col_num)
    print(His2A.data)
  }
  if (RIL4[1,i]== "His2B") {
    col_num <- i
    His2B.data <- list.append(His2B.data, col_num)
    print(col_num)
    print(His2B.data)}
  
  if (RIL4[1,i]== "His3") {
    col_num <- i
    His3.data <- list.append(His3.data, col_num)
    print(col_num)
    print(His3.data)}
  
  if (RIL4[1,i]== "His4") {
    col_num <- i
    His4.data <- list.append(His4.data, col_num)
    print(col_num)
    print(His4.data)}
}

RIL4$His2A <- NA
RIL4$His2B <- NA
RIL4$His3 <- NA
RIL4$His4 <- NA
His2A_list <- list(0) #another list for keeping the values of  histone 1
His2B_list <- list(0)
His3_list <- list(0)
His4_list <- list(0)
for(k in 2:dim(RIL4)[1])
{
  for(i in 2:length(His2A.data)){
    
    His2A_list <- list.append(His2A_list,RIL4[k, His2A.data[[i]]])
    #print(i)
  }
  
  for(i in 2:length(His2B.data)){
    
    His2B_list <- list.append(His2B_list,RIL4[k, His2B.data[[i]]])
    #print(i)
  }
  
  for(i in 2:length(His3.data)){
    
    His3_list <- list.append(His3_list,RIL4[k, His3.data[[i]]])
    #print(i)
  }
  
  for(i in 2:length(His4.data)){
    
    His4_list <- list.append(His4_list,RIL4[k, His4.data[[i]]])
    #print(i)
  }
  #print(His1_list)
  RIL4[k,115] <- mean(as.numeric(His2A_list[-1]))
  RIL4[k,116] <- mean(as.numeric(His2B_list[-1]))
  RIL4[k,117] <- mean(as.numeric(His3_list[-1]))
  RIL4[k,118] <- mean(as.numeric(His4_list[-1]))
  
  His2A_list <- list(0)
  His2B_list <- list(0)
  His3_list <- list(0)
  His4_list <- list(0)
}
#checking if the code above worked
#mean(as.numeric(RIL6[4,c(8,19,28,33,38,42,46,55,60,69,72,80,85,87,92,93,102,106,107,112)]))
Histone_mean <- subset(RIL4, select = c(founder, patRIL ,His1, His2A, His2B, His3, His4))
Histone_mean <- Histone_mean[-1,]
rownames(Histone_mean) <- c(1:563)

#boxplot

His1 <- subset(Histone_mean, select= c(founder,His1))

His1$phenotype[His1$founder == "BB1"] <- "tolerant"
His1$phenotype[His1$founder == "BB4"] <- "tolerant"
His1$phenotype[His1$founder == "BB8"] <- "tolerant"
His1$phenotype[His1$founder == "BB6"] <- "sensitive"
His1$phenotype[His1$founder == "BB2"] <- "intermediate"
His1$phenotype[His1$founder == "BB3"] <- "intermediate"
His1$phenotype[His1$founder == "BB7"] <- "intermediate"

colnames(His1) <- c("founder", "expression", "phenotype")

#To find if founder phenotype is a factor of histone expression variation: if there is variation between the groups 
#F= ratio(between group expression variation/ within group expression variation). Larger F --> more likely the group have different means.
resAnova_His1 <- aov(expression~founder, His1)
summary(resAnova_His1)

#To find the difference in B6 and B8 is significant

linearMod_His1 <- lm(expression~founder, His1)
summary(linearMod_His1)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (20), hjust = 0.5), 
                      legend.title = element_text(colour = "black",  face = "bold", family = "Helvetica", size = (18)), 
                      legend.text = element_text(face = "bold", colour="black",family = "Helvetica", size = (15)), 
                      axis.title = element_text(family = "Helvetica", size = (18), colour = "black",  face = "bold"),
                      axis.text = element_text(family = "Courier", colour = "black", size = (15),  face = "bold"))

ggplot(His1, aes(x=founder, y=expression,fill = phenotype)) + geom_hline(yintercept=0)+
  geom_boxplot() + scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) + geom_point()+
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  theme_bw() + mynamestheme + labs(y=" His1 expression") +
   annotate("text", x = 5, y= 3.2, label="***", size= 7) + annotate("text", x = 7, y= 1.8, label="**", size= 7)

ggplot(His1, aes(founder, expression))+ geom_hline(yintercept=0) + geom_violin(aes(founder, expression,fill = phenotype)) +
  scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) + geom_point() +
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  theme_bw() + mynamestheme + labs(y=" His1 expression") +
  annotate("text", x = 5, y= 3.2, label="***", size= 7) + annotate("text", x = 7, y= 1.8, label="**", size= 7)

#Histone 2A
His2A <- subset(Histone_mean, select= c(founder,His2A))
His2A <- subset(Histone_mean, select= c(founder,His2A))
His2A$phenotype[His2A$founder == "BB1"] <- "tolerant"
His2A$phenotype[His2A$founder == "BB4"] <- "tolerant"
His2A$phenotype[His2A$founder == "BB8"] <- "tolerant"
His2A$phenotype[His2A$founder == "BB6"] <- "sensitive"
His2A$phenotype[His2A$founder == "BB2"] <- "intermediate"
His2A$phenotype[His2A$founder == "BB3"] <- "intermediate"
His2A$phenotype[His2A$founder == "BB7"] <- "intermediate"

colnames(His2A) <- c("founder", "expression", "phenotype")
resAnova_His2A <- aov(expression~founder, His2A)
summary(resAnova_His2A)

linearMod_His2A <- lm(expression~founder, His2A)
summary(linearMod_His2A)



ggplot(His2A, aes(x=founder, y=expression,fill = phenotype)) + geom_hline(yintercept=0)+
  geom_boxplot() + scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) + geom_point()+
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  theme_bw() + mynamestheme + labs(y=" His2A expression") +
  annotate("text", x = 5, y= 3.4, label="***", size= 7) + annotate("text", x = 7, y= 1.8, label="*", size= 7)

ggplot(His2A, aes(founder, expression))+ geom_hline(yintercept=0) + geom_violin(aes(founder, expression,fill = phenotype)) +
  scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) + geom_point() +
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  theme_bw() + mynamestheme + labs(y=" His2A expression") +
  annotate("text", x = 5, y= 3.4, label="***", size= 7) + annotate("text", x = 7, y= 1.8, label="*", size= 7)

#Histone 2B
His2B <- subset(Histone_mean, select= c(founder,His2B))
His2B <- subset(Histone_mean, select= c(founder,His2B))
His2B$phenotype[His2B$founder == "BB1"] <- "tolerant"
His2B$phenotype[His2B$founder == "BB4"] <- "tolerant"
His2B$phenotype[His2B$founder == "BB8"] <- "tolerant"
His2B$phenotype[His2B$founder == "BB6"] <- "sensitive"
His2B$phenotype[His2B$founder == "BB2"] <- "intermediate"
His2B$phenotype[His2B$founder == "BB3"] <- "intermediate"
His2B$phenotype[His2B$founder == "BB7"] <- "intermediate"

resAnova_His2B <- aov(His2B~founder, His2B)
summary(resAnova_His2B)
colnames(His2B) <- c("founder" , "expression", "phenotype")

linearMod_His2B <- lm(expression~founder, His2B)
summary(linearMod_His2B)

ggplot(His2B, aes(x=founder, y=expression,fill = phenotype)) + geom_hline(yintercept=0)+
  geom_boxplot() + scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) + geom_point() +
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  theme_bw() + mynamestheme + labs(y=" His2B expression") +
  annotate("text", x = 5, y= 3, label="***", size= 7)

ggplot(His2B, aes(founder, expression)) + geom_hline(yintercept=0)+ geom_violin(aes(founder, expression,fill = phenotype)) +
  scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) + geom_point() +
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  theme_bw() + mynamestheme + labs(y=" His2B expression") +
  annotate("text", x = 5, y= 3, label="***", size= 7)

#Histone 3

His3 <- subset(Histone_mean, select= c(founder,His3))
His3$phenotype[His3$founder == "BB1"] <- "tolerant"
His3$phenotype[His3$founder == "BB4"] <- "tolerant"
His3$phenotype[His3$founder == "BB8"] <- "tolerant"
His3$phenotype[His3$founder == "BB6"] <- "sensitive"
His3$phenotype[His3$founder == "BB2"] <- "intermediate"
His3$phenotype[His3$founder == "BB3"] <- "intermediate"
His3$phenotype[His3$founder == "BB7"] <- "intermediate"

resAnova_His3 <- aov(His3~founder, His3)
summary(resAnova_His3)

colnames(His3) <- c("founder" , "expression", "phenotype")

linearMod_His3 <- lm(expression~founder, His3)
summary(linearMod_His3)

ggplot(His3, aes(x=founder, y=expression,fill = phenotype)) + geom_hline(yintercept=0)+
  geom_boxplot() + scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) +
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) + geom_point()+ theme_bw() + mynamestheme + labs(y= "His3 expression") +
  annotate("text", x = 5, y= 3.2, label="***", size= 7) + annotate("text", x = 7, y= 1.8, label="**", size= 7)


ggplot(His3, aes(founder, expression))+ geom_hline(yintercept=0) + geom_violin(aes(founder, expression,fill = phenotype)) +
  scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) + geom_point() +
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  theme_bw() + mynamestheme + labs(y=" His3 expression") +
  annotate("text", x = 5, y= 3.2, label="***", size= 7) + annotate("text", x = 7, y= 1.8, label="**", size= 7)

#Histone 4
His4 <- subset(Histone_mean, select= c(founder,His4))
His4$phenotype[His4$founder == "BB1"] <- "tolerant"
His4$phenotype[His4$founder == "BB4"] <- "tolerant"
His4$phenotype[His4$founder == "BB8"] <- "tolerant"
His4$phenotype[His4$founder == "BB6"] <- "sensitive"
His4$phenotype[His4$founder == "BB2"] <- "intermediate"
His4$phenotype[His4$founder == "BB3"] <- "intermediate"
His4$phenotype[His4$founder == "BB7"] <- "intermediate"

resAnova_His4 <- aov(His4~founder, His4)
summary(resAnova_His4)

colnames(His4) <- c("founder" , "expression", "phenotype")

linearMod_His4 <- lm(expression~founder, His4)
summary(linearMod_His4)

ggplot(His4, aes(x=founder, y=expression,fill = phenotype))+ geom_hline(yintercept=0)+
  geom_boxplot() + scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) +
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  geom_point()+theme_bw() + mynamestheme+ labs(y="His4 expression") +
  annotate("text", x = 5, y= 3, label="***", size= 7) +  annotate("text", x = 2, y= 3, label="*", size= 7)

ggplot(His4, aes(founder, expression)) + geom_hline(yintercept=0) + geom_violin(aes(founder, expression,fill = phenotype)) +
  scale_y_continuous(breaks = c(-4,-3, -2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4)) + geom_point() +
  scale_fill_manual(values=c("deeppink", "pink1", "darkred")) +  theme_bw() + mynamestheme +labs(y=" His4 expression") +
  annotate("text", x = 5, y= 3, label="***", size= 7) +  annotate("text", x = 2, y= 3, label="*", size= 7)


summary(resAnova_His1)
summary(resAnova_His2A)
summary(resAnova_His2B)
summary(resAnova_His3)
summary(resAnova_His4)

summary(linearMod_His1)
summary(linearMod_His2A)
summary(linearMod_His2B)
summary(linearMod_His3)
summary(linearMod_His4)
