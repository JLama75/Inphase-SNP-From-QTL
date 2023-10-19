setwd("~/Documents/RNAseq/Count_file/TE-count/unaligned-gene-counts/kallisto/abundance.tsv") 

library(tximport)
library(dplyr)
library(DESeq2)

dir <- "~/Documents/RNAseq/Count_file/TE-count/unaligned-gene-counts/kallisto/abundance.tsv"
samples <- c("21147-R1.aligned", "21147-R2.aligned", "21147-R3.aligned" , "21183-R1.aligned", "21183-R3.aligned", "21188-R1.aligned", "21188-R2.aligned", "21188-R3.aligned", "21213-R1.aligned", "21213-R2.aligned", "21213-R3.aligned", "21291-R1.aligned", "21291-R2.aligned", "21291-R3.aligned", "21346-R1.aligned", "21346-R2.aligned", "21346-R3.aligned")
rownames.file <-c("21147-R1", "21147-R2", "21147-R3" , "21183-R1", "21183-R2", "21183-R3", "21188-R1", "21188-R3", "21213-R1", "21213-R2", "21213-R3", "21291-R1", "21291-R2", "21291-R3", "21346-R1", "21346-R2", "21346-R3")

samples <- data.frame(samples)
colnames(samples) <- "run"
files <- file.path(dir, samples$run, "abundance.tsv")
names(files) <- rownames.file
tx2gene <- read.csv("mart_export_ensembleID.csv", sep=",", header = T)
tx2gene <- subset(tx2gene, select = c("Gene", "Gene.stable.ID"))
colnames(tx2gene) <- c("TXNAME", "GENEID")
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)

head(txi.kallisto.tsv$counts)
txi <- data.frame(txi.kallisto.tsv$counts)
txi[is.na(txi)] <- 0
txi2 <- round(txi) 

all.counts <- txi2
all.counts$FBgn <- rownames(all.counts) #14240

#convert FBGN to gene name.
B <- read.csv("coordi_Gene2.csv", header = T)
B<- data.frame(B$FBgn, B$gene_symbol)
colnames(B) <- c("FBgn", "gene_symbol")
AA <- merge(all.counts, B, by = "FBgn", no.dups = FALSE) #14222


missing <- anti_join(all.counts, AA, by= "FBgn")
#searching alternative FBgn IDs for missing genes and substituting it.
missing_gene <- read.csv("missingGenes_fromB.csv", header = T, sep = ",")
missing <- merge(missing, missing_gene, by ="FBgn")
missing <- missing[,-c(1)]
colnames(missing)[18] <- "FBgn"
missing <- missing[,c(18,1:17,19)]
AA <- rbind(AA, missing) #14240

####################@#############################################################################
#Preparing the count files before DESEQ...

#list of all histones (named as FBtr)
AA$gene_symbol <- gsub("His1:CG\\d\\d\\d\\d\\d", "His1", AA$gene_symbol)
AA$gene_symbol <- gsub("His2A:CG\\d\\d\\d\\d\\d", "His2A", AA$gene_symbol)
AA$gene_symbol <- gsub("His2B:CG\\d\\d\\d\\d\\d", "His2B", AA$gene_symbol)
AA$gene_symbol <- gsub("His3:CG\\d\\d\\d\\d\\d", "His3", AA$gene_symbol)
AA$gene_symbol <- gsub("His4:CG\\d\\d\\d\\d\\d", "His4", AA$gene_symbol)
AA$gene_symbol <- gsub("His-Psi:CR\\d\\d\\d\\d\\d", "His-psi", AA$gene_symbol)

#Merging all the copies of histones 
AA1 <- subset(AA, AA$gene_symbol == "His1" | AA$gene_symbol == "His2A" | AA$gene_symbol == "His2B" | AA$gene_symbol == "His3" | AA$gene_symbol == "His4" | AA$gene_symbol == "His-psi")
#write.csv(AA1, "All.histone.csv")
setwd("~/Documents/Paper/Supplementals")
AA1<- read.csv("All.histone.csv", sep=",", header = T)
all.Histones <- AA1$Gene
all.Histones <- AA1$FBgn %>%  data.frame()
colnames(all.Histones) <- "FBgn"
AA1 <- AA1[,1:18]

# deleting rows in all.counts that matches the values in all Histones. using the function %in% (includes in) and ! (not)
deleteRows <- all.counts[ !(all.counts$FBgn %in% all.Histones$FBgn), ] #4121+119 =14240

library(plyr)
#Once the original histones are removed, adding the merged histone list(summed up) to the all.count data frame

x <- ddply(AA, "gene_symbol", numcolwise(sum))
#after summing up selecting only the summed up histone genes
count_histone <-subset(x, x$gene_symbol == "His1" | x$gene_symbol == "His2A" | x$gene_symbol == "His2B" | x$gene_symbol == "His3" | x$gene_symbol == "His4" | x$gene_symbol == "His-psi")
#giving the rowname for gene symbol as Gene 
count_histone$Gene <- count_histone$gene_symbol
count_histone <- count_histone[,c(19,2:18)]
colnames(count_histone)[1] <- "FBgn"

all.counts_final <- merge(deleteRows, count_histone, all= T)

####################@#############################################################################

#all.counts_final <- AA[c(1:18)]
rownames(all.counts_final) <- all.counts_final$FBgn
#txi <- all.counts_final[,-1]

txi <- all.counts_final[,-18]
txi <- round(txi)

#makes an array indicating strain
RIL <- c(rep("21147",3),rep("21183",3),rep("21188",2),rep("21213",3),rep("21291",3),rep("21346",3))
#makes an array indicating genotype
genotype <- c(rep("B6",3),rep("B4",3),rep("B4",2),rep("B6",3),rep("B6",3),rep("B4",3))
#makes an array indicating which pair
pair <- c(rep("A",3),rep("B",3),rep("C",2),rep("B",3),rep("C",3),rep("A",3))
#makes a table combining strain, genotype and pair
exp.Design <- data.frame(genotype,pair,rownames=names(txi))
#exp.Design <- data.frame(genotype,pair,rownames=names(txi))

#run a DEseq analysis estimating abundance of different repeats in BB6 vs BB4
dds <- DESeqDataSetFromMatrix(as.matrix(txi), colData = exp.Design, design =~ pair + genotype)
dds <- DESeq(dds)
res <- results(dds, contrast=c("genotype","B6","B4")) #14127
res <- data.frame(res)
res$FBgn <- rownames(res)
res1 <- res #14127


###########################################################################################################################
#Looking at differential expression in QTL/pericentromere vs. euchromatin

setwd("~/Documents/RNAseq/Count_file/TE-count/unaligned-gene-counts/kallisto/abundance.tsv") 
#res <- res %>% filter(baseMean >= 100)
B <- read.csv("coordi_Gene2.csv", header = T)
res <- merge(res,B, by="FBgn") #14111
res2 <- res
#res[is.na(res)] <- 0

#Chromosome 2
#select for chromosome arm 2L and 2R
Count_2L <- subset(res, Chr_arm == "2L") #475 genes
Count_2R <- subset(res, Chr_arm == "2R") #488 genes
Count_2 <- rbind(Count_2L,Count_2R)

##22160000-5692495
#select for chromosome position: pericentromere
Count_QTL2L <- subset(Count_2L, start >= 22160000) #86 genes
Count_QTL2R <- subset(Count_2R, start <= 5692495) #29 genes
Count_QTL2 <- rbind(Count_QTL2L, Count_QTL2R)
Count_QTL2$arm <- "pericentromere"


#euchromatin
Count_euchromatin_2L <- subset(Count_2L, start <= 22160000) #26 ///38 clusters
Count_euchromatin_2L <- subset(Count_euchromatin_2L, start >=5041) #37 clusters
Count_euchromatin_2R <- subset(Count_2R, start >= 5692495) #30 genes
Count_euchromatin_2R <- subset(Count_euchromatin_2R, start <= 25258060)
Count_euchromatin_2 <- rbind(Count_euchromatin_2L, Count_euchromatin_2R)
Count_euchromatin_2$arm <- "euchromatin"


#merge pericentromere and euchromatin together to represent in a single graph
Count_merge2 <- rbind(Count_QTL2, Count_euchromatin_2)

#Telomere associated region (TAS)

Count_telomere_2L <- subset(Count_2L, start <= 5041) #0
Count_telomere_2R <- subset(Count_2R, start >= 25258060 & start <= 25261551)
Count_telomere_2 <- rbind(Count_telomere_2L, Count_telomere_2R)
#Count_telomere_2$arm <- "telomere"

#make a dataframe for the mean log2FC value for each arm
data <- c(mean(Count_QTL2$log2FoldChange), mean(Count_euchromatin_2$log2FoldChange)) %>% data.frame
colnames(data) <- "grp_mean" 
data$arm <- c("pericentromere", "euchromatin")
mu1 <- data
colnames(mu1) <- c("grp_mean", "arm")

#Generating density plot graph for differentially expressed gene density in pericentromere, euchromatin and telomere.

library(ggplot2)
mynamestheme <- theme(plot.title = element_text( face = "bold", size = (20), hjust = 0.5), 
                      legend.title = element_text(colour = "black",  face = "bold", size = (18)), 
                      legend.text = element_text(face = "bold", colour="black", size = (15)), 
                      axis.title = element_text(size = (20), colour = "black",  face = "bold"),
                      axis.text = element_text(colour = "black", size = (16),  face = "bold"))

#peri and euchromatin arm
ggplot(Count_merge2, aes(x=log2FoldChange, fill=arm)) + geom_density(alpha=.2)  +
  theme_bw() + geom_vline(data = mu1,  aes(xintercept= grp_mean,
                                           color=arm), linetype="dashed", size=1) + 
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)+ mynamestheme + xlim(-2,2) +
  labs(title="2nd Chromosome\n",x="\nLog2FoldChange (Tolerant/Sensitive)", y = "Density\n") +
  annotate("text", x = -0.1, y= 1, label= "ns", size= 7)  + annotate("text", x = 0.7, y= 1, label= "***", size= 7)


# 3rd chromosome peri/centromere
Count_3L <- subset(res, Chr_arm == "3L") 
Count_3R <- subset(res, Chr_arm == "3R")
Count_3 <- rbind(Count_3L,Count_3R)

Count_QTL3L <- subset(Count_3L, start >= 22926900) #19 genes
Count_QTL3R <- subset(Count_3R, start <= 4552934) #12 genes
Count_QTL3 <- rbind(Count_QTL3L, Count_QTL3R) #31
Count_QTL3$arm <- "pericentromere"

#euchromatin
Count_euchromatin_3L <- subset(Count_3L, start <= 22926900) #8///#38 clusters
Count_euchromatin_3L <- subset(Count_euchromatin_3L, start >= 19608) #37
Count_euchromatin_3R <- subset(Count_3R, start >= 4552934) #42
Count_euchromatin_3R <- subset(Count_euchromatin_3R, start <= 31173015) #40
Count_euchromatin_3 <- rbind(Count_euchromatin_3L, Count_euchromatin_3R)
Count_euchromatin_3$arm <- "euchromatin"
Count_arm3 <- rbind(Count_QTL3, Count_euchromatin_3)
#telomere associated clusters

Count_telomere_3L <- subset(Count_3L, start <= 19608)  #1
Count_telomere_3R <- subset(Count_3R, start >= 31173015 ) # 2 clusters in telomere #& start <= 31179331

Count_telomere_3 <- rbind(Count_telomere_3L, Count_telomere_3R)
Count_telomere_3$arm <- "telomere"
Count_arm3 <- rbind(Count_arm3, Count_telomere_3)

#X chromosome
Count_X <- subset(res, Chr_arm == "X") 
Count_QTLX <- subset(Count_X, start >= 22838164) #14 genes
Count_QTLX$arm <- "pericentromere"
#euchromatin
Count_euchromatin_X <- subset(Count_X, start <= 22838164 )
Count_euchromatin_X$arm <- "euchromatin"

Count_armX <- rbind(Count_QTLX, Count_euchromatin_X)

#4th chrom
Count_4 <- subset(res, Chr_arm == "4") #85 genes
Count_4$arm <- "4th chromosome" 

#colors -##00BFC4  #F8766D
mynamestheme <- theme(plot.title = element_text( size = (20), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_text( colour="black", size = (15)), 
                      axis.title = element_text(size = (18), colour = "black"),
                      axis.text = element_text(colour = "black", size = (13)))

B <- ggplot(data=Count_4, aes(x=log2FoldChange, fill=arm))  + scale_fill_manual( values=c("gray3"))  +
  geom_density(alpha=.2) +
  theme_bw() + geom_vline(aes(xintercept=mean(log2FoldChange)),
                          color="black", linetype="dashed", size=1) + 
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)  + theme(legend.position = "top") +
  mynamestheme + labs(title="",x="\nLog2FC (Sensitive/Tolerance)", y = "")+
  annotate("text", x = -0.4, y= 0.86, label="ns", size= 10) + xlim(-1.5,1.5)  

B

qqnorm(Count_4$log2FoldChange) 
shapiro.test(Count_4$log2FoldChange)
#W = 0.93638, p-value = 0.004515 Data is not  normal. 
SIGN.test(pull(Count_4,log2FoldChange))

t.test(Count_4$log2FoldChange, mu = 0) # ***
#t = 4.3848, df = 56, p-value = 5.169e-05
#mitochondrial gene
mitochon_gene <- subset(res, Chr_arm == "mitochondrion_genome") #85 genes
mitochon_gene$arm <- "mitochondrialGenes"

#all pericentromere genes
all_peri_eu <- rbind(Count_arm3, Count_armX, Count_merge2) #2375
setwd("~/Documents/Rotation1/project kelleher/project 2-piRNA alignment/R folder/newPiRNA-align-newRepbase")
#write.csv(all_peri_eu, "all_peri_eu")
all_peri <- subset(all_peri_eu, arm=="pericentromere") #150
all_eu <- subset(all_peri_eu, arm == "euchromatin") #2225
all_tel <- subset(all_peri_eu, arm == "telomere")
all_peri_eu1 <- rbind(all_peri, all_eu, all_tel, Count_4)
all_genes <- rbind(all_peri_eu1,mitochon_gene)

data <- c(mean(all_peri$log2FoldChange), mean(all_eu$log2FoldChange), mean(all_tel$log2FoldChange), mean(Count_4$log2FoldChange))
data <- data.frame(data)
colnames(data) <- "grp_mean"
data$arm <- c("pericentromere", "euchromatin", "telomere", "4th chromosome")
mu <- data
colnames(mu) <- c("grp_mean", "arm")


mynamestheme <- theme(plot.title = element_text( size = (20), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_text( colour="black", size = (24)), 
                      axis.title = element_text(size = (24), colour = "black"),
                      axis.text.x = element_text(colour = "black", size = (20)),
                      axis.text.y = element_text(colour = "black", size = (20))
)

#peri and euchromatin arm
library(ggplot2)
A <- ggplot(all_peri_eu1, aes(x=log2FoldChange, color=arm)) + geom_density(alpha=.2, size =1.1)  +
  theme_bw()+ geom_vline(data = mu,  aes(xintercept= grp_mean,
                                         color=arm), linetype="dashed", size=1) + 
  geom_vline( alpha=.5,aes(xintercept= 0), linetype="solid", size=1)+ mynamestheme + xlim(-1.5,2) +
  labs(title="",x="\nLog2FoldChange(Sensitive/Tolerant)", y = "Density\n") + theme(legend.position = "top", legend.title = element_blank()) +
  annotate("text", x = -0.15, y= 1.55, label= "ns", size= 6, color="chartreuse3")  +  
  annotate("text", x = 0.12, y= 1.5, label= "*", size= 10, color="blue")  + 
  annotate("text", x = 0.6, y= 1.5, label= "***", size= 10, color="lightcoral") +
  annotate("text", x = 0.18, y= 1.3, label= "ns", size= 6, color="gray") +
  scale_fill_manual(values =c("gray","deepskyblue", "lightcoral", "chartreuse3")) +
  scale_color_manual(values =c("gray","blue", "lightcoral", "chartreuse3")) 
A
t.test(all_eu$log2FoldChange,all_peri$log2FoldChange)
t.test(all_eu$log2FoldChange,Count_4$log2FoldChange)
#t = -9.3067, df = 141.03, p-value = 2.335e-16 ***
#t = -4.5627, df = 53.178, p-value = 3.014e-05 ***
t.test(all_eu$log2FoldChange,all_tel$log2FoldChange)
#t = -8.9887, df = 141.21, p-value = 1.468e-15
#t = -4.3057, df = 53.217, p-value = 7.186e-05

qqnorm(all_tel$log2FoldChange,datax=TRUE) + qqline(all_tel$log2FoldChange,datax=TRUE)
shapiro.test(all_tel$log2FoldChange)
#W = 0.96509, p-value = 0.171 Data is normal
t.test(all_tel$log2FoldChange, mu = 0)
#t = -0.22656, df = 49, p-value = 0.8217 m =-0.01045122 

library(ggpubr)
setwd("~/Documents/Paper/Supplementals")
supplemental <- all_genes[,c(1,11,8,9,10,2,3,7,12)]
supplementals <- supplemental

supplemental <- supplemental %>% filter(baseMean >= 100) 
supplemental <- supplemental %>% filter(padj <= 0.05)
supplemental <- supplemental %>% filter(log2FoldChange >= 0.586 | log2FoldChange <= -0.586 ) #530
#write.csv(supplemental, "DEG_outlier_histone.csv")
res <- res1 %>% filter(baseMean >= 100) 
res <- res %>% filter(padj <= 0.05) #1312
res <- res %>% filter(log2FoldChange >= 0.586 | log2FoldChange <= -0.586 ) #530 (histone paralogs added)

#select for chromosome arm 2L and 2R
Count_2L <- subset(supplemental, Chr_arm == "2L") #475 genes
Count_2R <- subset(supplemental, Chr_arm == "2R") #488 genes
Count_2 <- rbind(Count_2L,Count_2R)

##22160000-5692495
#select for chromosome position: pericentromere
#2L:20,820,000- 2R:6,942,495
#2L:19,010,000-20,000,000

Count_QTL2L <- subset(Count_2L, start >= 20820000) # genes
Count_QTL2R <- subset(Count_2R, start <= 6942495) # genes
Count_QTL2 <- rbind(Count_QTL2L, Count_QTL2R)#31

Count_QTL1 <- subset(Count_2L, start >= 19010000 & start <= 20000000 ) #13 genes
#13+31 = 44 out of 530 DEG within the QTL

##########################
#selecting all genes that encode proteins that are members of TIP60 complex
TIP60<- c("FBgn0026577","FBgn0267398","FBgn0020306","FBgn0040075","FBgn0040078","FBgn0025716","FBgn0000046","FBgn0031873","FBgn0034537","FBgn0032321","FBgn0039654","FBgn0053554","FBgn0026080","FBgn0027378","FBgn0035624","FBgn0033341","FBgn0000581","FBgn0030945")
TIP60 <- TIP60 %>% data.frame
colnames(TIP60) <- "FBgn" #18

TIP60_merge<- merge(res2,TIP60, by="FBgn")
TIP60_merge <- TIP60_merge[,c(1,11,8,9,10,2,3,7)]
setwd("~/Documents/Paper/Supplementals")

write.csv(TIP60_merge,"Tip60_Exp.csv")
