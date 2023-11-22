#you will need to change this to whatever your directory is.
#setwd("~/Documents/computational/R_working/association_study/PopB-Final/paper_analyses/QTL_phasing")
setwd("~/Documents/QTLphasing_responder_260")
library("ggplot2")
#this is a data frame of the phenotypes of every RIL we measured, after accounting for block, experimenter and age
load("conc.phenos.11-9-17")
#these are the HMM genotypes of the phenotyped RILS at the LOD peak close to ranGAP
load("conc.2L.19420000")
#these are the HMM genotypes of the phenotyped RILS at the LOD peak that spans the centromere
load("conc.2L.20770000")
#these are the HMM genotypes of the phenotyped RILS at the LOD peak that includes 42AB
load("conc.2R.2080000")


#here is where he "hard code" the genotypes, assigning RILs with a >95% chance of having the allele of 1 particular founder as carrying that allele
conc.2L.19420000$founder[conc.2L.19420000$BB1 > 0.95] <- "BB1"
conc.2L.19420000$founder[conc.2L.19420000$BB2 > 0.95] <- "BB2"
conc.2L.19420000$founder[conc.2L.19420000$BB3 > 0.95] <- "BB3"
conc.2L.19420000$founder[conc.2L.19420000$BB4 > 0.95] <- "BB4"
conc.2L.19420000$founder[conc.2L.19420000$BB5 > 0.95] <- "BB5"
conc.2L.19420000$founder[conc.2L.19420000$BB6 > 0.95] <- "BB6"
conc.2L.19420000$founder[conc.2L.19420000$BB7 > 0.95] <- "BB7"
conc.2L.19420000$founder[conc.2L.19420000$BB8 > 0.95] <- "BB8"

conc.2L.20770000$founder[conc.2L.20770000$BB1 > 0.95] <- "BB1"
conc.2L.20770000$founder[conc.2L.20770000$BB2 > 0.95] <- "BB2"
conc.2L.20770000$founder[conc.2L.20770000$BB3 > 0.95] <- "BB3"
conc.2L.20770000$founder[conc.2L.20770000$BB4 > 0.95] <- "BB4"
conc.2L.20770000$founder[conc.2L.20770000$BB5 > 0.95] <- "BB5"
conc.2L.20770000$founder[conc.2L.20770000$BB6 > 0.95] <- "BB6"
conc.2L.20770000$founder[conc.2L.20770000$BB7 > 0.95] <- "BB7"
conc.2L.20770000$founder[conc.2L.20770000$BB8 > 0.95] <- "BB8"

conc.2R.2080000$founder[conc.2R.2080000$BB1 > 0.95] <- "BB1"
conc.2R.2080000$founder[conc.2R.2080000$BB2 > 0.95] <- "BB2"
conc.2R.2080000$founder[conc.2R.2080000$BB3 > 0.95] <- "BB3"
conc.2R.2080000$founder[conc.2R.2080000$BB4 > 0.95] <- "BB4"
conc.2R.2080000$founder[conc.2R.2080000$BB5 > 0.95] <- "BB5"
conc.2R.2080000$founder[conc.2R.2080000$BB6 > 0.95] <- "BB6"
conc.2R.2080000$founder[conc.2R.2080000$BB7 > 0.95] <- "BB7"
conc.2R.2080000$founder[conc.2R.2080000$BB8 > 0.95] <- "BB8"

#HERE is where we do QTL phasing for the first LOD peak

conc.2L.19420000.dat <- cbind(conc.2L.19420000$founder,phenos)
names(conc.2L.19420000.dat)[1] <- "founder"
tapply(conc.2L.19420000.dat$founder,conc.2L.19420000.dat$founder,length)

conc.2L.19420000.full.lm <- lm(res~founder, conc.2L.19420000.dat) 
summary(conc.2L.19420000.full.lm)
library(ggplot2)
tapply(conc.2L.19420000.dat$res,conc.2L.19420000.dat$founder,mean)
#the Estimate in the summary above is the adjusted atrophy associated with each founder allele at the LOD peak
#we use the estimate to sort the founder alleles according to the odds of atrophy, starting with the lowest odds and moving to the highest
founder.order <- c("BB4","BB8","BB1","BB3","BB7","BB2","BB6")


#here we explore sequential partitions of the founders into two allelic classes "A" and "B", to see which group best fits the data
conc.2L.19420000.null.lm <- lm(res~ 1, na.exclude(conc.2L.19420000.dat)) 
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[1:6]] <- "A"
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[1:5]] <- "A"
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[6:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[1:4]] <- "A"
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[5:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[1:3]] <- "A"
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[4:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.19420000.dat))  
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[1:2]] <- "A"
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[3:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[1]] <- "A"
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[2:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)

#The best fit is 1:4, 5:7, which has the highest F stat and the lowest P-value
#we fix those two partitions here
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[1:4]] <- "A"
conc.2L.19420000.dat$partition[conc.2L.19420000.dat$founder %in% founder.order[5:7]] <- "B"
conc.2L.19420000.null.lm <- lm(res~ partition, na.exclude(conc.2L.19420000.dat)) 

#now we ask if the model can be improved by adding a 3rd partition
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[1:6]] <- "A"
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[1:5]] <- "A"
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[6:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[1:3]] <- "A"
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[4:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[1:2]] <- "A"
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[3:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[1]] <- "A"
conc.2L.19420000.dat$partition2[conc.2L.19420000.dat$founder %in% founder.order[2:7]] <- "B"
conc.2L.19420000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.19420000.dat)) 
anova(conc.2L.19420000.null.lm,conc.2L.19420000.alt.lm)
# we find marginal support for 1:5, 5:6, 7 (p = 0.0556) but it is way above 10-4 so we find no evidence for 2 QTL at this locus. 
#final paritioning is "BB4" "BB8" "BB1" "BB3" (A) "BB7" "BB2" "BB6" (B)

#plot founder means and SE
se <- function(x) {sd(x)/sqrt(length(x))}
conc.2L.19420000.conc.plot <- data.frame("founder" = c("BB1","BB2","BB3","BB4","BB6","BB7","BB8","BB5"), "mean.residual" = c(tapply(conc.2L.19420000.dat$res,conc.2L.19420000.dat$founder,mean),0),"se.residual" = c(tapply(conc.2L.19420000.dat$res,conc.2L.19420000.dat$founder,se),0))
conc.2L.19420000.conc.plot$groups <- c("A","C","A","A","C","C","A",NA) 
ggplot(conc.2L.19420000.conc.plot, aes(x=as.factor(founder),y=mean.residual,color=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none") + scale_colour_manual(values=c("pink1","red3","cadetblue3"))
ggplot(conc.2L.19420000.conc.plot, aes(x=as.factor(founder),y=mean.residual,fill=groups)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=1,width=0.2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none") + scale_fill_manual(values=c("red3","pink1","cadetblue3"))




#HERE is where we do QTL phasing for the second LOD peak

conc.2L.20770000.dat <- cbind(conc.2L.20770000$founder,phenos)
names(conc.2L.20770000.dat)[1] <- "founder"
tapply(conc.2L.20770000.dat$founder,conc.2L.20770000.dat$founder,length)

conc.2L.20770000.full.lm <- lm(res~founder, conc.2L.20770000.dat) 
summary(conc.2L.20770000.full.lm)
plot(res~founder, conc.2L.20770000.dat)
abline(conc.2L.20770000.full.lm)
#the Estimate in the summary above is the adjusted atrophy associated with each founder allele at the LOD peak
#we use the estimate to sort the founder alleles according to the odds of atrophy, starting with the lowest odds and moving to the highest
founder.order <- c("BB4","BB1","BB8","BB7","BB3","BB2","BB6")


#here we explore sequential partitions of the founders into two allelic classes "A" and "B", to see which group best fits the data
conc.2L.20770000.null.lm <- lm(res~ 1, na.exclude(conc.2L.20770000.dat)) 
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[1:6]] <- "A"
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[1:5]] <- "A"
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[6:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[1:4]] <- "A"
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[5:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[1:3]] <- "A"
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[4:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.20770000.dat))  
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[1:2]] <- "A"
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[3:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[1]] <- "A"
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[2:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)

#The best fit is 1:6, 7, which has the highest F stat and the lowest P-value
#we fix those two partitions here
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[1:6]] <- "A"
conc.2L.20770000.dat$partition[conc.2L.20770000.dat$founder %in% founder.order[7]] <- "B"
conc.2L.20770000.null.lm <- lm(res~ partition, na.exclude(conc.2L.20770000.dat)) 

#now we ask if the model can be improved by adding a 3rd partition
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[1:5]] <- "A"
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[6:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[1:4]] <- "A"
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[5:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[1:3]] <- "A"
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[4:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[1:2]] <- "A"
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[3:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[1]] <- "A"
conc.2L.20770000.dat$partition2[conc.2L.20770000.dat$founder %in% founder.order[2:7]] <- "B"
conc.2L.20770000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2L.20770000.dat)) 
anova(conc.2L.20770000.null.lm,conc.2L.20770000.alt.lm)
# we find significant support for 1:3, 4:6, 7 (0.007981) but it is above the typical 10-4 so its in a grey area
#final paritioning is "BB4","BB1","BB8" (A1) "BB7","BB3","BB2" (A2),"BB6" (B)

write.csv(conc.2L.20770000.dat, "conc.2L.20770000.dat.csv")
#plot founder means and SE
se <- function(x) {sd(x)/sqrt(length(x))}
conc.2L.20770000.conc.plot <- data.frame("founder" = c("BB1","BB2","BB3","BB4","BB6","BB7","BB8","BB5"), "mean.residual" = c(tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,mean),0),"se.residual" = c(tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,se),0))
conc.2L.20770000.conc.plot$groups <- c("A","A","A","A","C","A","A",NA) 


ggplot(conc.2L.20770000.conc.plot, aes(x=as.factor(founder),y=mean.residual,color=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none") + scale_colour_manual(values=c("pink1","red3"))
ggplot(conc.2L.20770000.conc.plot, aes(x=as.factor(founder),y=mean.residual,fill=groups)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=1,width=0.2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none") + scale_fill_manual(values=c("red3","cadetblue3","pink1"))
colnames(conc.2L.20770000.dat)
BB6 <- subset(conc.2L.20770000.dat, founder == "BB6")
BB8 <- subset(conc.2L.20770000.dat, founder == "BB8")

#HERE is where we do QTL phasing for the third LOD peak

conc.2R.2080000.dat <- cbind(conc.2R.2080000$founder,phenos)
names(conc.2R.2080000.dat)[1] <- "founder"
tapply(conc.2R.2080000.dat$founder,conc.2R.2080000.dat$founder,length)

conc.2R.2080000.full.lm <- lm(res~founder, conc.2R.2080000.dat) 
summary(conc.2R.2080000.full.lm)
tapply(conc.2R.2080000.dat$res,conc.2R.2080000.dat$founder,mean)
#the Estimate in the summary above is the adjusted atrophy associated with each founder allele at the LOD peak
#we use the estimate to sort the founder alleles according to the odds of atrophy, starting with the lowest odds and moving to the highest
founder.order <- c("BB8","BB1","BB7","BB3","BB4","BB2","BB6")


#here we explore sequential partitions of the founders into two allelic classes "A" and "B", to see which group best fits the data
conc.2R.2080000.null.lm <- lm(res~ 1, na.exclude(conc.2R.2080000.dat)) 
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[1:6]] <- "A"
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[1:5]] <- "A"
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[6:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[1:4]] <- "A"
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[5:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[1:3]] <- "A"
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[4:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition, na.exclude(conc.2R.2080000.dat))  
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[1:2]] <- "A"
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[3:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[1]] <- "A"
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[2:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)

#The best fit is 1:6, 7, which has the highest F stat and the lowest P-value
#we fix those two partitions here
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[1:6]] <- "A"
conc.2R.2080000.dat$partition[conc.2R.2080000.dat$founder %in% founder.order[7]] <- "B"
conc.2R.2080000.null.lm <- lm(res~ partition, na.exclude(conc.2R.2080000.dat)) 

#now we ask if the model can be improved by adding a 3rd partition
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[1:5]] <- "A"
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[6:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[1:4]] <- "A"
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[5:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[1:3]] <- "A"
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[4:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)

conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[1:2]] <- "A"
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[3:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[1]] <- "A"
conc.2R.2080000.dat$partition2[conc.2R.2080000.dat$founder %in% founder.order[2:7]] <- "B"
conc.2R.2080000.alt.lm <- lm(res~ partition + partition2, na.exclude(conc.2R.2080000.dat)) 
anova(conc.2R.2080000.null.lm,conc.2R.2080000.alt.lm)
# we find significant support for 1:3, 4:6, 7 (0.007981) but it is above the typical 10-4 so its in a grey area
#final paritioning is "BB8","BB1","BB7" (A1) "BB3","BB4","BB2",(A2),"BB6" (B)

se <- function(x) {sd(x)/sqrt(length(x))}
conc.2R.2080000.conc.plot <- data.frame("founder" = c("BB1","BB2","BB3","BB4","BB6","BB7","BB8","BB5"), "mean.residual" = c(tapply(conc.2R.2080000.dat$res,conc.2R.2080000.dat$founder,mean),0),"se.residual" = c(tapply(conc.2R.2080000.dat$res,conc.2R.2080000.dat$founder,se),0))
conc.2R.2080000.conc.plot$groups <- c("A","A","A","A","C","A","A",NA) 
ggplot(conc.2R.2080000.conc.plot, aes(x=as.factor(founder),y=mean.residual,color=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none") + scale_colour_manual(values=c("pink1","red3"))
ggplot(conc.2R.2080000.conc.plot, aes(x=as.factor(founder),y=mean.residual,fill=groups)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=1,width=0.2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none") + scale_fill_manual(values=c("red3","cadetblue3","pink1"))

#Bar graph for all 3 qtl peaks
three.qtl <- rbind(conc.2L.19420000.conc.plot,conc.2L.20770000.conc.plot,conc.2R.2080000.conc.plot)
three.qtl$qtl <- c(rep("QTL 1",8),rep("QTL 2",8),rep("QTL 3",8))
levels(three.qtl$founder) <- c("B1","B2","B3","B4","B5","B6","B7","B8")
#vertical facet
ggplot(three.qtl, aes(x=as.factor(founder),y=mean.residual,fill=groups)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=1,width=0.2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none",strip.text = element_text(size=18)) + scale_fill_manual(values=c("red3","deeppink","pink1")) + facet_grid(qtl~.)
ggsave("qtl_vertical.stack.pdf",width=5,height=10)
#horizontal facet
ggplot(three.qtl, aes(x=as.factor(founder),y=mean.residual,fill=groups)) + geom_bar(stat="identity") + geom_errorbar(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=1,width=0.2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none",strip.text = element_text(size=18)) + scale_fill_manual(values=c("red3","pink1")) + facet_grid(.~qtl)
ggsave("qtl_horizontal.stack.pdf",width=10,height=5)

ggplot(three.qtl, aes(x=as.factor(founder),y=mean.residual,color=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2) + xlab("") + ylab("Adjusted F1 atrophy\n") + theme_bw() + theme(axis.title=element_text(size=20),axis.text=element_text(size=16),legend.position="none",strip.text = element_text(size=18)) + scale_colour_manual(values=c("red3","pink1")) + facet_grid(.~qtl)
ggsave("qtl_horizontal.stackPOINT.pdf",width=10,height=5)



#install.packages("ggpubr")
library(ggpubr)

test <- ggarrange(Phasing,LODpeak)
#responder abundance plot
#you should update the responder abundance values, you can also modify this code to look at 260
#I suggest that you use estimated copy numbers here, which will be:
#(rsp left reads/length rsp left + rsp right reads/length rsp right + riX reads/ length riX) / (total aligned reads / genome size (nt))
#
resp.df <- data.frame("founder" = c("B1", "B2", "B3", "B4", "B6", "B7", "B8"), "responder" = c(780.5822008, 216.939064, 222.055646, 1038.4281720, 58.5214216, 261.44243719,1176.701673), "mean.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,mean),"se.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,se), "groups" = c("A","A","A","A","C","A","A") )
Rsp_graph <- ggplot(resp.df, aes(log(responder,2), mean.residual,colour=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2)  + scale_colour_manual(values=c("red3","pink1")) + xlab(paste("Log2 estimated Rsp copy number")) + ylab("Adjusted F1 atrophy") +theme_bw() + theme(axis.text=element_text(size=16, face = "bold"),axis.title=element_text(size=20, face= "bold"),legend.position="none") + geom_vline(xintercept = 8.525048887, colour="darkblue",linetype=2, size = 1.2) +
  annotate("text", x = 9.35, y= 0.15, label="Spearman's rho = -0.89", size= 7) + annotate("text", x = 10, y= 0.12, label=expression(paste(italic("p"), "= 0.01")), size= 7) +
  annotate("text", x = 5.9, y= 0.09, label=expression(bold("B6")), size= 7, color="pink1") + annotate("text", x = 10.2, y= -0.03, label=expression(bold("B8")), size= 7, color="red3") +
  annotate("text", x = 10, y= -0.06, label=expression(bold("B4")), size= 7, color="red3")+ annotate("text", x = 9.6, y= -0.06, label=expression(bold("B1")), size= 7, color="red3") +
  annotate("text", x = 8.05, y= 0.01, label=expression(bold("B7")), size= 7, color="red3") + annotate("text", x = 7.8, y= -0.07, label=expression(bold("B3")), size= 7, color="red3") + annotate("text", x = 7.75, y= 0.04, label=expression(bold("B2")), size= 7, color="red3") +
  annotate("text", x = 8.65, y= -0.02, label=expression(bold("Paternal Abundance")), size= 7, color="royalblue4", angle = 90)
Rsp_graph
ggsave("responder_by_atrophy.pdf",width=7.5,height=6)

#for grant
resp.df <- data.frame("founder" = c("B1", "B2", "B3", "B4", "B6", "B7", "B8"), "responder" = c(780.5822008, 216.939064, 222.055646, 1038.4281720, 58.5214216, 261.44243719,1176.701673), "mean.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,mean),"se.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,se), "groups" = c("A","A","A","A","C","A","A") )
   Rsp_graph <- ggplot(resp.df, aes(log(responder,2), mean.residual,colour=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=0.5)  + 
  scale_colour_manual(values=c("red3","pink1")) + xlab(paste("Log2 estimated Rsp copy number")) + ylab("Adjusted F1 atrophy") +theme_bw() + theme(axis.text.x=element_text(family ="Helvetica",size=9),
                                                                                                                                                  axis.text.y =element_text(family ="Helvetica",size=10),
                                                                                                                                                  axis.title=element_text(family ="Helvetica",size=9),
                                                                                                                                                  legend.position="none") +
  geom_vline(xintercept = 8.525048887, colour="darkblue",linetype=2, size = 0.3) +
  annotate("text", x = 8.5, y= 0.15, label="Spearman's rho = -0.89", size= 3) + annotate("text", x = 9.7, y= 0.12, label=expression(paste(italic("p"), "= 0.01")), size= 3) +
  annotate("text", x = 5.9, y= 0.09, label=expression(bold("B6")), size=3, color="pink1") + annotate("text", x = 10.2, y= -0.028, label=expression(bold("B8")), size=3, color="red3") +
  annotate("text", x = 10, y= -0.06, label=expression(bold("B4")), size=3, color="red3")+ annotate("text", x = 9.6, y= -0.06, label=expression(bold("B1")), size=3, color="red3") +
  annotate("text", x = 8.05, y= 0.02, label=expression(bold("B7")), size=3, color="red3") + annotate("text", x = 7.3, y= -0.046, label=expression(bold("B3")), size=3, color="red3") +
  annotate("text", x = 7.75, y= 0.05, label=expression(bold("B2")), size=3, color="red3") +
  annotate("text", x = 8.8, y= -0.02, label=expression("paternal abundance"), size=3, color="royalblue4", angle = 90)
Rsp_graph
setwd("~/Documents/RNAseq/Count_file/TE-count/unaligned-gene-counts/kallisto/abundance.tsv") 

setwd("~/Desktop/FIG/kell_grant")

ggsave("responder_by_atrophy2.pdf",width=3,height=2.5)
ggsave("responder_by_atrophyNew2.3X2.2.pdf",width=2.3,height=2.2)
ggsave("responder_by_atrophyNew2.5X2.2.pdf",width=2.5,height=2.2)

#resp.df <- data.frame("founder" = c("B1", "B2", "B3", "B4", "B6", "B7", "B8"), "responder" = c(496.5106701, 157.9767952, 160.2820984, 738.5434813, 46.5259152, 185.0740944, 892.3127056), "mean.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,mean),"se.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,se), "groups" = c("A","B","B","A","C","B","A") )
#ggplot(resp.df, aes(log(responder,2), mean.residual,colour=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2)  + scale_colour_manual(values=c("red3","violetred1","pink1")) + xlab("\nLog2 responder derived reads per million") + ylab("Adjusted F1 atrophy\n") +theme_bw() + theme(axis.text=element_text(size=20),axis.title=element_text(size=24),legend.position="none") + geom_vline(xintercept = 8.252590282, colour="darkblue",linetype=2)
#ggsave("responder_by_atrophy.pdf",width=5,height=5)


#correlation test, 260 abundance by atrophy

bp260.df <- data.frame("founder" = c("B1", "B2", "B3", "B4", "B6", "B7", "B8"), "bp260satellite" = c(826.936799, 540.0101403, 554.7590474, 551.7772497, 128.6357472, 513.6024090, 379.7321537), "mean.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,mean),"se.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,se), "groups" = c("A","A","A","A","B","A","A") )
ggplot(bp260.df, aes(log(bp260satellite,2), mean.residual,colour=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2)  + scale_colour_manual(values=c("red3","pink1")) + xlab("\nLog2 estimated 260 bp copy number") + ylab("Adjusted F1 atrophy") +theme_bw() + theme(axis.text=element_text(size=20, face = "bold"),axis.title=element_text(size=20, face = "bold"),legend.position="none") +
  annotate("text", x = 9.08, y= 0.15, label="Spearman's rho = -0.5", size= 7) +
  annotate("text", x = 9.4, y= 0.11, label=expression(paste(italic("p"), "= 0.2667")), size= 7) +
  annotate("text", x = 7, y= 0.09, label=expression(bold("B6")), size=7, color="pink1") +
  annotate("text", x = 8.3, y= -0.09, label=expression(bold("B8")), size=7, color="red3") +
  annotate("text", x = 9.5, y= -0.12, label=expression(bold("B4")), size=7, color="red3")+ 
  annotate("text", x = 9.69, y= -0.07, label=expression(bold("B1")), size=7, color="red3") +
  annotate("text", x = 8.85, y= -0.06, label=expression(bold("B7")), size=7, color="red3") + 
  annotate("text", x = 9.3, y= -0.045, label=expression(bold("B3")), size=7, color="red3") +
  annotate("text", x = 9.0, y= 0.05, label=expression(bold("B2")), size=7, color="red3") 
  
ggsave("260bp_by_atrophy.pdf",width=7.5,height=5)


cor.test(bp260.df$bp260satellite,bp260.df$mean.residual,method="spearman")
#-0.5 S = 84, p-value = 0.2667
#correlation test, responder abundance by atrophy
cor.test(bp260.df$bp260satellite,bp260.df$mean.residual,method="spearman")

colnames(resp.df) [2] <- "estimated_copy_number"
resp.df$satellite <- "Rsp"
colnames(bp260.df) [2] <- "estimated_copy_number"
bp260.df$satellite <- "260 bp"
Rsp_260.data <- rbind(resp.df,bp260.df)

p <-ggplot(Rsp_260.data, aes(log(estimated_copy_number,2), mean.residual,colour=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2)  +
  scale_colour_manual(values=c("darkcyan","chartreuse1", "chartreuse1")) + xlab("\nLog2 estimated satellite copy number") + ylab("Adjusted F1 atrophy\n") +theme_bw() + 
  theme(axis.text=element_text(size=22, face = "bold"),axis.title=element_text(size=22, face= "bold"),legend.position="none", strip.text = element_text(size=20)) +
  geom_vline(data=data.frame(xint=8.525048887,satellite="Rsp"),aes(xintercept=xint), colour="darkblue",linetype=2, size = 1.2) + facet_grid(.~satellite, scales="free_x")

p 
ann_text <- data.frame(estimated_copy_number = 9.13,mean.residual = 0.12,lab = "Spearman's rho = -0.89\n p = 0.01",satellite = factor("Rsp",levels = c("260 bp", "Rsp")))
ann_text2 <- data.frame(estimated_copy_number = 9,mean.residual = 0.12,lab = "Spearman's rho = -0.5\n  p = 0.26",
                        satellite = factor("260 bp",levels = c("260 bp", "Rsp")))
ann_text3 <- data.frame(estimated_copy_number = 5.9,mean.residual = 0.09,lab = "B6",
                        satellite = factor("Rsp",levels = c("260 bp", "Rsp")))
ann_text4 <- data.frame(estimated_copy_number = 7,mean.residual = 0.09,lab = "B6",
                        satellite = factor("260 bp",levels = c("260 bp", "Rsp")))
ann_text5 <- data.frame(estimated_copy_number = 10,mean.residual = -0.06,lab = "B4",
                        satellite = factor("Rsp",levels = c("260 bp", "Rsp")))
ann_text6 <- data.frame(estimated_copy_number = 10.2,mean.residual = -0.03,lab = "B8",
                        satellite = factor("Rsp",levels = c("260 bp", "Rsp")))
ann_text7 <- data.frame(estimated_copy_number = 9.2,mean.residual = -0.17,lab = "B4",
                        satellite = factor("260 bp",levels = c("260 bp", "Rsp")))
ann_text8 <- data.frame(estimated_copy_number = 8.6,mean.residual = -0.03,lab = "B8",
                        satellite = factor("260 bp",levels = c("260 bp", "Rsp")))
ann_text9 <- data.frame(estimated_copy_number = 8.65,mean.residual = -0.02,lab = "Paternal Abundance",
                        satellite = factor("Rsp",levels = c("260 bp", "Rsp"))) 

ann_text10 <- data.frame(estimated_copy_number = 9.6,mean.residual = -0.14,lab = "B1",
                        satellite = factor("Rsp",levels = c("260 bp", "Rsp")))
ann_text11 <- data.frame(estimated_copy_number = 7.76,mean.residual = 0.05,lab = "B2",
                        satellite = factor("Rsp",levels = c("260 bp", "Rsp")))
ann_text12 <- data.frame(estimated_copy_number = 7.70,mean.residual = -0.09,lab = "B3",
                        satellite = factor("Rsp",levels = c("260 bp", "Rsp")))
ann_text13 <- data.frame(estimated_copy_number = 8.03,mean.residual = 0.01,lab = "B7",
                        satellite = factor("Rsp",levels = c("260 bp", "Rsp")))

ann_text14 <- data.frame(estimated_copy_number = 9.69,mean.residual = -0.16,lab = "B1",
                         satellite = factor("260 bp",levels = c("260 bp", "Rsp")))
ann_text15 <- data.frame(estimated_copy_number = 9.07,mean.residual = 0.05,lab = "B2",
                         satellite = factor("260 bp",levels = c("260 bp", "Rsp")))
ann_text16 <- data.frame(estimated_copy_number = 9.3,mean.residual = -0.041,lab = "B3",
                         satellite = factor("260 bp",levels = c("260 bp", "Rsp")))
ann_text17 <- data.frame(estimated_copy_number = 8.85,mean.residual = -0.05,lab = "B7",
                         satellite = factor("260 bp",levels = c("260 bp", "Rsp")))




p + geom_text(aes(estimated_copy_number,mean.residual, label= lab),
              data = ann_text, size = 7, inherit.aes = FALSE ) + geom_text(aes(estimated_copy_number,mean.residual, label= lab),
              data = ann_text2, size = 7, inherit.aes = FALSE ) + geom_text(aes(estimated_copy_number,mean.residual, label= lab),
              data = ann_text3, size = 7, color = "chartreuse1" ) + geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text4, size = 7, color = "chartreuse1" ) +
  geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text5, size = 7, color = "darkcyan" ) +   geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text6, size = 7, color = "darkcyan" ) +
  geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text7, size = 7, color = "darkcyan" ) +
  geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text8, size = 7, color = "darkcyan" ) + geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text9, size = 7, color = "royalblue4", angle=90 )+
  geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text10, size = 7, color = "darkcyan" ) + geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text11, size = 7, color = "darkcyan" )+
  geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text12, size = 7, color = "darkcyan" ) +geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text13, size = 7, color = "darkcyan" ) +
  geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text14, size = 7, color = "darkcyan" ) +geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text15, size = 7, color = "darkcyan" ) +
  geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text16, size = 7, color = "darkcyan" ) +geom_text(aes(estimated_copy_number,mean.residual, label= lab),data = ann_text17, size = 7, color = "darkcyan" )


ggsave("responder_260_by_atrophy.pdf",width=13,height=6)

