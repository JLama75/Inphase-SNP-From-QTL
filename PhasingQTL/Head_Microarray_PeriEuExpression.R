setwd("~/Documents/HMM_R4_QTL/QTL2")
library(tidyverse)
#Pericentromeric genes
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
#ggplot(BB1Peri_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB2_rowCountPeri <- BB2_new %>% mutate(m=rowMeans(BB2_new[,5:579]))
BB2_rowCountPeri <- BB2_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB2Peri_plot <- data.frame("founder" = c("BB2"), "mean.expression" = mean(BB2_rowCountPeri$m),"se.expression" = se(BB2_rowCountPeri$m))
#ggplot(BB2Peri_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB3_rowCountPeri <- BB3_new %>% mutate(m=rowMeans(BB3_new[,5:579]))
BB3_rowCountPeri <- BB3_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB3Peri_plot <- data.frame("founder" = c("BB3"), "mean.expression" = mean(BB3_rowCountPeri$m),"se.expression" = se(BB3_rowCountPeri$m))
#ggplot(BB3Peri_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB4_rowCountPeri <- BB4_new %>% mutate(m=rowMeans(BB4_new[,5:579]))
BB4_rowCountPeri <- BB4_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB4Peri_plot <- data.frame("founder" = c("BB4"), "mean.expression" = mean(BB4_rowCountPeri$m),"se.expression" = se(BB4_rowCountPeri$m))
#ggplot(BB4Peri_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB6_rowCountPeri <- BB6_new %>% mutate(m=rowMeans(BB6_new[,5:579]))
BB6_rowCountPeri <- BB6_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB6Peri_plot <- data.frame("founder" = c("BB6"), "mean.expression" = mean(BB6_rowCountPeri$m),"se.expression" = se(BB6_rowCountPeri$m))
#ggplot(BB6Peri_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB7_rowCountPeri <- BB7_new %>% mutate(m=rowMeans(BB7_new[,5:579]))
BB7_rowCountPeri <- BB7_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB7Peri_plot <- data.frame("founder" = c("BB7"), "mean.expression" = mean(BB7_rowCountPeri$m),"se.expression" = se(BB7_rowCountPeri$m))
#ggplot(BB7Peri_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB8_rowCountPeri <- BB8_new %>% mutate(m=rowMeans(BB8_new[,5:579]))
BB8_rowCountPeri <- BB8_rowCountPeri[,c(2:4,580)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB8Peri_plot <- data.frame("founder" = c("BB8"), "mean.expression" = mean(BB8_rowCountPeri$m),"se.expression" = se(BB8_rowCountPeri$m))
#ggplot(BB8Peri_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

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

library(ggplot2)

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


#Euchromatic genes
RIL6_euchrom <- read.csv("RIL6_euchrom.csv", header = T, sep=",")
RIL6_euchrom_cp <- RIL6_euchrom[1:3,1:10]
RIL6_euchrom <- RIL6_euchrom[,-1]

BB1_euchrom <- subset(RIL6_euchrom, RIL6_euchrom$founder == "BB1")
BB2_euchrom <- subset(RIL6_euchrom, RIL6_euchrom$founder == "BB2")
BB3_euchrom <- subset(RIL6_euchrom, RIL6_euchrom$founder == "BB3")
BB4_euchrom <- subset(RIL6_euchrom, RIL6_euchrom$founder == "BB4")
BB5_euchrom <- subset(RIL6_euchrom, RIL6_euchrom$founder == "BB5")
BB6_euchrom <- subset(RIL6_euchrom, RIL6_euchrom$founder == "BB6")
BB7_euchrom <- subset(RIL6_euchrom, RIL6_euchrom$founder == "BB7")
BB8_euchrom <- subset(RIL6_euchrom, RIL6_euchrom$founder == "BB8")

#obtaining the mean of all pericentromeric gene expression for each RILs with BB1 QTL2 allele.
BB1_euchrom <- BB1_euchrom %>% mutate(m=rowMeans(BB1_euchrom[,4:10481]))
BB1_euchrom <- BB1_euchrom[,c(1:3,10482)]

#now taking the mean of all RIL expression and calculating SE.
se <- function(x) {sd(x)/sqrt(length(x))}
BB1EU_plot <- data.frame("founder" = c("BB1"), "mean.expression" = mean(BB1_euchrom$m),"se.expression" = se(BB1_euchrom$m))
#ggplot(BB1EU_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB2_euchrom <- BB2_euchrom %>% mutate(m=rowMeans(BB2_euchrom[,4:10481]))
BB2_euchrom <- BB2_euchrom[,c(1:3,10482)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB2EU_plot <- data.frame("founder" = c("BB2"), "mean.expression" = mean(BB2_euchrom$m),"se.expression" = se(BB2_euchrom$m))
#ggplot(BB2EU_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB3_euchrom <- BB3_euchrom %>% mutate(m=rowMeans(BB3_euchrom[,4:10481]))
BB3_euchrom <- BB3_euchrom[,c(1:3,10482)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB3EU_plot <- data.frame("founder" = c("BB3"), "mean.expression" = mean(BB3_euchrom$m),"se.expression" = se(BB3_euchrom$m))
#ggplot(BB3EU_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB4_euchrom <- BB4_euchrom %>% mutate(m=rowMeans(BB4_euchrom[,4:10481]))
BB4_euchrom <- BB4_euchrom[,c(1:3,10482)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB4EU_plot <- data.frame("founder" = c("BB4"), "mean.expression" = mean(BB4_euchrom$m),"se.expression" = se(BB4_euchrom$m))
#ggplot(BB4EU_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB6_euchrom <- BB6_euchrom %>% mutate(m=rowMeans(BB6_euchrom[,4:10481]))
BB6_euchrom <- BB6_euchrom[,c(1:3,10482)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB6EU_plot <- data.frame("founder" = c("BB6"), "mean.expression" = mean(BB6_euchrom$m),"se.expression" = se(BB6_euchrom$m))
#ggplot(BB6EU_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB7_euchrom <- BB7_euchrom %>% mutate(m=rowMeans(BB7_euchrom[,4:10481]))
BB7_euchrom <- BB7_euchrom[,c(1:3,10482)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB7EU_plot <- data.frame("founder" = c("BB7"), "mean.expression" = mean(BB7_euchrom$m),"se.expression" = se(BB7_euchrom$m))
#ggplot(BB7EU_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB8_euchrom <- BB8_euchrom %>% mutate(m=rowMeans(BB8_euchrom[,4:10481]))
BB8_euchrom <- BB8_euchrom[,c(1:3,10482)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB8EU_plot <- data.frame("founder" = c("BB8"), "mean.expression" = mean(BB8_euchrom$m),"se.expression" = se(BB8_euchrom$m))
#ggplot(BB8EU_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

EU_plot <- rbind(BB1EU_plot,BB2EU_plot,BB3EU_plot,BB4EU_plot,BB6EU_plot,BB7EU_plot,BB8EU_plot)
EU_plot$bais <- c("A","A","A","A","C","A","A")
EU_plot$founder <- c("B1", "B2","B3","B4","B6","B7","B8")


#Combining all the founders to perform anova test
rowCount_euchrom <- rbind(BB1_euchrom,BB2_euchrom,BB3_euchrom,BB4_euchrom,BB6_euchrom,BB7_euchrom,BB8_euchrom) 


resAnova_euchrom <- aov(m~founder, rowCount_euchrom)
summary(resAnova_euchrom)
#F = 0.208 P= 0.974 ns   
# Compare with alternative method by calculating mean RILs instead of mean genes: 
R2 <- 6/(6+494)
R2 #0.012
#0.2% variation in expression explained by founders
TukeyHSD(resAnova_euchrom,conf.level=0.95)

linearMod_euchrom<- lm(m~founder, rowCount_euchrom)
summary(linearMod_euchrom)
#anova(linearMod_euchrom)


g <- ggplot(EU_plot, aes(as.factor(founder),mean.expression, color=bais)) +
  geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2) + xlab("") +
  ylab("Adjusted Expression\n") + theme_bw() +  mynamestheme + theme(legend.position="none") + 
  scale_colour_manual(values=c("red3","pink1")) + geom_hline(aes(yintercept=0), alpha=.5, size=1) + ylim(-0.25, 0.2)

g


#4th chromosome
setwd("~/Documents/HMM_R4_QTL/QTL2")
RIL5_fourth <- read.csv("RIL5fourth.csv", sep=",", header = T) 

BB1_fourth <- subset(RIL5_fourth, RIL5_fourth$founder == "BB1")
BB2_fourth <- subset(RIL5_fourth, RIL5_fourth$founder == "BB2")
BB3_fourth <- subset(RIL5_fourth, RIL5_fourth$founder == "BB3")
BB4_fourth <- subset(RIL5_fourth, RIL5_fourth$founder == "BB4")
BB5_fourth <- subset(RIL5_fourth, RIL5_fourth$founder == "BB5")
BB6_fourth <- subset(RIL5_fourth, RIL5_fourth$founder == "BB6")
BB7_fourth <- subset(RIL5_fourth, RIL5_fourth$founder == "BB7")
BB8_fourth <- subset(RIL5_fourth, RIL5_fourth$founder == "BB8")


#obtaining the mean of all pericentromeric gene expression for each RILs with BB1 QTL2 allele.
BB1_rowCount4 <- BB1_fourth %>% mutate(m=rowMeans(BB1_fourth[,5:139]))
BB1_rowCount4 <- BB1_rowCount4[,c(2:4,140)]

#now taking the mean of all RIL expression and calculating SE.
se <- function(x) {sd(x)/sqrt(length(x))}
BB1fourth_plot <- data.frame("founder" = c("BB1"), "mean.expression" = mean(BB1_rowCount4$m),"se.expression" = se(BB1_rowCount4$m))
ggplot(BB1fourth_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB2_rowCount4 <- BB2_fourth %>% mutate(m=rowMeans(BB2_fourth[,5:139]))
BB2_rowCount4 <- BB2_rowCount4[,c(2:4,140)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB2fourth_plot <- data.frame("founder" = c("BB2"), "mean.expression" = mean(BB2_rowCount4$m),"se.expression" = se(BB2_rowCount4$m))
#ggplot(BB2fourth_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB3_rowCount4 <- BB3_fourth %>% mutate(m=rowMeans(BB3_fourth[,5:139]))
BB3_rowCount4 <- BB3_rowCount4[,c(2:4,140)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB3fourth_plot <- data.frame("founder" = c("BB3"), "mean.expression" = mean(BB3_rowCount4$m),"se.expression" = se(BB3_rowCount4$m))
#ggplot(BB3fourth_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB4_rowCount4 <- BB4_fourth %>% mutate(m=rowMeans(BB4_fourth[,5:139]))
BB4_rowCount4 <- BB4_rowCount4[,c(2:4,140)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB4fourth_plot <- data.frame("founder" = c("BB4"), "mean.expression" = mean(BB4_rowCount4$m),"se.expression" = se(BB4_rowCount4$m))
#ggplot(BB4fourth_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB6_rowCount4 <- BB6_fourth %>% mutate(m=rowMeans(BB6_fourth[,5:139]))
BB6_rowCount4 <- BB6_rowCount4[,c(2:4,140)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB6fourth_plot <- data.frame("founder" = c("BB6"), "mean.expression" = mean(BB6_rowCount4$m),"se.expression" = se(BB6_rowCount4$m))
#ggplot(BB6fourth_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB7_rowCount4 <- BB7_fourth %>% mutate(m=rowMeans(BB7_fourth[,5:139]))
BB7_rowCount4 <- BB7_rowCount4[,c(2:4,140)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB7fourth_plot <- data.frame("founder" = c("BB7"), "mean.expression" = mean(BB7_rowCount4$m),"se.expression" = se(BB7_rowCount4$m))
#ggplot(BB7fourth_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB8_rowCount4 <- BB8_fourth %>% mutate(m=rowMeans(BB8_fourth[,5:139]))
BB8_rowCount4 <- BB8_rowCount4[,c(2:4,140)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB8fourth_plot <- data.frame("founder" = c("BB8"), "mean.expression" = mean(BB8_rowCount4$m),"se.expression" = se(BB8_rowCount4$m))
#ggplot(BB8fourth_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

fourth_plot <- rbind(BB1fourth_plot,BB2fourth_plot,BB3fourth_plot,BB4fourth_plot,BB6fourth_plot,BB7fourth_plot,BB8fourth_plot)
fourth_plot$bais <- c("A","A","A","A","C","A","A")
fourth_plot$founder <- c("B1", "B2","B3","B4","B6","B7","B8")

ggplot(fourth_plot, aes(as.factor(founder),mean.expression, color=bais)) + 
  geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+
                        se.expression),size=2) + xlab("") + 
  ylab("Adjusted Expression\n") + theme_bw() +  mynamestheme + 
  theme(legend.position="none") + scale_colour_manual(values=c("red3","pink1")) +
  ylim(-0.25,0.2)+ geom_hline(aes(yintercept=0), alpha=.5, size=1)

#Combining all the founders
rowCount4 <- rbind(BB1_rowCount4,BB2_rowCount4,BB3_rowCount4,BB4_rowCount4,BB6_rowCount4,BB7_rowCount4,BB8_rowCount4) 


resAnova_fourth_rowCount4 <- aov(m~founder, rowCount4)
summary(resAnova_fourth_rowCount4)

#  F=1.047 P<0.394 ns
R2 <- 1.15/(1.15+90.66)
R2 #0.1038795
#1% variation in expression explained by founders
TukeyHSD(resAnova_fourth_rowCount4,conf.level=0.95)

linearMod_eu_fourth_rowCount4<- lm(m~founder, rowCount4)
summary(linearMod_eu_fourth_rowCount4)
anova(linearMod_eu_fourth_rowCount4)

#telomere
setwd("~/Documents/HMM_R4_QTL/QTL2")
RIL5_tel <- read_csv("RIL5_tel.csv") #69 GENES


BB1_tel <- subset(RIL5_tel, RIL5_tel$founder == "BB1")
BB2_tel <- subset(RIL5_tel, RIL5_tel$founder == "BB2")
BB3_tel <- subset(RIL5_tel, RIL5_tel$founder == "BB3")
BB4_tel <- subset(RIL5_tel, RIL5_tel$founder == "BB4")
BB5_tel <- subset(RIL5_tel, RIL5_tel$founder == "BB5")
BB6_tel <- subset(RIL5_tel, RIL5_tel$founder == "BB6")
BB7_tel <- subset(RIL5_tel, RIL5_tel$founder == "BB7")
BB8_tel <- subset(RIL5_tel, RIL5_tel$founder == "BB8")


#obtaining the mean of all pericentromeric gene expression for each RILs with BB1 QTL2 allele.
BB1_rowCount_tel <- BB1_tel %>% mutate(m=rowMeans(BB1_tel[,5:70]))
BB1_rowCount_tel <- BB1_rowCount_tel[,c(2:4,71)]

#now taking the mean of all RIL expression and calculating SE.
se <- function(x) {sd(x)/sqrt(length(x))}
BB1tel_plot <- data.frame("founder" = c("BB1"), "mean.expression" = mean(BB1_rowCount_tel$m),"se.expression" = se(BB1_rowCount_tel$m))
ggplot(BB1tel_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB2_rowCount_tel <- BB2_tel %>% mutate(m=rowMeans(BB2_tel[,5:70]))
BB2_rowCount_tel <- BB2_rowCount_tel[,c(2:4,71)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB2tel_plot <- data.frame("founder" = c("BB2"), "mean.expression" = mean(BB2_rowCount_tel$m),"se.expression" = se(BB2_rowCount_tel$m))
#ggplot(BB2tel_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB3_rowCount_tel <- BB3_tel %>% mutate(m=rowMeans(BB3_tel[,5:70]))
BB3_rowCount_tel <- BB3_rowCount_tel[,c(2:4,71)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB3tel_plot <- data.frame("founder" = c("BB3"), "mean.expression" = mean(BB3_rowCount_tel$m),"se.expression" = se(BB3_rowCount_tel$m))
#ggplot(BB3tel_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)


BB4_rowCount_tel <- BB4_tel %>% mutate(m=rowMeans(BB4_tel[,5:70]))
BB4_rowCount_tel <- BB4_rowCount_tel[,c(2:4,71)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB4tel_plot <- data.frame("founder" = c("BB4"), "mean.expression" = mean(BB4_rowCount_tel$m),"se.expression" = se(BB4_rowCount_tel$m))
#ggplot(BB4tel_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB6_rowCount_tel <- BB6_tel %>% mutate(m=rowMeans(BB6_tel[,5:70]))
BB6_rowCount_tel <- BB6_rowCount_tel[,c(2:4,71)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB6tel_plot <- data.frame("founder" = c("BB6"), "mean.expression" = mean(BB6_rowCount_tel$m),"se.expression" = se(BB6_rowCount_tel$m))
#ggplot(BB6tel_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB7_rowCount_tel <- BB7_tel %>% mutate(m=rowMeans(BB7_tel[,5:70]))
BB7_rowCount_tel <- BB7_rowCount_tel[,c(2:4,71)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB7tel_plot <- data.frame("founder" = c("BB7"), "mean.expression" = mean(BB7_rowCount_tel$m),"se.expression" = se(BB7_rowCount_tel$m))
#ggplot(BB7tel_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

BB8_rowCount_tel <- BB8_tel %>% mutate(m=rowMeans(BB8_tel[,5:70]))
BB8_rowCount_tel <- BB8_rowCount_tel[,c(2:4,71)]
se <- function(x) {sd(x)/sqrt(length(x))}
BB8tel_plot <- data.frame("founder" = c("BB8"), "mean.expression" = mean(BB8_rowCount_tel$m),"se.expression" = se(BB8_rowCount_tel$m))
#ggplot(BB8tel_plot, aes(as.factor(founder),mean.expression)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2)

tel_plot <- rbind(BB1tel_plot,BB2tel_plot,BB3tel_plot,BB4tel_plot,BB6tel_plot,BB7tel_plot,BB8tel_plot)
tel_plot$bais <- c("A","A","A","A","C","A","A")
tel_plot$founder <- c("B1", "B2","B3","B4","B6","B7","B8")

ggplot(tel_plot, aes(as.factor(founder),mean.expression, color=bais)) + 
  geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+
                        se.expression),size=2) + xlab("") + 
  ylab("Adjusted Expression\n") + theme_bw() +  mynamestheme + 
  theme(legend.position="none") + scale_colour_manual(values=c("red3","pink1")) +
  ylim(-0.25,0.2)+ geom_hline(aes(yintercept=0), alpha=.5, size=1)

#Combining all the founders
rowCount_tel <- rbind(BB1_rowCount_tel,BB2_rowCount_tel,BB3_rowCount_tel,BB4_rowCount_tel,BB6_rowCount_tel,BB7_rowCount_tel,BB8_rowCount_tel) 


resAnova_tel_rowCount_tel <- aov(m~founder, rowCount_tel)
summary(resAnova_tel_rowCount_tel)

#  F=0.641 P<0.698 ns
R2 <- 6/(6+494)
R2 #0.012
#1% variation in expression explained by founders
TukeyHSD(resAnova_tel_rowCount_tel,conf.level=0.95)

linearMod_eu_tel_rowCount_tel<- lm(m~founder, rowCount)
summary(linearMod_eu_tel_rowCount_tel)
anova(linearMod_eu)




#COMBINED PLOT

EU_plot$arm <- "euchromatin"
fourth_plot$arm <- "fourth chromsome"
Peri_plot$arm <- "pericentromere"
tel_plot$arm <- "telomere"

Plot_arm <- rbind(Peri_plot,EU_plot,tel_plot,fourth_plot)  


Plot_arm$arm <- factor(Plot_arm$arm, levels = c("pericentromere","euchromatin", "telomere","fourth chromsome"))
levels(as.factor(Plot_arm$arm)) 

plotRIL <- ggplot(Plot_arm, aes(x=as.factor(founder),y=mean.expression,color=bais)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2) +
  xlab("") + ylab("Adjusted Expression\n") + theme_bw() +  mynamestheme + theme(legend.position="none",strip.text = element_text(size=20)) + 
  scale_colour_manual(values=c("darkcyan","green3"))  + geom_hline(aes(yintercept=0), alpha=.5, size=1)  + facet_grid(.~arm, scales="free_x")
plotRIL





#correlation
setwd("~/Documents/QTLphasing_responder_260")
conc.2L.20770000.dat <- read.csv("conc.2L.20770000.dat.csv", sep=",", header = T)
correlation.df <- data.frame("founder" = c("B1", "B2", "B3", "B4", "B6", "B7", "B8"), "peri_expression" = Peri_plot$mean.expression, "mean.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,mean),"se.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,se), "groups" = c("A","A","A","A","C","A","A") )

ggplot(correlation.df, aes(peri_expression, mean.residual,colour=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2)  + scale_colour_manual(values=c("red3","pink1")) + xlab(paste("Adjusted head expression")) + ylab("Adjusted F1 atrophy") +theme_bw() + theme(axis.text=element_text(size=16, face = "bold"),axis.title=element_text(size=20, face= "bold"),legend.position="none") + 
annotate("text", x = -0.11, y= 0.15, label="Spearman's rho = 0.75", size= 7) +
  annotate("text", x = -0.15, y= 0.11, label=expression(paste(italic("p"), "= 0.06627")), size= 7) 

cor.test(correlation.df$peri_expression,correlation.df$mean.residual,method="spearman")  



#Combining all the founders
rowCount <- rbind(BB1_rowCount4,BB2_rowCount4,BB3_rowCount4,BB4_rowCount4,BB6_rowCount4,BB7_rowCount4,BB8_rowCount4) 


resAnova_fourth_rowCount4 <- aov(m~founder, rowCount)
summary(resAnova_fourth_rowCount4)

#  F=1.047 P<0.394 ns
R2 <- 1.15/(1.15+90.66)
R2 #0.1038795
#1% variation in expression explained by founders
TukeyHSD(resAnova_fourth_rowCount4,conf.level=0.95)

linearMod_eu_fourth_rowCount4<- lm(m~founder, rowCount)
summary(linearMod_eu_fourth_rowCount4)
anova(linearMod_eu)

#COMBINED PLOT

EU_plot$arm <- "euchromatin"
fourth_plot$arm <- "fourth chromsome"
Peri_plot$arm <- "pericentromere"
Plot_arm <- rbind(EU_plot,Peri_plot)  


Plot_arm$arm <- factor(Plot_arm$arm, levels = c("pericentromere","euchromatin"))
levels(as.factor(Plot_arm$arm)) 

ggplot(Plot_arm, aes(x=as.factor(founder),y=mean.expression,color=bais)) + geom_pointrange(aes(ymin=mean.expression-se.expression,ymax=mean.expression+se.expression),size=2) +
  xlab("") + ylab("Adjusted Expression\n") + theme_bw() +  mynamestheme + theme(legend.position="none",strip.text = element_text(size=20)) + 
  scale_colour_manual(values=c("darkcyan","green3"))  + geom_hline(aes(yintercept=0), alpha=.5, size=1)  + facet_grid(.~arm, scales="free_x")

#correlation
setwd("~/Documents/QTLphasing_responder_260")
conc.2L.20770000.dat <- read.csv("conc.2L.20770000.dat.csv", sep=",", header = T)
correlation.df <- data.frame("founder" = c("B1", "B2", "B3", "B4", "B6", "B7", "B8"), "peri_expression" = Peri_plot$mean.expression, "mean.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,mean),"se.residual" = tapply(conc.2L.20770000.dat$res,conc.2L.20770000.dat$founder,se), "groups" = c("A","A","A","A","C","A","A") )

ggplot(correlation.df, aes(peri_expression, mean.residual,colour=groups)) + geom_pointrange(aes(ymin=mean.residual-se.residual,ymax=mean.residual+se.residual),size=2)  + scale_colour_manual(values=c("red3","pink1")) + xlab(paste("Adjusted head expression")) + ylab("Adjusted F1 atrophy") +theme_bw() + theme(axis.text=element_text(size=16, face = "bold"),axis.title=element_text(size=20, face= "bold"),legend.position="none") + 
  annotate("text", x = -0.11, y= 0.15, label="Spearman's rho = 0.75", size= 7) +
  annotate("text", x = -0.15, y= 0.11, label=expression(paste(italic("p"), "= 0.06627")), size= 7) 

cor.test(correlation.df$peri_expression,correlation.df$mean.residual,method="spearman")  
