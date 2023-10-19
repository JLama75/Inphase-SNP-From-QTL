#Proportion of differentially enriched genes in QTL and rest of the genome
library(ggsignif)
library(ggplot2)

Total_QTL <- 572
TotalDiffExp <- 603
QTL_diffexp <- 88
QTL_notDiffexp <- 572-88
Diffexp_restgenome <- 603-88
Totalgenes <- 17773
NotDiffexp_restgenome <- 17773-572-515
QTL_enrich <- matrix(c(88,484,515,16686),nrow=2)
rownames(QTL_enrich) <- c("QTL" , "rest_of_genome")
colnames(QTL_enrich) <- c("diffExp", "Not_DiffExp")
#take a look
QTL_enrich
#we can do chisq.test(matrix,correct=FALSE) 
#correct = FALSE indicates no continuity correction, which can be overly conservative Yates continuity correction-- false negative
chisq.test(QTL_enrich)
chisq.test(QTL_enrich, correct = FALSE)

mynamestheme <- theme(plot.title = element_text(family = "Helvetica", size = (20), hjust = 0.5), 
                      legend.title = element_blank(), 
                      legend.text = element_text( colour="black",family = "Helvetica", size = (18)), 
                      axis.title = element_text(family = "Helvetica", size = (24), colour = "black"),
                      axis.text = element_text(family = "Helvetica", colour = "black", size = (20)))


QTL_enrich_graph <- data.frame("position" = c("QTL", "rest of genome"), "proportion"= c(((88/572)), ((515/17201))))

ggplot(QTL_enrich_graph, aes(position,proportion)) + geom_bar(stat = "identity", width = 0.5) +
  geom_hline(aes(yintercept=((603/17773))), size=1, linetype= "dashed") + theme_bw() +
  scale_fill_manual(values =c("gray", "gray"))  +  mynamestheme +  ylab("proportion of DEG\n") + xlab("")+

  geom_signif(y_position = c(0.16), xmin = c(1), xmax = c(2), annotation = c("***"), textsize =c(10), size = c(0.8), vjust = c(0.6)) +
  ylab("proportion of differentially\n expressed genes\n") + xlab("")

  #(88/572)*100=15.38462 ((515/17201))*100)=2.994012
  #FC = 15.38462/2.994012 = 5.1384
