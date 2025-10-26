library(ggplot2)
library(plyr)

setwd("~/Downloads")

load("methylation_data.rda") 

#This includes meth.df, the data frame with the input data which is used to construct clocks
#head(meth.df)
#  Chromosome Position Methylation_level Sample Age       Locus
#      Chr_7    23000                 0   AL11   5 Chr_7_23000
#      Chr_7    23007                 0   AL11   5 Chr_7_23007
#      Chr_7    23014                 0   AL11   5 Chr_7_23014
#      Chr_7    23020                 0   AL11   5 Chr_7_23020
#      Chr_7    23035                 0   AL11   5 Chr_7_23035
#      Chr_7    23039                 0   AL11   5 Chr_7_23039

mean_per_sample<-plyr::ddply(meth.df, c("Age", "Sample"), function(x) mean(x$Methylation_level, na.rm=T) )

colnames(mean_per_sample) <- c("Age", "Sample", "Avg")
mean_per_sample$Age <- ifelse(mean_per_sample$Age == "meconium", "ME", mean_per_sample$Age)
mean_per_sample$Age <- ifelse(mean_per_sample$Age == "long-distance", "LD", mean_per_sample$Age)
mean_per_sample$Age <- factor(mean_per_sample$Age , levels=c("0", "5", "10", "15", "20", "25", "30", "ME", "LD"))

mean_per_sample$Experiment <-ifelse(grepl("AL", mean_per_sample$Sample) | grepl("LI", mean_per_sample$Sample), "E2", "E1")
mean_per_sample$Experiment <-ifelse(grepl("ME", mean_per_sample$Age) | grepl("LD", mean_per_sample$Age), "Wild", mean_per_sample$Experiment)

N<-as.data.frame(table(meth.df$Sample))
mean_per_sample$margins<-qt(0.975,df=N$Freq-1)*(mean_per_sample$Avg*100)/sqrt(N$Freq)

mean_per_sample$Sample_type <- ifelse(mean_per_sample$Experiment == "Wild", "Wild", "Experimental")

ggplot(mean_per_sample, aes(x=Age, y=Avg*100, col=Experiment))+geom_point(size=4, position=position_dodge2(width=0.1), alpha=0.6)+
  geom_linerange(data=mean_per_sample, aes(y=Avg*100, ymin=(Avg*100)-margins, ymax=(Avg*100)+margins), size=1, alpha=1, col="black",position=position_dodge2(width=0.1))+
  theme_bw()+
  ylab("CDS Methylation (%)")+
  xlab("Age (Days)")+
  scale_color_manual(values=c("#44AA99", "#882255", "orange"), name="")+
  facet_grid(~Sample_type,  scales = "free_x", space = "free")+
  theme(strip.text = element_blank(), element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))

