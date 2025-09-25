library(ggplot2)
library(lmodel2)
library(plyr)
library(ggpubr)
setwd("~/Downloads")


# Wild samples ####
varseed.df <- read.table(file="df.out.1-1000.bigdata.comb.10.noW.MELD_DMRs_Only")
#varseed.df <- read.table(file="df.out.varying_seed_dinuc_data.comb.10.noW.MELD_DMRs_Only")
#varseed.df <- read.table(file="df.out.1-1000.data.comb.10.noW.MELD_DMRs_Only.Full_Lifespan")
#varseed.df <- read.table(file="df.out.1-1000.varying_seed_GT_CDS_dinuc.Full_Lifespan")


colnames(varseed.df) <- c("Sample", "Pred_age", "Age", "Seed", "eltest.p", "eltest.intercept", "eltest.slope", "eltest.adj.r.squared", "eltest.lm.p")

ggplot(varseed.df, aes(x=eltest.p, y=Pred_age, col=as.factor(Age)))+geom_point()

sub.df<-varseed.df[varseed.df$eltest.p <= 0.05,]
test<-t.test(sub.df[sub.df$Age == "LD",]$Pred_age, sub.df[sub.df$Age == "ME",]$Pred_age)

summary(lm(Pred_age~Age+Seed, data=varseed.df ))

diff_per_seed<-plyr::ddply(varseed.df, c("Seed"), function(x) mean(x[x$Age == "LD",]$Pred_age)-mean(x[x$Age == "ME",]$Pred_age))
colnames(diff_per_seed)<-c("Seed", "Diff_of_means_pred_age")

w.testolas<-function(x){ 
  test=wilcox.test(x[x$Age == "LD",]$Pred_age, x[x$Age == "ME",]$Pred_age) 
  return(test$p.value)
}

w.test_per_seed<-plyr::ddply(varseed.df, c("Seed"), w.testolas)

diff_per_seed$w.test.p <- w.test_per_seed$V1

ggplot(diff_per_seed, aes(x=Diff_of_means_pred_age))+geom_histogram(colour="red", fill="orange", alpha=0.3)+
  geom_histogram(data=diff_per_seed[diff_per_seed$w.test.p < 0.05,], aes(x=Diff_of_means_pred_age), colour="red", fill="orange")+
  geom_vline(xintercept=0)+
  theme_bw()+
  ylab("Count")+
  xlab("Difference")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))

dim(diff_per_seed[diff_per_seed$w.test.p < 0.05,])
dim(diff_per_seed[diff_per_seed$w.test.p < 0.05 & diff_per_seed$Diff_of_means_pred_age>0,])

sub.df<-varseed.df[varseed.df$eltest.p <1,]

mpa<-plyr::ddply(varseed.df, c("Seed"), function(x) anova(lm(Pred_age~Age, x) )$`Pr(>F)`)

table(mpa$V1<0.05)

varseed.df$anova.p <- NA
for(i in 1:length(mpa$Seed)){
  matches <- mpa$Seed[i] == varseed.df$Seed
  m_list <- (1:length(matches))[matches]
  varseed.df$anova.p[m_list] <- mpa$V1[i]
}


mps<-plyr::ddply(sub.df, c("Sample", "Age"), function(x) mean(x$Pred_age, na.rm=T) )
colnames(mps) <- c("Sample", "Age", "Pred_age")

N<-1000

margins<-qt(0.975,df=N-1)*(mps$Pred_age)/sqrt(N)


mps$margins<-margins

sample.map<-data.frame(num=c(54:63), name=c("GTcoll18B720", "GTcoll18B721", "GTcoll18B773", "GTcoll19H116", "GTcoll19H128", "GTcoll19M359", "GTcoll19M448", "GTcoll19M450", "GTcoll19M460", "GTcoll19M471"))

mps$Sample<-sample.map$name
mps$Lower95CI<-mps$Pred_age-mps$margins
mps$Upper95CI<-mps$Pred_age+mps$margins



ggplot(mps, aes(x=rev(as.factor(Sample)), y=Pred_age, col=Age )) +
  geom_pointrange(mapping=aes(y=Pred_age, ymin=Lower95CI, ymax=Upper95CI))+
  ylab("W.clock age (Days)")+
  xlab("Sample")+
  theme_bw()+
  ylim(0,17)+
  scale_colour_manual(values=c("#1E88E5", "#D81B60"), name="Group")+
 # theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text.y=element_text(size=16, colour="black"), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text.y=element_text(size=16, colour="black"), axis.text.x = element_text(angle = 90, vjust = 0, hjust=0), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))

 # ylim(0,25)



ggplot(varseed.df, aes(x=Sample, y=Pred_age, group=Seed))+geom_line(alpha=0.1, col="blue")+
  geom_line(data=varseed.df[varseed.df$eltest.lm.p < 0.05 & abs(varseed.df$eltest.intercept) < 0.2 & varseed.df$eltest.slope > 0.9 & varseed.df$eltest.slope < 1.1,], aes(x=Sample, y=Pred_age, group=Seed), alpha=0.7, col="red" )#+
#  geom_line(data=varseed.df[varseed.df$Seed == 238,], aes(x=Sample, y=Pred_age, group=Seed), alpha=1, col="red", lwd=2 )


varseed.df2 <- read.table(file="df.out.1-1000.bigdata.comb.10.noW.MELD_DMRs_Only")
colnames(varseed.df2) <- c("Sample", "Pred_age", "Age", "Seed", "eltest.p", "eltest.intercept", "eltest.slope", "eltest.adj.r.squared", "eltest.lm.p")




varseed.df$Sample<-sample.map$name
varseed.df2$Sample<-sample.map$name

ggplot(mps, aes(x=Sample, y=Pred_age, col=Age )) +
  geom_line(data=varseed.df, aes(x=Sample, y=Pred_age, group=Seed), alpha=0.05, col="blue")+
  geom_line(data=varseed.df2, aes(x=Sample, y=Pred_age, group=Seed), alpha=0.05, col="turquoise2")+
  
  #geom_pointrange(mapping=aes(y=Pred_age, ymin=Pred_age-margins, ymax=Pred_age+margins))+
  theme_bw()+

  ylab("Predicted age")+
  scale_colour_manual(values=c("#FFC107", "#D81B60"), name="Group")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text.y=element_text(size=16, colour="black"), axis.text.x = element_text(angle = 90, vjust = 0, hjust=0), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))



#Significantly "better" clocks for those with significant ANOVA between ME and LD
t.test(varseed.df[varseed.df$anova.p < 0.05,]$eltest.p, varseed.df[varseed.df$anova.p >= 0.05,]$eltest.p)

t.test(varseed.df[varseed.df$anova.p < 0.05 & varseed.df$Age == "LD",]$Pred_age, varseed.df[varseed.df$anova.p < 0.05 & varseed.df$Age == "ME",]$Pred_age, paired=T)


good.lab.clock <- varseed.df[varseed.df$eltest.lm.p < 0.05 & abs(varseed.df$eltest.intercept) < 0.2 & varseed.df$eltest.slope > 0.9 & varseed.df$eltest.slope < 1.1,]

t.test(good.lab.clock[good.lab.clock$Age == "LD",]$Pred_age, good.lab.clock[good.lab.clock$Age == "ME",]$Pred_age)



#### Migration distance versus age ####



w.clock.dist.df<-data.frame(samples=unique(mps$Sample),
           dist.min=c(2790,4250,2450,2570,3310, rep(0, 5)), dist.max=c(7320, 7160, 6520, 6330, 6530, rep(150, 5)),
           age.lower95CI=mps$Lower95CI, age.upper95CI=mps$Upper95CI, age.average=mps$Pred_age, SG=mps$Age)

#Estimated flight distance per day (km)
#Lower bound 
w.clock.dist.df$dist.min/w.clock.dist.df$age.lower95CI

#Upper bound
w.clock.dist.df$dist.max/w.clock.dist.df$age.upper95CI

#Flight velocity - 24 hours of flying per day
#Lower bound
((w.clock.dist.df$dist.min/w.clock.dist.df$age.lower95CI)*1000)/(24*60*60)

#Upper bound
((w.clock.dist.df$dist.max/w.clock.dist.df$age.upper95CI)*1000)/(24*60*60)

#Flight velocity - 12 hours of flying per day
#Lower bound
((w.clock.dist.df$dist.min/w.clock.dist.df$age.lower95CI)*1000)/(12*60*60)

#Upper bound
((w.clock.dist.df$dist.max/w.clock.dist.df$age.upper95CI)*1000)/(12*60*60)





#Distance vs age
ggplot(w.clock.dist.df, aes(x=age.average, y=(dist.max-dist.min)/2, col=SG))+geom_point()+
  geom_pointrange(aes(ymax=dist.max, ymin=dist.min), alpha=0.8)+
  geom_errorbarh(aes(xmax = age.average+age.upper95CI, xmin = age.average-age.lower95CI, height = 0), alpha=0.8)+
  theme_bw()+
  ylab("Predicted flight distance")+
  xlab("Predicted age")+
  scale_colour_manual(values=c("#1E88E5", "#D81B60"), name="Group")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"),  axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))







#### Predictions ####
pred.df <- read.table(file=file.choose())
pred.df <- read.table(file="all.pred.CDS.dinuc.df")
pred.df <- read.table(file="all.pred.varying_seed_dinuc_bigdata.comb.10.noW.MELD_DMRs_Only")

colnames(pred.df) <- c("Age", "Pred_age", "Model", "Seed")

set.seed(1000)

lmod.df.int_and_slope <- plyr::ddply(pred.df[pred.df$Model == "Elastic",], c("Seed"), function(x) lmodel2(jitter(Age, amount=2)~Pred_age, data=x)$regression.results[1, 2:3])
#lmod.df.CIs <- plyr::ddply(pred.df[pred.df$Model == "Elastic",], c("Seed"), function(x) lmodel2(Age~Pred_age, data=x)$confidence.intervals[2, 2:5] )

#lm_eqn <- function(df){
#  m <- lm(df[,1] ~ df[,2], df);
#  eq <- substitute(italic(y) == a + b %.% italic(x),
#                        list(a = format(unname(coef(m)[1]), digits = 2),
#                        b = format(unname(coef(m)[2]), digits = 2)))
#  as.character(as.expression(eq));
#}

avg_of_slopes<-as.character(as.expression(substitute(italic(y) == a + b %.% italic(x),
           list(a = format(mean(lmod.df.int_and_slope$Intercept), digits = 2),
                b = format(mean(lmod.df.int_and_slope$Slope), digits = 2)))))




ggplot(lmod.df.int_and_slope)+
  geom_abline(linetype=2, size=1)+
  geom_abline(aes(slope = Slope, intercept = Intercept), alpha=0.05, col="blue") +
#  geom_text(x = 4, y = 25, label = lm_eqn(pred.df), parse = TRUE, col="red")+
  geom_text(x = 4, y = 28, label = avg_of_slopes, parse = TRUE, col="red", size=5)+
  ylab("Age (Days)")+
  theme_bw()+
  ylim(0,30)+
  xlim(0,30)+
  xlab("E.clock age (Days)")+
  geom_abline(data=NULL,aes(slope = mean(lmod.df.int_and_slope$Slope), intercept=mean(lmod.df.int_and_slope$Intercept)), col="red", size=1)+

  
  #geom_line(data=pred.df[pred.df$Model == "Elastic",], aes(x=Pred_age, y=Age),stat="smooth",method="lm", alpha=1, size=1, col="red")+ 
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))


###

ggplot(pred.df[pred.df$Model == "Elastic",], aes(x=Pred_age, y=Age))+
  geom_abline(linetype=1, size=1)+
  geom_smooth(method="lm", alpha=0.5, size=0,col="blue", fill="turquoise")+
  geom_line(stat="smooth",method="lm", alpha=0.5, size=1,col="blue", fill="turquoise")+
  theme_bw()+
  ylim(0,30)+
  xlim(0,30)+
  #geom_abline(linetype=1, size=2)+
  theme(aspect.ratio=1)

ggplot(pred.df[pred.df$Model == "Elastic",], aes(x=Pred_age, y=Age, group=Seed))+
  geom_abline(linetype=1, size=1)+
  geom_smooth(method="lm", alpha=0.05, size=0,col="blue", fill="turquoise")+
  geom_line(stat="smooth",method="lm", alpha=0.1, size=0.4,col="blue")+
  theme_bw()+
  #geom_abline(linetype=1, size=2)+
  theme(aspect.ratio=1)+
  theme( element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))



#RMSE and Rsquared from predictions per seed ####

library(caret) 
seeds<-unique(pred.df$Seed)

PR_func<- function(x) {postResample(pred.df[pred.df$Model == "Elastic" & pred.df$Seed == x ,]$Pred_age, pred.df[pred.df$Model == "Elastic"  & pred.df$Seed == x,]$Age)}

pred.stats<-lapply(seeds, PR_func)

RMSEs<-sapply(pred.stats, function(x) if(length(x) >= 1) x[[1]] else NA)
R2s<-sapply(pred.stats, function(x) if(length(x) >= 1) x[[2]] else NA)

mean(RMSEs)
mean(R2s)

ggplot(data=NULL, aes(x=RMSEs, y=R2s))+geom_point(alpha=0.4)+
  theme_bw()+
  ylim(0,1)+
  xlim(0,11)+
  ylab("Variance explained")+
  xlab("RMSE (days)")+
  theme(aspect.ratio=1)+
  theme( element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))



#How many non-zero coefficients?####

nzc.df<-read.table("varying_seed_dinuc_CDS_elastic_nonZero.coefs")
nzc.df<-as.vector(nzc.df$V1)

N<-1000
min(nzc.df)
max(nzc.df)
mean(nzc.df)
margins<-qt(0.975,df=N-1)*(mean(nzc.df))/sqrt(N)
mean(nzc.df)+margins
mean(nzc.df)-margins


#How many loci predict positive/negative relationships with age? ####

df.nonzero<-read.table("all.meth.df.elastic_nonZero.CDS.dinuc")

colnames(df.nonzero) <- c("Chromosome", "Position", "Methylation_level", "Sample", "Age", "Locus", "Experiment", "Seed")


mean_per_sample_nonzero<-plyr::ddply(df.nonzero, c("Age", "Sample","Experiment"), function(x) mean(x$Methylation_level, na.rm=T) )
colnames(mean_per_sample_nonzero) <- c("Age", "Sample", "Experiment", "Avg")
mean_per_sample_nonzero$Age <- ifelse(mean_per_sample_nonzero$Age == "ME", "SD", mean_per_sample_nonzero$Age)
mean_per_sample_nonzero$Age <- factor(mean_per_sample_nonzero$Age , levels=c("0", "5", "10", "15", "20", "25", "30", "SD", "LD"))


ggplot(mean_per_sample_nonzero, aes(x=Age, y=Avg*100, col=Experiment))+geom_point(size=2, position=position_dodge2(width=0.1))+
  theme_bw()+
  ylab("Methylation (%)")+
  scale_color_manual(values=c("#44AA99", "#882255"))+
  theme(aspect.ratio = 1)+
  theme( element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))




lmod<-lmodel2(Methylation_level~as.integer(Age), df.nonzero[df.nonzero$Seed == 1 & df.nonzero$Locus == "Chr_1_6118963" & df.nonzero$Age != "ME" & df.nonzero$Age != "LD",])

locus.coefs.df<-plyr::ddply(df.nonzero[df.nonzero$Age != "ME" & df.nonzero$Age != "LD",], c("Seed", "Locus"), function(x) lmodel2(Methylation_level~as.integer(Age), data=x)$regression.results[1, 2:3])


ggplot(locus.coefs.df, aes(x=Slope))+geom_histogram(colour="red", fill="orange")+geom_vline(xintercept=0)+
  theme_bw()+
  ylab("Count")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))


prop.table(table(sign(locus.coefs.df$Slope)))

locus.count.df<-as.data.frame(table(locus.coefs.df$Locus))
colnames(locus.count.df)<-c("Locus", "Count")

locus.count.df<-locus.count.df[rev(order(locus.count.df$Count)),]

ggplot(locus.count.df, aes(x=Count/10))+geom_histogram(colour="red", fill="orange")+
  theme_bw()+
  ylab("Count")+
  xlab("Included in percent of models")+
  geom_bracket(xmin = 10, xmax = 99, y.position = 300, tip.length=0.2,label.size=6, label = "160 genes")+
  geom_bracket(xmin = 50, xmax = 99, y.position = 41, tip.length=0.1, label.size=6, label = "24 genes")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))

