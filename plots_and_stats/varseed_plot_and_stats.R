library(ggplot2)
library(lmodel2)
library(plyr)
library(ggpubr)
setwd("~/Downloads")


#### Wild samples ####
varseed.df <- read.table(file="df.out.1-1000.bigdata.comb.10.noW.MELD_DMRs_Only")


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



#Plot wild age mean predictions ±95% CI
mps$Sample<-as.factor(mps$Sample)
mps$Sample<-factor(mps$Sample, levels = c("GTcoll19M359", "GTcoll19M448", "GTcoll19M450", "GTcoll19M460", "GTcoll19M471", "GTcoll18B720", "GTcoll18B721", "GTcoll18B773", "GTcoll19H116", "GTcoll19H128"))

ggplot(mps, aes(x=Sample, y=Pred_age, col=Age )) +
  geom_pointrange(mapping=aes(y=Pred_age, ymin=Lower95CI, ymax=Upper95CI))+
  ylab("W.clock age (Days)")+
  xlab("Sample")+
  theme_bw()+
  ylim(0,17)+
  scale_colour_manual(values=c("#1E88E5", "#D81B60"), name="Group")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text.y=element_text(size=16, colour="black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))





#### Migration distance versus age ####



w.clock.dist.df<-data.frame(samples=unique(mps$Sample),
           dist.min=c(2790,4250,2450,2570,3310, rep(0, 5)), dist.max=c(7320, 7160, 6520, 6330, 6530, rep(150, 5)),
           age.lower95CI=mps$Lower95CI, age.upper95CI=mps$Upper95CI, age.average=mps$Pred_age, SG=mps$Age)

#Estimated flight distance per day (km)
#Lower bound 
round(w.clock.dist.df$dist.min/w.clock.dist.df$age.lower95CI, digits=0)

#Upper bound
round(w.clock.dist.df$dist.max/w.clock.dist.df$age.upper95CI, digits=0)

#Flight velocity - 24 hours of flying per day
#Lower bound
round(((w.clock.dist.df$dist.min/w.clock.dist.df$age.lower95CI)*1000)/(24*60*60), digits=1)

#Upper bound
round(((w.clock.dist.df$dist.max/w.clock.dist.df$age.upper95CI)*1000)/(24*60*60), digits=1)

#Flight velocity - 12 hours of flying per day
#Lower bound
round(((w.clock.dist.df$dist.min/w.clock.dist.df$age.lower95CI)*1000)/(12*60*60), digits=1)

#Upper bound
round(((w.clock.dist.df$dist.max/w.clock.dist.df$age.upper95CI)*1000)/(12*60*60), digits=1)


#Sahara calculation
1800/(w.clock.dist.df$dist.min/w.clock.dist.df$age.lower95CI)
1800/(w.clock.dist.df$dist.max/w.clock.dist.df$age.upper95CI)

w.clock.dist.df$Wing_wear <-c(2.5, 2, 5, 2, 4, 1.5, 2, 1.5, 2, 2)

#Distance vs age
ggplot(w.clock.dist.df, aes(x=age.average, y=dist.min+(dist.max-dist.min)/2, col=SG))+geom_point()+
  geom_pointrange(aes(ymax=dist.max, ymin=dist.min), alpha=0.8)+
  geom_pointrange(aes(xmax = age.upper95CI, xmin = age.lower95CI), alpha=0.8)+
  theme_bw()+
  ylab("Predicted flight distance")+
  xlab("W.clock age")+
  scale_colour_manual(values=c("#1E88E5", "#D81B60"), name="Group")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"),  axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))

#Wing wear vs age
ggplot(w.clock.dist.df, aes(x=age.average, y=Wing_wear, col=SG))+geom_point()+
  geom_pointrange(aes(xmax = age.upper95CI, xmin = age.lower95CI), alpha=0.8)+
  theme_bw()+
  ylim(1,5)+
  ylab("Wing wear")+
  xlab("W.clock age")+
  scale_colour_manual(values=c("#1E88E5", "#D81B60"), name="Group")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"),  axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))

cor.test(w.clock.dist.df$age.average, w.clock.dist.df$Wing_wear)






#### Predictions ####
pred.df <- read.table(file=file.choose())
pred.df <- read.table(file="all.pred.CDS.dinuc.df")
#pred.df <- read.table(file="all.pred.varying_seed_dinuc_bigdata.comb.10.noW.MELD_DMRs_Only")

colnames(pred.df) <- c("Age", "Pred_age", "Model", "Seed")

set.seed(1000)

lmod.df.int_and_slope <- plyr::ddply(pred.df[pred.df$Model == "Elastic",], c("Seed"), function(x) lmodel2(jitter(Age, amount=2)~Pred_age, data=x)$regression.results[1, 2:3])
#lmod.df.CIs <- plyr::ddply(pred.df[pred.df$Model == "Elastic",], c("Seed"), function(x) lmodel2(Age~Pred_age, data=x)$confidence.intervals[2, 2:5] )


avg_of_slopes<-as.character(as.expression(substitute(italic(y) == a + b %.% italic(x),
           list(a = format(mean(lmod.df.int_and_slope$Intercept), digits = 2),
                b = format(mean(lmod.df.int_and_slope$Slope), digits = 2)))))

mad_of_slopes <- as.expression(
  substitute(
    atop(MAD[I] == a, MAD[S] == b),
    list(
      a = format(mad(lmod.df.int_and_slope$Intercept), digits = 2),
      b = format(mad(lmod.df.int_and_slope$Slope), digits = 2)
    )
  )
)



ggplot(lmod.df.int_and_slope)+
  geom_abline(linetype=2, size=1)+
  geom_abline(aes(slope = Slope, intercept = Intercept), alpha=0.05, col="blue") +
  #  geom_text(x = 4, y = 25, label = lm_eqn(pred.df), parse = TRUE, col="red")+
  #geom_text(x = 6, y = 25, label = mad_of_slopes, parse = TRUE, col="red", size=5)+
  ylab("Age (Days)")+
  theme_bw()+
  ylim(0,30)+
  xlim(0,30)+
  xlab("E.clock age (Days)")+
  geom_abline(data=NULL,aes(slope = mean(lmod.df.int_and_slope$Slope), intercept=mean(lmod.df.int_and_slope$Intercept)), col="red", size=1)+
  #geom_line(data=pred.df[pred.df$Model == "Elastic",], aes(x=Pred_age, y=Age),stat="smooth",method="lm", alpha=1, size=1, col="red")+ 
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))



#Full dataset fig
ggplot(lmod.df.int_and_slope)+
  geom_abline(linetype=2, size=1)+
  geom_abline(aes(slope = Slope, intercept = Intercept), alpha=0.05, col="blue") +
#  geom_text(x = 4, y = 25, label = lm_eqn(pred.df), parse = TRUE, col="red")+
  geom_text(x = 4, y = 28, label = avg_of_slopes, parse = TRUE, col="red", size=5)+
  geom_text(x = 4, y = 23, label = mad_of_slopes, parse = TRUE, col="red", size=5)+
  
  ylab("Age (Days)")+
  theme_bw()+
  ylim(0,30)+
  xlim(0,30)+
  xlab("E.clock age (Days)")+
  geom_abline(data=NULL,aes(slope = mean(lmod.df.int_and_slope$Slope), intercept=mean(lmod.df.int_and_slope$Intercept)), col="red", size=1)+


#Plot to show a specific seed
ggplot(pred.df[pred.df$Model == "Elastic" & pred.df$Seed == 1000,], aes(x=Pred_age, y=Age))+
  geom_abline(linetype=2, size=1, alpha=0.3)+
  geom_point()+
  ylab("Age (Days)")+
  theme_bw()+
  ylim(-5,30)+
  xlim(-5,30)+
  xlab("E.clock age (Days)")+
  theme( element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))


#RMSE, MAE and Rsquared from predictions per seed ####

library(caret) 
PR_func<- function(x) {postResample(pred.df[pred.df$Model == "Elastic" & pred.df$Seed == x ,]$Pred_age, pred.df[pred.df$Model == "Elastic"  & pred.df$Seed == x,]$Age)}


library(patchwork)

loopVec<-c(20,25,30,35,40,45)
loopVec<-c("all.pred.males_CDS_dinuc_fvsc5", "all.pred.females_CDS_dinuc_fvsc5", "all.pred.GT_CDS_dinuc_fvsc5", "all.pred.GT_no30_CDS_dinuc_fvsc5")

p <- list()
for(i in loopVec){
  pred.df <- read.table(paste0("subsampleN.", i, ".pred.df"))
  pred.df <- read.table(i)
  colnames(pred.df) <- c("Age", "Pred_age", "Model", "Seed")
  
  seeds<-unique(pred.df$Seed)
  
  pred.stats<-lapply(seeds, PR_func)
  
  RMSEs<-sapply(pred.stats, function(x) if(length(x) >= 1) x[[1]] else NA)
  R2s<-sapply(pred.stats, function(x) if(length(x) >= 1) x[[2]] else NA)
  MAEs<-sapply(pred.stats, function(x) if(length(x) >= 1) x[[3]] else NA)
  
  df_rmse <- data.frame(RMSE = RMSEs, R2 = R2s)
  df_mae <- data.frame(MAE = MAEs)

  print(paste("Sample size", i))
  print(mean(RMSEs))
  print(mean(R2s))
  print(mean(MAEs))
  
  p1<-ggplot(data=df_rmse, aes(x=RMSE, y=R2))+geom_point(alpha=0.4)+
    theme_bw()+
    ylim(0,1)+
    xlim(0,20)+
    ylab("Variance explained")+
    xlab("RMSE (days)")+
    theme(aspect.ratio=1)+
    theme( element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))
  
  p2<-ggplot(data=df_mae, aes(x=MAE))+geom_histogram(alpha=0.7, bins=50)+
    geom_vline(xintercept=mean(MAEs), lty=1, col="red", lwd=1.8)+
    theme_bw()+
    ylab("Counts")+
    xlab("MAE (Days)")+
    theme(aspect.ratio=1)+
    xlim(0,20)+
    theme( element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))
  
  
  if(i != 45){
    p1 <- p1 +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank())
    
    p2 <- p2 +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank())
  }
  

  p[[as.character(i)]] <- p1 + p2
  
  
  
}


wrap_plots(p, ncol = 1) &
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12)
  )

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

#Original
locus.coefs.df<-plyr::ddply(df.nonzero[df.nonzero$Age != "ME" & df.nonzero$Age != "LD",], c("Seed", "Locus"), function(x) lmodel2(Methylation_level~as.integer(Age), data=x)$regression.results[1, 2:3])

#Flatten
locus.coefs.df.flat <- locus.coefs.df[!duplicated(locus.coefs.df[c("Locus", "Intercept", "Slope")]), ]

ggplot(locus.coefs.df.flat, aes(x=Slope))+geom_histogram(colour="red", fill="orange")+geom_vline(xintercept=0)+
  theme_bw()+
  ylab("Count")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))


prop.table(table(sign(locus.coefs.df.flat$Slope)))

locus.count.df<-as.data.frame(table(locus.coefs.df$Locus))
colnames(locus.count.df)<-c("Locus", "Count")

locus.count.df<-locus.count.df[rev(order(locus.count.df$Count)),]

ggplot(locus.count.df, aes(x=Count/10))+geom_histogram(colour="red", fill="orange")+
  theme_bw()+
  ylab("Count")+
  xlab("Included in percent of models")+
  geom_bracket(xmin = 10, xmax = 99, y.position = 300, tip.length=0.2,label.size=6, label = "197 genes")+
  geom_bracket(xmin = 50, xmax = 99, y.position = 41, tip.length=0.1, label.size=6, label = "64 genes")+
  theme(aspect.ratio=1, element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text=element_text(size=16, colour="black"), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))

