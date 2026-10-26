#!/usr/bin/env Rscript
#### A script to train an epigenetic aging clock
#### Usage: Rscript fancier_varying_seed_clock5_F.R $i $outdir $inputdata $threshold_age
#### I run the script in the "varying seed clock approach" using the bash script master.clocksmith.sh, which in turn runs smither.sh
#### Jesper Boman - 2026-09-18

library(reshape2)
library(tidymodels) 
library(caret)
library(dplyr)
library(doParallel)
library(zoo)

#### CHANGELOG ####
#fancier_varying_seed_clock3
#Builds on fancier_varying_seed_clock

#Changed this below:
#preProcValues <- preProcess(age_test, method=c("center","scale"))
#testTransformed <- predict(preProcValues, age_test)

#To
#testTransformed <- predict(preProcValues, age_test)

#To avoid "data leakage" issue

#fancier_varying_seed_clock4_F
#Builds on fancier_varying_seed_clock3, adds a piecewise age transformation to account for faster aging at young ages 


#fancier_varying_seed_clock5_F
#Second data leakage issue fix
#preProcValues <- preProcess(meth.final.impute, method = c("center", "scale")) 

#And adding adult age as an input parameter

#Fixed code for checking for highly correlated variables

###################


#This script is inspired Anastasiadi and Piferrer (2023): https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2023.1096909/full


args = commandArgs(trailingOnly = TRUE)

#Arguments are
#1) seed (a number to be used for setting the seed)
#2) output directory
#3) input data (which is the output of i) comb_data.sh only or after also running ii) add_E2_data.sh


#Read input data
meth.df <- read.table(args[3])

colnames(meth.df) <- c("Chromosome", "Position", "Methylation_level", "Sample", "Age")

#First we optimize a clock based on training and test samples

#We exclude wild samples
meth.df.sub <- meth.df[meth.df$Age != "meconium",]
meth.df.sub <- meth.df.sub[meth.df.sub$Age != "long-distance",]

meth.df.sub$Age <- as.numeric(meth.df.sub$Age)
meth.df.sub$Locus <- paste(meth.df.sub$Chromosome, meth.df.sub$Position, sep="_")

#For CDS only, this is to remove duplicated positions with overlapping CDS annotations
meth.df.sub <- meth.df.sub[!duplicated(meth.df.sub),]
#

#We make the dataset wide
meth.df.wide <- dcast(meth.df.sub, Age+Sample~Locus, value.var = "Methylation_level")

#Imputation using the mean across all samples for sites with missing data
meth.df.impute <- meth.df.wide[ , colSums(is.na(meth.df.wide))<3,]
meth.df.impute<- na.aggregate(meth.df.impute[,3:ncol(meth.df.impute)])
meth.df.impute<-meth.df.impute[,-1]
meth.df.impute$Age<-meth.df.wide$Age



set.seed(args[1])

print(paste("Seed", args[1]))

#Here we split the data into training+test splits. This is a random process, with the randomness controlled by the seed
splits <- initial_split(meth.df.impute, prop=0.75, strata = Age) 
age_other <- training(splits) 
age_test <- testing(splits) 


print("Training set proportions by age class")
age_other %>% count(Age) %>% mutate(prop = n/sum(n)) 

print("Test set proportions by age class")
age_test %>% count(Age) %>% mutate(prop = n/sum(n))


#Excluding features with zero or near-zero variance among groups
nzv.cpg <- nearZeroVar(age_other, saveMetrics= TRUE, names=TRUE, freqCut = 85/15, uniqueCut = 50)

## Detect features, exclude them and save the object 
nzv.cpg.list <- nearZeroVar(age_other, freqCut = 85/15, uniqueCut = 50) 

filteredDescr <- age_other[, -nzv.cpg.list] 

print("Dimensions filteredDescr")
dim(filteredDescr)


#Exclude highly correlated variables
filteredDescr$Age <- as.numeric(filteredDescr$Age)

cor.mat <- cor(
  filteredDescr[, !(names(filteredDescr) %in% "Age")],
  use = "pairwise.complete.obs"
)

highlyCorDescr <- findCorrelation(
  cor.mat,
  cutoff = 0.8
)
print("Number of highly correlated CpGs:")
length(highlyCorDescr)

filteredDescr.cor <- filteredDescr[, -highlyCorDescr]

#For some reason have to skip previous step for top50 DMRs
#filteredDescr.cor <- filteredDescr

print("Number of CpGs in final input set:")
dim(filteredDescr.cor)

print("Transformation via preProcess data")
preProcValues <- preProcess(
  filteredDescr.cor[, !(names(filteredDescr.cor) %in% "Age")],
  method = c("center","scale")
)
trainTransformed <- predict(preProcValues, filteredDescr.cor)



print("Training")
fitControl <- trainControl(method = "repeatedcv" ,number=10, repeats=10)
lambda_vals <- 10^seq(-3, 3, length = 10)

adult_age<-as.numeric(args[4])

F_age <- function(age, adult.age = adult_age) {
  ifelse(
    age <= adult.age,
    log(age + 1) - log(adult.age + 1),
    (age - adult.age) / (adult.age + 1)
  )
}

Finv <- function(y, adult.age = adult_age) {
  ifelse(
    y <= 0,
    (adult.age + 1) * exp(y) - 1,
    adult.age + y * (adult.age + 1)
  )
}


trainTransformed$Age <- F_age(age_other$Age)


cl <- makePSOCKcluster(5)
registerDoParallel(cl)

options(expressions = 5e5)
#https://stackoverflow.com/questions/32826906/how-to-solve-protection-stack-overflow-issue-in-r-studio

ridge_model <- train(Age ~., data = trainTransformed, method ="glmnet", trControl = fitControl, tuneGrid = data.frame(alpha = 0, lambda = lambda_vals), tuneLength = 10)

lasso_model <- train(Age ~., data = trainTransformed, method ="glmnet", trControl = fitControl, tuneGrid = data.frame(alpha = 1, lambda = lambda_vals), tuneLength = 10)

elastic_model <- train(Age ~ ., data = trainTransformed, method ="glmnet", trControl = fitControl, tuneLength = 10)

elastic_model.05 <- train(Age ~., data = trainTransformed, method = "glmnet", trControl = fitControl, tuneGrid = data.frame(alpha = 0.5, lambda = lambda_vals), tuneLength = 10)


print("Compare metrics of the models")
models_compare <- resamples(list(R=ridge_model, LM=lasso_model, EM=elastic_model, EM05=elastic_model.05)) 
summary(models_compare) 


print("Count features (CpGs) kept by each model.")
sum(coef(ridge_model$finalModel, ridge_model$bestTune $lambda)!=0) 
sum(coef(lasso_model$finalModel, lasso_model$bestTune $lambda)!=0) 
sum(coef(elastic_model$finalModel, elastic_model$bestTune $lambda)!=0) 


print("Compare metrics in the training datasets")
#This is of lesser interest since you predict age based on the samples you've built the clock on. Predictability should be very good.

predicted.age <-  Finv(predict.train(ridge_model, trainTransformed))
postResample(pred = predicted.age, age_other$Age) 
try(cor.test(predicted.age, age_other$Age))

predicted.age <-  Finv(predict.train(lasso_model, trainTransformed))
postResample(pred = predicted.age, age_other$Age) 
try(cor.test(predicted.age, age_other$Age))

predicted.age <-  Finv(predict.train(elastic_model, trainTransformed))
postResample(pred = predicted.age, age_other$Age) 
try(cor.test(predicted.age, age_other$Age))







print("Test using test dataset")
#This is the real test of the clock where you test using your test dataset

testTransformed <- predict(preProcValues, age_test)
testTransformed$Age <- as.numeric(age_test$Age)


print("Ridge")
predict.ridge.test <- Finv(predict(ridge_model, testTransformed))
postResample(pred = predict.ridge.test, testTransformed$Age) 
try(cor.test(predict.ridge.test, testTransformed$Age))

print("Lasso")
predict.lasso.test <- Finv(predict(lasso_model, testTransformed))
postResample(pred = predict.lasso.test, testTransformed$Age) 
try(cor.test(predict.lasso.test, testTransformed$Age))

print("Elastic")
predict.elastic.test <- Finv(predict(elastic_model, testTransformed))
postResample(pred = predict.elastic.test, testTransformed$Age) 
try(cor.test(predict.elastic.test, testTransformed$Age))

print("Elastic.05")
predict.elastic.05.test <- Finv(predict(elastic_model.05, testTransformed))
postResample(pred = predict.elastic.05.test, testTransformed$Age) 
try(cor.test(predict.elastic.05.test, testTransformed$Age))






#Save predictions into a dataframe
test_predictions_df<-as.data.frame(cbind(Age=age_test$Age, Pred_age=predict.ridge.test, Model=rep("Ridge", length(predict.ridge.test))))
test_predictions_df<-rbind(test_predictions_df, cbind(Age=age_test$Age, Pred_age=predict.lasso.test, Model=rep("Lasso", length(predict.lasso.test))))
test_predictions_df<-rbind(test_predictions_df, cbind(Age=age_test$Age, Pred_age=predict.elastic.test, Model=rep("Elastic", length(predict.elastic.test))))
test_predictions_df<-rbind(test_predictions_df, cbind(Age=age_test$Age, Pred_age=predict.elastic.05.test, Model=rep("Elastic0.5", length(predict.elastic.05.test))))

test_predictions_df$Age <- as.numeric(test_predictions_df$Age)
test_predictions_df$Pred_age <- as.numeric(test_predictions_df$Pred_age)

#Here we just create a dataset of mean methylation levels per sample, useful for plotting.
#This can be skipped but good for exploratory analysis
mean_per_sample<-plyr::ddply(meth.df, c("Age", "Sample"), function(x) mean(x$Methylation_level, na.rm=T) )


colnames(mean_per_sample) <- c("Age", "Sample", "Avg")
mean_per_sample$Age <- ifelse(mean_per_sample$Age == "meconium", "SD", mean_per_sample$Age)
mean_per_sample$Age <- ifelse(mean_per_sample$Age == "long-distance", "LD", mean_per_sample$Age)
mean_per_sample$Age <- factor(mean_per_sample$Age , levels=c("0", "5", "10", "15", "20", "25", "30", "SD", "LD"))

mean_per_sample$Experiment <-ifelse(grepl("AL", mean_per_sample$Sample) | grepl("LI", mean_per_sample$Sample), "E2", "E1")



#Predict age of wild samples (i.e. in our case neither training nor test samples)
meth.df$Locus <- paste(meth.df$Chromosome, meth.df$Position, sep="_")

meth.final<-meth.df[ meth.df$Locus %in% colnames(trainTransformed), ]

meth.final.wide<-dcast(meth.final, Age+Sample~Locus, value.var = "Methylation_level")
meth.final.impute<- na.aggregate(meth.final.wide[,3:ncol(meth.final.wide)])
meth.final.impute$Age <- meth.final.wide$Age



meth.final.predictors <- meth.final.impute[, setdiff(colnames(filteredDescr.cor), "Age")]


lwTransformed <- predict(preProcValues, meth.final.predictors)

lwTransformed$Age <- meth.final.impute$Age

lwTransformed.LD<-lwTransformed[lwTransformed$Age == "long-distance",]
lwTransformed.ME<-lwTransformed[lwTransformed$Age == "meconium",]



LD.df <- data.frame(
  Prediction = Finv(
    predict(elastic_model, lwTransformed.LD)
  )
)

LD.df$Age <- "LD"


ME.df <- data.frame(
  Prediction = Finv(
    predict(elastic_model, lwTransformed.ME)
  )
)

ME.df$Age <- "ME"

#Here we create a data frame to output with relevant information
out.df<-rbind(LD.df, ME.df)

out.df$Seed<-args[1]

eltest<-cor.test(predict.elastic.test, testTransformed$Age)
out.df$eltest.p <- eltest$p.value

lmod<-lm(Age~Pred_age, test_predictions_df[test_predictions_df$Model == "Elastic",])

out.df$lm.intercept<-coef(lmod)[1]
out.df$lm.slope<-coef(lmod)[2]

sumlmod<-summary(lmod)

out.df$lm.adj.r.squared <- sumlmod$adj.r.squared

#http://gettinggeneticsdone.blogspot.com/2011/01/rstats-function-for-extracting-f-test-p.html
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

out.df$lm.p <- lmp(lmod)

test_predictions_df$Seed<-args[1]




elastic_cl<-coef(elastic_model$finalModel, elastic_model$bestTune $lambda)

elastic_nonZero_coefs<-elastic_cl@Dimnames [[1]] [elastic_cl@i+1]

#Here we create a data frame of sites included in the final model of the elastic net regression
meth.df.elastic_nonZero <- meth.df[meth.df$Locus %in% elastic_nonZero_coefs , ]


meth.df.elastic_nonZero$Experiment <-ifelse(grepl("AL", meth.df.elastic_nonZero$Sample) | grepl("LI", meth.df.elastic_nonZero$Sample), "E2", "E1")
meth.df.elastic_nonZero$Age <- ifelse(meth.df.elastic_nonZero$Age == "meconium", "ME", meth.df.elastic_nonZero$Age)
meth.df.elastic_nonZero$Age <- ifelse(meth.df.elastic_nonZero$Age == "long-distance", "LD", meth.df.elastic_nonZero$Age)
meth.df.elastic_nonZero$Age <- factor(meth.df.elastic_nonZero$Age , levels=c("0", "5", "10", "15", "20", "25", "30", "ME", "LD"))
meth.df.elastic_nonZero$Seed <- args[1]




#save(elastic_model, file=paste(args[2], "/elastic.model.", args[1], ".rda", sep=""))
write.table(out.df, file = paste(args[2], "/out.df.", args[1], sep=""), quote = F, sep = "\t", row.names = T, col.names = F)
write.table(test_predictions_df, file = paste(args[2], "/pred.df.", args[1], sep=""), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(meth.df.elastic_nonZero, file = paste(args[2], "/meth.df.elastic_nonZero.df.", args[1], sep=""), quote = F, sep = "\t", row.names = F, col.names = F)
