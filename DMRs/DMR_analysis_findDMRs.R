#!/usr/bin/env Rscript
#### A script to perform an analysis of differentially methylated regions  (DMRs) using BSmooth
#### Usage: Rscript DMR_analysis_findDMRs.R 
#### Jesper Boman - 2024-10-13



library("bsseq")
library("BiocParallel")

args = commandArgs(trailingOnly = TRUE)

load(args[1])
load(args[2])


BS.cov <- getCoverage(BS_data.fit)

#Assuming sample sizes of 5 and 5
keepLoci<- which(rowSums(BS.cov[, 1:5] >= 10) >= 2 &
                            rowSums(BS.cov[, 6:10] >= 10) >= 2)
length(keepLoci)

comp.fit.filt <- BS_data.fit[keepLoci,]


comp.fit.filt.tstat <- BSmooth.tstat(comp.fit.filt, 
                                          group1 = colnames(BS_data)[1:5],
                                          group2 = colnames(BS_data)[6:10], 
                                          estimate.var = "same",
                                          local.correct = TRUE,
                                          verbose = TRUE)

png(filename = "t_plot.png")
plot(comp.fit.filt.tstat )
dev.off()

#Quantile cutoff
dmrs0 <- dmrFinder(comp.fit.filt.tstat , qcutoff = c(0.01, 0.99))

#Number of CpGs and abs mean diff cutoff
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)



save(comp.fit.filt, file="comp.fit.cov10.rda")
save(dmrs, file="dmrs.cov10.rda")
write.table(dmrs, file="dmrs.cov10.txt", sep="\t", quote =F, row.names=F)
dmrs$start <- dmrs$start - 1 
write.table(dmrs, file="dmrs.cov10.bed", sep="\t", quote =F, row.names=F)
