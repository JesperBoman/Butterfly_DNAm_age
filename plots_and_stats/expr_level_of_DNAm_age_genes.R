library(data.table)
library(ggplot2)


setwd("~/Downloads")

df.tpm<-read.table(file="salmon.merged.gene_counts_length_scaled.tsv", header=T)

head(df.tpm)

GOIset<-c("Vcard_DToL01744", "Vcard_DToL01808", "Vcard_DToL01881", "Vcard_DToL01904", "Vcard_DToL01978","Vcard_DToL02146", "Vcard_DToL02147", 
          "Vcard_DToL02373", "Vcard_DToL02435", "Vcard_DToL08452", "Vcard_DToL08638","Vcard_DToL08671", "Vcard_DToL08716", "Vcard_DToL08824", 
          "Vcard_DToL08891", "Vcard_DToL09872", "Vcard_DToL09949", "Vcard_DToL10243", "Vcard_DToL11074", "Vcard_DToL17456","Vcard_DToL17620", 
          "Vcard_DToL17718", "Vcard_DToL17719", "Vcard_DToL17720")


df.tpm <- as.data.frame(melt(setDT(df.tpm[,-2]), id.vars = c("gene_id"), variable.name = "Sample"))

df.tpm[df.tpm$gene_id %in% GOIset,]$Type <- "GOI"

df.tpm$Type <- NA
for(i in 1:length(GOIset)){
  matches <-GOIset[i] == df.tpm$gene_id
  m_list <- (1:length(matches))[matches]
  df.tpm$Type[m_list] <- "500+ seeds"
}

df.tpm$Type <- ifelse(is.na(df.tpm$Type), "Background", df.tpm$Type)


ggplot(df.tpm, aes(x=Sample, y=log(value), col=Type))+geom_boxplot()


plyr::ddply(df.tpm, c("Type", "Sample"), function(x) median(x$value, na.rm=T) )


vs_geneID_and_loci_counts<-read.table("vs_geneID_and_loci_counts.list")
colnames(vs_geneID_and_loci_counts) <- c("Seeds", "Gene", "Chromosome", "Position")

vs_geneID_and_loci_counts_not500 <- vs_geneID_and_loci_counts[ ! (vs_geneID_and_loci_counts$Gene %in% GOIset),]

vs_geneID_and_loci_counts_lessthan100 <- vs_geneID_and_loci_counts_not500[vs_geneID_and_loci_counts_not500$Seeds < 100,]

for(i in 1:length(vs_geneID_and_loci_counts_lessthan100$Gene)){
  matches <-vs_geneID_and_loci_counts_lessthan100$Gene[i] == df.tpm$gene_id
  m_list <- (1:length(matches))[matches]
  df.tpm$Type[m_list] <- "1-99 seeds"
}

fit<-lm(log(value+1)~Type+Sample, df.tpm[df.tpm$value<3e4,])
hist(residuals(fit))
plot(fit)

anova(fit)
TukeyHSD(aov(fit))$Type

ggplot(df.tpm[,], aes(x=Sample, y=log(value), col=Type))+geom_boxplot()+
  theme_bw()+
  ylim(-1,15)+
  theme(aspect.ratio=1)+
  ylab("Log TPM")+
  theme( element_text(face = "bold", hjust = 0.5),  panel.border = element_rect(colour = "black", fill=NA, linewidth=1), axis.text.y=element_text(size=16, colour="black"), axis.text.x=element_text(size=10, colour="black", angle = 90), axis.title=element_text(size=18), legend.text=element_text(size=14),  legend.title=element_text(size=16))


length(unique(vs_geneID_and_loci_counts_lessthan100$Gene))
