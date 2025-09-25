
BiocManager::install("topGO")
library(topGO)
library(ggplot2)

setwd("~/Downloads")
geneID2GO <- readMappings(file="vanCar_GO_genes.list", sep = "\t", IDsep = ",")

str(geneID2GO)
geneNames <- names(geneID2GO)



#read the list of genes of interest
candidate_genes <- read.table("vs_geneID_atLeast500seeds.list", header = F)
colnames(candidate_genes) <- c("geneid")


head(candidate_genes)
candidate_genes <- as.character(candidate_genes$geneid)
geneList <- factor(as.integer(geneNames %in% candidate_genes))
names(geneList) <- geneNames
head(geneList)
str(geneList)
#create the topGOdata object for each ontology
GO_data_BP <- new("topGOdata", 
                  ontology="BP", 
                  allGenes=geneList, 
                  annot=annFUN.gene2GO, 
                  gene2GO=geneID2GO, 
                  nodeSize=1)
GO_data_BP
GO_data_CC <- new("topGOdata", 
                  ontology="CC", 
                  allGenes=geneList, 
                  annot=annFUN.gene2GO, 
                  gene2GO=geneID2GO, 
                  nodeSize=1)
GO_data_CC
GO_data_MF <- new("topGOdata", 
                  ontology="MF", 
                  allGenes=geneList, 
                  annot=annFUN.gene2GO, 
                  gene2GO=geneID2GO, 
                  nodeSize=1)
GO_data_MF


#"when only a list of interesting genes is provided, 
#the user can use only tests statistics that are based on gene counts, 
#like Fisher’s exact test, Z score and alike."

##Algorithms:
#classic uses all GO terms separately
#elim eliminates genes from parentGOterms if child is more specific (bottom up analysis)
#weight trying to determin if GOterm is better representing the list of interesting genes (is more enriched) than any other node from its neighborhood
#parentChild 

#biological process
BP_resultFisher_weight01 <- runTest(GO_data_BP, algorithm = "weight01", statistic = "fisher")
BP_resultFisher_weight01

BP_resultFisher_classic <- runTest(GO_data_BP, algorithm = "classic", statistic = "fisher")
BP_resultFisher_classic

BP_resultFisher_parentChild <- runTest(GO_data_BP, algorithm = "parentChild", statistic = "fisher")
BP_resultFisher_parentChild

BP_resultFisher_elim <- runTest(GO_data_BP, algorithm = "elim", statistic = "fisher")
BP_resultFisher_elim



#allGO = usedGO(object = GOdata)
#topNodes = length(allGO) in GenTable() to get a table with all GO:s to do fdr
allGO = usedGO(object = GO_data_BP)
allRes_BP <- GenTable(GO_data_BP, 
                      weight01Fisher = BP_resultFisher_weight01, 
                      classicFisher = BP_resultFisher_classic,
                      parentChFisher=BP_resultFisher_parentChild,
                      elimFisher=BP_resultFisher_elim,
                      orderBy = "weight01Fisher", 
                      ranksOf = "weight01Fisher", 
                      topNodes = length(allGO),
                      numChar=1000)
allRes_BP$weight01Fisher <- as.numeric(allRes_BP$weight01Fisher)
allRes_BP[allRes_BP$weight01Fisher<=0.05, ]
allRes_BP$parentChFisher <- as.numeric(allRes_BP$parentChFisher)
allRes_BP[allRes_BP$parentChFisher<=0.05, ]
allRes_BP$elimFisher <- as.numeric(allRes_BP$elimFisher)
allRes_BP[allRes_BP$elimFisher<=0.05, ]

write

#multiple test correction, method Benjamini-Hochberg
#TopGO manual: https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf
#For the methods that account for the GO topology like elim and weight, the problem of multiple
#testing is even more complicated. Here one computes the p-value of a GO term conditioned on the
#neighbouring terms. The tests are therefore not independent and the multiple testing theory does not
#directly apply. We like to interpret the p-values returned by these methods as corrected or not affected
#by multiple testing.
allRes_BP$p_adj <- p.adjust(allRes_BP$weight01Fisher, method = "BH")
allRes_BP$p_adj_pc <- p.adjust(allRes_BP$parentChFisher, method = "BH")
allRes_BP$p_adj_elim <- p.adjust(allRes_BP$elimFisher, method = "BH")

#adj p-value below 0.1, NOTE IN THE TOPGO MANUAL THEY RECOMMEND TO NOT ADJUST P-VALUES
allRes_BP[allRes_BP$p_adj<0.1, ]
allRes_BP[allRes_BP$p_adj_pc<0.1, ]
allRes_BP[allRes_BP$p_adj_elim<0.1, ]
allRes_BP$GO_class <- "BP"

allRes_BP$Sig_over_exp<-allRes_BP$Significant/allRes_BP$Expected
write.table(allRes_BP[allRes_BP$weight01Fisher<=0.05, c(1:6, dim(allRes_BP)[2])], file="EpiClock_GO_BP_results_vs_geneID_atLeast500seeds.txt", sep="\t", row.names=F, quote=F)

#cellular compartment
CC_resultFisher_weight01 <- runTest(GO_data_CC, algorithm = "weight01", statistic = "fisher")
CC_resultFisher_weight01
CC_resultFisher_classic <- runTest(GO_data_CC, algorithm = "classic", statistic = "fisher")
CC_resultFisher_classic
CC_resultFisher_parentChild <- runTest(GO_data_CC, algorithm = "parentChild", statistic = "fisher")
CC_resultFisher_parentChild
CC_resultFisher_elim <- runTest(GO_data_CC, algorithm = "elim", statistic = "fisher")
CC_resultFisher_elim
#allGO = usedGO(object = GOdata)
#topNodes = length(allGO) in GenTable() to get a table with all GO:s to do fdr
allGO = usedGO(object = GO_data_CC)
allRes_CC <- GenTable(GO_data_CC, 
                      weight01Fisher = CC_resultFisher_weight01, 
                      classicFisher = CC_resultFisher_classic,
                      parentChFisher=CC_resultFisher_parentChild,
                      elimFisher=CC_resultFisher_elim,
                      orderBy = "weight01Fisher", 
                      ranksOf = "weight01Fisher", 
                      topNodes = length(allGO),
                      numChar=1000)
allRes_CC$weight01Fisher <- as.numeric(allRes_CC$weight01Fisher)
allRes_CC[allRes_CC$weight01Fisher<=0.05, ]
allRes_CC$parentChFisher <- as.numeric(allRes_CC$parentChFisher)
allRes_CC[allRes_CC$parentChFisher<=0.05, ]
allRes_CC$elimFisher <- as.numeric(allRes_CC$elimFisher)
allRes_CC[allRes_CC$elimFisher<=0.05, ]
#multiple test correction, method Benjamini-Hochberg
allRes_CC$p_adj <- p.adjust(allRes_CC$weight01Fisher, method = "BH")
allRes_CC$p_adj_pc <- p.adjust(allRes_CC$parentChFisher, method = "BH")
allRes_CC$p_adj_elim <- p.adjust(allRes_CC$elimFisher, method = "BH")
#adj p-value below 0.1, NOTE IN THE TOPGO MANUAL THEY RECOMMEND TO NOT ADJUST P-VALUES
allRes_CC[allRes_CC$p_adj<0.1, ]
allRes_CC[allRes_CC$p_adj_pc<0.1, ]
allRes_CC[allRes_CC$p_adj_elim<0.1, ]
allRes_CC$GO_class <- "CC"

#molecular function
MF_resultFisher_weight01 <- runTest(GO_data_MF, algorithm = "weight01", statistic = "fisher")
MF_resultFisher_weight01
MF_resultFisher_classic <- runTest(GO_data_MF, algorithm = "classic", statistic = "fisher")
MF_resultFisher_classic
MF_resultFisher_parentChild <- runTest(GO_data_MF, algorithm = "parentChild", statistic = "fisher")
MF_resultFisher_parentChild
MF_resultFisher_elim <- runTest(GO_data_MF, algorithm = "elim", statistic = "fisher")
MF_resultFisher_elim
#allGO = usedGO(object = GOdata)
#topNodes = length(allGO) in GenTable() to get a table with all GO:s to do fdr
allGO = usedGO(object = GO_data_MF)
allRes_MF <- GenTable(GO_data_MF, 
                      weight01Fisher = MF_resultFisher_weight01, 
                      classicFisher = MF_resultFisher_classic,
                      parentChFisher=MF_resultFisher_parentChild,
                      elimFisher=MF_resultFisher_elim,
                      orderBy = "weight01Fisher", 
                      ranksOf = "weight01Fisher", 
                      topNodes = length(allGO),
                      numChar=1000)
allRes_MF$weight01Fisher <- as.numeric(allRes_MF$weight01Fisher)
allRes_MF[allRes_MF$weight01Fisher<=0.05, ]
allRes_MF$parentChFisher <- as.numeric(allRes_MF$parentChFisher)
allRes_MF[allRes_MF$parentChFisher<=0.05, ]
allRes_MF$elimFisher <- as.numeric(allRes_MF$elimFisher)
allRes_MF[allRes_MF$elimFisher<=0.05, ]
#multiple test correction, method Benjamini-Hochberg
allRes_MF$p_adj <- p.adjust(allRes_MF$weight01Fisher, method = "BH")
allRes_MF$p_adj_pc <- p.adjust(allRes_MF$parentChFisher, method = "BH")
allRes_MF$p_adj_elim <- p.adjust(allRes_MF$elimFisher, method = "BH")

#adj p-value below 0.1, NOTE IN THE TOPGO MANUAL THEY RECOMMEND TO NOT ADJUST P-VALUES
allRes_MF[allRes_MF$p_adj<0.1, ]
allRes_MF[allRes_MF$p_adj_pc<0.1, ]
allRes_MF[allRes_MF$p_adj_elim<0.1, ]
allRes_MF$GO_class <- "MF"
#add gene names per GO (candidate genes)
allRes_BP$genes <- sapply(allRes_BP$GO.ID,
                          function(x)
                          {
                            genes_sel <- genesInTerm(GO_data_BP, x)
                            genes_sel[[1]][genes_sel[[1]] %in% candidate_genes]
                          })
allRes_BP$genes <- vapply(allRes_BP$genes, paste, collapse =",", character(1L))
allRes_CC$genes <- sapply(allRes_CC$GO.ID,
                          function(x)
                          {
                            genes_sel <- genesInTerm(GO_data_CC, x)
                            genes_sel[[1]][genes_sel[[1]] %in% candidate_genes]
                          })
allRes_CC$genes <- vapply(allRes_CC$genes, paste, collapse =",", character(1L))
allRes_MF$genes <- sapply(allRes_MF$GO.ID,
                          function(x)
                          {
                            genes_sel <- genesInTerm(GO_data_MF, x)
                            genes_sel[[1]][genes_sel[[1]] %in% candidate_genes]
                          })
allRes_MF$genes <- vapply(allRes_MF$genes, paste, collapse =",", character(1L))



allRes_all_comb <- rbind(allRes_BP[allRes_BP$weight01Fisher<=0.01, ], allRes_CC[allRes_CC$weight01Fisher<=0.01, ], allRes_MF[allRes_MF$weight01Fisher<=0.01, ])
allRes_all_comb$GO_descr <- paste(allRes_all_comb$Term," (",allRes_all_comb$GO.ID,")", sep = "")

#plot go-terms and number of significant genes
ggplot(allRes_all_comb, aes(reorder(GO_descr, -Significant), Significant)) +
  geom_bar(stat = "identity", aes(fill=GO_class)) +
  facet_grid(~GO_class, scales = "free_x", space = "free_x") +
  ylab("Number of genes (p-value < 0.01)") +
  scale_fill_discrete(guide = guide_legend(label.position = "top")) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.line = element_line(size = 0.2, colour = "grey"),
        axis.text.x = element_text(size = 2, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12, angle = 90, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, angle = 90, vjust = 0.5),
        legend.position = "top"
  )



