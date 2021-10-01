#Lior's version of Generating PCA plot. Need to make adjustment to working
#directories. Also can remove first 3 lines if DESeq2 is installed.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(dplyr)
library(DESeq2)
setwd("Bioinformatics-Group/GSE143791_condensed")
temp = list.files(pattern="*.count")
myfiles = lapply(temp, read.csv)
matrix <- bind_cols(myfiles, .id="column_label")
expressionMatrix <- select(matrix, -c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65))
expressionMatrix <- distinct(expressionMatrix, X...1, .keep_all=TRUE)
rownames(expressionMatrix) <- expressionMatrix[,1]
expressionMatrix[,1] <- NULL
View(expressionMatrix)

####Generating PCA plot######

setwd("./..")


groupClassification <- read.csv("group_classification.csv",header = T,row.names=1,
                                stringsAsFactors=T)

dds <- DESeqDataSetFromMatrix(countData = expressionMatrix,
                              colData = groupClassification,
                              design = ~ dex)
#dds <- dds[rowSums(count(dds)) > 1,]

vst <- vst(dds, blind = FALSE)

plotPCA(vst, intgroup=c("dex"))

#create the tsne plot
library(Rtsne)
tsneResults <- Rtsne(expressionMatrix, perplexity=1, check_duplicates = FALSE)
plot(tsneResults$Y,col=groupClassification$dex)