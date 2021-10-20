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


####clustProfiler#####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(clusterProfiler)
library(DOSE)
##BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

getwd()
d = read.csv("diff_expression_results.txt")
View(d)

d2 <- d[, c(1, 3, 2, 4, 5, 6, 7)]
View(d2)

## assume 1st column is ID
## 2nd column is FC

## feature 1: numeric vector
myGeneList = d2[,2]
## feature 2: named vector
names(myGeneList) = as.character(d2[,1])
## feature 3: decreasing orde
myGeneList = sort(myGeneList, decreasing = TRUE)

Name <- as.character(d2[,1])
FoldChange <- d2[,2]
myGeneListDf <- data.frame(Name, FoldChange)
#myGeneListDf <- myGeneListDf[order(myGeneListDf$Name),]
View(myGeneListDf)

myGenes2 <- data.frame(matrix(mapIds(as.character(unlist(myGeneListDf[,1]), org.Hs.eg.db), "ENTREZID","SYMBOL"), 2))
View(myGenes2)
geneNames = as.character(myGeneList[,1])

##convert to entrez id

#eg = bitr(geneNames, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#myGeneList[,1] <- eg[,2]
#head(eg)

View(myGeneList)
View(geneList)

eg = bitr(geneNames, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
View(eg)

#data(myGeneList)
gene <- names(myGeneList)[abs(myGeneList) > 1.5]
head(eg)
geneIds <- eg[,2]
head(geneIds)
head(gene)

data <- myGeneList

head(x)

dgn <- enrichDGN(geneIds)
head(dgn)
