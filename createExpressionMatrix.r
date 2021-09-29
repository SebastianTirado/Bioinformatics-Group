#This script creates a matrix of the counts of genes for each sample
library(dplyr)
setwd("GSE143791_condensed")
temp = list.files(pattern="*.count")
myfiles = lapply(temp, read.csv)
matrix <- bind_cols(myfiles, .id="column_label")
expressionMatrix <- select(matrix, -c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65))
View(expressionMatrix)
