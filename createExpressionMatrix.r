#This script creates a matrix of the counts of genes for each sample
library(dplyr)
setwd("GSE143791_condensed")
temp = list.files(pattern="*.count")
myfiles = lapply(temp, read.csv)
matrix <- bind_cols(myfiles, .id="column_label")
expressionMatrix <- select(matrix, -c(3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65))
expressionMatrix <- distinct(expressionMatrix, X...1, .keep_all=TRUE)
rownames(expressionMatrix) <- expressionMatrix[,1]
expressionMatrix[,1] <- NULL
View(expressionMatrix)

# find range for each gene expression
range_df <- data.frame("gene", "range")
for(i in 1:nrow(expressionMatrix)){     # for each gene in the matrix
  gene_name <- expressionMatrix[i,1]
  num_col <- ncol(expressionMatrix)
  numeric_row <- as.numeric(expressionMatrix[i, 2:num_col])
  
  # finds the difference between the largest and smallest value
  range_val <- max(numeric_row)- min(numeric_row)
  range_df[nrow(range_df) + 1,] = list(gene_name, as.numeric(range_val))
}
view(range_df)
