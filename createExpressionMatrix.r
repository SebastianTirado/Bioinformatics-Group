#This script creates a matrix of the counts of genes for each sample
library(dplyr)
setwd("GSE143791_condensed")
temp = list.files(pattern="*.count")
myfiles = lapply(temp, read.csv)
matrix <- bind_cols(myfiles, id="column_label")

# remove duplicate label columns
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

### Generating DESeq2 object ###
setwd("./..")


group_classification <- read.csv("group_classification.csv",header = T,row.names=1,
                                stringsAsFactors=T)
filt_expression_matrix <- expressionMatrix %>% dplyr::filter(rowSums(.) >= 10)

dds <- DESeqDataSetFromMatrix(countData = filt_expression_matrix,
                              colData = group_classification,
                              design = ~ dex)

deseq_object <- DESeq(dds)
deseq_results <- results(deseq_object)

deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)

deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

# example showing difference in counts for tumor vs non-tumor
# plotCounts(dds, gene = "AP006222.2", intgroup = "dex")

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
