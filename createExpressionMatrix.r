#This script creates a matrix of the counts of genes for each sample
install.packages("tidyverse")
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
<<<<<<< HEAD
view(range_df)
=======
View(range_df)
>>>>>>> 2b5c11aab3e5b3c88b906219eab3b98771d1ea6a

### Generating DESeq2 object ###
setwd("./..")

# holds list of sample names a factor which puts them into tumor or non-tumor group
group_classification <- read.csv("group_classification.csv",header = T,row.names=1,
                                stringsAsFactors=T)

# filter out genes that are not highly expressed
filt_expression_matrix <- expressionMatrix %>% dplyr::filter(rowSums(.) >= 10)

<<<<<<< HEAD
=======
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")  # in case package is not installed
BiocManager::install("apeglm")
library("DESeq2")
library("apeglm")

>>>>>>> 2b5c11aab3e5b3c88b906219eab3b98771d1ea6a
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

# normalize values
dds_norm <- rlog(dds)

# get upper 25% of high-variance genes
variances <- apply(assay(dds_norm), 1, var)
upper_var <- quantile(variances, 0.75)

# put these genes into dataframe
df_by_var <- data.frame(assay(dds_norm)) %>%
  dplyr::filter(variances > upper_var)

heatmap <- pheatmap(
  df_by_var,
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Non-Annotated Heatmap",
  colorRampPalette(c(
    "blue",
    "white",
    "red"
  ))(25
  ),
  scale = "row" # Scale values in the direction of genes (rows)
)

# create list of gene names for enrichment analysis
gene_names <- range_df$X.gene.[-1]

# use gprofiler3 to find relationships
gostres <- gost(query = gene_names, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
<<<<<<< HEAD
=======

install.packages("matrixStats") # only need to install once and then ignore these lines
library(matrixStats)

# create new column with variance values
expressionMatrix$variance = rowVars(as.matrix(expressionMatrix[,c(-1)]))
# order data by variance decreasing
expressionMatrix = expressionMatrix[order(expressionMatrix[,33],decreasing=TRUE),]

expression_var <- expressionMatrix
expressionMatrix <- expressionMatrix[-33]

# get subset of X amount of top values based on variance - convert to df
greatest5000VarGenes.df = head(as.data.frame(expressionMatrix),5000)
greatest10VarGenes.df = head(as.data.frame(expressionMatrix), 10)
greatest100VarGenes.df = head(as.data.frame(expressionMatrix), 100)
greatest1000VarGenes.df = head(as.data.frame(expressionMatrix), 1000)
greatest10000VarGenes.df = head(as.data.frame(expressionMatrix), 10000)
# RUN INDIVIDUAL CLUSTERING ALGORITHM - PUSH TO GITHUB
test_df <- greatest5000VarGenes.df
test_df <- t(test_df)
# necessary for hclust
dist_matrix = dist(test_df)

hc <- hclust(dist_matrix, method="ward.D2")
plot(hc, labels = NULL)

# alluvial diagram for hclust
cluster_cut <- cutree(hc, 4)

# alluvial_df <- data.frame(table(cluster_cut))
alluvial_df <- rbind(alluvial_df, data.frame(table(cluster_cut)))
View(alluvial_df)
group <- c('10', '10', '10', '10', '10', '100', '100', '100', '100', '1000',
           '1000','1000', '1000','5000', '5000', '5000', '5000', '10000', 
           '10000', '10000', '10000')
alluvial_df$num_genes <- group
mod_alluvial_df <- alluvial_df[-c(1,2,3,4,5,6), ]  # get a subset for better scales

ggplot(mod_alluvial_df,
               aes(y = Freq, axis1 = num_genes, axis2 = cluster_cut)) +
       geom_alluvium(aes(fill = cluster_cut), width = 1/12) +
       geom_stratum(width = 1/12, fill = "black", color = "grey") +geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
       scale_x_discrete(limits = c("Number of Clustered Genes", "Cluster"), expand = c(.05, .05)) +
       scale_fill_brewer(type = "qual", palette = "Set1") +
       ggtitle("Breakdown of Clustered Samples")

data_cluster <- data.frame(cutree(hc, 3))

copy_5000 <- greatest5000VarGenes.df
copy_5000$gene <- row.names(copy_5000)
copy_5000$cluster <- cluster_cut[copy_5000$gene]

copy_5000 <- copy_5000[-33]
library(pheatmap)
heatmap <- pheatmap(
  copy_5000,
  cluster_rows = TRUE, # Cluster the rows of the heatmap (genes in this case)
  cluster_cols = TRUE, # Cluster the columns of the heatmap (samples),
  show_rownames = FALSE, # There are too many genes to clearly show the labels
  main = "Non-Annotated Heatmap",
  colorRampPalette(c(
    "blue",
    "white",
    "red"
  ))(25
  ),
  scale = "row" # Scale values in the direction of genes (rows)
)

cancer_samps <- c(6, 1, 2, 0)
noncancer_sampls <- c(7, 1, 9, 6)

chi_table = rbind(cancer_samps, noncancer_sampls)
View(chi_table)
colnames(chi_table) <- c('cluster_1', 'cluster_2', 'cluster_3', 'cluster_4')
