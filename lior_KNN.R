#This script creates a matrix of the counts of genes for each sample
library(dplyr)
getwd()
setwd("Bioinformatics-Group/GSE143791_condensed")
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
View(range_df)

### Generating DESeq2 object ###
setwd("./..")

# holds list of sample names a factor which puts them into tumor or non-tumor group
group_classification <- read.csv("group_classification.csv",header = T,row.names=1,
                                 stringsAsFactors=T)

# filter out genes that are not highly expressed
filt_expression_matrix <- expressionMatrix %>% dplyr::filter(rowSums(.) >= 10)

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = filt_expression_matrix,
                              colData = group_classification,
                              design = ~ dex)

deseq_object <- DESeq(dds)
deseq_results <- results(deseq_object)

BiocManager::install("apeglm")

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

BiocManager::install("EnhancedVolcano")

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

BiocManager::install("pheatmap")
library(pheatmap)

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

BiocManager::install("gprofiler2")
library(gprofiler2)

gostres <- gost(query = gene_names, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

# STARTING HERE - USE EXPRESSION MATRIX TO GENERATE VARIANCE AND MAKE NEW DATA FRAME
# OF TOP X AMOUNT OF GENES WITH HIGHEST VARIANCE

install.packages("matrixStats") # only need to install once and then ignore these lines
library(matrixStats)

# create new column with variance values
expressionMatrix$variance = rowVars(as.matrix(expressionMatrix[,c(-1)]))
# order data by variance decreasing
expressionMatrix = expressionMatrix[order(expressionMatrix[,33],decreasing=TRUE),]
# get subset of X amount of top values based on variance - convert to df
greatestVarGenes.df = head(as.data.frame(expressionMatrix),5000)
# RUN INDIVIDUAL CLUSTERING ALGORITHM - PUSH TO GITHUB
greatestVarGenes.df <- t(greatestVarGenes.df)
greatestVarGenes.df <- head(greatestVarGenes.df, - 1) 
View(greatestVarGenes.df)

BiocManager::install("Seurat")



# library(class)
# normalize <- function(x) {
#   return ((x - min(x)) / (max(x) - min(x))) }
# #gvg_n <- as.data.frame(lapply(greatestVarGenes.df, normalize))
# #View(gvg_n)
# 
# gvg_train <- greatestVarGenes.df[1:20,]
# gvg_test <- greatestVarGenes.df[21:32,]
# View(group_classification)
# gvg_train_labels <- group_classification[1:20, 1]
# gvg_test_labels <- group_classification[21:32, 1]  
# head(gvg_train_labels)
# gvg_test_pred <- knn(train = gvg_train, test = gvg_test ,cl = gvg_train_labels, k=2)
# gvg_test_pred

BiocManager::install("ClusterR")
library(gtools)
library(ClusterR)
# opt_gmm = Optimal_Clusters_GMM(greatestVarGenes.df, max_clusters = 10, criterion = "BIC", 
#                                
#                                dist_mode = "maha_dist", seed_mode = "random_subset",
#                                
#                                km_iter = 10, em_iter = 10, var_floor = 1e-10, 
#                                
#                                plot_data = T)
# 
# gmm = GMM(greatestVarGenes.df, 6, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 10,
#           
#           em_iter = 10, verbose = F)          
# 
# # predict centroids, covariance matrix and weights
# #View(gmm)
# pr = predict_GMM(greatestVarGenes.df, gmm$centroids, gmm$covariance_matrices, gmm$weights) 
# View(pr)

opt = Optimal_Clusters_KMeans(greatestVarGenes.df, max_clusters = 10, plot_clusters = T,
                              
                              criterion = 'distortion_fK', fK_threshold = 0.85,
                              
                              initializer = 'optimal_init', tol_optimal_init = 0.2)

kmeans_clust = KMeans_rcpp(greatestVarGenes.df, clusters = 2, num_init = 5, max_iters = 100, 
                    
                    initializer = 'optimal_init', verbose = F)
kmeans_clust[["clusters"]]
greatestVarGenes.df.clustered <- greatestVarGenes.df[c(4,6,7,8,9,13,14,15,18,19,22,26,27,29,30,31,32,1,2,3,5,10,11,12,16,17,20,21,23,24,25,28),]
View(greatestVarGenes.df.clustered)
View(greatestVarGenes.df)
greatestVarGenes.df.clustered <- t(greatestVarGenes.df.clustered)
View(greatestVarGenes.df.clustered)
library(pheatmap)
heatmap <- pheatmap(
  greatestVarGenes.df.clustered,
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


greatest10VarGenes.df = head(as.data.frame(expressionMatrix), 10)
greatest100VarGenes.df = head(as.data.frame(expressionMatrix), 100)
greatest1000VarGenes.df = head(as.data.frame(expressionMatrix), 1000)
greatest10000VarGenes.df = head(as.data.frame(expressionMatrix), 10000)
greatest10VarGenes.df <- t(greatest10VarGenes.df)
greatest10VarGenes.df <- head(greatest10VarGenes.df, -1)
greatest100VarGenes.df <- t(greatest100VarGenes.df)
greatest100VarGenes.df <- head(greatest100VarGenes.df, -1)
greatest1000VarGenes.df <- t(greatest1000VarGenes.df)
greatest1000VarGenes.df <- head(greatest1000VarGenes.df, -1)
greatest10000VarGenes.df <- t(greatest10000VarGenes.df)
greatest10000VarGenes.df <- head(greatest10000VarGenes.df, -1)
View(greatest10VarGenes.df)
opt10 = Optimal_Clusters_KMeans(greatest10VarGenes.df, max_clusters = 10, plot_clusters = T,
                              
                              criterion = 'distortion_fK', fK_threshold = 0.85,
                              
                              initializer = 'optimal_init', tol_optimal_init = 0.2)
k10 <- 7
opt100 = Optimal_Clusters_KMeans(greatest100VarGenes.df, max_clusters = 10, plot_clusters = T,
                                       
                                       criterion = 'distortion_fK', fK_threshold = 0.85,
                                       
                                       initializer = 'optimal_init', tol_optimal_init = 0.2)
k100 <- 7
opt1000 = Optimal_Clusters_KMeans(greatest1000VarGenes.df, max_clusters = 10, plot_clusters = T,
                              
                              criterion = 'distortion_fK', fK_threshold = 0.85,
                              
                              initializer = 'optimal_init', tol_optimal_init = 0.2)
k1000 <- 2

opt10000 = Optimal_Clusters_KMeans(greatest10000VarGenes.df, max_clusters = 10, plot_clusters = T,
                              
                              criterion = 'distortion_fK', fK_threshold = 0.85,
                              
                              initializer = 'optimal_init', tol_optimal_init = 0.2)

k10000 <- 2

numGenes <- c(10,100,1000,10000)
numClusters <- c(7,7,2,2)
genesToClusters <- data.frame(numGenes, numClusters)
View(genesToClusters)

BiocManager::install("ggalluvial")
library(ggplot2)
library(alluvial)
library(ggalluvial)

ggplot(genesToClusters,
       aes(y = numClusters, axis1 = numGenes, axis2 = numClusters)) +
  geom_alluvium(aes(fill = factor(numClusters)), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Number of Clustered Genes", "Cluster"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Breakdown of Clustered Samples")

View(kmeans_clust)

cluster_results <- kmeans_clust[["clusters"]] #2 is non tumor, 1 is tumor
true_labels <- c(1,2,2,2,1,1,2,2,1,2,2,1,2,2,1,2,2,2,2,1,1,2,2,1,2,2,2,2,2,2,2,2)
results <- data.frame(true_labels, cluster_results)
View(results)

test <- chisq.test(table(results$true_labels, results$cluster_results))
test
