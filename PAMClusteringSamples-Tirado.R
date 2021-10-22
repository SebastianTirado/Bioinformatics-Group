library(matrixStats)

# create new column with variance values
expressionMatrix$variance = rowVars(as.matrix(expressionMatrix[,c(-1)]))
# order data by variance decreasing
expressionMatrix = expressionMatrix[order(expressionMatrix[,33],decreasing=TRUE),]
# get subset of X amount of top values based on variance - convert to df
fiveThousandSamples.df = t(head(as.data.frame(expressionMatrix),5000)) # k = 3 optimized

tenSamples.df = t(head(as.data.frame(expressionMatrix),10)) # k = 
oneHundredSamples.df = t(head(as.data.frame(expressionMatrix),100)) # k = 
oneThousandSamples.df = t(head(as.data.frame(expressionMatrix),1000)) # k = 
tenThousandSamples.df = t(head(as.data.frame(expressionMatrix),10000)) # k = 

# PAM Clustering unsupervised analysis
library(cluster)
fiveThousandSamples.df = head(fiveThousandSamples.df,-1)
pamClust <- pam(dist(fiveThousandSamples.df), 3)
plot(pamClust, ask = TRUE)

# Alluvial Time
cluster_cut <- cutree(pamClust, 4)

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

data_cluster <- data.frame(cutree(pamClust, 3))

copy_5000 <- fiveThousandSamples.df
copy_5000$gene <- row.names(copy_5000)
copy_5000$cluster <- cluster_cut[copy_5000$gene]

copy_5000 <- copy_5000[-32]
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

cancer_samps <- c(6, 1, 2)
noncancer_sampls <- c(7, 1, 9)

chi_table = rbind(cancer_samps, noncancer_sampls)
View(chi_table)
colnames(chi_table) <- c('cluster_1', 'cluster_2', 'cluster_3')
chisq.test(chi_table)