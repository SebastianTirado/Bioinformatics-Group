library(ConsensusClusterPlus)

#Organizing data for algorithm usage
CCP5000 = greatest5000VarGenes.df[-33]
CCP10 = greatest10VarGenes.df[-33]
CCP100 = greatest100VarGenes.df[-33]
CCP1000 = greatest1000VarGenes.df[-33]
CCP10000 = greatest10000VarGenes.df[-33]
CCP5000 = data.matrix(CCP5000)
CCP10 = data.matrix(CCP10)
CCP100 = data.matrix(CCP100)
CCP1000 = data.matrix(CCP1000)
CCP10000 = data.matrix(CCP10000)

#Running CCP algorithm
results5000 = ConsensusClusterPlus(CCP5000,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
results10 = ConsensusClusterPlus(CCP10,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
results100 = ConsensusClusterPlus(CCP100,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
results1000 = ConsensusClusterPlus(CCP1000,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
results10000 = ConsensusClusterPlus(CCP10000,maxK=6,reps=50,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

#create new data frame for alluvial plot
x = data.frame("num_genes" = c(10, 10, 10, 100, 100, 100, 1000, 1000, 1000, 5000, 5000, 5000, 10000, 10000, 10000), "cluster" = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3), "count" = c(1, 2, 30, 1, 1, 31, 1, 8, 24, 1, 8, 24, 1, 8, 24))

#create alluvial plot
ggplot(x,
       aes(y=count, axis1=num_genes, axis2=cluster)) +
  geom_alluvium(aes(fill = factor(cluster)), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Number of Clustered Genes", "Cluster"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Breakdown of Clustered Samples - ConsensusClusterPlus")

#heatmap
library(pheatmap)
heatmap <- pheatmap(
  CCP5000,
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

cancer_samps <- c(6, 3, 0)
noncancer_sampls <- c(17, 5, 5)
chi_table = rbind(cancer_samps, noncancer_sampls)
colnames(chi_table) <- c('cluster_1', 'cluster_2', 'cluster_3')
chisq.test(chi_table)
