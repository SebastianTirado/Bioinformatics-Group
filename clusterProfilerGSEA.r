library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(tidyverse)
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

df = read_tsv("diff_expression_results.tsv", col_names=TRUE)
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$Gene
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing=TRUE)

gse <- gseGO(geneList = gene_list, 
             ont = "ALL", 
             keyType = "GENENAME", 
             nPerm = 1000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

summary(gse)
