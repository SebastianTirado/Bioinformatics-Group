library(ConsensusClusterPlus)
title = tempdir()

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
