setwd("C:/Users/willi/Downloads/#1skolaht18/1_thesis/code")
getwd()

# load

dataRaw = read.csv(file="TCGA-PANCAN-HiSeq-801x20531/data.csv",head=TRUE,sep=",")
dataClean = subset( dataRaw, select = -X ) # remove first col, "sample nr"

# normalize

data = data.frame(scale(dataClean))
data[is.na(data)] = 0

data.som = scale(dataClean)
data.som[is.na(data.som)] = 0

# set random

set.seed(426902794)
metrics.kmeans = c("correlation","euclidean", "manhattan" )
methods = c("single", "complete", "average") #"ward" does not work
metrics = c("euclidean", "manhattan")


# data desc
#summary(data)
#library(Hmisc)
#describe(data)

############
# PCA
options(max.print=1000000)
library("ggplot2") # install
library("factoextra")
library("FactoMineR")
lib = library(cluster)
library(kohonen)
library(clValid)

library(mclust) # EM
library(rgl)# 3D plot

packageVersion(lib)
sessionInfo()

facto.pca <- PCA(data, graph = FALSE)

eigenvalues <- facto.pca$eig
head(eigenvalues[, 1:2])

# plot - elbow method
#fviz_screeplot(facto.pca, ncp=120) 
#fviz_screeplot(facto.pca, ncp=31)
# take 15 & 30

pca = prcomp(data)
#print(pca)
#plot(pca)
#summary(pca)


pcas = estim_ncp(data, scale=F) # estmiate pc - return 227
print(pcas)


# plot first 15 largest pc - very dense
comp15 = data.frame(pca$x[,1:15])
#plot(comp15, pch=15, col=rgb(0,0,0,0.5))

comp227 = data.frame(pca$x[,1:227])



#############
# SOM
############
som10 = som(data.som, grid = somgrid(10, 10, "hexagonal"))
groups = 12
som10.hc = cutree(hclust(dist(som10$codes[[1]])), groups) # use dist to calc distance
plot(som20, type="codes", bgcol=rainbow(groups)[som20.hc])
add.cluster.boundaries(som20, som20.hc)

carte <- som(Z,grid=somgrid(15,10,"hexagonal"))
nb <- table(carte$unit.classif)
hclust(dc,method="ward.D2",members=nb)



###########
# correct SOM


# try 40 x 20

som28 = som(data.som, grid = somgrid(28, 28, "hexagonal")) # sqrt 801
som28.codes = getCodes(som28)  

for ( metric in metrics.kmeans) {
  kmeans = clValid(
    som28.codes, 2:12, 
    clMethods = c("kmeans","fanny"), 
    validation = "internal", 
    maxitems = 801, 
    metric = metric,
  )
  print("som 28")
  print(metric)
  summary(kmeans)
}

for ( method in methods) {
  for ( metric in metrics) {
    hier = clValid(
      som28.codes, 2:12, 
      clMethods = "hierarchical", 
      validation = "internal", 
      maxitems = 801, 
      metric = metric,
      method = method
    )
    print("som 28")
    print( method) 
    print(metric)
    summary(hier)
  }
}

## 
# end


# about measures
# Connectivity - minmize to 0
# dunn - maximize from 0
# silhouette - close to 1

###########
# change sink
###########

con <- file("results/results.txt")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")
source("test1.R", echo=TRUE, max.deparse.length=10000)

# Restore output to console

sink() 
sink(type="message")


############# 
# K-Means
############


for ( metric in metrics.kmeans) {
  kmeans = clValid(
    data, 2:12, 
    clMethods = c("kmeans","fanny"), 
    validation = "internal", 
    maxitems = 801, 
    metric = metric,
  )
  print("baseline")
  print(metric)
  summary(kmeans)
}
for ( metric in metrics.kmeans) {
  kmeans = clValid(
    comp227, 2:12, 
    clMethods = c("kmeans","fanny"), 
    validation = "internal", 
    maxitems = 801, 
    metric = metric,
  )
  print("pc 227")
  print(metric)
  summary(kmeans)
}
for ( metric in metrics.kmeans) {
  kmeans = clValid(
    comp15, 2:12, 
    clMethods = c("kmeans","fanny"), 
    validation = "internal", 
    maxitems = 801, 
    metric = metric,
  )
  print("pc 15")
  print(metric)
  summary(kmeans)
}





############
# Hierarchical agglomerative
###########

for ( method in methods) {  
  for ( metric in metrics) {
    hier = clValid(
    data, 2:12, 
    clMethods = "hierarchical", 
    validation = "internal", 
    maxitems = 801, 
    metric = metric,
    method = method
   )
  print("baseline")
  print( method) 
  print(metric)
  summary(hier)
  }
}
for ( method in methods) {
  for ( metric in metrics) {
    hier = clValid(
      comp227, 2:12, 
      clMethods = "hierarchical", 
      validation = "internal", 
      maxitems = 801, 
      metric = metric,
      method = method
    )
    print("pc 227")
    print( method) 
    print(metric)
    summary(hier)
  }
}
for ( method in methods) {
  for ( metric in metrics) {
    hier = clValid(
    comp15, 2:12, 
    clMethods = "hierarchical", 
    validation = "internal", 
    maxitems = 801, 
    metric = metric,
    method = method
   )
  print("pc 15")
  print( method) 
  print(metric)
  summary(hier)
  }
}


