#getwd()
setwd("C:/Users/willi/Downloads/#1skolaht18/1_thesis/code")

# load
dataRaw = read.csv(file="TCGA-PANCAN-HiSeq-801x20531/data.csv",head=TRUE,sep=",")
dataClean = subset( dataRaw, select = -X )


# normalize

data = data.frame(scale(dataClean))
data[is.na(data)] = 0

data.short = dataClean[1:100,]# select fist 100 rows

data.short.som = scale(data.short)
data.short.som[is.na(data.short.som)] = 0

desc = data[1:5,1:5]


# set random
set.seed(426902794)


# data desc
#summary(data)
#library(Hmisc)
#describe(data)

# PCA
options(max.print=1000000)
library("factoextra")
library("FactoMineR")
library(cluster)
library(kohonen)
library(clValid)

library(mclust) # EM

facto.pca <- PCA(data.short, graph = FALSE)

eigenvalues <- facto.pca$eig
head(eigenvalues[, 1:2])

# plot - elbow method
#fviz_screeplot(facto.pca, ncp=120) 
#fviz_screeplot(facto.pca, ncp=31)
# take 15 & 30

pca = prcomp(data.short)
#print(pca)
#plot(pca)
#summary(pca)

pcas = estim_ncp(data.short, scale=F) # estmiate pc - return 28
print(pcas)




# plot first 15 largest pc - very dense
comp15 = data.frame(pca$x[,1:15])
plot(comp15, pch=15, col=rgb(0,0,0,0.5))

# plot first 30 largest pc - very dense
comp28 = data.frame(pca$x[,1:28])
#plot(comp28, pch=28, col=rgb(0,0,0,0.5))


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
# SOM
############


som10 = som(data.short.som, grid = somgrid(10, 10, "hexagonal")) # sqrt 100
som10.codes = getCodes(som10)  

for ( metric in metrics.kmeans) {
  kmeans = clValid(
    som10.codes, 2:12, 
    clMethods = c("kmeans","fanny"), 
    validation = "internal", 
    maxitems = 801, 
    metric = metric,
  )
  print("som 10")
  print(metric)
  summary(kmeans)
}

for ( method in methods) {
  for ( metric in metrics) {
    hier = clValid(
      som10.codes, 2:12, 
      clMethods = "hierarchical", 
      validation = "internal", 
      maxitems = 801, 
      metric = metric,
      method = method
    )
    print("som 10")
    print( method) 
    print(metric)
    summary(hier)
  }
}


# try 
method = ward.D2


#############
# 
#############



############# 
# K-Means
############
metrics.kmeans = c("correlation","euclidean", "manhattan" )

# add fanny

for ( metric in metrics.kmeans) {
  kmeans = clValid(
    data.short, 2:12, 
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
    comp28, 2:12, 
    clMethods = c("kmeans","fanny"), 
    validation = "internal", 
    maxitems = 801, 
    metric = metric,
  )
  print("pc 28")
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
methods = c("single", "complete", "average") #"ward" does not work
metrics = c("euclidean", "manhattan")


for ( method in methods) {
  for ( metric in metrics) {
    hier = clValid(
      data.short, 2:12, 
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
      comp28, 2:12, 
      clMethods = "hierarchical", 
      validation = "internal", 
      maxitems = 801, 
      metric = metric,
      method = method
    )
    print("pc 28")
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



