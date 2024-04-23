# Packages
library(readr)           # Read package for txt file
library(writexl)         # Write excel file .xlsx
library(WriteXLS)        # Write excel file .xlsx
library(readxl)          # Read package for xlsx file
library(limma)           # Use TMM normalization on raw counts
library(edgeR)           # Use TMM normalization on raw counts
library(plyr)
library(dplyr)
library(data.table)      # fread function
library(purrr)
library(rsample) 
library(amap)
library(ggplot2)
library(ggrepel)
library(factoextra)      # fviz
library(cluster)
library(dtwclust)
library(proxy)
library(kohonen)         # SOM package
library(pheatmap)       # Heatmap

## Clustering ##
par(mfrow=c(1,2))
set.seed(2023)

# Elbow Method
## Drug cholestrol toxin
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_drug, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters",
     ylab="Total Within sum of square",
     main="K-means",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)
### SOM
som_grid <- somgrid(xdim = 376, ydim = 10, topo = "hexagonal")   # depend with the size of the dataset (3760)
gene_som <- som(log_drug, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight1<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight1, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters",
     ylab="Total Within sum of square",
     main="SOM",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)

## Clark London
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_clark, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters",
     ylab="Total within sum of squares",
     main="k-means",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)
### SOM
som_grid <- somgrid(xdim = 3827, ydim = 1, topo = "hexagonal")
gene_som <- som(log_clark, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight2<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight2, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters",
     ylab="Total within sum of squares",
     main="SOM",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)

## ESX WT Mutant
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_ESX, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within sum of squares",
     main="K-means",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)
### SOM
som_grid <- somgrid(xdim = 1849, ydim = 2, topo = "hexagonal")
gene_som <- som(log_ESX, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight3<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight3, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within sum of squares",
     main="SOM",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)

## Inaki counts
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_Inaki, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within sum of squares",
     main="K-means",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 4, col = "red", lwd = 2)
### SOM
som_grid <- somgrid(xdim = 1976, ydim = 2, topo = "hexagonal")
gene_som <- som(log_Inaki, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight4<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight4, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within sum of squares",
     main="SOM",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 4, col = "red", lwd = 2)


## Primary drug
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_primary, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within sum of squares",
     main="k-means",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)
### SOM
som_grid <- somgrid(xdim = 1957, ydim = 2, topo = "hexagonal")
gene_som <- som(log_primary, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight5<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight5, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within sum of squares",
     main="SOM",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 4, col = "red", lwd = 2)

## timecourse nitricacid
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_timecourse, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within sum of squares",
     main="k-means",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)
### SOM
som_grid <- somgrid(xdim =3947, ydim = 1, topo = "hexagonal")
gene_som <- som(log_timecourse, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight6 <- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight6, k, nstart=50,iter.max = 40)$tot.withinss})
plot(2:kmax, wss, type="b", pch =10, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within sum of squares",
     main="SOM",
     cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
abline(v = 3, col = "red", lwd = 2)



## Cluster Size ##
# drug
## k-means
kmc1 <- kmeans(log_drug, centers = 3, iter.max=40,nstart=50)
kmc1$size
## SOM
SOM <- kmeans(weight1, centers = 3, iter.max=40,nstart=50)
SOM$size

# clark
## K-means
kmc2 <- kmeans(log_clark,centers = 3, iter.max=40,nstart=50)
kmc2$size
## SOM
SOM2 <- kmeans(weight2, centers = 3, iter.max=40,nstart=50)
SOM2$size

# ESX
## K-means
kmc3 <- kmeans(log_ESX, centers = 3, iter.max=40,nstart=50)
kmc3$size
## SOM
SOM3 <- kmeans(weight3, centers = 3, iter.max=40,nstart=50)
SOM3$size

# Inaki
## K-means
kmc4 <- kmeans(log_Inaki, centers = 4, iter.max=40,nstart=50)
kmc4$size
## SOM
SOM4 <- kmeans(weight4, centers = 4, iter.max=40,nstart=50)
SOM4$size

# Primary
## K-means
kmc5 <- kmeans(log_primary,centers = 3, iter.max=40,nstart=50)
kmc5$size
## SOM
SOM5 <- kmeans(weight5, centers = 4, iter.max=40,nstart=50)
SOM5$size

# timecourse
## K-means
kmc6 <- kmeans(log_timecourse,centers = 3, iter.max=40,nstart=50)
kmc6$size
## SOM
SOM6 <- kmeans(weight6, centers = 3, iter.max=40,nstart=50)
SOM6$size


## Cluster distribution ##
par(mfrow=c(1,1))
# K-means
## drug
fviz_cluster(kmc1, log_drug, ellipse.type = 'euclid', main="drug_cholestrol_toxin", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## clark
fviz_cluster(kmc2, log_clark, ellipse.type = 'euclid', main="Clark_london dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## ESX
fviz_cluster(kmc3, log_ESX, ellipse.type = 'euclid', main="ESX WT Mutant dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## Inaki
fviz_cluster(kmc4, log_Inaki, ellipse.type = 'euclid', main="Inaki_counts dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## primary
fviz_cluster(kmc5, log_primary, ellipse.type = 'euclid', main="primary_drug dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## timecourse
fviz_cluster(kmc6, log_timecourse, ellipse.type = 'euclid', main="timecourse_nitricacid dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)


# SOM
## drug
fviz_cluster(SOM, weight1, ellipse.type = 'euclid', main="drug_cholestrol_toxin", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## Clark
fviz_cluster(SOM2, weight2, ellipse.type = 'euclid', main="Clark_london dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## ESX
fviz_cluster(SOM3, weight3, ellipse.type = 'euclid', main="ESX WT Mutant dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## Inaki
fviz_cluster(SOM4, weight4, ellipse.type = 'euclid', main="Inaki_counts dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## Primary
fviz_cluster(SOM5, weight5, ellipse.type = 'euclid', main="primary_drug dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
## timecourse
fviz_cluster(SOM6, weight6, ellipse.type = 'euclid', main="timecourse_nitricacid dataset", xlab="Dimension 1",
             ylab="Dimension 2", cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)






























