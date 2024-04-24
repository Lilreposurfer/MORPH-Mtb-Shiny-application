## Normalization, Filtering and histograms after normalization

log <- function(data){
  a <- DGEList(data, group = NULL) # create DGEList object
  b <- calcNormFactors(a, method = "TMM") # perform TMM normalization
  norm_data <- cpm(b, log=FALSE) # retrieve normalized counts
  sd_expr <- apply(norm_data, 1, sd)   # SD for each gene
  threshold <- 1    # threshold
  norm_data_filtered <- norm_data[sd_expr >= threshold, ] # remove gene with sd<1
  log_data <- log2(norm_data_filtered + 1)
  log_data
}


wsskmeans <- function(data){
# Elbow Method
## K-means
wssk <- sapply(2:kmax, function(k){
  kmeans(data, k, nstart=50,iter.max = 40)$tot.withinss})
wssk
}



wsssom <- function(data){
# Elbow Method
## SOM
som_grid <- somgrid(xdim = 376, ydim = 10, topo = "hexagonal")   # depend with the size of the dataset (3760)
gene_som <- som(data, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight1<- getCodes(gene_som)
#kmax <- 10
wsss <- sapply(2:kmax, function(k){
  kmeans(weight1, k, nstart=50,iter.max = 40)$tot.withinss})
wsss
}

