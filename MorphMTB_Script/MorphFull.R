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
library(pheatmap)        # Heatmap
library(caroline) 

# Seeds for reproducibility
set.seed(2023)

## Pre-proccessing ##
# Loading data
#setwd("Path to the files")
setwd("C:/Users/elise/Documents/Schuul/Stage/Project/MorphMTB_Script/data")
drug <- read.table("drug_cholestrol_toxin.txt", header = TRUE, row.names = 1)
clark <- read.table("Clark_london.txt", header = TRUE, row.names = 1)
ESX <- read.table("ESX_WT_mutant.txt", header = TRUE, row.names = 1)
Inaki <- read.table("Inaki_counts.txt", header = TRUE, row.names = 1)
primary <- read.table("primary_data.txt", header = TRUE, row.names = 1)
timecourse <- read.table("timecourse_nitricacid.txt", header = TRUE, row.names = 1)

# drug_cholestrol_toxin
## Normalization, Filtering and histograms after normalization
a <- DGEList(counts = drug, group = NULL) # create DGEList object
a <- calcNormFactors(a, method = "TMM") # perform TMM normalization
norm_drug <- cpm(a, log=FALSE) # retrieve normalized counts
sd_expr <- apply(norm_drug, 1, sd)   # SD for each gene
threshold <- 1    # threshold
norm_drug_filtered <- norm_drug[sd_expr >= threshold, ] # remove gene with sd<1
log_drug <- log2(norm_drug_filtered + 1)

# Clark_london
##  Normalization, Filtering and histograms after normalization
b <- DGEList(counts = clark, group = NULL) # create DGEList object
b <- calcNormFactors(b, method = "TMM") # perform TMM normalization
norm_clark <- cpm(b, log=FALSE) # retrieve normalized counts
sd_expr <- apply(norm_clark, 1, sd)   # SD for each gene
threshold <- 1    # threshold
norm_clark_filtered <- norm_clark[sd_expr >= threshold, ] # remove gene with sd<1
log_clark <- log2(norm_clark_filtered + 1)

# ESX_WT_Mutant
##  Normalization, Filtering and histograms after normalization 
c <- DGEList(counts = ESX, group = NULL) # create DGEList object
c <- calcNormFactors(c, method = "TMM") # perform TMM normalization
norm_ESX <- cpm(c, log=FALSE) # retrieve normalized counts
sd_expr <- apply(norm_ESX, 1, sd)   # SD for each gene
threshold <- 1    # threshold
norm_ESX_filtered <- norm_ESX[sd_expr >= threshold, ] # remove gene with sd<1
log_ESX <- log2(norm_ESX_filtered + 1)

# Inaki_counts
## Normalization, Filtering and histograms after normalization
d <- DGEList(counts = Inaki, group = NULL) # create DGEList object
d <- calcNormFactors(d, method = "TMM") # perform TMM normalization
norm_Inaki <- cpm(d, log=FALSE) # retrieve normalized counts
sd_expr <- apply(norm_Inaki, 1, sd)   # SD for each gene
threshold <- 1    # threshold
norm_Inaki_filtered <- norm_Inaki[sd_expr >= threshold, ] # remove gene with sd<1
log_Inaki <- log2(norm_Inaki_filtered + 1)

# Primary_drug
## Normalization, Filtering and histograms after normalization
e <- DGEList(counts = primary, group = NULL) # create DGEList object
e <- calcNormFactors(e, method = "TMM") # perform TMM normalization
norm_primary <- cpm(e, log=FALSE) # retrieve normalized counts
sd_expr <- apply(norm_primary, 1, sd)   # SD for each gene
threshold <- 1    # threshold
norm_primary_filtered <- norm_primary[sd_expr >= threshold, ] # remove gene with sd<1
log_primary <- log2(norm_primary_filtered + 1)

# timecourse_nitricacid
## Normalization, Filtering and histograms after normalization
f <- DGEList(counts = timecourse, group = NULL) # create DGEList object
f <- calcNormFactors(f, method = "TMM") # perform TMM normalization
norm_timecourse <- cpm(f, log=FALSE) # retrieve normalized counts
sd_expr <- apply(norm_timecourse, 1, sd)   # SD for each gene
threshold <- 1    # threshold
norm_timecourse_filtered <- norm_timecourse[sd_expr >= threshold, ] # remove gene with sd<1
log_timecourse <- log2(norm_timecourse_filtered + 1)

################################################################################################"

## Clustering ##

# Elbow Method
## Drug cholestrol toxin
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_drug, k, nstart=50,iter.max = 40)$tot.withinss})
### SOM
som_grid <- somgrid(xdim = 376, ydim = 10, topo = "hexagonal")   # depend with the size of the dataset (3760)
gene_som <- som(log_drug, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight1<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight1, k, nstart=50,iter.max = 40)$tot.withinss})

## Clark London
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_clark, k, nstart=50,iter.max = 40)$tot.withinss})
### SOM
som_grid <- somgrid(xdim = 3827, ydim = 1, topo = "hexagonal")
gene_som <- som(log_clark, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight2<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight2, k, nstart=50,iter.max = 40)$tot.withinss})

## ESX WT Mutant
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_ESX, k, nstart=50,iter.max = 40)$tot.withinss})
### SOM
som_grid <- somgrid(xdim = 1849, ydim = 2, topo = "hexagonal")
gene_som <- som(log_ESX, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight3<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight3, k, nstart=50,iter.max = 40)$tot.withinss})

## Inaki counts
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_Inaki, k, nstart=50,iter.max = 40)$tot.withinss})
### SOM
som_grid <- somgrid(xdim = 1976, ydim = 2, topo = "hexagonal")
gene_som <- som(log_Inaki, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight4<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight4, k, nstart=50,iter.max = 40)$tot.withinss})

## Primary drug
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_primary, k, nstart=50,iter.max = 40)$tot.withinss})
### SOM
som_grid <- somgrid(xdim = 1957, ydim = 2, topo = "hexagonal")
gene_som <- som(log_primary, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight5<- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight5, k, nstart=50,iter.max = 40)$tot.withinss})

## timecourse nitricacid
### K-means
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(log_timecourse, k, nstart=50,iter.max = 40)$tot.withinss})
### SOM
som_grid <- somgrid(xdim =3947, ydim = 1, topo = "hexagonal")
gene_som <- som(log_timecourse, grid = som_grid, rlen = 100, alpha = c(0.05, 0.01), keep.data = TRUE)
weight6 <- getCodes(gene_som)
kmax <- 10
wss <- sapply(2:kmax, function(k){kmeans(weight6, k, nstart=50,iter.max = 40)$tot.withinss})


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


#################################################################################################

## Data Preparation for MORPH ##
# Loading Pathways
Pathway1 <- read_excel("Pathways_virulence.xlsx",  sheet = "Pathway1_allVirulaence", col_names = FALSE)
Pathway2 <- read_excel("Pathways_virulence.xlsx",  sheet = "Pathway2", col_names = FALSE)
Pathway3 <- read_excel("Pathways_virulence.xlsx",  sheet = "Pathway3", col_names = FALSE)
Pathway4 <- read_excel("Pathways_virulence.xlsx",  sheet = "pathway4", col_names = FALSE)
Pathway5 <- read_excel("Pathways_virulence.xlsx",  sheet = "Pathway5", col_names = FALSE)
# Loading gene expression
drug <- read.table("drug_cholestrol_toxin.txt", header = TRUE)
clark <- read.table("Clark_london.txt", header = TRUE)
ESX <- read.table("ESX_WT_mutant.txt", header = TRUE)
Inaki <- read.table("Inaki_counts.txt", header = TRUE)
primary <- read.table("primary_data.txt", header = TRUE)
timecourse <- read.table("timecourse_nitricacid.txt", header = TRUE)


# Preparing 
## Gene Expression
write.delim(drug, "drug.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.delim(clark, "clark.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.delim(Inaki, "inaki.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.delim(ESX, "ESX.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.delim(primary, "primary.txt", sep="\t", col.names=FALSE, row.names=FALSE)
write.delim(timecourse, "timecourse.txt", sep="\t", col.names=FALSE, row.names=FALSE)

## Config
A=c("clark.txt","clark.txt", "drug.txt", "drug.txt", "ESX.txt", "ESX.txt", "inaki.txt", 
    "inaki.txt", "primary.txt", "primary.txt", "timecourse.txt", "timecourse.txt")
B=c("kmeansclark.txt", "somclark.txt", "kmeansdrug.txt", "somdrug.txt", 
    "kmeansESX.txt", "somESX.txt", "kmeansinaki.txt", "sominaki.txt",
    "kmeansprimary.txt", "somprimary.txt", "kmeanstimecourse.txt", "somtimecourse.txt")
Configs=data.frame(A,B)
write.delim(Configs, "Configs.txt", sep="\t", col.names=FALSE, row.names=FALSE)

## Pathways
Pathway1=Pathway1$...1
Pathway2=Pathway2$...1
Pathway3=Pathway3$...1
Pathway4=Pathway4$...1
Pathway5=Pathway5$...1
writeLines(Pathway1, "Pathway1.txt", sep="\t")
writeLines(Pathway2, "Pathway2.txt", sep="\t")
writeLines(Pathway3, "Pathway3.txt", sep="\t")
writeLines(Pathway4, "Pathway4.txt", sep="\t")
writeLines(Pathway5, "Pathway5.txt", sep="\t")

## K-means clusters
### drug
genes1 <- rownames(log_drug)
cluster1 <- kmc1$cluster
kmeans_cluster1 <- data.frame(genes1, cluster1)
write.delim(kmeans_cluster1, "kmeansdrug.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### clark
genes2 <- rownames(log_clark)
cluster2 <- kmc2$cluster
kmeans_cluster2 <- data.frame(genes2, cluster2)
write.delim(kmeans_cluster2, "kmeansclark.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### ESX
genes3 <- rownames(log_ESX)
cluster3 <- kmc3$cluster
kmeans_cluster3 <- data.frame(genes3, cluster3)
write.delim(kmeans_cluster3, "kmeansESX.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### Inaki
genes4 <- rownames(log_Inaki)
cluster4 <- kmc4$cluster
kmeans_cluster4 <- data.frame(genes4, cluster4)
write.delim(kmeans_cluster4, "kmeansinaki.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### primary
genes5 <- rownames(log_primary)
cluster5 <- kmc5$cluster
kmeans_cluster5 <- data.frame(genes5, cluster5)
write.delim(kmeans_cluster5, "kmeansprimary.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### Timecourse
genes6 <- rownames(log_timecourse)
cluster6 <- kmc6$cluster
kmeans_cluster6 <- data.frame(genes6, cluster6)
write.delim(kmeans_cluster6, "kmeanstimecourse.txt", sep="\t", col.names=FALSE, row.names=FALSE)


## SOM clusters
### drug
gene_names1 <- rownames(log_drug)
cluster_assignments1 <- SOM$cluster
gene_clusters1 <- data.frame(gene = gene_names1, cluster = cluster_assignments1)
write.delim(gene_clusters1, "somdrug.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### clark
gene_names2 <- rownames(log_clark)
cluster_assignments2 <- SOM2$cluster
gene_clusters2 <- data.frame(gene = gene_names2, cluster = cluster_assignments2)
write.delim(gene_clusters2, "somclark.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### ESX
gene_names3 <- rownames(log_ESX)
cluster_assignments3 <- SOM3$cluster
gene_clusters3 <- data.frame(gene = gene_names3, cluster = cluster_assignments3)
write.delim(gene_clusters3, "somESX.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### Inaki
gene_names4 <- rownames(log_Inaki)
cluster_assignments4 <- SOM4$cluster
gene_clusters4 <- data.frame(gene = gene_names4, cluster = cluster_assignments4)
write.delim(gene_clusters4, "sominaki.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### primary
gene_names5 <- rownames(log_primary)
cluster_assignments5 <- SOM5$cluster
gene_clusters5 <- data.frame(gene = gene_names5, cluster = cluster_assignments5)
write.delim(gene_clusters5, "somprimary.txt", sep="\t", col.names=FALSE, row.names=FALSE)
### timecourse
gene_names6 <- rownames(log_timecourse)
cluster_assignments6 <- SOM6$cluster
gene_clusters6 <- data.frame(gene = gene_names6, cluster = cluster_assignments6)
write.delim(gene_clusters6, "somtimecourse.txt", sep="\t", col.names=FALSE, row.names=FALSE)

## All genes for sampling
g = drug$Gene_ID

###############################################################################################

## MORPH ALGORITHM ##
# Functions
## Reading data functions
### [1.1] getClusteringInformation
getClusteringInformation = function(ClusterFile)
{
  Cluster = read.delim(ClusterFile, sep="\t", header=FALSE, row.names = 1,stringsAsFactors=F) #Read tab-delimited file.
  colnames(Cluster) = ClusterFile #Change cluster column name to file-name for future identification.
  return(Cluster)
}

### [1.2] Get Pathways
GetGOIs = function(GOI_File)
{
  GOI = readLines(GOI_File, n=1) #Reads the first (and only) line from the pathway-genes file.
  GOI = strsplit(GOI,"\t")[[1]] #Splits names into a vector by tabs.
  return(GOI)
}


### [1.3] getGeneExpression
getGeneExpression = function(InputGE)
{
  GE = (read.delim(InputGE, sep="\t", header=TRUE,row.names=1)) #Reads the gene-expression data from the tab-delimited file.
}


## MORPH Algorithm
### [2.1] prepare Morph Object From Files {Input;Configuration file (data and cluster solution), Output:MORPH object}
prepareMorphObjectFromFiles = function(InputConfig= NULL,InputGOI = NULL) {
  Config = read.delim(InputConfig, sep = "\t", header=FALSE) #Reads the configs.txt file.
  List_GE = as.character(Config[,1]) #Reads the first column (containing paths to gene-expression data files)
  List_C = as.character(Config[,2]) #Reads the second column (containing paths to clustering solution files)
  G = c() #Initialize the vector to contain names of pathway-genes.
  if (! is.null(InputGOI)){ #Is the path to the pathway-genes file defined?
    G = GetGOIs(InputGOI) #Reads file into variable.
  }
  ge_datasets = list() #Initialize list to contain all gene-expression matrices.
  clustering_solutions = list() #Initialize list to contain all clustering solutions.
  #Go over all the datasets and solutions pairs, and insert them to the data structures.
  for (I in (1:length(List_C))) {
    cl = List_C[[I]]
    ge = List_GE[[I]]
    #Read the dataset only if it still doesn't exists in ge_datasets.
    #As each dataset can appear many time in the config file, we insert it only once to the data
    #structure (preformance reasons)
    if (is.null(ge_datasets[[ge]])) {
      ge_data = getGeneExpression(ge) #Read data into variable.
      ge_datasets[[ge]] = ge_data #Add gene-expression dataset to the list.
    }
    #Insert the current solution to the clustering solutions data stucture.
    #The processing of the solution is against the current dataset.
    #Clustering solution should be unique to each dataset (recommended).
    raw_cl =  getClusteringInformation(cl) #Read data into variable.
    genes_to_keep = intersect(rownames(raw_cl),rownames(ge_data)) #Keep only the genes that appear both in the gene-expression data AND the clustering solution.
    raw_cl = as.matrix(raw_cl[genes_to_keep,],ncol=1) #Extract the kept-genes' lines as a matrix.
    rownames(raw_cl) = genes_to_keep #Rename matrix's rows accordingly.
    raw_cl_genes = rownames(raw_cl) #Extract names genes in clustering solution.
    raw_cl_clusters = raw_cl[,1] #Extract parent clusters of genes in clustering solution.
    unclustered = setdiff(rownames(ge_data),raw_cl_genes) #Detect unclustered genes (those who appear in gene-expression but NOT in clustering solution).
    if (length(unclustered) == 0){ #Have we met NO unclustered genes?
      clustering_solutions[[cl]] = raw_cl #Add clustering solution data to the list.
      next #Move to next clustering solution.
    }
    #Have we met unclustered genes?
    uncl = rep ("unclustered",length(unclustered)) #Set their parent clusters as "unclustered".
    clustering = matrix(c(raw_cl_clusters,uncl),ncol=1) #Create a matrix with the unclustered genes at the bottom.
    rownames(clustering) = c(raw_cl_genes,unclustered) #Rename the rows according to gene names.
    clustering_solutions[[cl]] = clustering #Add clustering solution data to the list.
  }									                          
  #Create and return a MORPH object (see function description).
  morph_obj = list(config = Config, clustering_solutions = clustering_solutions, ge_datasets = ge_datasets, pathway_genes = G)
  class(morph_obj)<-"morph_data" # equvallent to Java class
  return (morph_obj)}

### [2.2] rankGenes {Input;G-pathway, C-cluster solution, GE-dataset, Output; Scores}
rankGenes = function(G,C,GE,corrs_mat = NULL) #Ranks all candidates-vs-GOIs for a "data-couple" (GE data & Clustering solution)
{if (is.null(corrs_mat)){ #Have we NOT received a coexpression (covariance) matrix (is it defined)?
  corrs_mat = cor(t(GE),t(GE[G,])) #Calculate the covariance matrix from scratch.
  print(corrs_mat)
}
  print(corrs_mat)
  G_Clusters = matrix(C[G,1], length(G),1) #Detect the parent cluster for all pathway-genes.
  rownames(G_Clusters) = G #Change the row names according to the pathway-genes.
  RelevantClusters =na.omit(unique(G_Clusters[,1]))  #Keep only the clusters that contain pathway-genes.
  Scores = list() #Initialize list to contain results.
  Scores[["Scored"]] = c() #Initialize the "Scored" field that will contain all valid calculated scores.
  for (Cluster in RelevantClusters) #Go through each cluster within those who contain pathway-genes.
  { #Acquire cluster-specific ("current") pathway-genes and candidates.
    CurrentG = na.omit(names(G_Clusters[(G_Clusters[,1]==Cluster),]))  #Acquire current pathway-genes
    CurrentB = rownames(C)[which(C[,1]==Cluster)] #Acquire current background genes.
    Candidates = setdiff(CurrentB,CurrentG) #setdiff function filters Candidates are genes which are present at the background but NOT in our pathway-genes group.
    #print(CurrentG)
    #print(CurrentB)
    #print(Candidates)
    
    #Create an ordered vector containing the normalized covariance scores for all candidates vs. pathway-genes.
    NormalizedCorrelations = getNormalizedCorrelations(corrs_mat, CurrentG, Candidates) #Calculate normalized scores.
    NormalizedCorrelations = NormalizedCorrelations[order(NormalizedCorrelations,decreasing=TRUE)] #Sort list in descending order.
    Scores[["Scored"]] = c(Scores[["Scored"]],NormalizedCorrelations) #Add normalized scores to the "Scored" field of our list.
  }
  # All genes that were not clustered with any pathway gene are placed at the bottom of the ranking (normalized score = NA)
  NonRelevantClusters = setdiff(unique(C[,1]),RelevantClusters) #Detect clusters that do NOT contain pathway-genes (non-relevant clusters).
  NonCandidates = character(0) #Initialize variable to contain names of genes who are not clustered with pathway-genes (non-candidates).
  for (NRC in NonRelevantClusters){ #Go through all non-relevant clusters.
    NC = rownames(C)[which(C[,1]==NRC)] #Extract names of all genes in current cluster.
    NonCandidates = c(NonCandidates,NC) #Add the names to the non-candidates list.
  }
  LastScores = rep(NA,length(NonCandidates)) #Set scores of all non-candidates as NA.
  names(LastScores) = NonCandidates #Rename vector according to non-candidates' names.
  Scores[["Rejected"]] = LastScores #Add non-candidates to the "Rejected" field of our list.
  #Sort all scores in a descending order.
  scores_vector = Scores[["Scored"]]
  scores_vector = scores_vector[order(scores_vector,decreasing=T)]
  Scores[["Scored"]] = scores_vector
  return(Scores)}


### [2.3] Normalized correlation {Input; correlation matrix, pathway, genes, Outputl;Ordered vector}

#Acquire cluster-specific ("current") pathway-genes and candidates.
getNormalizedCorrelations = function(corrs_mat, CurrentG, Candidates)
{AverageCorrelations = c() #Initialize variable to contain average correlations.
if (length(CurrentG)==1){ #Is there only ONE pathway gene in this cluster?
  AverageCorrelations = corrs_mat[Candidates,CurrentG] #Set average correlation as it's coexpression (covariance) score from our matrix.
}
else{ #Are there MORE than one pathway genes in this cluster?
  AverageCorrelations = apply(corrs_mat[Candidates,CurrentG],1,mean) #Calculate the mean of their coexpression (covariance) scores.
}
#Normalize and return the scores.
names(AverageCorrelations)<-Candidates #Rename the average correlations according to the appropriate candidates' names.
Mean = mean(AverageCorrelations) #Calculate the mean of all average scores.
SD = sd(AverageCorrelations) #Calculate the standard deviation of all average scores.
NormalizedCorrelations = ((AverageCorrelations - Mean) / SD) #Normalize scores according to mean and standard deviation.
return(NormalizedCorrelations)}


### [2.4] LOOCV (Leave-One-Out Cross Validation) {Input; G, C, GE, Output;Score}
LOOCV = function(G,C,GE,K,corrs_mat = NULL, NameGE = NULL, NameC = NULL)
{if (is.null(corrs_mat)){ #Have we NOT received a coexpression (covariance) matrix (is it defined)?
  corrs_mat = cor(t(GE),t(GE[G,])) #Calculate the covariance matrix from scratch.
}
  SelfRanks = numeric(0) #Initialize vector to contain self-ranks.
  for (V in G) #Go through all pathway-genes and leave one out at a time.
  { S = setdiff(G,V) #Reduce pathway genes group to contain all pathway genes EXCEPT the one left out.
  Scores = rankGenes(S,C,GE,corrs_mat) #Run rankGenes function with the reduced pathway-genes.
  Rank = which(names(Scores$Scored)==V) #Acquire the location of the left-out gene among the others.
  if (length(Rank)==0 || is.nan(Rank)){Rank = (2*K)} #Was the left-out genes NOT detected at all? Set score to 2K as penalty.
  SelfRanks=c(SelfRanks,Rank) #Add self-rank to our vector.
  }
  return(AUSR(SelfRanks,K,NameGE,NameC))}

### [2.5] AUSR (Area Under the curve of the Self-Ranked genes) {Input; Self-ranks, K-threshold, Output;Score}
AUSR = function(SelfRanks, K, NameGE = NULL, NameC = NULL)
{nGOIs = length(SelfRanks) #Count the number of self-ranks in the vector.
Fractions = numeric(0) #Initialize vector to contain fractions (i.e. fractions of pathway genes with self rank below threshold)
AUC = 0 #Initialize accumulative variable to contain Area-Under-the-Curve score.
PlotFlag = 0 #Initialize flag to indicate whether plotting should be executed.
if ((!is.null(NameGE)) & (!is.null(NameC)))	{ #Have we received both the names of gene-expression data file AND the clustering solution file?
  PlotFlag = 1 #User requested for plotting, set flag to 1.
}
for (Index in 1:K) { #Go through threshold calculation-range (1 to K by 1)
  Fraction = ((length(SelfRanks[SelfRanks<=Index])) / nGOIs) #Calculate the fraction of genes that have a self-rank smaller or equal to K of current iteration.
  Fractions = c(Fractions,Fraction) #Add the fraction score to our vector.
  AUC = AUC + Fraction #Accumulate fraction score as AUC score.
}
if (PlotFlag == 1) { #Is the plotting flag "up"?
  #Set X-axis to calculation-range and plot the data.
  X = (1:K)
  plot(X,Fractions, type="l", col="blue",xlab="Self-rank threshold", ylab="Fraction of pathway genes with self rank below threshold",main=paste(NameGE,"with",NameC))
}
#Calculate and return the AUSR score (i.e. AUC divided by threshold value).
return((AUC / K))}


### [2.6] MORPH {Input; MORPH object, k, output;Ranks, Score, AUSR, C, GE, }
MORPH = function(morph_obj, K=1000, view = FALSE)
{#Acquire data from input MORPH object
  G = morph_obj$pathway_genes #Pathway-genes.
  List_C = as.character(morph_obj$config[,2]) #Paths to clustering solution files.
  List_GE = as.character(morph_obj$config[,1]) #Paths to gene-expression data files.
  ge_datasets = morph_obj$ge_datasets #Gene-expression matrices.
  cl_solutions = morph_obj$clustering_solutions #Clustering solutions.
  #Initialize global variables to contain - 
  Gene_Scores = list() #Normalized gene scores.
  AUSR_Scores = numeric(0) #AUSR scores.
  Final_Scores = list() #MORPH Results object to contain all approriate information (see function description for more details).
  #Calculate the coexpression (covariance) matrices for all gene-expression matrices.
  corr_matrices = list() #Initialize list to contain all covariance matrices.
  ges = unique(List_GE) #Eliminate duplicates from paths to gene-expression data files.
  for (ge in ges){ #Go through all gene-expression matrices
    GE = ge_datasets[[ge]] #Extract gene-expression matrix.
    P = intersect(G,rownames(GE)) #Acquire relevant genes (those which appear both in the pathway-genes AND the gene-expression data).
    corr_matrices[[ge]] = cor(t(GE),t(GE[P,])) #Calculate and save coexpression (covariance) matrix into out list.
  }
  if (view == TRUE){ #Is the "view" parameter set to TRUE?
    pdf("Threshold_Plots.pdf", onefile=TRUE) #Open PDF formatted file to hold all self-ranking plots.
  }
  for (I in (1:length(List_C))){ #Go through all configurations (clustering solution vs. gene-expression dataset)
    #Acquire necessary information
    C = cl_solutions[[(List_C[I])]] #Clustering solution.
    curr_corrs = corr_matrices[[List_GE[I]]] #Coexpression (covariance) matrix.
    Pathway = removeAbscentGOIs(G,C,rownames(curr_corrs)) #Eliminate all pathway genes that do not appear in this configuration.
    
    Ranking = rankGenes(Pathway,C,GE=NULL,curr_corrs) #Rank candidates vs. pathway-genes.
    AUC = numeric(0) #Initialize variable to contain AUSR score.
    
    if (view == 0){ #Is the "view" parameter set to FALSE?
      AUC = LOOCV(Pathway,C,GE=NULL,K,curr_corrs) #Proceed to execute LOOCV as usual.
    } else { #Is the "view" parameter set to TRUE?
      AUC = LOOCV(Pathway,C,GE=NULL,K,curr_corrs, List_GE[I],List_C[I]) #Pass names of configuration files for plotting purposes.
    }
    #Add approriate information to the list (see function description for more details).
    Score = list(C=List_C[I],GE=List_GE[I],AUSR = AUC, Ranking = Ranking) #Arrange information into a single data block.
    Final_Scores[[I]] = Score #Add data block to the list.
  }
  if (view == TRUE){ #Is the "view" parameter set to TRUE?
    dev.off() #Turn off console presentation in order to write into output file and not to screen.
  }
  #Set final scores class as "morph_results", and return it.
  class(Final_Scores)<-"morph_results"
  return(Final_Scores)}

### [2.7] remove abscent GOIs {Input; G, C, Output; Pathway genes appeared in both}
removeAbscentGOIs = function(G,C,GENames)
{int = intersect(G,rownames(C)) #Keep only pathway-genes that appear in clustering solution.
int = intersect(int,GENames) #Keep only pathway-genes that appear in gene-expression data.
return (int)}

### [2.8] LOOCV_MORPH {Input; Morph object, k, output; AUSR score}
LOOCV_MORPH = function(morph_obj,K=1000)
{SelfRanks = numeric(0) #Initialize variable to contain all self-ranks.
for (V in morph_obj$pathway_genes) #Go through all pathway-genes and leave one out at a time.
{# Create a new MORPH input object but EXCLUDE the tested gene.
  new_input = list(config = morph_obj$config,
                   clustering_solutions = morph_obj$clustering_solutions,
                   ge_datasets = morph_obj$ge_datasets, 
                   pathway_genes = setdiff(morph_obj$pathway_genes,V))
  class(new_input)<-"morph_data"
  morph_run_results = MORPH(morph_obj=new_input,K) #Run the MORPH algorithm on the current setting.
  Scores = getMorphPredictions(morph_run_results) #Get the normalized correlation scores of current run.
  Rank = which(names(Scores)==V) #Acquire the location of the left-out gene among the others.
  if (length(Rank)==0 || is.nan(Rank)){Rank = (2*K)} #Was the left-out genes NOT detected at all? Set score to 2K as penalty.
  SelfRanks=c(SelfRanks,Rank) #Add self-rank to our vector.
}
#Calculate and return AUSR score for the MORPH process.
return(AUSR(SelfRanks,K))}

### [2.9] getMorphResultBestConfig {Input; MORPH-result object, Output; Ranks, c, GE, AUSR, Score}
getMorphResultBestConfig = function(morph_res_obj){ # list of lists [like PERL hush table]
  #Initialize variables to contain - 
  best_ausr = -1 #Best AUSR score.
  best_ind = -1 #Index of best-AUSR configuration.
  clustering_solution = "" #Name of clustering solution file in best configuration.
  gene_expression = "" #Name of gene-expression data file in best configuration.
  for (i in 1:length(morph_res_obj)){ #Go through all configurations in MORPH Results object.
    curr_res = morph_res_obj[[i]] #Extract current results data block.
    ausr = curr_res$AUSR #Acquire AUSR score.
    if (ausr>best_ausr){ #Is the AUSR score better than the last recorded best?
      best_ind = i; best_ausr = ausr #Save both the score and the configuration index.
      clustering_solution = morph_res_obj[[i]]$C #Acquire the clustering solution file name.
      gene_expression = morph_res_obj[[i]]$GE #Acquire the gene-expression data file name.
    }
  }
  # Create and return a list containing information regarding best configuration (for more information see function description).
  morph_best_res = list(AUSR = best_ausr, Ranking = morph_res_obj[[best_ind]]$Ranking, C = clustering_solution, GE = gene_expression)
  return (morph_best_res)}

### [2.10] get Morph Predictions {Input; Morph-result object, output;Ranking score}
getMorphPredictions = function(morph_res_obj){
  return (getMorphResultBestConfig(morph_res_obj)$Ranking$Scored) #Acquire and return scores of best configuration.}
  
  ### [2.11] get Scores Distribution Plots  {Input:Scores,path for saving,  output;plot}
  getScoresDistributionPlots = function(Scores, OutputFile="ScoresDistribution.pdf", Color="blue", Type = "b")
  {pdf(OutputFile, onefile=TRUE) #Open the output file for writing in PDF format.
    for (i in 1:length(Scores))	{ #Go through all MORPH Results object data blocks
      Ranks = Scores[[i]] #Acquire current data block.
      Y = as.numeric((Ranks$Ranking)$Scored) #Extract normalized correlation scores into a numeric vector.
      X = 1:(length(Y)) #Create an X-axis in the appropriate length.
      plot(X,Y,col=Color,type=Type, xlim=c(0,1), xlab="Self-rank threshold", ylab="Fraction of pathway genes with self rank below threshold",main=paste("Clustering:",Ranks$C,"; Gene expression:",Ranks$GE)) #Plot data.
    }
    dev.off() #Turn off console presentation in order to write into output file and not to screen.
  }}


# Solutions
#setwd("D:/UHasselt-Master of statistics and data science/Biostat/Second year/MORPH_OG/morph")
setwd("C:/Users/elise/Documents/Schuul/Stage/Project/MorphMTB_Script/data")
## MORPH Input
InputConfig = "Configs.txt"
InputGOI = "Pathway2.txt"
morph_input = prepareMorphObjectFromFiles(InputConfig,InputGOI)
Scores = MORPH(morph_input, view=TRUE)
## Removing Absent genes
G = morph_input$pathway_genes #Get pathway genes
C = (morph_input$clustering_solution)[[1]] #Get clustering solution
GE = (morph_input$ge_datasets)[[1]] #Get gene expression dataset
GENames = rownames(GE) #Get names of genes in dataset
Intersection = removeAbscentGOIs(G,C,GENames)
print(Intersection)
## Validation
InputConfig = "Configs.txt"
InputGOI = "Pathway2.txt"
morph_input = prepareMorphObjectFromFiles(InputConfig,InputGOI)
LOOCVc = LOOCV_MORPH(morph_input)
print(LOOCVc) #0.8283571
## AUSR
Scores2 = LOOCV_MORPH(morph_input) # --> gives 1 score?? 0.8283571
BestConfig <- getMorphResultBestConfig(Scores) #Error in curr_res$AUSR : $ operator is invalid for atomic vectors
print(names(BestConfig))
print(BestConfig$AUSR) #0.8283571
## MORPH gene scores
Predictions <- getMorphPredictions(Scores)
print(head(Predictions))
## Threshold Plot 
Scores <- MORPH(morph_input, view = TRUE)
getScoresDistributionPlots(Scores) #RStudioGD  --> gives empty plots?
#        2 

### Solutions for each dataset for a specific pathway
Scores[[1]]$AUSR   ## Numbers are correspondings to where the clustering solution is in Configs.txt
head(Scores[[1]]$Ranking$Scored, 5)

# Sampling for Real pathway recognition
uniform_sample <- function(vector) {
  index <- sample(length(vector), 14)
  return(vector[index])}
Score<-c()
B <-100
for (i in 1:B) {
  G <- uniform_sample(g) #Error: object 'g' not found
  writeLines(G, "random.txt", sep = "\t")
  InputConfig <- "Configs.txt"
  InputGOI <- "random.txt"
  morph_input <- prepareMorphObjectFromFiles(InputConfig, InputGOI)
  print(names(morph_input)) 
  G <- morph_input$pathway_genes
  NameC <- names(morph_input$clustering_solutions)[1] 
  C <- (morph_input$clustering_solutions)[[NameC]] #Get first clustering solution
  NameGE <- names(morph_input$ge_datasets)[1] 
  GE <- (morph_input$ge_datasets)[[NameGE]] #Get first gene expression dataset
  InputConfig <- "Configs.txt"
  InputGOI <- "random.txt"
  morph_input <- prepareMorphObjectFromFiles(InputConfig,InputGOI)
  Scores <- MORPH(morph_input, view = TRUE)
  Score[i]<-Scores}

# Assuming your list is named "my_list"
scores_AUSR <- vector("numeric", length(Score))
for (i in seq_along(Score)) {
  scores_AUSR[i] <- Score[[i]]$AUSR}
scores_AUSR <- sapply(Score, function(x) x$AUSR)
Pathways2<-scores_AUSR

## Boxplot for comparability
## Data preparation
dataset1<-data.frame(Pathways1, Pathways2, Pathways3, Pathways4, Pathways5) #Error: object 'Pathways1' not found
x<-c("Pathways1","Pathways2", "Pathways3", "Pathways4", "Pathways5")
y<-c(0.2448, 0.8481,0.6927, 0.5087,0.6713)
dataset2<-data.frame(x,y)

##### Boxplot
boxplot(dataset1, ylim=c(0,1))
points(x = match(dataset2$x, names(dataset1)), y = dataset2$y, col = "red", pch = 16)
text(x = match(dataset2$x, names(dataset1)), y = dataset2$y, labels = dataset2$y, pos = 3)

