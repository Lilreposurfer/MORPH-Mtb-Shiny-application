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

## Data Preparation for MORPH ##
# Loading Pathways
#setwd("D:/UHasselt-Master of statistics and data science/Biostat/Second year/MORPH_OG/Data")
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
#setwd("D:/UHasselt-Master of statistics and data science/Biostat/Second year/MORPH_OG/morph1")
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
g <- drug$Gene_ID


























