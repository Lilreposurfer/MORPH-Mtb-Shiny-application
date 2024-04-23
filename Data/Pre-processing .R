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

# Seeds for reproducibility
set.seed(2023)

## Data Preparation ##
# Loading data
#setwd("Path to the files")
setwd("C:/Users/elise/Documents/Schuul/Stage/Project/Data")
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



























