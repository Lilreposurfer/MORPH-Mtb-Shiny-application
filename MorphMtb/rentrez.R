# Character vector of installed packages
pkg <- installed.packages()[, "Package"]
# Check if rentrez package is installed
'rentrez' %in% pkg
# If package is not installed --> install!
if(!('rentrez' %in% pkg)) {install.packages("rentrez")}
library("rentrez")


Organism <- function(id){
  id = id
  geneOrg = c()
  for (i in id){
    res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN])"))
    res$count
    res$ids
    esums <- entrez_summary(db="gene", id=res$ids)
    geneOrg <- c(geneOrg, esums$organism$scientificname)
  }
  return(geneOrg)
}

GeneID <- function(names){
  id = id
  geneIDs = c()
  for (i in id){
    res <- entrez_search(db="gene", term=paste("(", id, "[GENE] AND Mycobacterium tuberculosis H37Rv [ORGN])"))
    res$count
    res$ids
    esums <- entrez_summary(db="gene", id=res$ids)
    geneIDs <- c(geneIDs, esums$otheraliases)
}
  return(geneIDs)
}

GeneName <- function(id){
  id = id
  geneName = c()
    for (i in id){
      res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN])"))
      res$count
      res$ids
      esums <- entrez_summary(db="gene", id=res$ids)
      geneName <- c(geneName, esums$name)
    }
  return(geneName)
}


Description <- function(id){
  id = id
  desc = c()
    for (i in id){
      res <- entrez_search(db="gene", term=paste("(", i, "[GENE] AND Mycobacterium tuberculosis[ORGN])"))
      res$count
      res$ids
      esums <- entrez_summary(db="gene", id=res$ids)
      extracteddesc <- extract_from_esummary(esums, "description")
      desc <- c(desc, extracteddesc)
    }
  return(desc)
}

