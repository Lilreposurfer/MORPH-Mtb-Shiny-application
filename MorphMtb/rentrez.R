# Character vector of installed packages
pkg <- installed.packages()[, "Package"]
# Check if rentrez package is installed
'rentrez' %in% pkg
# If package is not installed --> install!
if(!('rentrez' %in% pkg)) {install.packages("rentrez")}
library("rentrez")

#https://www.rdocumentation.org/packages/rentrez/versions/1.2.3/topics/entrez_summary
#https://www.rdocumentation.org/packages/rentrez/versions/1.2.3/topics/extract_from_esummary
#https://www.rdocumentation.org/packages/rentrez/versions/1.2.3/topics/entrez_search

#https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
#https://statsandr.com/blog/web-scraping-in-r/
#https://bioconnector.github.io/workshops/r-ncbi.html
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://journal.r-project.org/archive/2017/RJ-2017-058/RJ-2017-058.pdf


#entrez_dbs()
#entrez_db_searchable(db="gene")
#ALL 	 All terms from all searchable fields 
#UID 	 Unique number assigned to a gene record 
#FILT 	 Limits the records 
#TITL 	 gene or protein name 
#WORD 	 Free text associated with record 
#ORGN 	 scientific and common names of organism 
#MDAT 	 The last date on which the record was updated 
#CHR 	 Chromosome number or numbers; also 'mitochondrial', 'unknown' properties 
#MV 	 Chromosomal map location as displayed in MapViewer 
#GENE 	 Symbol or symbols of the gene 
#ECNO 	 EC number for enzyme or CAS registry number 
#MIM 	 MIM number from OMIM 
#DIS 	 Name(s) of diseases associated with this gene. When available, OMIM name will be used 
#ACCN 	 Nucleotide or protein accession(s) associated with this gene 
#UGEN 	 UniGene cluster number for this gene 
#PROP 	 Properties of Gene record 
#CDAT 	 The date on which this record first appeared 
#NCAC 	 nucleotide accessions of sequences 
#NUID 	 nucleotide uids of sequences 
#PACC 	 protein accessions 
#PUID 	 protein uids 
#PMID 	 PubMed ids of accessions linked to the record 
#TID 	 taxonomy id 
#GO 	 Gene Ontology 
#DOM 	 Domain Name 
#DDAT 	 The date on which the record was discontinued 
#CPOS 	 Chromosome base position 
#GFN 	 Gene full name 
#PFN 	 Protein full name 
#GL 	 Gene length 
#XC 	 Exon count 
#GRP 	 Relationships for this gene 
#PREF 	 Preferred symbol of the gene 
#AACC 	 Assembly accession 
#ASM 	 Assembly name 
#EXPR 	 Gene expression 

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

