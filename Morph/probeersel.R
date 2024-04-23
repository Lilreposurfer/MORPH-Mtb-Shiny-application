
### [2.8] LOOCV_MORPH {Input; Morph object, k, output; AUSR score}
LOOCV_MORPH <- function(morph_input,K=1000, view=TRUE)
{SelfRanks <- numeric(0) #Initialize variable to contain all self-ranks.
for (V in morph_input$pathway_genes) #Go through all pathway-genes and leave one out at a time.
{ # Create a new MORPH input object but EXCLUDE the tested gene.
  new_input <- list(config = morph_input$config,
                    clustering_solutions = morph_input$clustering_solutions,
                    ge_datasets = morph_input$ge_datasets, 
                    pathway_genes = setdiff(morph_input$pathway_genes,V))
  class(new_input)<-"morph_data"
  morph_input<-new_input
  morph_run_results <- MORPH(morph_input,1000) #Run the MORPH algorithm on the current setting.
  #Scores <- getMorphPredictions(morph_run_results) #Get the normalized correlation scores of current run.
  Scores = getMorphResultBestConfig(morph_run_results)[[2]]$Scored
  Rank = which(names(Scores)==V) #Acquire the location of the left-out gene among the others.
  if (length(Rank)==0 || is.nan(Rank)){Rank = (2*1000)} #Was the left-out genes NOT detected at all? Set score to 2K as penalty.
  SelfRanks<-c(SelfRanks,Rank) #Add self-rank to our vector.
}
#Calculate and return AUSR score for the MORPH process.
return(AUSR(SelfRanks,1000))}


getMorphPredictions <- function(morph_run_results){
  return (getMorphResultBestConfig(morph_run_results)$Scored) #Acquire and return scores of best configuration.
}
bestconfigg <- getMorphResultBestConfig(morph_run_results)
bestconfigg[[2]]$Scored

LOOCV_MORPH(morph_input, view=TRUE)
Scores = AUSR(SelfRanks,K)
