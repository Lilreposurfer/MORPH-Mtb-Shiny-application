
# MORPH-Mtb Shiny web application

This Shiny web application was made based on the MORPH algorithm, originally designed for gene discovery in plants. MORPH-Mtb aims to detect and rank peviously unidentified genes linked to a particular pathway in *Mycobacterium tuberculosis*.  
  

In this repository, two folders where attached named Data and MorphMtb.

## Data folder
This data file constist out of all files that are needed as input to the MORPH-Mtb algorithm, including 6 pre-processed data sets used in the algorithm:
- ESX.txt
- clark.txt
- drug.txt
- inaki.txt
- primary.txt
- timecourse.txt
Together with their clustering solution files from K-means and SOM.

## MorphMtb folder
In this folder, the 4 scripts needed for the Shiny web application are stored:  
1. **MorphMtb.R**  
    This R script contains the code for the Shiny app itself. From here, you can start the web application.
2. **functions.R**  
    This R script has stored all of the functions used for the calculations of for example the AUSR-score.
3. **rentrez.R**  
    In this R script, a connection with NCBI is made to retrieve information about the candidate genes.
4. **styles.css**  
    In this css script, the GUI of the Shiny app is changed to make it more pretty.

## Using the Shiny web application