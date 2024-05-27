
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

To make use of the Shiny web application, R and RStudio need to be installed first. 

### Installing R on Windows

Go to the website using [this link](https://cran.r-project.org/bin/windows/base/old/4.3.2/) to install version 4.3.2 for Windows. 

Click on the "R-4.3.2-win.exe" link. A download should be started.  
1. Open the downloaded file
2. Choose English as installation language
3. Start the setup
4. Do not customize the startup options
5. Finish setup

### Installing RStudio on Windows
Go to the website using [this link](https://dailies.rstudio.com/version/2023.12.1+402/) to install version 2023.12.1+402.  

Click on the "RStudio-2023.12.1-402.exe" for Windows. A download should be started.  
1. Open the downloaded file
2. Click on install
3. Finish setup  
4. You can now open the RStudio application

### The MORPH-Mtb Shiny web application
To make use of the MORPH-Mtb Shiny web application, you first have to open the MorphMTB.R script with RStudio (Right-click on the file and choose open with RStudio).   
You will see a button with a green arrow and "Run App". Click on this button and the Shiny app will be started. You can choose to optionally open this in a web browser by clicking the "Open in Browser" button. 

The web application consists out of 3 pages:
1. **Gene centric query**  
    You can use this page if you only want to use the implemented 6 datasets.  
2. **Implement own expression** data    
    You can use this page if you want to upload your own expression data and use this together with the 6 datasets for the analysis.  
3. **About**  
    This is an informative page about the MORPH-Mtb algorithm. 