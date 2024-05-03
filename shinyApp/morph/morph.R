# Load packages ----
library(shiny)
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

# RV numbers!!

# Source ----
source("clustering.R")

# Define UI for application that draws a histogram
ui <- fluidPage(
  #Title of web application
  titlePanel("MORPH"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file","Upload the file", multiple = TRUE), # fileinput() function is used to get the file upload control option
      helpText("Default max. file size is 5MB"),
      helpText("Select the table parameters below"),
      checkboxInput(inputId = 'header', label = 'Header', value = TRUE), #checkbox to select if file has a header
      radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ''), #radiobuttons to specify separator
      numericInput("kmax", "Maximum K-clusters:", 10, min = 1, max = 100),
      numericInput("candidate", "Max of candidate genes to display:", 30, min = 1, max = 1000),
      uiOutput("selectfile")
    ),
    mainPanel(
      uiOutput("tb"),
      #fluidRow(
      #  splitLayout(cellWidths = c("50%","50%"), plotOutput("kmeans"), plotOutput("som"))
      #)
      #plotOutput("kmeans")
      
      
    )
    
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  ## set seed for reproducibility
  global_seed = reactive(input$chosen_seed)
  set.seed(2023)
  
  ## input$file is a data frame and contains the details around the name, size and temp location of the files uploaded
  # this reactive output display the content of the input$file dataframe
  output$filedf <- renderTable({
    if(is.null(input$file)){return ()}
    input$file # the file input data frame object that contains the file attributes
  })
  
  # Extract the file path for file
  output$filedf2 <- renderTable({
    if(is.null(input$file)){return ()}
    input$file$datapath # the file input data frame object that contains the file attributes
  })
  
  ## Below code to display the structure of the input file object
  output$fileob <- renderPrint({
    if(is.null(input$file)){return ()}
    str(input$file)
  })
  
  ## Side bar select input widget coming through renderUI()
  # Following code displays the select input widget with the list of file loaded by the user
  output$selectfile <- renderUI({
    if(is.null(input$file)) {return()}
    list(hr(), 
         helpText("Select the files for which you need to see data and summary stats"),
         selectInput("Select", "Select", choices=input$file$name)
    )
    
  })
  
  ## Summary Stats code ##
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$summ <- renderPrint({
    if(is.null(input$file)){return()}
    summary(read.table(file=input$file$datapath[input$file$name==input$Select], 
                       sep=input$sep, 
                       header = input$header))})
  
  ## Dataset code ##
  # This reactive output contains the dataset and display the dataset in table format
  output$table <- renderTable({ 
    if(is.null(input$file)){return()}
    read.table(file=input$file$datapath[input$file$name==input$Select], 
               sep=input$sep, 
               header = input$header)},
    striped=TRUE,
    caption="Uploaded table")
  
  ## Clustering code ##
  # This reactive output contains the K-means and the SOM plot of clustering next to each other
  file <- reactive({
    read.table(file=input$file$datapath[input$file$name==input$Select], 
               sep=input$sep, 
               header = input$header)
  })
  logs <- reactive({
    log(file())
  })
  kmax <- reactive({
    input$kmax
  })
  output$wssk <- reactive({
    wsskmeans(logs(),kmax())
  })
  output$wsss <- reactive({
    wsssom(logs(), kmax())
  })
  output$x <- reactive({
    2:kmax
  })
  
  output$kmeans = renderPlot({
        plot(kmax(), wssk(), 
                   type="b", pch =10, frame = FALSE, 
                   xlab="Number of clusters",
                   ylab="Total Within sum of square",
                   main="K-means",
                   cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
       #abline(v = 3, col = "red", lwd = 2)
    })
  #output$som = renderPlot({
        #plot(kmax, wsss, type="b", pch =10, frame = FALSE, 
     #            xlab="Number of clusters",
     #            ylab="Total Within sum of square",
     #            main="SOM",
     #            cex=2, cex.main=1.5, cex.lab=1.5, cex.axis=2)
     #abline(v = 3, col = "red", lwd = 2)
    #})
  
  
  ## MainPanel tabset renderUI code ##
  # the following renderUI is used to dynamically generate the tabsets when the file is loaded. 
  # Until the file is loaded, app will not show the tabset.
  output$tb <- renderUI({
    if(is.null(input$file)) {return()}
    else
      tabsetPanel(
        tabPanel("Dataset", tableOutput("table")),
        tabPanel("Summary", verbatimTextOutput("summ")),
        tabPanel("Input File Object DF ", tableOutput("filedf"), tableOutput("filedf2")),
        tabPanel("Structure", verbatimTextOutput("fileob")),
        #tabPanel("Pre-processing", textOutput("wssk"), textOutput("wsss")),
        tabPanel("Clustering", plotOutput("kmeans"))
      )
  })
}


# Run the application 
shinyApp(ui = ui, server = server)


###############################################################################################
output$contents <- eventReactive(input$button, {
  output$contents <- reactive({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    #if (input$input == "Upload file"){
    ## Side bar select input widget coming through renderUI()
    # Following code displays the select input widget with the list of file loaded by the user
    output$selectfile <- renderUI({
      if(is.null(input$file)) {return()}
      list(hr(), 
           helpText("Select the files for which you need to see data and summary stats"),
           selectInput("Select", "Select", choices=input$file$name))
    })
    output$contents <- renderTable({
      #inFile <- input$file
      #gene <- unlist(read.csv(inFile$datapath, header=FALSE))
      #return(data.frame(Genes=gene))})
      read.table(file=input$file$datapath[input$file$name==input$Select], 
                 sep='', 
                 header = TRUE)},
      striped=TRUE)
    
    
    #} else {
    output$pathway <- renderTable({
      gene <- unlist(strsplit(input$genes, "\n"))
      return(data.frame(Genes=gene))},
      striped=TRUE)
    #} 
    
    
    ## Clustering code ##
    # This reactive output contains the K-means and the SOM plot of clustering next to each other
    #file <- reactive({
    #  read.table(file=input$file$datapath[input$file$name==input$Select], 
    #             sep='', 
    #             header = TRUE)
    #})
    #logs <- reactive({
    #  log(file())
    #})
    #kmax <- reactive({
    #  kmax <- 10
    #})
    #wssk <- reactive({
    #  wsskmeans(logs(),kmax())
    #})
    #wsss <- reactive({
    #  wsssom(logs(), kmax())
    #})
    #x <- reactive({
    #  2:kmax
    #})
    #output$Kmeans <- reactive({
    #  renderPlot({
    #    clusterK(x, wssk())
    #  })
    #})
    #output$SOM <- reactive({
    #  renderPlot({
    #    clusterS(x, wsss())
    #  })
    #})
    
    
    
    
    #ScoreRandom <- c()
    #B <- 100  
    #for (i in 1:B) {
    #  InputGOIRandom <- reactive({
    #    uniform_sample(g())
    #  })
    #  morphinputRandom <- reactive({
    #    prepareMorphObjectFromFiles(InputGOIRandom())
    #  })
    #  #output$morphinputNames <- reactive({
    #  #  names(morphinputRandom())
    #  #})
    #  output$GRandom <- reactive({
    #    morphinputRandom()$pathway_genes
    #  }) 
    #  NameC <- reactive({
    #    names(morphinputRandom()$clustering_solutions)[1]
    #  })
    #  output$C <- reactive({
    #    (morphinputRandom()$clustering_solutions)[[NameC()]] #Get first clustering solution
    #  })
    #  NameGE <- reactive({
    #    names(morphinputRandom()$ge_datasets)[1]
    #  })
    #  output$GE <- reactive({
    #    (morphinputRandom()$ge_datasets)[[NameGE()]] #Get first gene expression dataset
    #  })
    #  ScoresRandom <- reactive({
    #    MORPH(morphinputRandom(), view = TRUE)
    #  })
    #  ScoreRandom()[i] <- reactive({
    #    ScoresRandom()
    #  })
    #}
