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

# Seeds for reproducibility
#set.seed(2023)

# Define UI for application that draws a histogram
ui <- fluidPage(

  titlePanel("MORPH"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file","Upload the file", multiple = TRUE), # fileinput() function is used to get the file upload control option
      helpText("Default max. file size is 5MB"),
      helpText("Select the read.table parameters below"),
      checkboxInput(inputId = 'header', label = 'Header', value = TRUE),
      radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = '\t'),
      uiOutput("selectfile")
    ),
    mainPanel(
      uiOutput("tb")
      
    )
    
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
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
                       header = input$header, 
                       stringsAsFactors = input$stringAsFactors))})
  
  ## Dataset code ##
  # This reactive output contains the dataset and display the dataset in table format
  output$table <- renderTable({ 
    if(is.null(input$file)){return()}
    read.table(file=input$file$datapath[input$file$name==input$Select], 
               sep=input$sep, 
               header = input$header)},
    striped=TRUE,
    caption="Uploaded table:")
  
  
  ## MainPanel tabset renderUI code ##
  # the following renderUI is used to dynamically generate the tabsets when the file is loaded. 
  # Until the file is loaded, app will not show the tabset.
  output$tb <- renderUI({
    if(is.null(input$file)) {return()}
    else
      tabsetPanel(
        tabPanel("Input File Object DF ", tableOutput("filedf"), tableOutput("filedf2")),
        tabPanel("Input File Object Structure", verbatimTextOutput("fileob")),
        tabPanel("Dataset", tableOutput("table")),
        tabPanel("Summary Stats", verbatimTextOutput("summ")))
  })
}


# Run the application 
shinyApp(ui = ui, server = server)
