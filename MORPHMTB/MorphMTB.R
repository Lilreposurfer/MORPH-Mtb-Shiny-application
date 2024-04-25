# Load packages ----
library(shiny)
library(shinyjs)

# Source ----
source("functions.R")

#Define the js method that resets the page
jsResetCode <- "shinyjs.reset = function() {history.go(0)}"

shinyApp(
  # Define UI for application
  shinyUI(
    #Make page with multiple panels
    navbarPage("MorphMTB", #title
               tabPanel("Gene centric query", uiOutput('page1')),
               tabPanel("Functional category query",  uiOutput('page2')),
               tabPanel("About", uiOutput('page3')),
   
    
   #   useShinyjs(),                                           # Include shinyjs in the UI
    #  extendShinyjs(text = jsResetCode, functions = "reset"), # Add the js code to the page
    #  actionButton("reset_button", "Reset Page")
    )
  ),
  # Define server logic
  shinyServer(function(input, output, session) {
    ## set seed for reproducibility
    global_seed = reactive(input$chosen_seed)
    set.seed(2023)
    
    
    #what happens on page1 (gene centric query)
    output$page1 <- renderUI({
      sidebarLayout(
        sidebarPanel(
          helpText("Here you can submit your genes of interest and get the MORPH prediction results in a table format. \n Species currently available in MorphMTB: Mycobacterium tuberculosis"),
          tags$hr(),
          numericInput("candidate", "Max of candidate genes to display:", 30, min = 1, max = 1000),
          tags$hr(),
          radioButtons("input","Select input:", 
                       choices=c("Genes","Upload file"),
                       selected="Genes"),
          textAreaInput("genes", "Please enter at least 10 genes as Input pathway (each on next row)."),
          fileInput("file","Upload file"), # fileinput() function is used to get the file upload control option
          tags$hr(),
          actionButton("button","Start"),
          actionButton("reset_inputs","Reset inputs")
        ),
        mainPanel(
          uiOutput("tb")
          
        )
      )
    })
    
    output$page3 <- renderUI({
            mainPanel(
              h4("Summary"),
              verbatimTextOutput("summary")
            )})
    
    output$contents <- eventReactive(input$button, {
      output$contents <- reactive({
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, it will be a data frame with 'name',
        # 'size', 'type', and 'datapath' columns. The 'datapath'
        # column will contain the local filenames where the data can
        # be found.
        if (input$input == "Upload file"){
          output$contents <- renderTable({
            inFile <- input$file
            gene <- unlist(read.csv(inFile$datapath, header=FALSE))
            return(data.frame(Genes=gene))})
        } else {
          output$contents <- renderTable({
            gene <- unlist(strsplit(input$genes, "\n"))
            return(data.frame(Genes=gene))})
        } 
      })
     
    })
    
  
    
    #output$contents <- eventReactive(input$restart, {
    #  return(NULL)
    #})

      
    output$tb <- renderUI({
      tabsetPanel(
        tabPanel("Pathway", tableOutput("contents")),
        tabPanel("Result", verbatimTextOutput("result"))
        )
    })
      

  })

)
