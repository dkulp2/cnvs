#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Copy Number Variants"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
      sidebarPanel(
          tabsetPanel(type="tabs",
                      tabPanel("Selection",
                               tags$table(tags$tr(tags$td(selectInput("block","Block Runs",c(TEST1='/humgen/cnp04/bobh/projects/segmentation/pipeline_test_1',DKTest1='/home/unix/dkulp/data/pipeline_test_1'))),
                                                  tags$td(textOutput("block.descr"))), width='100%'),
                               tags$table(tags$tr(tags$td(numericInput("bin.start","Start Bin", NA)),
                                                  tags$td(numericInput("bin.end", "End Bin", NA)),
                                                  tags$td(actionButton("block.go", "Update")))),
                               tags$table(tags$tr(tags$td(selectizeInput("samples", "Samples", choices=list(), multiple=TRUE))),width='100%'),
                               tags$table(tags$tr(tags$td(checkboxInput("show_binMap_gc","Show G+C", value=FALSE)),
                                                  tags$td(checkboxInput("show_binData_ratio","Show Profiles (Obs/Exp)", value=FALSE)),
                                                  tags$td(checkboxInput("show_binData_coverage","Show Coverage", value=FALSE))),
                                          tags$tr(tags$td(checkboxInput("show_winGeno_cn","Show Windowed ProfileGenotyper Calls", value=FALSE)),
                                                  tags$td(checkboxInput("show_winGeno_cnq","Show Windowed ProfileGenotyper CNQ Calls", value=FALSE)),
                                                  tags$td(checkboxInput("show_winGeno_ratio","Show Windowed Profiles", value=FALSE))
                                                  )
                                          )
                               ),
                      tabPanel("UI Settings",
                               sliderInput("sample.bg.n", "Number of Background Samples:", min = 0, max = 400, value = 10, step=1),
                               sliderInput("padding.bins", "Number of Extra Flanking Bins Displayed:", min=0, max=10, value=4, step=1),
                               sliderInput("simplify.resolution", "Maximum Number of Bins for Full Display:", min=200, max=5000, value=500, step=100),
                               hr(),
                               checkboxInput("show.wings", "Display Bars Showing Window Range", value=FALSE),
                               checkboxInput("use.steps", "Display Profiles as Steps", value=TRUE),
                               checkboxInput("genomic.scale", "Display in Genomic Scale (NOT IMPLEMENTED)", value=FALSE)
                               )
                      ) # tabsetPanel
      ), # sidebarPanel
    
    # Show a plot of the generated distribution
      mainPanel(
          tabsetPanel(type="tabs",
                      tabPanel("Main Display",
                               conditionalPanel("input.show_binMap_gc", plotOutput("gc.plot")),
                               conditionalPanel("input.show_binData_ratio",plotOutput("profile.plot")),
                               conditionalPanel("input.show_binData_coverage",plotOutput("coverage.plot")),
                               conditionalPanel("input.show_winGeno_cn",plotOutput("win.genotype.plot")),
                               conditionalPanel("input.show_winGeno_cnq",plotOutput("win.CNQ.plot")),
                               conditionalPanel("input.show_winGeno_ratio",plotOutput("win.profile.plot"))
                               )
                      )
      ) # mainPanel
  ) # sidebarLayout
) # fluidPage
) # shinyUI
