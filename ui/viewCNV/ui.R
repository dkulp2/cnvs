#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

textareaInput <- function(id, label, value, rows=20, cols=35, class="form-control"){
  tags$div(
    class="form-group shiny-input-container",
    tags$label('for'=id,label),
    tags$textarea(id=id,class=class,rows=rows,cols=cols,value))
}

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  singleton(
    tags$head(tags$script(src = "message-handler.js"))
  ),
  
  # Application title
  titlePanel("Copy Number Variants"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      tags$table(tags$td(selectInput("saved.site","Saved Sites",list(),width='350px')),
                 tags$td(actionButton("saved.prev","Prev")),
                 tags$td(actionButton("saved.next","Next"))),
      actionButton("reset","Reset"),actionButton("delete","Delete"),
      actionButton("save.button","Save Site"),
      tags$table(tags$td(selectInput("candidate","Known Deletions", list(), width='350px')),
                 tags$td(actionButton("candidate.prev","Prev")),
                 tags$td(actionButton("candidate.next","Next"))),
      tags$table(tags$td(selectInput("predicted","Predicted Deletions", list(), width='350px')),
                 tags$td(actionButton("predicted.prev","Prev")),
                 tags$td(actionButton("predicted.next","Next"))),
      tags$table(tags$tr(tags$td(selectInput("pred.order","Order", c('Chromosome'='chrom','Largest First'='desc','Smallest First'='asc'))),
                         tags$td(selectInput("pred.cn","Predict CN", c('Het Del'=1, 'Hom Del'=0, 'Dels'=2, 'Non-WT'=3,'Dups'=4)))),width='100%'),
      hr(),
      tags$table(tags$td(textInput("seg.chr","Chromosome",NA)),
                 tags$td(numericInput("seg.start","Start Position", NA)),
                 tags$td(numericInput("seg.end", "End Position", NA)),
                 tags$td(textOutput("seg.size")),
                 tags$td(actionButton("zoom.out", "Zoom Out"))),
      selectizeInput("seg.sample", "Samples", choices=list(), multiple=TRUE),
      textareaInput("note","Note (first line is title)", "", rows = 3),
      h3("Show"),
      selectInput("truth_data", "Truth Data Set", choices=list('GStrip Sn DEL Data', 'CNV Pipeline Sp Data'), selected='Gstrip Sn Data'),
      checkboxInput("show_CN", "Color Copy Number", value=TRUE),
      tags$table(tags$tr(tags$td(checkboxInput("show_cnv","CNV Extents", value=TRUE)),
                 tags$td(checkboxInput("show_winGeno","Windowed CN Genotypes", value=TRUE)),
                 tags$td(checkboxInput("show_frag","Fragment Profiles", value=TRUE)),
                 tags$td(checkboxInput("show_readPairs","Read Pair Counts", value=FALSE)),
                 tags$td(checkboxInput("show_expected","Expected Depth", value=FALSE)))),
      conditionalPanel("input.show_winGeno",
                       checkboxInput("show_winGeno_2","Show Jitter, Transparency, Window", value=TRUE),
                       checkboxInput("show_hires","Show Smaller Windowed Genotype Calls", value=FALSE),
                       sliderInput("winGeno_nudge_L","Nudge Windowed CN Genotypes (inches)", min=0, max=3, value=0.5, step=0.25)),
      conditionalPanel("input.show_frag",
                       sliderInput("win.size", "Number of Bins in Profile Window:", min = 1, max = 30, value = 12, step=1)),
      numericInput("pad","Left/Right Padding (bases)", min=0, max=500000, value=5000, step=1000),
      conditionalPanel("input.show_cnv",
                       checkboxInput("show_target_only","Show CNVs of Selected Samples Only", value=FALSE),
                       checkboxInput("show_wildtype","Show WildType (CN=2) Extents", value=FALSE),
                       checkboxInput("show_extended_cnvs","Show Extended Bounds", value=FALSE),
                       checkboxInput("show_extended_ML", "Show Maximum Likelihood Bounds", value=FALSE),
                       sliderInput("min.cnv.len", "Minimum Length of Displayed CNV Extents", min=0, value=100, max=5000, step=100)),
      sliderInput("paired.reads","Minimum Number of Paired Reads", min=0, max=10, value=0, step=1)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type="tabs",
                  tabPanel("Main Display",
                           plotOutput("allThree", height="750px"),
                           conditionalPanel("input.show_expected", 
                                            plotOutput("expPlot",height="200px"))),
                  tabPanel("Per Sample CNQs",
                           plotOutput("cnqs", height="800px")))
    )
  )
))
