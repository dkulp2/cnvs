#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

source("../../conf/config.R")
read.conf("env.txt")  # pre-computed environment vars. Rerun setup.R to regenerate.

if.knowns <- function(x) {
  if (check.conf("USE_KNOWNS")) { x } else { tags$br() }
}

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
      actionButton("save.button","Save Site"),actionButton("dump","Save Debug State"),
      if.knowns(tags$table(tags$td(selectInput("candidate","Known Deletions", list(), width='350px')),
                 tags$td(actionButton("candidate.prev","Prev")),
                 tags$td(actionButton("candidate.next","Next")))),
      tags$table(tags$td(selectInput("predicted","Predicted Deletions", list(), width='350px')),
                 tags$td(actionButton("predicted.prev","Prev")),
                 tags$td(actionButton("predicted.next","Next"))),
      sliderInput("min.cnv.len", "Minimum Prediction Length", min=0, value=8, max=50, step=1),
      tags$table(tags$tr(tags$td(selectInput("pred.order","Order", c('Chromosome'='chrom','Largest First'='desc','Smallest First'='asc'))),
                         tags$td(selectInput("pred.cn","Predict CN", c('Het Del'=1, 'Hom Del'=0, 'Dels'=2, 'Non-WT'=3,'Dups'=4),selected=4))),width='100%'),
      hr(),
      tags$table(tags$tr(tags$td(textInput("seg.chr","Chromosome",NA)),
                         tags$td(numericInput("seg.start","Start Position", NA)),
                         tags$td(numericInput("seg.end", "End Position", NA)),
                         tags$td(textOutput("seg.size")),
                         tags$td(actionButton("zoom.out", "Zoom Out"))),
                 tags$tr(tags$td(numericInput("bin.start","Start Bin", NA)),
                         tags$td(numericInput("bin.end","End Bin", NA)),
                         tags$td(actionButton("bin.go", "Go")))),
      tags$table(tags$tr(tags$td(selectizeInput("seg.sample", "Samples", choices=list(), multiple=TRUE)),
                         tags$td(actionButton("add_quartet","Add Family Quartets"))),width="100%"),
      textareaInput("note","Note (first line is title)", "", rows = 3),
      if.knowns(selectInput("truth_data", "Truth Data Set", choices=list('GStrip Sn DEL Data', 'GStrip Sn Best DEL Data', 'CNV Pipeline Sp Data'), selected='Gstrip Sn DEL Data')), 
      tags$table(tags$tr(tags$td(checkboxInput("show_cnv","CNV Extents", value=FALSE)),
                         tags$td(checkboxInput("show_CI","Conf Intervals", value=TRUE)),
                         tags$td(checkboxInput("show_frag","Fragment Profiles", value=TRUE)),
#                         tags$td(checkboxInput("show_readPairs","Read Pair Counts", value=FALSE)),
                         tags$td(checkboxInput("show_expected","Expected Depth", value=FALSE))),
                 tags$tr(#tags$td(checkboxInput("show_prior","Loss/Gain Pseudo-Prior", value=FALSE)),
                         tags$td(checkboxInput("show_posterior","Likelihood & Posterior", value=FALSE)),
                         tags$td(checkboxInput("show_bayes_prior","Mean Prior", value=FALSE)),
                         tags$td(checkboxInput("show_each_bkpt","Loss/Gain By Sample", value=FALSE)))),
      conditionalPanel("input.show_frag",
                       selectInput("color_frag", "Color Fragments", choices=list('Selected','By Genotype', 'By Sample','By Family'), selected='Selected'),
                       sliderInput("win.size", "Number of Bins in Profile Window:", min = 1, max = 30, value = 12, step=1),
                       sliderInput("sample.bg.n", "Number of Background Samples:", min = 0, max = 400, value = 10, step=1)),
      conditionalPanel("input.show_prior || input.show_posterior || input.show_bayes_prior",
                       selectInput("bkpt_prior_samples", "Samples Prior", choices=list('All','Selected','Exclude Selected'), selected='All'),
                       selectInput("bkpt_prior_CN", "CN Prior", choices=list('Any','0 or 1','0','1'), selected='0 or 1')),
      numericInput("pad","Left/Right Padding (bases)", min=0, max=500000, value=5000, step=1000),
      conditionalPanel("input.show_cnv | input.show_CI",
                       checkboxInput("show_target_only","Show CNVs of Selected Samples Only", value=TRUE),
                       checkboxInput("show_extended_ML", "Show Maximum Likelihood Bounds", value=FALSE),
                       checkboxInput("show_bayes", "Show Bayes Adjusted Bounds", value=TRUE),
		       conditionalPanel("input.show_cnv", checkboxInput("show_wildtype","Show WildType (CN=2) Extents", value=FALSE),
                                                          checkboxInput("show_basic","Show Basic Extents", value=TRUE))
                       ),
      checkboxInput("show_table","Show Profile Data as Table", value=FALSE)
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type="tabs",
                  tabPanel("Main Display",
                           plotOutput("allThree", height="1000px"), #750px"),
                           plotOutput("geno.frag"),
                           conditionalPanel("input.show_expected", 
                                            plotOutput("expPlot",height="200px"))),
                  tabPanel("Per Sample CNQs",
                           fluidRow(class='CNQ',plotOutput("cnqs", height="800px")),
                           conditionalPanel("input.show_table",
                                            fluidRow(class='ProfileTable',
                                                     DT::dataTableOutput("profile.table")))),
                  tabPanel("Per Sample Breakpoints",
                           plotOutput("bkpts", height="800px"),
                           radioButtons("show_probs","Display value:", c("Prob"="probs", "Log"="log",'Normalized'="Z","No Loss/Gain"="no"), 'log')),
                  tabPanel("Per Sample Odds",
                           plotOutput("bkpt.odds", height="800px"),
                           radioButtons("show_bayes_odds","Display value:", c("Log Odds of Bayesian Estimates"="bayes", "Log Odds of Normalized Bayesian Estimates"="bayesZ","Log Odds of Likelihood Estimates"="likelihood")))
      )
    )
  )
))
