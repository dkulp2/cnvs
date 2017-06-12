# The shiny server log files don't seem to be reliably writing to 
# /home/unix/dkulp/bobh/shiny/var/log/shiny-server
# so I just redirect all output to my own file!
zz <- file("/tmp/cnv_seg.shiny.log", open = "at")
sink(zz)
sink(zz, type = "message")

library(shiny)
library(plyr)
library(dplyr)
library(ggplot2)
library(stats)
library(foreach)
library(zoo)
library(gridExtra)
library(tidyr)
library(reshape2)
library(metaTbl)

# uncomment to debug within RStudio
#reactive <- function(expr, env=parent.frame()) exprToFunction(expr, env=env)

# write either warning or message to console with a time stamp
log <- function(..., warn=FALSE) {
  level.func <- ifelse(warn, warning, message)
  level.func(Sys.time(),": ",...)
}

selected.colors <- c(selected='#FF0000',other='#222222')
selected.alpha <- c(selected=1,other=0.2)

# common colors for copy number
cn.colors <- scale_color_manual(values=c('0'="#7fc97f", '1'="#c51b7d", '2'="#fdc086", '3'="#ffff99", '4'="#386cb0", '5+'="#f0027f",'Disc'="#999999",'NA'='#333333'), 
                                name="CN")

# max out CN at 5. Set fixed levels and labels for fixed legend
cn.label <- function(cn) {
  factor(ifelse(cn==99, cn, ifelse(cn > 5, 5, cn)), levels=c('0','1','2','3','4','5','99'), labels=c('0','1','2','3','4','5+','Disc'))
}

# returns a nt position corresponding to bin.
# pos='start' : returns left side of bin
# pos='end'   : returns right side of bin
# pos='center': returns midpoint
bin2nt <- function(bin.map, bin, pos='start') {
  r <- filter(bin.map, binIndex==bin)
  if (pos=='start') {
    return(r[['startpos']])
  } else if (pos=='end') {
    return(r[['endpos']])
  } else if (pos=='center') {
    return(r[['startpos']] + (r[['endpos']]+1-r[['startpos']]) %/% 2)
  }
  stop("Bad pos: ",pos,". Must be start, end or center")
}

# if block is a path (starts with a '/'), then read the file.
# otherwise, read the table from the schema named block
readTblDb <- function(block, table) {
  # FIXME: security. UI can change file path
  fn <- sprintf("%s/%s.tbl.gz",block, table)
  log("Reading ",fn)
  tbl <- suppressWarnings(readTbl(fn))
  log("Read ",nrow(tbl)," rows")
  return(tbl)
}

updateLogNumericInput <- function(session, label, value) {
  log("Setting ",label," to ",value)
  updateNumericInput(session, label, value=value)
}

shinyServer(function(input, output, session) {

  ################################################################################
  # Reactive tables

  # region block{TEXT}, refId{TEXT}, chrom{TEXT}, startpos{INT}, endpos{INT};
  regions <- reactive({
    readTblDb(input$block,'regions')
  })
  
  # mapping of bin to nt, effective length, G+C
  # region block{TEXT}, binIndex{INT}, chrom{TEXT}, startpos{INT}, endpos{INT}, effLength{FLOAT}, gcnum{INT}, gcdenom{INT};
  bin.map <- reactive({
    bm <- readTblDb(input$block,'binMap')
    return(bm)
  })

  # data from original bins (observed, expected)
  # region block{TEXT}, binIndex{INT}, sample{TEXT}, oberved{FLOAT}, expected{FLOAT}
  bin.data <- reactive({
    readTblDb(input$block,'binData')
  })

  # winMap: block{TEXT}, winIndex{INT}, startBin{INT}, endBin{INT}, targetBin{INT}; [block, winIndex]
  win.map <- reactive({
    readTblDb(input$block,'winMap')
  })
  
  # region block{TEXT}, winIndex{INT}, sample{TEXT}, oberved{FLOAT}, expected{FLOAT}, cn{INT}, cnf{FLOAT}, cnq{FLOAT}; [block] (block, winIndex, sample)
  win.geno <- reactive({
    readTblDb(input$block,'winGeno')
  })


  ################################################################################
  # Reactive variables
  
  # all available samples in this block.
  # choose by arbitrarily querying the first bin, assuming that data is present
  # for all samples at all positions.
  all.samples <- reactive({
    unique(filter(bin.data(), binIndex==bin.min())[['sample']])
  })

  # returns the chromosome for the current block
  block.chr <- reactive({
    regions()$chrom
  })

  # bin bounds for this block
  bin.min <- reactive({
    min(bin.map()$binIndex)
  })
  bin.max <- reactive({
    max(bin.map()$binIndex)
  })

  ################################################################################
  # Initialize client inputs
  
  # return an initial coordinate before user chooses one
  # currently returns 0..500
  init.start <- reactive({
    updateLogNumericInput(session, 'bin.start', value=bin.min())
    bin.min()
  })
  init.end <- reactive({
    ie <- min(init.start()+500, bin.max())
    updateLogNumericInput(session, 'bin.end', value=ie)
    ie
  })

  observe({
    updateSelectizeInput(session, 'samples', choices=all.samples())
  })

  ################################################################################
  # Selected Position

  # use bin.start() and bin.end() instead of input$bin.start and input$bin.end
  # if you want to be updated only when the user presses the update button. otherwise
  # you run the risk of executing an expensive redraw as the user is typing
  bin.start <- reactive({
    input$block.go
    b.s <- isolate(input$bin.start)
    if (is.na(b.s)) { b.s <- init.start() }
    log("bin.start: ", b.s)
    b.s
  })
  bin.end <- reactive({
    input$block.go
    b.e <- isolate(input$bin.end)
    if (is.na(b.e)) { b.e <- init.end() }
    log("bin.end: ", b.e)
    b.e
  })
  
  nt.start <- reactive({
    bin2nt(bin.map(),bin.start(),'start')
  })

  nt.end <- reactive({
    bin2nt(bin.map(),bin.end(),'end')
  })

  region.descr <- reactive({
    sprintf("Chr %s:B%s-B%s [%s-%s] (%s bins, %skb nt)", block.chr(), bin.start(), bin.end(),
            nt.start(), nt.end(),
            bin.end()-bin.start()+1, format( (nt.end() + 1 - nt.start())/1000, digits=4, nsmall=2))
  })

  ################################################################################
  # Reactive variables and convenience functions with reactive elements

  # filters the df for only rows within the selected region including padding.
  # The df must have a column 'binIndex'.
  regionFilter <- function(df) {
    filter(df, binIndex >= max(bin.start()-input$padding.bins, bin.min()) & binIndex <= min(bin.end()+input$padding.bins, bin.max()))
  }

  # returns a subset of win.map that is within the selected region including padding
  # that can be joined with tables that are indexed by winIndex.
  regionWinMapFilter <- reactive({
    filter(win.map(), startBin >= max(bin.start()-input$padding.bins, bin.min()) & endBin <= min(bin.end()+input$padding.bins, bin.max()))
  })
  
  # filter the df so that all selected samples are included and a subset of input$sample.bh.n additional samples
  # are also included. The df must have a column 'sample'.
  sampleFilter <- function(df) {
    log("sampleFilter: sample.bg.n=",input$sample.bg.n)
    filter(df, sample %in% c('', input$samples, sample(all.samples(), input$sample.bg.n)))
  }

  # returns a collected tibble of win.geno() in the selected region
  genoWinMap <- reactive({
    df <- sampleFilter(inner_join(win.geno(), regionWinMapFilter(), by=c('block','winIndex'))) %>% collect
    df$sample.type <- ifelse(df$sample %in% c('',input$samples),'selected','other')

    log("genoWinMap: ",nrow(df)," rows. ",nrow(regionWinMapFilter())," win.map rows.")
    df
  })

  ################################################################################
  # Common Plot Attributes

  # For plots of win.geno, plot dots, wings and/or steps depending on UI params
  addGeomFeatures <- function(p) {
    if (input$use.steps) {
      p <- p + geom_step(aes(x=targetBin))
    } else {
      p <- p + geom_point(aes(x=targetBin))
    }
    if (input$show.wings) {
      p <- p + geom_segment(aes(x=startBin, xend=endBin, yend=obsexp))
    }
    p
  }

  # when plotting samples, highlight selected and show rest lightly
  selected.theme <- function(p) {
    p + theme(axis.title.x = element_blank()) + 
      scale_color_manual(name = "sample.type",values = selected.colors) +
      scale_alpha_manual(name = "sample.type",values = selected.alpha) + guides(color=FALSE,alpha=FALSE)
  }


  ################################################################################
  # Outputs

  output$block.descr <- renderText({
    nt.min <- bin2nt(bin.map(), bin.min(),'start')
    nt.max <- bin2nt(bin.map(), bin.max(), 'end')
    sprintf("Chr %s:B%s-B%s [%s-%s] (%skb nt)", block.chr(), bin.min(), bin.max(),
            nt.min, nt.max,
            format( (nt.max+1-nt.min)/1000, digits=4, nsmall=2))
  })

  output$region.descr <- renderText({
    region.descr()
  })

  # plot G+C per bin
  output$gc.plot <- renderPlot({
    log("gc.plot: ",bin.start(),'-',bin.end())
    df <- regionFilter(bin.map()) %>% collect

    ggplot(df, aes(x=binIndex, y=gcnum/gcdenom)) + geom_point() + ggtitle(sprintf("G+C\n%s", region.descr())) + xlab("Bin") + ylab("Fraction G+C")
  })

  # plot obs/exp
  output$profile.plot <- renderPlot({
    log("profile.plot: ",bin.start(),'-',bin.end())
    df <- sampleFilter(regionFilter(bin.data())) %>% collect

    df$obsexp <- df$observed / df$expected
    df$sample.type <- ifelse(df$sample %in% c('',input$samples),'selected','other')

    (ggplot(df, aes(x=binIndex, group=sample, color=sample.type, alpha=sample.type, y=obsexp)) + geom_line() +
      xlab("Bin") + ylab("Observed / Expected") +
      ggtitle(sprintf("Observed/Expected\n%s", region.descr()))) %>% selected.theme
  })

  # plot observed
  output$coverage.plot <- renderPlot({
    log("coverage.plot: ",bin.start(),'-',bin.end())
    df <- sampleFilter(regionFilter(bin.data())) %>% collect %>% melt(measure.vars=c('observed','expected'))
    df$sample.type <- ifelse(df$sample %in% c('',input$samples),'selected','other')

    (ggplot(df, aes(x=binIndex, group=sample, color=sample.type, alpha=sample.type, y=value)) + geom_line() +
       facet_grid(variable~., scales="free_y") + xlab("Bin") +
       ggtitle(sprintf("Coverage\n%s", region.descr()))) %>% selected.theme
  })

  # plot obs/exp from windowed
  output$win.profile.plot <- renderPlot({
    log("win.profile.plot: ",bin.start(),'-',bin.end())
    df <- genoWinMap()

    df$obsexp <- df$observed/df$expected

    p <- ggplot(df, aes(color=sample.type, alpha=sample.type, group=sample, y=obsexp)) +
      xlab("Bin") + ylab("Observed / Expected") +
      ggtitle(sprintf("Windowed Observed/Expected\n%s", region.descr()))

    p %>% addGeomFeatures %>% selected.theme
  })

  # plot predicted genotypes in windows
  output$win.genotype.plot <- renderPlot({
    log("win.genotype.plot: ",bin.start(),'-',bin.end())
    df <- genoWinMap()

    p <- ggplot(df, aes(color=sample.type, alpha=sample.type, group=sample, y=cn)) +
      xlab("Bin") + ylab("Copy Number") + 
      ggtitle(sprintf("Windowed ProfileGenotyper Calls\n%s", region.descr()))

    p %>% addGeomFeatures %>% selected.theme
  })

  # plot predicted genotypes in windows, color CN and plot CNQ on Y
  output$win.CNQ.plot <- renderPlot({
    log("win.CNQ.plot: ",bin.start(),'-',bin.end())
    df <- genoWinMap()
    df$cn <- cn.label(df$cn)
    save(df,file="/tmp/df.Rdata")
    ggplot(filter(df, sample.type=='selected'), aes(y=cnq, color=cn)) +
      xlab("Bin") + ylab("Copy Number") + facet_grid(sample~.) +
      ggtitle(sprintf("Windowed ProfileGenotyper CNQ Calls\n%s", region.descr())) + cn.colors +
      geom_point(aes(x=targetBin))

  })

})
