library(shiny)
library(plyr)
library(dplyr)
library(ggplot2)
library(stats)
library(foreach)
library(zoo)
library(gridExtra)
library(tidyr)
library(RPostgreSQL)
library(reshape2)

# uncomment to debug within RStudio
#reactive <- function(expr, env=parent.frame()) exprToFunction(expr, env=env)

# Some hard-coded constants
pad <- 5000  # add pad bases to L and R of target region
too.big <- 500000  # don't retrieve big data when window is larger than too.big

source("../../conf/config.R")
read.conf("env.txt")  # pre-computed environment vars. Rerun setup.R to regenerate.

# Some file locations
tmp.dir <- Sys.getenv(("TMPDIR"))
input.dump.fn <- Sys.getenv("INPUT_DUMP")

# connection determined by environment variables
db <- src_postgres()

# Flags to enable different components
USE_IRS <- check.conf("USE_IRS")
USE_KNOWNS <- check.conf("USE_KNOWNS")
USE_QUARTETS <- check.conf("USE_QUARTETS")

if (USE_IRS) {
  # indata.dir is where the original inputs live. 
  indata.dir <- paste0(data.dir,"/../../gpc_wave2_batch1")  # FIXME
  
  irs.fn <- paste0(indata.dir,"/cnv_segs.irs")
  probe.fn <- paste0(indata.dir,"/probes.txt")

  irs.orig <- read.table(irs.fn,header=T,sep="\t", as.is=T)
  colnames(irs.orig) <- c('seg','chr','start','end','Pval','nprobes','# Samples','Lower Pval','Lower # Samples','Higher Pval', 'Higher # Samples')
  cat(sprintf("Loaded %s IRS rows from %s\n", nrow(irs.orig),irs.fn))
  
  probes.orig <- read.table(probe.fn)
  colnames(probes.orig) <- c('seg','chr','start.map','end.map')
  cat(sprintf("Loaded %s probes from %s\n",nrow(probes.orig),probe.fn))
}

# KNOWNS
if (USE_KNOWNS) {
  # DEL pipeline only for Sn
  gs_dels.fn <- paste0(data.dir,"/../../gpc_wave2_batch1/gs_dels_flt.genotypes.txt") # flattened, filtered
  
  # DEL pipeline — only sampleseg with most readpair evidence
  gs_dels_best.fn <- paste0(data.dir,"/../../gpc_wave2_batch1/gs_dels_best.genotypes.txt") # filtered, highest read pair
  
  # More generous set for Sp
  gs_cnvdels_flat.fn <- paste0(data.dir,"/../../gpc_wave2/gs_cnv_del_flt.genotypes.txt")
}

# NOTES - contains "bookmarks" of loci with corresponding notes
if (dbExistsTable(db$con, "notes")) {
  notes <- tbl(db, 'notes') %>% collect
} else {
  notes <- data.frame()
}

# common colors for copy number
cn.colors <- scale_color_manual(values=c('0'="#7fc97f", '1'="#c51b7d", '2'="#fdc086", '3'="#ffff99", '4'="#386cb0", '5+'="#f0027f",'Disc'="#999999"), 
                                name="CN")

cnv.mrg <- tbl(db, "cnv_mrg") %>% collect %>% mutate(method='Basic')
cnv.mle <- tbl(db, "cnv_mle") %>% collect %>% mutate(method='MLE')
cnv.post <- tbl(db, "cnv_post") %>% collect %>% mutate(method='Bayes')

geno <- tbl(db, "geno")
profile.segments <- tbl(db, "profile_segment")
profile.counts <- tbl(db, "profile_counts")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  pred.deletions <- reactive({
    if (input$pred.order == 'chrom') {
      cn.segs.basic <- arrange(cnv.mrg, start.map)
    } else if (input$pred.order == 'desc') {
      cn.segs.basic <- arrange(cnv.mrg, -(end.map-start.map))
    } else if (input$pred.order == 'asc') {
      cn.segs.basic <- arrange(cnv.mrg, (end.map-start.map))
    }
    cn.segs.basic$labels <- sprintf("%s: %s (%s)", seq(1,nrow(cn.segs.basic)), cn.segs.basic$seg, cn.segs.basic$end.map-cn.segs.basic$start.map)
    if (input$pred.cn==0) {
      x <- subset(cn.segs.basic, cn==0)
    } else if (input$pred.cn==1) {
      x <- subset(cn.segs.basic, cn==1)
    } else if (input$pred.cn==2) {
      x <- subset(cn.segs.basic, cn<2)
    } else if (input$pred.cn==3) {
      x <- subset(cn.segs.basic, cn!=2)
    } else {
      x <- subset(cn.segs.basic, cn>2)
    }
    segs <- x$seg
    names(segs) <- x$labels
    return(segs)
  })

  # # load the selection options for predicted sites
  observe({
    updateSelectInput(session, 'predicted', choices=c("Choose Predicted Deletion"="",pred.deletions()))
  })
  
  # load the selection options for saved sites
  updateSelectInput(session, 'saved.site', choices=c("Choose Preset"='',notes$name))
  
  if (USE_KNOWNS) {
    gs.dels.orig <- reactive({
      if (input$truth_data == 'GStrip Sn DEL Data') {
        gdo.fn <- gs_dels.fn
      } else if (input$truth_data == 'GStrip Sn Best DEL Data') {
        gdo.fn <- gs_dels_best.fn
      } else {
        gdo.fn <- gs_cnvdels_flat.fn
      }
      gdo <- read.table(gdo.fn, header=TRUE, sep="\t", stringsAsFactors = FALSE)
      cat(sprintf("Loaded %s known deletions from %s\n", nrow(gdo),gdo.fn))
      
      colnames(gdo) <- c('.id','seg','chr','start.map','end.map','cn','cq','paired.reads')
      gdo$evidence <- ' Known'
      gdo$copy.number <- addNA(as.factor(gdo$cn))
      return(gdo)
    })
  
    # load the selection options for known sites
    observe({
      updateSelectInput(session, 'candidate', choices=c("Choose Known Deletion"="",unique(gs.dels.orig()$seg)))
    })

    # if the candidate site changes, then set the samples, chromosome, start and end positions
    observe({
      if (input$candidate != "") {
        isolate(gdo <- gs.dels.orig())
        gs.row <- gdo[gdo$seg == input$candidate,]
        seg.sample <- gs.row[,'.id']
        seg.chr <- gs.row[1,'chr']
        seg.start <- gs.row[1,'start.map']
        seg.end <- gs.row[1,'end.map']
        updateTextInput(session,'seg.chr', value=seg.chr)
        updateNumericInput(session,'seg.start', value=seg.start)
        updateNumericInput(session,'seg.end', value=seg.end)
        updateSelectizeInput(session, 'seg.sample', selected=seg.sample)
      }
    })
    
    # prev candidate
    observe({
      input$candidate.prev
      candidates <- unique(gs.dels.orig()$seg)
      idx <- which(candidates == isolate(input$candidate))
      if (length(idx)>0 && idx > 1) { updateSelectInput(session,'candidate',selected=candidates[idx-1]) }
    })
    
    # next candidate
    observe({
      input$candidate.next
      candidates <- unique(gs.dels.orig()$seg)
      idx <- which(candidates == isolate(input$candidate))
      if (length(idx)>0 && idx < length(candidates)) { updateSelectInput(session,'candidate',selected=candidates[idx+1]) }
    })
  }
  
  observe({
    # load the values for all sample names
    updateSelectizeInput(session, 'seg.sample', choices=unique(cnv.mrg$.id))
  })
  
  # if the user selects a saved site, load those coords
  observe({
    input$reset
    input$saved.site
    isolate({
      
      if (input$saved.site != "") {
        note <- notes[notes$name == input$saved.site,]
        updateTextInput(session,'seg.chr', value=note$chr)
        updateNumericInput(session,'seg.start', value=note$start.map)
        updateNumericInput(session,'seg.end', value=note$end.map)
        updateSelectizeInput(session, 'seg.sample', selected=unlist(strsplit(note$sample,',',fixed=TRUE)))
        
        # must call special added hook to set value of textarea because there is no builtin R func
        session$sendCustomMessage(type = 'setnote', message = list(note = note$note))
      }
      
    })    
  })
  
  # prev saved
  observe({
    input$saved.prev
    idx <- which(notes$name == isolate(input$saved.site))
    if (length(idx)>0 && idx > 1) { updateSelectInput(session,'saved.site',selected=notes$name[idx-1]) }
  })
  
  # next saved
  observe({
    input$saved.next
    idx <- which(notes$name == isolate(input$saved.site))
    if (length(idx)>0 && idx < nrow(notes)) { updateSelectInput(session,'saved.site',selected=notes$name[idx+1]) }
  })
  
  # if the predicted site changes, then set the samples, chromosome, start and end positions
  observe({
    
    if (input$predicted != "") {
      gs.row <- filter(cnv.mrg, seg == input$predicted & !is.na(cn))
      cat(sprintf("Found %s rows for prediction %s\n", nrow(gs.row), input$predicted))
      print(gs.row)
      seg.sample <- gs.row$.id[gs.row$cn != 2]   # hack. very rarely the seg name is the same for the WT
      seg.chr <- gs.row$chr[1]
      seg.start <- gs.row$start.map[1]
      seg.end <- gs.row$end.map[1]
      updateTextInput(session,'seg.chr', value=seg.chr)
      updateNumericInput(session,'seg.start', value=seg.start)
      updateNumericInput(session,'seg.end', value=seg.end)
      updateSelectizeInput(session, 'seg.sample', selected=seg.sample)
    }
  })
  
  # prev predicted
  observe({
    input$predicted.prev
    isolate({
      idx <- which(pred.deletions() == input$predicted)
      if (length(idx)>0 && idx > 1) { updateSelectInput(session,'predicted',selected=pred.deletions()[idx-1]) }
    })
  })
  
  # next predicted
  observe({
    input$predicted.next
    isolate({
      idx <- which(pred.deletions() == input$predicted)
      if (length(idx)>0 && idx < length(pred.deletions())) { updateSelectInput(session,'predicted',selected=pred.deletions()[idx+1]) }
    })
  })
  
  
  # save input data structure to disk.
  # reload with:
  #   source(input.dump.fn); input <- input.dump
  observe({
    input$dump
    # save the latest selected site for interactive debug
    isolate({cat("Dump\n"); input.dump <- reactiveValuesToList(input); dump('input.dump',file=input.dump.fn)})
  })

  # zoom out 1.5x
  observe({
    input$zoom.out
    cat("zoom out\n")
    isolate({
      incr <- isolate(floor(.25 * (input$seg.end-input$seg.start)))
      updateNumericInput(session,'seg.start', value=max(0,isolate(input$seg.start - incr )))
      updateNumericInput(session,'seg.end', value=isolate(input$seg.end + incr))
    })
  })
  
  observe({
    input$delete
    cat("Delete Note\n")
    isolate({
      if (nrow(notes)>0) {
        notes <<- notes[notes$name != input$saved.site,]
        dbWriteTable(db$con, 'notes', notes, overwrite=TRUE, row.names=FALSE)
        notes <<- tbl(db, 'notes') %>% collect
        updateSelectInput(session, 'saved.site', choices=c("Choose Preset"='',notes$name))
      }
    })
  })
  
  # prepend site with notes to log file.
  observe({
    input$save.button
    isolate({
      if (input$note != "" && length(input$seg.sample) > 0) {
        cat("save\n")
        note.name <- sprintf("chr%s %sMb (%skb) %s",
                             input$seg.chr, format(input$seg.start/1000000,digits=4,nsmall=2),
                             format((input$seg.end-input$seg.start)/1000,digits=4,nsmall=2),
                             regmatches(input$note, regexpr("[^\n]{1,50}",input$note,perl=TRUE)))
        new.note <- tibble(name=note.name, chr=input$seg.chr, start.map=input$seg.start, end.map=input$seg.end,
                           sample=paste(input$seg.sample,collapse=','),note=input$note)
        if (nrow(notes)==0) { 
          notes <<- new.note 
        } else { 
          if (any(notes$name==note.name)) {
            notes <<- rbind(new.note, notes[-which(notes$name==note.name),])
          } else {
            notes <<- rbind(new.note, notes)  
          }
        }
        notes <<- unique(notes)
        dbWriteTable(db$con, 'notes', notes, overwrite=TRUE, row.names=FALSE)
        notes <<- tbl(db, 'notes') %>% collect
        updateSelectInput(session, 'saved.site', choices=c("Choose Preset"='',notes$name), selected=note.name)
      }
    })
  })
  
  region.filter <- function(cnv) {
    filter(cnv, (chr==input$seg.chr & end.map > input$seg.start-input$pad & start.map < input$seg.end+input$pad &
                   (cn != 2 | input$show_wildtype) & end.map-start.map > input$min.cnv.len & !is.na(cn) &
                   (!input$show_target_only | .id %in% input$seg.sample)))
  }
  
  # All the predicted extents. rbind the rows from different prediction stages 
  csm <- reactive({
    cat("Loading cn.segs.merged...\n")
    csm.all <- data.frame()
    
    if (input$show_basic) {
      csm.all <- rbind(csm.all, select(region.filter(cnv.mrg), .id, cn, chr, start.map, end.map, label=seg, evidence=method))
    }
    
    if (input$show_extended_ML) {
      csm.all <- rbind(csm.all, select(region.filter(cnv.mle), .id, cn, chr, start.map, end.map, label, evidence=method))
    }
    
    if (input$show_bayes) {
      csm.all <- rbind(csm.all, select(region.filter(cnv.post), .id, cn, chr, start.map, end.map, label, evidence=method))
    }
    csm.all$len <- csm.all$end.map - csm.all$start.map
    csm.all$target <- csm.all$.id %in% input$seg.sample
    cat("csm.all\n")
    print(head(csm.all))
    save(csm.all,file=sprintf("%s/csm.all.Rdata",tmp.dir))
    return(csm.all)
  })
  
  
  if (USE_IRS) {
    irs <- reactive({
      subset(irs.orig, chr==input$seg.chr & end > input$seg.start-input$pad & start < input$seg.end+input$pad & nprobes > 0)
    })
    
    probes <- reactive({
      subset(probes.orig, chr==input$seg.chr & end.map > input$seg.start-input$pad & start.map < input$seg.end+input$pad)
    })
  }
  
  if (USE_KNOWNS) {
    gs.dels <- reactive({
      x <- subset(gs.dels.orig(), chr==input$seg.chr & end.map > input$seg.start-input$pad & start.map < input$seg.end+input$pad & paired.reads >= input$paired.reads)
      x$target <- x$.id %in% input$seg.sample
      cat("gs.dels\n")
      print(head(x))
      save(x,file=sprintf('%s/gs.dels.Rdata',tmp.dir))
      return(x)  
    })
  }
  
  # format gs.dels, cn.segs, probes, and irs into a single data.frame for ggplot
  cnv.disp <- reactive({
    disp.cols <- c('.id','cn','chr','start.map','end.map','label', 'target', 'evidence')
    disp <- data.frame()

    if (USE_KNOWNS) {
      gs.dels <- gs.dels()
      gs.del.sel <- gs.dels[,c(disp.cols,'paired.reads')] 
      disp <- rbind(disp, gs.del.sel)
    }
    
    cn.segs <- csm()
    if (nrow(cn.segs) > 0) {
      cn.segs.sel <- cn.segs[,disp.cols]
      cn.segs.sel$paired.reads <- NA      
      disp <- rbind(disp, cn.segs.sel)
    } 
    
    if (USE_IRS) {
      irs <- irs()
      probes <- probes()
      if (nrow(probes) > 0 && FALSE) {  # remove probes 
        disp <- rbind(disp, 
                      data.frame(.id='Probe',cn=NA,copy.number=NA,chr=probes$chr,start.map=probes$start.map,end.map=probes$end.map,seg=probes$seg,target=FALSE,evidence="Array",paired.reads=NA))
      } 
      
      disp <- join(disp, subset(irs, select=-c(chr,start,end)), by="seg", type="left")
      disp$Pval <- ifelse(is.na(disp$Pval),1,disp$Pval)
      disp$Phred <- ifelse(disp$Pval==0, 100, -log10(disp$Pval)*10)
    }
    
    # max out CN at 5. Set fixed levels and labels for fixed legend
    disp$cn.disp <- factor(ifelse(disp$cn==99, disp$cn, ifelse(disp$cn > 5, 5, disp$cn)), levels=c('0','1','2','3','4','5','99'), labels=c('0','1','2','3','4','5+','Disc'))
    
    save(disp,file=sprintf('%s/disp.Rdata',tmp.dir))
    return(disp)
  })
  
  get.profile.segments <- reactive({
    x <- filter(profile.segments, chrom == input$seg.chr & start_pos > 
                  input$seg.start - input$pad & end_pos < input$seg.end + 
                  input$pad) %>% select(bin, chrom, start_pos, end_pos) %>% collect
  })
  
  frags <- reactive({
    if (input$seg.end - input$seg.start > too.big) {
      cat(sprintf("Skipping profile - region too big %s > %s", 
                  input$seg.end - input$seg.start, too.big))
      f <- data.frame()
    }
    else {
      start.bin <- min(get.profile.segments()$bin)
      end.bin <- max(get.profile.segments()$bin)
      
      pc <- filter(profile.counts, chrom == input$seg.chr & 
                     bin >= start.bin & bin <= end.bin) 
      
      if (input$show_target_frag_only) {
        pc <- filter(pc, sample %in% c("", input$seg.sample))
      }

      f <- inner_join(pc %>% collect, get.profile.segments(), by = c("chrom", "bin"))
      
      f$sample.type <- factor(ifelse(f$sample %in% 
                                       input$seg.sample, "target", "other"))
    }
    f
  })
  
  get.geno <- reactive({
    if (input$seg.end - input$seg.start > too.big) {
      cat(sprintf("Skipping profile - region too big %s > %s", 
                  input$seg.end - input$seg.start, too.big))
      return(data.frame())
    }
    else {
      start.bin <- min(get.profile.segments()$bin)
      end.bin <- max(get.profile.segments()$bin)
      g <- mutate(filter(geno, chr == input$seg.chr & bin >= start.bin & bin <= end.bin), bin = bin + 6)
      if (input$show_target_frag_only) {
        return(g %>% filter(sample %in% c("", input$seg.sample)))
      }
      else {
        return(g)
      }
    }
  })
  
  # breakpoint log likelihoods
  bkpts.ll <- reactive({
    if (input$seg.end-input$seg.start > too.big) {
      cat(sprintf("Skipping profile - region too big %s > %s", input$seg.end-input$seg.start, too.big))
      bkpts <- data.frame(sample=character(0), start_pos=integer(0), end_pos=integer(0), bkpt_ll=numeric(0), no_bkpt_ll=numeric(0))
    } else {
      # query each sample separately
      bkpts <-
        ldply(input$seg.sample, function(sample) {
          cat(sample,"\n")
          dbGetQuery(db$con, sprintf("SELECT b.sample, ps.start_pos, ps.end_pos, loss_ll, gain_ll, any_ll, no_bkpt_ll FROM bkpt b, profile_segment ps WHERE b.sample='%s' AND ps.chrom='%s' AND ps.start_pos < %s AND ps.end_pos > %s AND b.chr = ps.chrom AND b.bkpt_bin = ps.bin", 
                                     sample, input$seg.chr, input$seg.end+input$pad, input$seg.start-input$pad))
      })
      
      # # db values are p(x|gain,loss,nc). Turn into p(gain,loss,nc|x) = p(x|gain,loss,nc)p(gain,loss,nc) and normalize
      # # p(gain) = p(loss) = 6/16. p(nc) = 4/16. Assuming 4 levels of CN from staircase.R.
      # bkpts <- mutate(mutate(bkpts, loss = 6/16. * 10^-loss_ll, gain= 6/16. * 10^-gain_ll, any=12/16.*10^-any_ll,
      #                        no_bkpt = 4/16. * 10^-no_bkpt_ll,
      #                        tot=loss+gain+no_bkpt),
      #                 loss=loss/tot, gain=gain/tot, any=any/tot, no_bkpt=no_bkpt/tot,
      #                 Lloss=-log10(loss), Lgain=-log10(gain), Lany=-log10(any), Lno_bkpt=-log10(no_bkpt),
      #                 Lno_loss=-log10(gain+no_bkpt), Lno_gain=-log10(loss+no_bkpt))
      
      bkpts <- mutate(bkpts,
                      loss=10^-loss_ll, 
                      gain=10^-gain_ll, 
                      any=10^-any_ll, 
                      no_bkpt=10^-no_bkpt_ll,
                      nc=no_bkpt,
                      lossZ=loss/(loss+gain+no_bkpt),
                      gainZ=gain/(loss+gain+no_bkpt),
                      no_bkptZ=no_bkpt/(loss+gain+no_bkpt),
                      LlossZ=-log10(lossZ),
                      LgainZ=-log10(gainZ),
                      Lno_loss=-log10(1-loss),
                      Lno_gain=-log10(1-gain),
                      Lno_lossZ=-log10(1-lossZ),
                      Lno_gainZ=-log10(1-gainZ))
    }
    save(bkpts,file=sprintf('%s/bkpts.Rdata',tmp.dir))
    return(bkpts)
  })
  
  
  # observed/expected ratio
  oer <- reactive({
    f <- frags()
    oer <-
      mutate(ddply(f, .(sample), mutate, 
                   obs=c(rep(NA,input$win.size/2), rollapply(observed, width=input$win.size, sum), rep(NA,(input$win.size/2-1))),
                   exp=c(rep(NA,input$win.size/2), rollapply(expected, width=input$win.size, sum), rep(NA,(input$win.size/2-1)))),
             obs_exp_ratio=obs/exp,
             sample.type=ifelse(sample %in% input$seg.sample, 'target','other'))
             
    oer.median <- mutate(ddply(oer, .(bin, start_pos, end_pos), summarize, obs_exp_ratio=median(obs_exp_ratio)),
                         chrom=as.character(input$seg.chr), sample='site.median',
                         sample.type='median')
    
    oer <- rbind(select(oer, chrom, bin, start_pos, end_pos, sample, obs_exp_ratio, sample.type),oer.median)
    save(oer,file=sprintf('%s/oer.Rdata',tmp.dir))
    return(oer)  
  })
  
  # helper functions for plotting
  xbounds <- reactive({
    xmin = min(get.profile.segments()$start_pos)
    xmax = max(get.profile.segments()$end_pos)
    cat(sprintf("xmin=%s xmax=%s\n", xmin, xmax))
    coord_cartesian(xlim = c(xmin, xmax))
  })
  ylims <- reactive({
    ylim(min(get.profile.segments()$start_pos), max(get.profile.segments()$end_pos))
  })
  flip.bounds <- function(coord) {
    coord_cartesian(xlim=coord$limits$y, ylim=coord$limits$x)
  }
  
  # display the rectangle segment plot of known, probes, array, and different predictions
  cnvPlot <- reactive({
    cat("cnvPlot...\n")
    
    cnv.disp <- cnv.disp()
    if (nrow(cnv.disp)==0) { return(NULL) }
    if (input$show_target_only && nrow(cnv.disp)>0) { cnv.disp <- subset(cnv.disp, target==TRUE | evidence=='Array') }
    cat("cnv.disp\n")
    print(head(cnv.disp))
    save(cnv.disp,file=sprintf('%s/cnv.disp.Rdata',tmp.dir))
    
    if (!input$show_readPairs) {
      plt <- ggplot(cnv.disp, aes(x=start.map, xend=end.map, y=.id, yend=.id, color=cn.disp, size=target)) + cn.colors 
    } else {
      plt <- ggplot(cnv.disp, aes(x=start.map, xend=end.map, y=.id, yend=.id, color=paired.reads, size=target)) + scale_color_gradient2() 
    }
    
    if (nrow(subset(cnv.disp,evidence=='Array'))>0) {
      plt <- plt + geom_point(data=subset(cnv.disp, evidence=='Array'), size=5, color='black')
    }
    plt + geom_segment() + facet_grid(evidence ~ ., scales="free", space="free") + guides(size="none") + theme(axis.title.x = element_blank(), axis.title.y=element_blank()) + xbounds()
  })
  
  cnvCIPlot <- reactive({
    cat("cnvCIPlot...\n")
    segs <- select(region.filter(cnv.mle), chr, end.map, start.map, 
                   cn, .id, start.map.L, start.map.R, end.map.L, end.map.R, 
                   method)
    if (input$show_bayes) {
      segs <- rbind(segs, select(region.filter(cnv.post), chr, 
                                 end.map, start.map, cn, .id, start.map.L = start.CI.L, 
                                 start.map.R = start.CI.R, end.map.L = end.CI.L, end.map.R = end.CI.R, 
                                 method))
    }
    segs$target <- segs$.id %in% input$seg.sample
    print(head(segs))
    ggplot(segs, aes(x = .id, y = start.map, ymin = start.map.L, 
                     ymax = start.map.R, color = target)) + geom_pointrange(size = 1.5) + 
      geom_pointrange(aes(y = end.map, ymin = end.map.L, ymax = end.map.R), 
                      size = 1.5) + geom_segment(aes(y = start.map, yend = end.map, 
                                                     xend = .id), size = 1.5, linetype = 3) + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(), plot.margin = unit(c(0, 0.75, 0, 0.1), "inches")) + guides(color = "none") + 
      ylims() + coord_flip() + facet_grid(method ~ .)
  })
  
  fragPlot <- reactive({
    cat("fragPlot...\n")
    df <- oer()
    ggplot(df, aes(x=start_pos, y=obs_exp_ratio, group=sample)) + geom_step(color='grey') +
      geom_step(data=df[df$sample.type=='median',], color='black') + 
      geom_step(data=df[df$sample.type=='target',], color='red') +
      theme(axis.title.x = element_blank(), plot.margin=unit(c(0,1,0,0.5),'inches')) +
      xbounds()
  })
  
  bkptPriors <- reactive({
    if (input$seg.end - input$seg.start > too.big) {
      cat(sprintf("Skipping profile - region too big %s > %s", 
                  input$seg.end - input$seg.start, too.big))
      return(data.frame())
    }
    else {
      bkpts <- dbGetQuery(db$con, sprintf("SELECT b.sample, ps.start_pos, ps.end_pos, loss_ll, gain_ll, no_bkpt_ll FROM bkpt b, profile_segment ps WHERE ps.chrom ='%s' AND ps.start_pos < %s AND ps.start_pos > %s AND b.bkpt_bin = ps.bin", 
                                          input$seg.chr, input$seg.end + input$pad, input$seg.start - 
                                            input$pad))
      if (input$bkpt_prior_samples == "Selected") {
        bkpts <- filter(bkpts, sample %in% input$seg.sample)
      }
      else if (input$bkpt_prior_samples == "Exclude Selected") {
        bkpts <- filter(bkpts, !sample %in% input$seg.sample)
      }
      if (input$bkpt_prior_CN != "Any") {
        stopifnot(USE_KNOWNS)
        if (input$bkpt_prior_CN == "0 or 1") {
          CN.match = c(0, 1)
        }
        else {
          CN.match = as.integer(input$bkpt_prior_CN)
        }
        gdo <- read.table(gs_dels.fn, header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE)
        colnames(gdo) <- c(".id", "seg", "chr", "start.map", 
                           "end.map", "cn", "cq", "paired.reads")
        x <- subset(gdo, chr == input$seg.chr & end.map > 
                      input$seg.start - input$pad & start.map < input$seg.end + 
                      input$pad & cn %in% CN.match)
        bkpts <- filter(bkpts, sample %in% x$.id)
      }
      cat(sprintf("bkpts: %s rows\n", nrow(bkpts)))
      cat(sprintf("samples: %s\n", length(unique(bkpts$sample))))
      bkpts <- mutate(bkpts, loss = 10^-loss_ll, gain = 10^-gain_ll, 
                      nc = 10^-no_bkpt_ll, tot = loss + gain + nc, lossZ = loss/tot, 
                      gainZ = gain/tot, ncZ = nc/tot, Lloss = -log10(lossZ), 
                      Lgain = -log10(gainZ), Lnc = -log10(ncZ), Lno_loss = -log10(1 - 
                                                                                    lossZ), Lno_gain = -log10(1 - gainZ))
      bkpt.prior <- mutate(ddply(bkpts, .(start_pos), summarize, 
                                 Lall_no_loss = sum(Lno_loss), Lall_no_gain = sum(Lno_gain), 
                                 Lall_nc = sum(Lnc), some_loss = 1 - 10^(-Lall_no_loss), 
                                 some_gain = 1 - 10^(-Lall_no_gain), some_nc = 1 - 
                                   10^(-Lall_nc), Lsome_loss = -log10(some_loss), 
                                 Lsome_gain = -log10(some_gain), Lsome_nc = -log10(some_nc)))
      save(bkpt.prior, file = sprintf("%s/bkpt_prior.Rdata", 
                                      tmp.dir))
      return(bkpt.prior)
    }
  })
  
  bayesPriorPlot <- reactive({
    cat("bayesPriorPlot...\n")
    priors <- dbGetQuery(db$con, sprintf("SELECT p.region_id, p.start_pos, p.\"loss.u\", p.\"gain.u\", pr.side FROM prior p, prior_region pr, profile_segment ps1, profile_segment ps2 WHERE p.region_id=pr.id AND pr.chr=ps1.chrom AND pr.binL>=ps1.bin AND ps1.chrom='%s' AND ps1.start_pos <= %s AND ps1.end_pos >= %s AND pr.chr=ps2.chrom AND pr.binR<=ps2.bin AND ps2.chrom='%s' AND ps2.start_pos <= %s AND ps2.end_pos >= %s", 
                                         input$seg.chr, input$seg.start - input$pad, input$seg.start - 
                                           input$pad, input$seg.chr, input$seg.end + input$pad, 
                                         input$seg.end + input$pad))
    if (nrow(priors) > 0) {
      priors.m <- melt(priors, c("start_pos", "side"), c("loss.u", "gain.u"))
      ggplot(priors.m, aes(x = start_pos, y = value, color = variable)) + 
        geom_point() + geom_line() + facet_grid(side ~ .) + 
        theme(axis.title.x = element_blank(), plot.margin = unit(c(0, 0.75, 0, 0.4), "in"), legend.position = "bottom") + 
        ylab("Prob") + xbounds()
    }
    else {
      NULL
    }
  })
  
  posteriorPlot <- reactive({
    cat("posteriorPlot...\n")
    posterior <- dbGetQuery(db$con, sprintf("SELECT po.start_pos, po.sample, po.\"loss.u\", po.\"gain.u\", po.loss, po.gain, po.bayes_loss, po.bayes_gain FROM posterior_dist po, profile_segment ps1, profile_segment ps2 WHERE po.chr=ps1.chrom AND po.bin>=ps1.bin AND ps1.chrom='%s' AND ps1.start_pos <= %s AND ps1.end_pos >= %s AND po.chr=ps2.chrom AND po.bin<=ps2.bin AND ps2.chrom='%s' AND ps2.start_pos <= %s AND ps2.end_pos >= %s AND po.sample IN ('%s')", 
                                            input$seg.chr, input$seg.start - input$pad, input$seg.start - 
                                              input$pad, input$seg.chr, input$seg.end + input$pad, 
                                            input$seg.end + input$pad, paste(input$seg.sample, collapse = "','")))
    if (nrow(posterior) > 0) {
      posterior.m <- mutate(melt(posterior, c("start_pos", 
                                              "sample"), c("loss", "gain", "bayes_loss", "bayes_gain")), 
                            kind = ifelse(grepl("bayes_", variable), "Bayes", 
                                          "Likelihood"))
      ggplot(posterior.m, aes(x = start_pos, y = value, color = variable)) + 
        geom_point() + geom_line() + facet_grid(kind + sample ~ ., scales = "free_y") + 
        theme(axis.title.x = element_blank(), plot.margin = unit(c(0, 0.6, 0, 0), "in"), legend.position = "bottom") + 
        ylab("Prob") + xbounds()
    }
    else {
      NULL
    }
  })
  
  bkptPlot <- reactive({
    cat("bkptPlot...\n")
    bkpt.prior <- bkptPriors()
    bkpt.prior.m <- mutate(melt(bkpt.prior, "start_pos", c("Lall_no_loss", 
                                                           "Lall_no_gain")), kind = ifelse(grepl("_ratio", variable), 
                                                                                           "Log Ratio", ifelse(grepl("_ll", variable), "Broken Log Likelihood", 
                                                                                                               "Log No Loss/Gain")))
    ggplot(bkpt.prior.m, aes(x = start_pos, y = value, color = variable)) + 
      geom_point() + geom_line() + facet_grid(kind ~ ., scales = "free_y") + 
      theme(axis.title.x = element_blank(), plot.margin = unit(c(0, 0.75, 0, 0.4), "in"), legend.position = "bottom") + 
      ylab("Log Prob") + xbounds()  
  })
  
  eachBkptPlot <- reactive({
    cat("eachBkptPlot...\n")
    bkpts <- bkpts.ll()
    if (nrow(bkpts)>0) {
      
      bkpts <- melt(bkpts, c('sample','start_pos','end_pos'),c('loss','gain')) 
      ggplot(bkpts, aes(x=start_pos, y=value, color=variable)) + geom_point() + geom_line() + facet_grid(sample~.) + theme(axis.title.x = element_blank(),plot.margin=unit(c(0,0.75,0,.25),'in'),legend.position="bottom") + ylab("Prob") + xbounds()
    } else {
      ggplot() + geom_blank()
    }
  })
  
  # merge frag and geno segments and color fragment by geno label
  # requires the windows to be the same size
  genoFragPlot <- reactive({
    oer.wg <- left_join(filter(oer(), sample.type != "median"), 
                        collect(get.geno()), by = c("bin", "sample"))
    oer.wg$cn <- as.numeric(as.character(oer.wg$cn))
    oer.wg$cn <- ifelse(is.na(oer.wg$cn), 99, oer.wg$cn)
    oer.wg$sample.type.size <- ifelse(oer.wg$sample.type == "target", 
                                      2, 1)
    oer.wg$cn.disp <- factor(ifelse(oer.wg$cn == 99, oer.wg$cn, 
                                    ifelse(oer.wg$cn > 5, 5, oer.wg$cn)), 
                             levels = c("0", "1", "2", "3", "4", "5", "99"), 
                             labels = c("0", "1", "2", "3", "4", "5+", "Disc"))
    ggplot(oer.wg, aes(x = start_pos, y = obs_exp_ratio, yend = obs_exp_ratio, 
                       color = cn.disp, group = sample)) + geom_step() + cn.colors + 
      theme(axis.title.x = element_blank(), plot.margin = unit(c(0, 1, 0, 0.5), "inches")) + 
      xbounds() + guides(color = FALSE)
    })

  output$expPlot <- renderPlot({
    expected <- frags()
    colors <- c(target='#FF0000',other='#222222')
    alpha <- c(target=1,other=0.2)
    ggplot(expected, aes(x=start_pos,y=expected,group=sample,color=sample.type,alpha=sample.type)) + geom_line() +
      theme(axis.title.x = element_blank(), plot.margin=unit(c(0,1,0,0),'inches')) + 
      scale_color_manual(name = "sample.type",values = colors) +
      scale_alpha_manual(name = "sample.type",values = alpha) +
      xbounds() + guides(color=FALSE,alpha=FALSE)
  })
  
  output$seg.size <- renderText({ paste0(format((input$seg.end-input$seg.start)/1000,digits=4,nsmall=2),"kb") })
  
  output$allThree <- renderPlot({
    plots <- list()
    if (input$show_cnv) { plots$cnv = cnvPlot() }
    if (input$show_CI) { plots$CI = cnvCIPlot() }
    if (input$show_frag) { if (input$color_frag) { plots$frag=genoFragPlot() } else { plots$frag = fragPlot() } }
    if (input$show_each_bkpt) { plots$bkpt = eachBkptPlot() }
    if (input$show_bayes_prior) { p <- bayesPriorPlot(); if (!is.null(p)) { plots$bayes.prior=p } }
    if (input$show_posterior) { p <- posteriorPlot(); if (!is.null(p)) { plots$posterior=p} }
    grid.arrange(do.call(arrangeGrob,c(plots,list(ncol=1))))
  })
  
  output$cnqs <- renderPlot({
    wg <- filter(get.geno(), sample %in% c('',input$seg.sample))
    if (USE_KNOWNS) {
      gs.dels <- subset(gs.dels(), .id %in% input$seg.sample)
    } else {
      gs.dels <- data.frame(.id=NULL)
    }

    plots <- lapply(input$seg.sample, function(this.sample) {
      p <- ggplot(collect(filter(wg, sample==this.sample)), aes(x=start,y=cnq, color=as.factor(cn))) + geom_point(size=3) + cn.colors + ggtitle(this.sample)
      gs.dels.id <- filter(gs.dels, .id==this.sample) 
      if (nrow(gs.dels.id)==0) {
        p
      } else {
        p + geom_segment(data=gs.dels.id, mapping=aes(x=start.map, xend=end.map, y=110, yend=110, 
                                                                         color=as.factor(copy.number)))
      }     
    })
    grid.arrange(do.call(arrangeGrob,c(plots,list(ncol=1))))
    
  })
  
  bkpts.bayes <- reactive({
    bkpt.prior <- bkptPriors()
    bkpts <- bkpts.ll()
    bkpt.mrg <- mutate(merge(bkpt.prior, bkpts), 
                       bayes_loss = loss * some_loss, 
                       bayes_gain = gain * some_gain, 
                       bayes_nc = nc * some_nc, 
                       bayes_tot = bayes_loss + bayes_gain + bayes_nc, 
                       bayes_lossZ = bayes_loss/bayes_tot, 
                       bayes_gainZ = bayes_gain/bayes_tot, 
                       bayes_no_loss = 1 - bayes_loss, 
                       bayes_no_gain = 1 - bayes_gain, 
                       bayes_no_loss2 = (1 - loss) * (1 - some_loss), 
                       bayes_no_gain2 = (1 - gain) * (1 - some_gain), 
                       bayes_no_lossZ = 1 - bayes_lossZ, 
                       bayes_no_gainZ = 1 - bayes_gainZ, 
                       Lbayes_loss = -log10(bayes_loss), 
                       Lbayes_gain = -log10(bayes_gain), 
                       Lbayes_lossZ = -log10(bayes_lossZ), 
                       Lbayes_gainZ = -log10(bayes_gainZ), 
                       Lbayes_no_loss = -log10(bayes_no_loss), 
                       Lbayes_no_gain = -log10(bayes_no_gain), 
                       Lbayes_no_loss2 = -log10(bayes_no_loss2), 
                       Lbayes_no_gain2 = -log10(bayes_no_gain2), 
                       loss.ratio = log10(loss/(no_bkpt + gain)), 
                       gain.ratio = log10(gain/(no_bkpt + loss)), 
                       loss.ratioZ = log10(lossZ/(no_bkptZ + gainZ)), 
                       gain.ratioZ = log10(gainZ/(no_bkptZ + lossZ)), 
                       bayes.loss.ratio = log10(bayes_loss/bayes_no_loss), 
                       bayes.gain.ratio = log10(bayes_gain/bayes_no_gain), 
                       bayes.loss.ratioZ = log10(bayes_lossZ/bayes_no_lossZ), 
                       bayes.gain.ratioZ = log10(bayes_gainZ/bayes_no_gainZ))  
  })
  
  output$bkpts <- renderPlot({
    bkpt.mrg <- bkpts.bayes()
        
    if (input$show_probs == 'probs') {
      show.vars <- c('bayes_loss','some_loss','loss','bayes_gain','some_gain','gain','no_bkpt')
    } else if (input$show_probs == 'no') {
      show.vars <- c('Lbayes_no_loss2', 'Lall_no_loss', 'Lno_loss', 'Lbayes_no_gain2','Lall_no_gain','Lno_gain')
    } else if (input$show_probs == 'Z') {
      show.vars <- c('Lbayes_lossZ','Lsome_loss','LlossZ','Lbayes_gainZ','Lsome_gain','LgainZ')
    } else {
      show.vars <- c('Lbayes_loss','Lsome_loss','loss_ll','Lbayes_gain','Lsome_gain','gain_ll')
    }
    bkpts.plt <- melt(bkpt.mrg, c('sample','start_pos','end_pos'), show.vars)
    bkpts.plt <- mutate(bkpts.plt, kind=ifelse(grepl('loss',variable),'Loss',ifelse(grepl('gain',variable),'Gain','NC')))
        
    last.sample <- input$seg.sample[length(input$seg.sample)]
    plots <- lapply(input$seg.sample, function(this.sample) {
      p <- ggplot(subset(bkpts.plt, sample==this.sample), aes(x=start_pos, y=value, color=variable)) + geom_point() + geom_line() + facet_grid(kind~.) + theme(axis.title.x = element_blank(),plot.margin=unit(c(0,0.35,0,.35),'in'),legend.position="bottom") + ggtitle(this.sample) 
      if (this.sample != last.sample) {
        return(p + guides(color='none') + theme(axis.text.x=element_blank()))
      } else {
        return(p)
      }
    })
    grid.arrange(do.call(arrangeGrob, c(plots,list(ncol=1))))
  })

  output$bkpt.odds <- renderPlot({
    bkpt.mrg <- bkpts.bayes()
    
    if (input$show_bayes_odds == 'bayesZ') {
      show.vars <- c('bayes.loss.ratioZ','bayes.gain.ratioZ')
    } else if (input$show_bayes_odds == 'bayes') {
      show.vars <- c('bayes.loss.ratio','bayes.gain.ratio')
    } else {
      show.vars <- c('loss.ratio','gain.ratio')
    }

    bkpts.plt <- melt(bkpt.mrg, c('sample','start_pos','end_pos'), show.vars)
    bkpts.plt <- mutate(bkpts.plt, kind=ifelse(grepl('loss',variable),'Loss',ifelse(grepl('gain',variable),'Gain','NC')),
                        pos=value>0)
    
    last.sample <- input$seg.sample[length(input$seg.sample)]
    plots <- lapply(input$seg.sample, function(this.sample) {
      p <- ggplot(subset(bkpts.plt, sample==this.sample), aes(x=start_pos, y=value, color=pos, group=variable)) + geom_point() + geom_line() + geom_hline(yintercept=0) + facet_grid(kind~.) + theme(axis.title.x = element_blank(),plot.margin=unit(c(0,0.35,0,.35),'in'),legend.position="bottom") + ggtitle(this.sample) 
      return(p + guides(color='none') + theme(axis.text.x=element_blank()))
    })
    grid.arrange(do.call(arrangeGrob, c(plots,list(ncol=1))))
  })
  
})
