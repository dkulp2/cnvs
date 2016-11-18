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

source("config.R")

db.conn.params <- as.list(unlist(strsplit(db.conn.str,":")))
names(db.conn.params) <- c('user','host','port','dbname')
db <- do.call(src_postgres, db.conn.params)

# input <-
#   structure(list(candidate = "DEL_P0563_29", seg.chr = 20L, seg.end = 36477348L, 
#                  seg.start = 36474476L), .Names = c("candidate", "seg.chr", 
#                                                     "seg.end", "seg.start"))

irs.fn <- paste0(data.dir,"/cnv_segs.irs")
probe.fn <- paste0(data.dir,"/probes.txt")

# KNOWNS
#gs_dels.fn <- paste0(data.dir,"/gs_dels.txt")
gs_dels.fn <- paste0(data.dir,"/../../gpc_wave2_batch1/gs_dels_flt.genotypes.txt") # flattened, filtered
# NOT USED: gs_cnvs.fn <- paste0(data.dir,"../gpc_wave2/gs_cnv.genotypes.txt")
gs_cnvdels_flat.fn <- paste0(data.dir,"/../../gpc_wave2/gs_cnv_del_flt.genotypes.txt")

cnv.geno.fn <- paste0(data.dir,"/sites_cnv_segs.txt.cnvgeno.srt.gz")
cnv.hires.geno.fn <- paste0(data.dir,"/hires_sites_cnv_segs.txt.cnvgeno.srt.gz")

log.fn <- paste0(data.dir,"/notes.Rdata")
if (file.exists(log.fn)) {
  load(log.fn)
  cat(sprintf("Loaded %s notes from %s",nrow(notes),log.fn))
} else {
  notes <- data.frame()
}

pad <- 5000  # add pad bases to L and R of target region
too.big <- 500000  # don't retrieve big data when window is larger than too.big

# common colors for copy number
cn.colors <- scale_color_manual(values=c('0'="#7fc97f", '1'="#c51b7d", '2'="#fdc086", '3'="#ffff99", '4'="#386cb0", '5+'="#f0027f",'Disc'="#999999"), 
                                name="CN")

frags.cmd <- sprintf("%s 'zcat %s | head -1'", shell, profile.fn)
cat(frags.cmd,"\n",file=stderr())
frags.header <- unlist(strsplit(system(frags.cmd, intern=T),"\t"))

genoWin.cmd <- sprintf("%s 'zcat %s | head -1'", shell, cnv.geno.fn)
cat(genoWin.cmd,"\n",file=stderr())
genoWin.header <- unlist(strsplit(system(genoWin.cmd, intern=T),"[\r\t]",perl=TRUE))

irs.orig <- read.table(irs.fn,header=T,sep="\t", as.is=T)
colnames(irs.orig) <- c('seg','chr','start','end','Pval','nprobes','# Samples','Lower Pval','Lower # Samples','Higher Pval', 'Higher # Samples')
cat(sprintf("Loaded %s IRS rows from %s\n", nrow(irs.orig),irs.fn))

probes.orig <- read.table(probe.fn)
colnames(probes.orig) <- c('seg','chr','start.map','end.map')
cat(sprintf("Loaded %s probes from %s\n",nrow(probes.orig),probe.fn))

cn.segs.merged.fn <- paste0(data.dir,"/sites_cnv_segs.txt.smlcsm.Rdata")
load(cn.segs.merged.fn)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
  pred.deletions <- reactive({
    if (input$pred.order == 'chrom') {
      cn.segs.basic <- arrange(cn.segs.merged, start.map)
    } else if (input$pred.order == 'desc') {
      cn.segs.basic <- arrange(cn.segs.merged, -(end.map-start.map))
    } else if (input$pred.order == 'asc') {
      cn.segs.basic <- arrange(cn.segs.merged, (end.map-start.map))
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
  
  gs.dels.orig <- reactive({
    if (input$truth_data == 'GStrip Sn DEL Data') {
      gdo.fn <- gs_dels.fn
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
  
  observe({
    # load the values for all sample names
    updateSelectizeInput(session, 'seg.sample', choices=unique(gs.dels.orig()$.id))
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
      gs.row <- cn.segs.merged[cn.segs.merged$seg == input$predicted & !is.na(cn.segs.merged$cn),]
      cat(sprintf("Found %s rows for prediction %s\n", nrow(gs.row), input$predicted))
      print(gs.row)
      seg.sample <- gs.row[gs.row$cn != 2,'.id']   # hack. very rarely the seg name is the same for the WT
      seg.chr <- gs.row[1,'chr']
      seg.start <- gs.row[1,'start.map']
      seg.end <- gs.row[1,'end.map']
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
    isolate({
      notes <<- notes[notes$name != input$saved.site,]
      save(notes, file=log.fn)
      updateSelectInput(session, 'saved.site', choices=c("Choose Preset"='',notes$name))
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
        new.note <- data.frame(name=note.name, chr=input$seg.chr, start.map=input$seg.start, end.map=input$seg.end,
                               sample=paste(input$seg.sample,collapse=','),note=input$note,
                               stringsAsFactors = FALSE)
        if (nrow(notes)==0) { notes <<- new.note } else { 
          if (any(notes$name==note.name)) {
            notes <<- rbind(new.note, notes[-which(notes$name==note.name),])
          } else {
            notes <<- rbind(new.note, notes)  
          }
        }
        notes <<- unique(notes)
        save(notes, file=log.fn)
        updateSelectInput(session, 'saved.site', choices=c("Choose Preset"='',notes$name), selected=note.name)
      }
    })
  })
  
  # All the predicted extents, i.e. merged (contiguous) copy number segments, aka cn.segs.merged 
  csm <- reactive({
    cat("Loading cn.segs.merged...\n")
    # if (input$show_extended_cnvs) {
    #   cn.segs.merged.fn <- paste0(data.dir,"/sites_cnv_segs.txt.ncsm.Rdata")
    # } else if (input$show_extended_ML) {
    #   cn.segs.merged.fn <- paste0(data.dir,"/sites_cnv_segs.txt.mlcsm.Rdata")
    # } else {
    #   cn.segs.merged.fn <- paste0(data.dir,"/sites_cnv_segs.txt.csm.Rdata")
    # }
    csm.all <- data.frame()
    
    if (input$show_basic) {
      cn.segs.merged.fn <- paste0(data.dir,"/sites_cnv_segs.txt.csm.Rdata")
      load(cn.segs.merged.fn)
      cn.segs.merged <- select(cn.segs.merged, .id, cn, chr, start.map, end.map, copy.number, len, seg)
      cn.segs.merged$label <- sprintf("%s_%s", cn.segs.merged$seg, cn.segs.merged$.id)
      cn.segs.merged$evidence <- 'Basic'
      cat(sprintf("Loaded %s CNVs from %s\n",nrow(cn.segs.merged), cn.segs.merged.fn))
      csm.all <- rbind(csm.all, cn.segs.merged)
    }
    
    if (input$show_extended_ML) {
      cn.segs.merged.fn <- paste0(data.dir,"/sites_cnv_segs.txt.smlcsm.Rdata")
      load(cn.segs.merged.fn)
      cn.segs.merged <- select(cn.segs.merged, .id, cn, chr, start.map, end.map, copy.number, len, seg)
      cn.segs.merged$label <- sprintf("%s_%s", cn.segs.merged$seg, cn.segs.merged$.id)
      cn.segs.merged$evidence <- 'MLE'
      cat(sprintf("Loaded %s CNVs from %s\n",nrow(cn.segs.merged), cn.segs.merged.fn))
      csm.all <- rbind(csm.all, cn.segs.merged)
    }
    
    return(csm.all)
  })
  
  
  # a subset of the extents
  cn.segs <- reactive({
    csm <- csm()
    x <- csm[csm$chr == input$seg.chr & 
               csm$end.map > input$seg.start-input$pad & 
               csm$start.map < input$seg.end+input$pad & 
               (csm$cn != 2 | input$show_wildtype) & 
               csm$len > input$min.cnv.len & 
               !is.na(csm$cn),]
    x$target <- x$.id %in% input$seg.sample
    cat("cn.segs\n")
    print(head(x))
    save(x,file=sprintf("%s/cn.segs.Rdata",tmp.dir))
    return(x)
  })
  
  irs <- reactive({
    subset(irs.orig, chr==input$seg.chr & end > input$seg.start-input$pad & start < input$seg.end+input$pad & nprobes > 0)
  })
  
  probes <- reactive({
    subset(probes.orig, chr==input$seg.chr & end.map > input$seg.start-input$pad & start.map < input$seg.end+input$pad)
  })
  
  gs.dels <- reactive({
    x <- subset(gs.dels.orig(), chr==input$seg.chr & end.map > input$seg.start-input$pad & start.map < input$seg.end+input$pad & paired.reads >= input$paired.reads)
    x$target <- x$.id %in% input$seg.sample
    cat("gs.dels\n")
    print(head(x))
    save(x,file=sprintf('%s/gs.dels.Rdata',tmp.dir))
    return(x)  
  })
  
  # format gs.dels, cn.segs, probes, and irs into a single data.frame for ggplot
  cnv.disp <- reactive({
    gs.dels <- gs.dels()
    cn.segs <- cn.segs()
    irs <- irs()
    probes <- probes()
    
    disp.cols <- c('.id','cn','copy.number','chr','start.map','end.map','seg', 'target', 'evidence')
    gs.del.sel <- gs.dels[,c(disp.cols,'paired.reads')] 
    cn.segs.sel <- cn.segs[,disp.cols]
    if (nrow(cn.segs.sel)>0) { cn.segs.sel$paired.reads <- NA }
    disp <- rbind(gs.del.sel, cn.segs.sel)
    
    if (nrow(probes) > 0 && FALSE) {  # remove probes 
      disp2 <- rbind(disp, 
                     data.frame(.id='Probe',cn=NA,copy.number=NA,chr=probes$chr,start.map=probes$start.map,end.map=probes$end.map,seg=probes$seg,target=FALSE,evidence="Array",paired.reads=NA))
    } else { disp2 <- disp }
    
    disp3 <- join(disp2, subset(irs, select=-c(chr,start,end)), by="seg", type="left")
    disp3$Pval <- ifelse(is.na(disp3$Pval),1,disp3$Pval)
    disp3$Phred <- ifelse(disp3$Pval==0, 100, -log10(disp3$Pval)*10)
    
    # max out CN at 5. Set fixed levels and labels for fixed legend
    disp3$cn.disp <- factor(ifelse(disp3$cn==99, disp3$cn, ifelse(disp3$cn > 5, 5, disp3$cn)), levels=c('0','1','2','3','4','5','99'), labels=c('0','1','2','3','4','5+','Disc'))
    
    save(disp3,file=sprintf('%s/disp3.Rdata',tmp.dir))
    return(disp3)
  })
  
  winGeno <- reactive({
    # load windowed genotypes
    wg.fn <- basename(tempfile("wg"))
    wg.in.fn <- ifelse(input$show_hires, cnv.hires.geno.fn, cnv.geno.fn) 
    if (input$seg.end-input$seg.start > too.big) {
      cat(sprintf("Skipping genotypes - region too big %s > %s", input$seg.end-input$seg.start, too.big))
      wg <- data.frame(chr=NULL,start=NULL,end=NULL,sample=NULL,cn=NULL,cnq=NULL,pos=NULL)
    } else {
      tabix <- sprintf("tabix -h %s %s:%s-%s > /tmp/%s", wg.in.fn, input$seg.chr, input$seg.start-input$pad, input$seg.end+input$pad, wg.fn)
      cat(tabix,"\n",file=stderr())
      system(sprintf("%s '%s'", shell, tabix))
      wg <- read.table(sprintf("%s/%s",tmp.dir,wg.fn), header=FALSE, sep="\t", as.is=TRUE, check.names=FALSE)
      unlink(wg.fn)
      colnames(wg) <- genoWin.header
      wg$cn <- as.factor(wg$cn)
      wg$pos <- wg$start + (wg$end - wg$start)/2
      cat("wg\n")
      print(head(wg))
      save(wg,file=sprintf('%s/wg.Rdata',tmp.dir))
    }
    return(wg)
  })
  
  frags <- reactive({
    # load fragment data
    frag.fn <- basename(tempfile("cnv"))
    if (input$seg.end-input$seg.start > too.big) {
      cat(sprintf("Skipping profile - region too big %s > %s", input$seg.end-input$seg.start, too.big))
      wg <- data.frame()
    } else {
      tabix <- sprintf("tabix -h %s %s:%s-%s > /tmp/%s", profile.fn, input$seg.chr, input$seg.start-input$pad, input$seg.end+input$pad, frag.fn)
      cat(tabix,"\n",file=stderr())
      system(sprintf("%s '%s'", shell, tabix))
      frags <- read.table(sprintf("%s/%s",tmp.dir,frag.fn), header=FALSE, sep="\t", as.is=TRUE, check.names=FALSE)
      colnames(frags) <- frags.header
    }
    return(frags)
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
  
  # returns a data.frame for plotting that is a row per sample and expected counts
  expected <- reactive({
    frags <- frags()
    sample.count <- ncol(frags) - 6
    exp <- foreach(i=1:sample.count) %do% {
      obs_exp <- strsplit(frags[,6+i], ',')
      as.numeric(unlist(lapply(obs_exp, function(x) { x[2] })))
    }
    names(exp) <- colnames(frags[1,7:(6+sample.count)])
    exp.df <- as.data.frame(exp)
    colnames(exp.df) <- names(exp)
    exp.df$pos <- as.integer(1:nrow(exp.df))
    exp2 <- gather(exp.df, sample, exp, -pos)
    exp2$start.map <- frags[exp2$pos,'START']
    exp2$sample.type <- factor(ifelse(exp2$sample %in% input$seg.sample, 'target','other'))
    return(exp2)
  })
  
  # observed/expected ratio
  oer <- reactive({
    frags <- frags()
    sample.count <- ncol(frags) - 6
    # returns a list of ratios for each window
    obs_exp_ratio <- foreach(i=1:sample.count) %do% {
      # each column is comma-separated observed,expected
      obs_exp <- strsplit(frags[,6+i], ',')
      obs <- as.integer(unlist(lapply(obs_exp, function(x) { x[1] })))
      exp <- as.numeric(unlist(lapply(obs_exp, function(x) { x[2] })))
      
      # window obs and exp
      obs_mean <- rollapply(obs, width=input$win.size, sum)
      exp_mean <- rollapply(exp, width=input$win.size, sum)
      obs_mean / exp_mean
    }
    
    names(obs_exp_ratio) <- colnames(frags[1,7:(6+sample.count)])
    obs_exp_ratio.df <- as.data.frame(obs_exp_ratio)
    colnames(obs_exp_ratio.df) <- names(obs_exp_ratio)
    
    obs_exp_ratio.df$site.median = apply(obs_exp_ratio.df, 1, median, na.rm=TRUE)
    obs_exp_ratio.df$pos <- as.integer(1:nrow(obs_exp_ratio.df))
    
    oer <- gather(obs_exp_ratio.df, sample, ratio, -pos)
    oer$left <- frags[oer$pos,'START']
    oer$right <- frags[oer$pos + input$win.size - 1,'END']
    oer$start.map <- frags[oer$pos,'START'] + (frags[oer$pos + input$win.size - 1,'END'] - frags[oer$pos,'START'])/2
    oer$sample.type <- ifelse(oer$sample == 'site.median','median', ifelse(oer$sample %in% input$seg.sample, 'target','other'))
    if (all(oer$sample.type != 'target')) {
      cat("Can't find any fragments with selected samples: \n"); print(input$seg.sample)
    }
    save(oer,file=sprintf('%s/oer.Rdata',tmp.dir))
    return(oer)  
  })
  
  xbounds <- reactive({
#    xlim(min(oer()$start.map,cnv.disp()$start.map),max(oer()$start.map,cnv.disp()$end.map))
    xmin=min(oer()$start.map)
    xmax=max(oer()$start.map)
    cat(sprintf("xmin=%s xmax=%s\n",xmin,xmax))
    coord_cartesian(xlim=c(xmin,xmax))
  })
  ylims <- reactive({
    ylim(min(oer()$start.map), max(oer()$start.map))
  })
  flip.bounds <- function(coord) {
    coord_cartesian(xlim=coord$limits$y, ylim=coord$limits$x)
  }
  
  cnvPlot <- reactive({
    cat("cnvPlot...\n")
    
    cnv.disp <- cnv.disp()
    if (input$show_target_only && nrow(cnv.disp)>0) { cnv.disp <- subset(cnv.disp, target==TRUE | evidence=='Array') }
    cat("cnv.disp\n")
    print(head(cnv.disp))
    save(cnv.disp,file=sprintf('%s/cnv.disp.Rdata',tmp.dir))
    
    if (!input$show_readPairs) {
      plt <- ggplot(cnv.disp, aes(x=start.map, xend=end.map, y=.id, yend=.id, color=cn.disp, size=target)) + cn.colors 
    } else {
      plt <- ggplot(cnv.disp, aes(x=start.map, xend=end.map, y=.id, yend=.id, color=paired.reads, size=target)) + scale_color_gradient2() 
    }# NOT USED else {
    #  plt <- ggplot(cnv.disp, aes(x=start.map, xend=end.map, y=.id, yend=.id, color=Phred, size=target)) + scale_color_gradient(limits=c(0,100))
    #}
    if (nrow(subset(cnv.disp,evidence=='Array'))>0) {
      plt <- plt + geom_point(data=subset(cnv.disp, evidence=='Array'), size=5, color='black')
    }
    return(plt + geom_segment() + facet_grid(evidence ~ ., scales="free", space="free") + guides(size="none") + theme(axis.title.x = element_blank(), axis.title.y=element_blank()) + xbounds())
  })
  
  cnvCIPlot <- reactive({
    cat("cnvCIPlot...\n")  
    cn.segs.merged.fn <- paste0(data.dir,"/sites_cnv_segs.txt.smlcsm.Rdata")
    load(cn.segs.merged.fn)
    x <- filter(cn.segs.merged, 
                chr == input$seg.chr & 
                  end.map > input$seg.start-input$pad & 
                  start.map < input$seg.end+input$pad & 
                  (cn != 2 | input$show_wildtype) & 
                  len > input$min.cnv.len & 
                  !is.na(cn))
    if (input$show_target_only) { x <- filter(x, .id %in% input$seg.sample) }
    x$target <- x$.id %in% input$seg.sample
    print(x)
    ggplot(x, aes(x=.id, y=start.map, ymin=start.map.L, ymax=start.map.R, color=target))+geom_pointrange(size=1.5) + 
      geom_pointrange(aes(y=end.map, ymin=end.map.L, ymax=end.map.R), size=1.5) +
      geom_segment(aes(y=start.map, yend=end.map, xend=.id), size=1.5, linetype=3) +
      theme(axis.title.x = element_blank(), axis.title.y=element_blank(), plot.margin=unit(c(0,1,0,0),'inches')) + guides(color='none') + ylims() + coord_flip()
  })
  
  winGenoPlot <- reactive({
    cat("winGenoPlot...\n")
    # if (length(gs.dels()$.id) > 0) {
    #   sample.set <- intersect(gs.dels()$.id,input$seg.sample)
    # } else {
      sample.set <- input$seg.sample
    # }

    df <- subset(winGeno(), sample %in% sample.set)
    df <- df[sample(nrow(df)),]  # randomly permute so that overlay is not always left to right
    
    if (input$show_winGeno_2) {
      aes.val <- aes(x=start, xend=end, group=sample, color=cnq, y=cn, yend=cn, alpha=0.25, size=5)
      shape.val <- geom_segment(position=position_jitter(height=0.4)) 
    } else {
      aes.val <- aes(x=pos, group=sample, color=cnq, y=cn)
      shape.val <- geom_point(size=5, shape=15)
    }
    return(ggplot(df) + aes.val + shape.val + facet_grid(sample ~ ., scale="free") + guides(size="none",alpha="none") + scale_color_gradient(limits=c(0,100)) + 
             xbounds() + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,0,0,as.numeric(input$winGeno_nudge_L)),'in')))
  })
  
  fragPlot <- reactive({
    cat("fragPlot...\n")
    df <- oer()
    ggplot(df, aes(x=start.map, y=ratio, group=sample)) + geom_step(color='grey') +
      geom_step(data=df[df$sample.type=='median',], color='black') + 
      geom_step(data=df[df$sample.type=='target',], color='red') +
      theme(axis.title.x = element_blank(), plot.margin=unit(c(0,1,0,0.5),'inches')) +
      xbounds()
  })
  
  bkptPriors <- reactive({
    if (input$seg.end-input$seg.start > too.big) {
      cat(sprintf("Skipping profile - region too big %s > %s", input$seg.end-input$seg.start, too.big))
      return(data.frame())
    } else {
      bkpts <- dbGetQuery(db$con, sprintf("SELECT b.sample, ps.start_pos, ps.end_pos, loss_ll, gain_ll, no_bkpt_ll FROM bkpt b, profile_segment ps WHERE ps.chrom ='%s' AND ps.start_pos < %s AND ps.start_pos > %s AND b.bkpt_bin = ps.bin", 
                                          input$seg.chr, input$seg.end+input$pad, input$seg.start-input$pad))
      
      if (input$bkpt_prior_samples == 'Selected') {
        bkpts <- filter(bkpts, sample %in% input$seg.sample)
      } else if (input$bkpt_prior_samples == 'Exclude Selected') {
        bkpts <- filter(bkpts, ! sample %in% input$seg.sample)
      }
      
      if (input$bkpt_prior_CN != 'Any') {
        if (input$bkpt_prior_CN == '0 or 1') {
          CN.match = c(0,1)
        } else {
          CN.match = as.integer(input$bkpt_prior_CN)
        }
        
        gdo <- read.table(gs_dels.fn, header=TRUE, sep="\t", stringsAsFactors = FALSE)
        colnames(gdo) <- c('.id','seg','chr','start.map','end.map','cn','cq','paired.reads')
        x <- subset(gdo, chr==input$seg.chr & end.map > input$seg.start-input$pad & start.map < input$seg.end+input$pad & cn %in% CN.match)
        bkpts <- filter(bkpts, sample %in% x$.id)
      }
      
      cat(sprintf("bkpts: %s rows\n",nrow(bkpts)))
      cat(sprintf("samples: %s\n",length(unique(bkpts$sample))))
      
      # first normalize, then compute NOT probs
      bkpts <- mutate(bkpts,
                      loss = 10^-loss_ll,
                      gain = 10^-gain_ll,
                      nc = 10^-no_bkpt_ll,
                      tot=loss+gain+nc,
                      lossZ = loss/tot,
                      gainZ = gain/tot,
                      ncZ = nc/tot,
                      Lloss = -log10(lossZ),
                      Lgain = -log10(gainZ),
                      Lnc = -log10(ncZ),
                      Lno_loss=-log10(1-lossZ),
                      Lno_gain=-log10(1-gainZ)
      )
      
      bkpt.prior <- mutate(ddply(bkpts, .(start_pos), summarize,
                                 Lall_no_loss=sum(Lno_loss),
                                 Lall_no_gain=sum(Lno_gain),
                                 Lall_nc=sum(Lnc),
                                 some_loss=1-10^(-Lall_no_loss), 
                                 some_gain=1-10^(-Lall_no_gain),
                                 some_nc=1-10^(-Lall_nc),
                                 Lsome_loss=-log10(some_loss), 
                                 Lsome_gain=-log10(some_gain),
                                 Lsome_nc=-log10(some_nc)))
      
      save(bkpt.prior,file=sprintf('%s/bkpt_prior.Rdata',tmp.dir))
      return(bkpt.prior)
    }
  })
  
  bkptPlot <- reactive({
    cat("bkptPlot...\n")
    bkpt.prior <- bkptPriors()
    
    bkpt.prior.m <- mutate(melt(bkpt.prior, 'start_pos', c('Lall_no_loss','Lall_no_gain')),
                           kind=ifelse(grepl('_ratio', variable), 'Log Ratio', ifelse(grepl('_ll',variable),'Broken Log Likelihood','Log No Loss/Gain')))

    ggplot(bkpt.prior.m, 
           aes(x=start_pos, y=value, color=variable))+geom_point()+geom_line() + facet_grid(kind~.,scales='free_y')+ theme(axis.title.x = element_blank(),plot.margin=unit(c(0,.75,0,.6),'in'),legend.position="bottom") + ylab("Log Prob") + xbounds()  
  })
  
  eachBkptPlot <- reactive({
    cat("eachBkptPlot...\n")
    bkpts <- bkpts.ll()
    if (nrow(bkpts)>0) {
      
      bkpts <- melt(bkpts, c('sample','start_pos','end_pos'),c('loss_ll','gain_ll','Lno_lossZ','Lno_gainZ')) 
      bkpts <- mutate(bkpts, kind=ifelse(grepl('Z',variable),'Z',ifelse(grepl('no_',variable),'No','')))
      ggplot(bkpts, aes(x=start_pos, y=value, color=variable)) + geom_point() + geom_line() + facet_grid(sample+kind~.,scales="free_y") + theme(axis.title.x = element_blank(),plot.margin=unit(c(0,0.75,0,.6),'in'),legend.position="bottom") + ylab("Log Prob") + xbounds()
    } else {
      ggplot() + geom_blank()
    }
  })
  
  # merge frag and geno segments and color fragment by geno label
  # requires the windows to be the same size
  genoFragPlot <- reactive({
    oer.df <- oer()
    wg.df <- winGeno()
    
    # fix to use common column names
    oer.df <- rename(oer.df, start=left, end=right)
    
    oer.wg <- left_join(subset(oer.df, sample.type != 'median'), wg.df, by=c('start','end','sample'))
    oer.wg$cn <- as.numeric(as.character(oer.wg$cn))
    oer.wg$cn <- ifelse(is.na(oer.wg$cn),99, oer.wg$cn)
    oer.wg$sample.type.size <- ifelse(oer.wg$sample.type=='target',2,1) # this doesn't work! I want 1px or 2px
    oer.wg$cn.disp <- factor(ifelse(oer.wg$cn==99, oer.wg$cn, ifelse(oer.wg$cn > 5, 5, oer.wg$cn)), levels=c('0','1','2','3','4','5','99'), labels=c('0','1','2','3','4','5+','Disc'))
    
    if (input$show_target_frag_only) {
      oer.wg <- filter(oer.wg, sample %in% input$seg.sample)
    }
    
    ggplot(oer.wg, aes(x=start.map, y=ratio, yend=ratio, color=cn.disp, group=sample)) +geom_step() + cn.colors + theme(axis.title.x = element_blank(), plot.margin=unit(c(0,1,0,0.5),'inches')) +
      xbounds() + guides(color=FALSE)
  })

  output$expPlot <- renderPlot({
    expected <- expected()
    colors <- c(target='#FF0000',other='#222222')
    alpha <- c(target=1,other=0.2)
    ggplot(expected, aes(x=start.map,y=exp,group=sample,color=sample.type,alpha=sample.type)) + geom_line() +
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
    if (input$show_winGeno) { plots$winGeno = winGenoPlot() }
    if (input$show_prior) { plots$prior = bkptPlot() }
    if (input$show_each_bkpt) { plots$bkpt = eachBkptPlot() }
    grid.arrange(do.call(arrangeGrob,c(plots,list(ncol=1))))
  })
  
  output$cnqs <- renderPlot({
    wg <- subset(winGeno(), sample %in% input$seg.sample)
    gs.dels <- subset(gs.dels(), .id %in% input$seg.sample)
    plots <- lapply(input$seg.sample, function(this.sample) {
      p <- ggplot(subset(wg, sample==this.sample), aes(x=pos,y=cnq, color=cn)) + geom_point() + ggtitle(this.sample)
      gs.dels.id <- subset(gs.dels, .id==this.sample) 
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
                       bayes_loss=loss*some_loss,
                       bayes_gain=gain*some_gain,
                       bayes_nc=nc*some_nc,
                       bayes_tot=bayes_loss+bayes_gain+bayes_nc,
                       bayes_lossZ=bayes_loss/bayes_tot,
                       bayes_gainZ=bayes_gain/bayes_tot,
                       bayes_no_loss=1-bayes_loss,
                       bayes_no_gain=1-bayes_gain,
                       bayes_no_loss2=(1-loss)*(1-some_loss),
                       bayes_no_gain2=(1-gain)*(1-some_gain),
                       bayes_no_lossZ=1-bayes_lossZ,
                       bayes_no_gainZ=1-bayes_gainZ,
                       Lbayes_loss=-log10(bayes_loss),
                       Lbayes_gain=-log10(bayes_gain),
                       Lbayes_lossZ=-log10(bayes_lossZ),
                       Lbayes_gainZ=-log10(bayes_gainZ),
                       Lbayes_no_loss=-log10(bayes_no_loss),
                       Lbayes_no_gain=-log10(bayes_no_gain),
                       Lbayes_no_loss2=-log10(bayes_no_loss2),
                       Lbayes_no_gain2=-log10(bayes_no_gain2),
                       loss.ratio=log10(loss/(no_bkpt+gain)),
                       gain.ratio=log10(gain/(no_bkpt+loss)),
                       loss.ratioZ=log10(lossZ/(no_bkptZ+gainZ)),
                       gain.ratioZ=log10(gainZ/(no_bkptZ+lossZ)),
                       bayes.loss.ratio=log10(bayes_loss/bayes_no_loss),
                       bayes.gain.ratio=log10(bayes_gain/bayes_no_gain),
                       bayes.loss.ratioZ=log10(bayes_lossZ/bayes_no_lossZ),
                       bayes.gain.ratioZ=log10(bayes_gainZ/bayes_no_gainZ))
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
