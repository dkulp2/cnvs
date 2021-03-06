library(plyr)
library(dplyr)
library(ggplot2)

db.conn.params <- as.list(unlist(strsplit(Sys.getenv('DBCONN'),":")))
names(db.conn.params) <- c('user','host','port','dbname')
db <- do.call(src_postgres, db.conn.params)

dbSendQuery(db$con, "SET search_path TO eval")

thresh.intvls <- c(0,0.00001,seq(0.1,0.9,0.1),0.99999,1)
less.or.equal <- function(a,b) { a <= b }
greater.or.equal <- function(a,b) { a >= b }
intvl.bin <- function(a,b) {
  idx <- which(thresh.intvls==b)
  if (idx < length(thresh.intvls)) { idx <- idx + 1 }
  b2 <- thresh.intvls[idx]
  return(a >= b & a <= b2)
}

snsp.stats <- function(kind, geno.expr, thresh.expr, id.fld, comp.func, filter.expr, intvls) {
  ovlp <- db %>% tbl('overlap') %>% filter((k_kind==kind || p_kind==kind)) %>% filter(filter.expr) %>% collect
  res <- ddply(ovlp, geno.expr, function(df) { 
    ldply(intvls, function(thresh) {
      mutate(with(df, data.frame(overlap=thresh,
                                 success=length(unique(df[[id.fld]][comp.func(eval(thresh.expr),thresh)])),
                                 total=(length(unique(df[[id.fld]]))))),
             frac=success/total)
    })
  })
  return(list(ovlp=ovlp, res=res))
}

snsp.plot <- function(df, title, xlabel, ylabel, with.lines=TRUE) {
  y50 <- filter(df, overlap==0.5) %>% select(frac, facet) 
  p <- ggplot(na.omit(df), aes(x=overlap, y=frac)) + geom_point() +
    facet_grid(~facet) +
    ggtitle(title) + 
    ylab(ylabel) + xlab(xlabel) + ylim(0,1) + geom_label(aes(label=sprintf("%.0f%%",frac*100), y=frac+0.05), y50, x=0.5)
  if (with.lines) { p <- p + geom_line() }
  return(p)
}

snsp.comp <- function(comp.func, comp.func.name, intvls, with.lines=TRUE) {
  ret <- snsp.stats('SITE', quote(k_geno), quote(overlap_len/known_len), 'k_id', comp.func,
                   quote(!is.na(k_geno) & k_geno %in% c('1','0') & known_len > 1200), intvls)
  ret$res$facet <- sprintf("%s (n=%s)", ret$res$k_geno, ret$res$total)
  print(snsp.plot(ret$res, 
                  sprintf("Site Sensitivity:\nSites Found / Known Sites (CN={0,1}, >1200 nt)\n at Fraction Overlap %s", comp.func.name), 
                  'Fraction of known DEL overlapped by prediction', 'Sn', with.lines))
  
  ret <- snsp.stats('SITESP', quote(p_geno), quote(overlap_len/pred_len), 'p_id', comp.func,
                   quote(!is.na(p_geno) & p_geno %in% c('1','0') & pred_len > 1200), intvls)
  ret$res$facet <- sprintf("%s (n=%s)", ret$res$p_geno, ret$res$total)
  print(snsp.plot(ret$res, 
                  sprintf("Site Specificity:\nTrue Sites / Predicted Sites (CN={0,1}, >1200 nt)\n at Fraction Overlap %s", comp.func.name), 
                  'Fraction of predicted DEL overlapped by known', 'Sp', with.lines))
  
  
  ret <- snsp.stats('SAMPLESEG',  
                    .(k_geno,k_id), quote(overlap_len/known_len), 'k_sample', comp.func,
                    quote(!is.na(k_geno) & k_geno %in% c('1','0') & known_len > 1200), intvls)

  x <- ret$res
  x <- ddply(x, .(k_geno, overlap), mutate, rank=rank(-frac, ties.method="first"))
  print(ggplot(filter(x, overlap==0.5), aes(x=rank, y=frac, color=k_geno)) + geom_line() + geom_point() + facet_grid(~k_geno, scales="free_x") + ylab("Sn") + xlab("Rank of Known Sites") + ggtitle("Ranked SampleSeg Sensitivity\n(at 50% overlap)"))

  x2 <- ldply(c(0.0001,0.1,0.25,0.5,0.75,0.9,1), function(sample.frac) {
    ddply(x, .(k_geno, overlap), summarize, sample.frac, successes=sum(success/total>=sample.frac), sites=length(success), frac=successes/sites)
  })
  x2$samples.chr <- ifelse(x2$sample.frac < .1, '>0%', sprintf("%.0f%%", x2$sample.frac*100))
  x2$samples <- factor(x2$samples.chr, ordered = TRUE, levels=unique(x2$samples.chr))
  x2$facet <- sprintf("%s (n=%s)", x2$k_geno, x2$sites)
  print(ggplot(na.omit(x2), aes(x=overlap, y=frac, color=samples)) + geom_line() + facet_grid(.~facet) + ylab("Sn") + ggtitle("Site Sn given required fraction of samples per site") + ylim(0,1))

  x$Frequency <- cut(x$total, c(0,1,2,3,6,15,40,100))
  x2 <- ddply(x, .(Frequency,overlap), summarize, success=sum(success), total=sum(total), frac=success/total)
  print(ggplot(x2, aes(x=overlap, y=frac, color=Frequency)) + geom_line() + ylab("Sn") + ggtitle("Site Sensitivity Partitioned by Frequency"))

  print(ggplot(filter(x2, overlap==0.5), aes(x=Frequency, y=frac))+geom_point() + ylab("Sn") + ggtitle("Site Sensitivity By Site Frequency (at 50% overlap)"))
    
  res <- ddply(ret$res, .(k_geno, overlap), summarize, success=sum(success), total=sum(total), frac=success/total)
  res$facet <- sprintf("%s (n=%s)", res$k_geno, res$total)
  print(snsp.plot(res, 
                  sprintf("SampleSeg Sensitivity:\nSampleSeg Found / Known SampleSegs (CN={0,1}, >1200 nt)\n at Fraction Overlap %s", comp.func.name), 
                  'Fraction of known SampleSeg overlapped by prediction', 'Sn', with.lines))
  
  
  ret <- snsp.stats('SAMPLESEGSP',  
                    .(p_geno,p_id), quote(overlap_len/pred_len), 'p_sample', comp.func,
                    quote(!is.na(p_geno) & p_geno %in% c('1','0') & pred_len > 1200), intvls)
  res <- ddply(ret$res, .(p_geno, overlap), summarize, success=sum(success), total=sum(total), frac=success/total)
  res$facet <- sprintf("%s (n=%s)", res$p_geno, res$total)
  print(snsp.plot(res, 
                  sprintf("SampleSeg Specificity:\nTrue SampleSeg / Predicted SampleSegs (CN={0,1}, >1200 nt)\n at Fraction Overlap %s", comp.func.name),
                  'Fraction of predicted SampleSeg overlapped by known', 'Sp', with.lines))

  x <- ret$res
  x <- ddply(x, .(p_geno, overlap), mutate, rank=rank(-frac, ties.method = 'first'))
  print(ggplot(filter(x, overlap==0.5), aes(x=rank, y=frac, color=p_geno)) + geom_line() + geom_point() + facet_grid(~p_geno, scales="free_x") + ylab("Sp") + xlab("Rank of Predicted Sites") + ggtitle("Ranked SampleSeg Specificity\n(at 50% overlap)"))
  
}

snsp.comp(greater.or.equal,'or Greater', thresh.intvls[2:length(thresh.intvls)])
#snsp.comp(less.or.equal,'or Lesser', thresh.intvls[1:(length(thresh.intvls)-1)])
#snsp.comp(intvl.bin, 'Bin', thresh.intvls[2:(length(thresh.intvls)-2)], with.lines=FALSE)


# BASE STATS - stratify by length
ovlp.sn <- db %>% tbl('overlap') %>% filter((k_kind=='SAMPLESEG' || p_kind=='SAMPLESEG')) %>% filter(!is.na(k_geno) & k_geno %in% c('1','0')) %>% collect

ovlp.sn.srt <- ddply(arrange(ovlp.sn, known_len), .(k_geno),
                     mutate,
                     tot_known_len=cumsum(known_len),
                     tot_overlap_len=cumsum(overlap_len),
                     sn=tot_overlap_len/tot_known_len)
base.sn <- ddply(ovlp.sn.srt, .(k_geno), summarize, sn=sum(overlap_len)/sum(known_len))

ggplot(ovlp.sn.srt, aes(x=known_len, y=sn)) + geom_point() + 
  xlab("Length of Known DEL") + facet_grid(.~k_geno) + ylim(0,1) +
  ggtitle("Base-Level Sensitivity\nup to Length") + xlim(0,10000) +
  geom_hline(aes(yintercept=sn), base.sn) + geom_label(aes(label=sprintf("%.0f%%",sn*100), y=sn), base.sn, x=5000) # max(ovlp.sn.srt$known_len)/2)

ovlp.sp <- db %>% tbl('overlap') %>% filter((k_kind=='SAMPLESEG' || p_kind=='SAMPLESEG')) %>% filter(!is.na(p_geno) & p_geno %in% c('1','0')) %>% collect

ovlp.sp.srt <- ddply(arrange(ovlp.sp, pred_len), .(p_geno),
                   mutate,
                   tot_pred_len=cumsum(pred_len),
                   tot_overlap_len=cumsum(overlap_len),
                   sp=tot_overlap_len/tot_pred_len)
base.sp <- ddply(ovlp.sp.srt, .(p_geno), summarize, sp=sum(overlap_len)/sum(pred_len))
ggplot(ovlp.sp.srt, aes(x=pred_len, y=sp)) + geom_point() + xlab("Length of predicted DEL") + 
  ggtitle("Base-Level Specificity\nup to Length") + facet_grid(.~p_geno) + ylim(0,1) +
geom_hline(aes(yintercept=sp), base.sp) + geom_label(aes(label=sprintf("%.0f%%",sp*100), y=sp), base.sp, x=max(ovlp.sp.srt$pred_len)/2)


# Edge stats - only for knowns that overlap predictions
ovlp3 <- db %>% tbl('overlap') %>% filter(k_kind=='SAMPLESEG' && p_kind=='SAMPLESEG') %>% collect
ggplot(ovlp3, aes(x=start_offset, y=end_offset,color=known_len)) + geom_point(alpha=0.20) + xlab('Offset from True Start') + ylab('Offset from True End') + ggtitle("Boundary Offsets\nPredictions that Overlap Conserved GStrip Deletions") + geom_vline(xintercept=0) + geom_hline(yintercept=0) + xlim(-1000,1000) + ylim(-1000,1000) + geom_density_2d(color='gray40') + guides(color=FALSE)

# What fraction of knowns are within the 95% confidence interval?
cil <-
ldply(c(0.1, 0.25, 0.5, 0.75, 0.9), function(ov) {
  k.ss <- db %>% tbl('overlap') %>% filter(k_kind=='SAMPLESEG' & as.numeric(overlap_len) / known_len > ov)
  mutate(na.omit(ddply(collect(k.ss), .(start_in_cil, end_in_cil), summarize, in_cil=length(k_sample))), in_cil_frac=in_cil/sum(in_cil), overlap=ov)
  
})
cil$inCI <- with(cil, ifelse(start_in_cil, ifelse(end_in_cil, 'Both', 'Only Start'), ifelse(end_in_cil, 'Only End', 'Neither')))
ggplot(cil, aes(x=overlap, y=in_cil_frac, color=inCI)) + geom_line() + ylab('Fraction of known DELs') + xlab("Minimum overlap")
print(cil)

k.ss <- db %>% tbl('overlap') %>% filter(k_kind=='SAMPLESEG' & as.numeric(overlap_len) / known_len > 0.75) %>% collect
k.ss$inCI <- with(k.ss, ifelse(start_in_cil, ifelse(end_in_cil, 'Both', 'Only Start'), ifelse(end_in_cil, 'Only End', 'Neither')))
ggplot(na.omit(k.ss), aes(x=start_offset, y=end_offset, color=inCI)) + geom_point() + facet_wrap(~inCI) + xlab('Offset from True Start') + ylab('Offset from True End') + ggtitle('Boundary Offsets and Confidence Intervals\nfor predictions with 75% overlap to knowns') + lims(x=c(-1000,1000),y=c(-1000,1000)) + geom_hline(yintercept=0) + geom_vline(xintercept = 0)
ggplot(filter(na.omit(k.ss),inCI!='Both'), aes(y=-start_offset, x=-end_offset, color=inCI)) + geom_point() + xlab('Offset from True Start') + ylab('Offset from True End') + ggtitle('Boundary Offsets and Confidence Intervals\nfor predictions with 75% overlap to knowns') + lims(x=c(-1000,1000),y=c(-1000,1000)) + geom_hline(yintercept=0) + geom_vline(xintercept = 0)


# Group the edge stats by site
k.sites <- db %>% tbl('all_known_cnvs') %>% filter(kind=='SITE' & source=='GS_DELS') %>% collect

# Split out known sites with a separate row per ID
k.sites.ids <- ddply(k.sites, .(id), function(df) return(data.frame(geno=df$geno, ids=unique(unlist(strsplit(df$id, ';'))),stringsAsFactors = FALSE)))
k.sites2 <- inner_join(k.sites, k.sites.ids, by=c('id','geno'))

# Use a single ID for prediction
ovlp4 <- ddply(ovlp3, .(k_id), mutate, k_id.first=first(unlist(strsplit(k_id, ';'))))  # choose first of N DEL IDs to join with ids in site

# one row for each site and sampleseg pair 
site.join.ovlp <- mutate(inner_join(k.sites2, ovlp4, by=c("ids"="k_id.first","geno"="k_geno")),
                         start.offset=start_pos.y-start_pos.x,
                         end.offset=end_pos.y-end_pos.x)

# for each site, compute stats on offset
site.offset.stats <- ddply(site.join.ovlp, .(id), summarize, 
                           start.offset.med=median(start.offset), 
                           start.offset.100=sum(abs(start.offset-start.offset.med)<=100),
                           start.offset.200=sum(abs(start.offset-start.offset.med)<=200),
                           n=length(id),
                           so.100.frac=start.offset.100/n,
                           so.200.frac=start.offset.200/n,
                           so.100.out=n-start.offset.100,
                           so.200.out=n-start.offset.200,
                           spread=as.integer(max(start.offset)-min(start.offset)),
                           multi.mode.candidate=spread>300)

hist(site.offset.stats$so.100.out, breaks=seq(0,15),xlab='Samples',ylab='Sites',main=sprintf('Histogram of number of samples\nwith breakpoint further than 1 bin\nfrom median breakpoint\n(site count=%s)',nrow(site.offset.stats)))
hist(site.offset.stats$so.200.out, breaks=seq(0,10),xlab='Samples',ylab='Sites',main=sprintf('Histogram of number of samples\nwith breakpoint further than 2 bin\nfrom median breakpoint\n(site count=%s)',nrow(site.offset.stats)))

plot(site.offset.stats$so.100.frac)
plot(density(site.offset.stats$so.100.frac),main="Density of Fraction of Samples within 1 Bin")

# If the bounds among different samples is far apart, then it's a good candidate for alternative starts
ggplot(site.offset.stats, aes(x=n, y=spread, color=multi.mode.candidate)) + geom_jitter() + ylim(0,5000) + xlab("Allele Frequency out of 100") + ggtitle("Min/Max Spread Among Start Bounds")
