library(plyr)
library(dplyr)
library(ggplot2)

db <- src_postgres(dbname='seq', host='localhost', port=5432, user='postgres')

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
  p <- ggplot(na.omit(df), aes(x=overlap, y=frac)) + geom_point() +
    facet_grid(~facet) +
    ggtitle(title) + 
    ylab(ylabel) + xlab(xlabel)
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
}

snsp.comp(greater.or.equal,'or Greater', thresh.intvls[2:length(thresh.intvls)])
snsp.comp(less.or.equal,'or Lesser', thresh.intvls[1:(length(thresh.intvls)-1)])
snsp.comp(intvl.bin, 'Bin', thresh.intvls[2:(length(thresh.intvls)-2)], with.lines=FALSE)



ovlp.srt <- ddply(arrange(subset(ovlp, !is.na(k_geno)), known_len), .(k_geno),
                  mutate,
                  tot_known_len=cumsum(known_len),
                  tot_overlap_len=cumsum(overlap_len),
                  sn=tot_overlap_len/tot_known_len)

sn.geno <- ddply(ovlp.srt, .(k_geno), summarize, sn=sum(overlap_len)/sum(known_len))




ggplot(ovlp.srt, aes(x=known_len, y=sn)) + geom_point() + 
  xlab("Length of Known DEL") + facet_grid(.~k_geno) +
  ggtitle("Base-Level Sensitivity\nup to Length") + # xlim(0,5000) +
  geom_hline(aes(yintercept=sn), sn.geno) + geom_label(aes(label=sprintf("%.0f%%",sn*100), y=sn), sn.geno, x=max(ovlp.srt$known_len)/2)


ovlp2.srt <- ddply(arrange(subset(ovlp2, !is.na(p_geno)), pred_len), .(p_geno),
                   mutate,
                   tot_pred_len=cumsum(pred_len),
                   tot_overlap_len=cumsum(overlap_len),
                   sp=tot_overlap_len/tot_pred_len)

sp.geno <- ddply(ovlp2.srt, .(p_geno), summarize, sp=sum(overlap_len)/sum(pred_len))
ggplot(ovlp2.srt, aes(x=pred_len, y=sp)) + geom_point() + xlab("Length of predicted DEL") + 
  ggtitle("Base-Level Specificity\nup to Length") + facet_grid(.~p_geno) +
geom_hline(aes(yintercept=sp), sp.geno) + geom_label(aes(label=sprintf("%.0f%%",sp*100), y=sp), sp.geno, x=max(ovlp2.srt$pred_len)/2)



ovlp3 <- db %>% tbl('overlap') %>% filter(k_kind=='SAMPLESEG' && p_kind=='SAMPLESEG') %>% collect
ggplot(ovlp3, aes(x=start_offset, y=end_offset,size=known_len)) + geom_point(alpha=0.20) + xlab('Offset from True Start') + ylab('Offset from True End') + ggtitle("Boundary Offsets\nPredictions that Overlap Conserved GStrip Deletions") + geom_vline(xintercept=0) + geom_hline(yintercept=0) + xlim(-1000,1000) + ylim(-1000,1000) # + geom_density_2d()

ovlp4 <- db %>% tbl('overlap') %>% filter(k_kind=='SAMPLESEGSP' && p_kind=='SAMPLESEGSP') %>% collect
ggplot(ovlp4, aes(x=start_offset, y=end_offset,size=known_len)) + geom_point(alpha=0.20) + xlab('Offset from True Start') + ylab('Offset from True End') + ggtitle("Boundary Offsets\nPredictions that Overlap Any GStrip CNV") + geom_vline(xintercept=0) + geom_hline(yintercept=0) + xlim(-10000,10000) + ylim(-10000,10000) # + geom_density_2d()

