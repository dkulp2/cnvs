library(plyr)
library(dplyr)
library(ggplot2)

db <- src_postgres(dbname='seq', host='localhost', port=5432, user='postgres')

plot.site <- function(kind, title, xlabel, ylabel, geno.expr, tresh.expr, id.fld, facet.expr, filter.expr=expression(TRUE)) {
  ovlp <- db %>% tbl('overlap') %>% filter((k_kind==kind || p_kind==kind)) %>% filter(filter.expr) %>% collect
  site <- ddply(ovlp, geno.expr, function(df) { 
    ldply(seq(0.1,1,0.1), function(thresh) {
      with(df, data.frame(overlap=thresh,
                          frac=length(unique(df[[id.fld]][eval(thresh.expr) >= thresh]))/(length(unique(df[[id.fld]])))))
    })
  })
  ggplot(na.omit(site), aes(x=overlap, y=frac)) + geom_point() + geom_line() + 
    facet_grid(facet.expr) +
    ggtitle(title) + 
    ylab(ylabel) + xlab(xlabel)
  return(ovlp)
}

plot.site('SITE', "Site Sensitivity:\nSites Found / Known Sites (CN={0,1}, >1200 nt)\n at Fraction Overlap or Greater", 
          'Fraction of known DEL overlapped by prediction', 'Sn', 
          quote(k_geno), quote(overlap_len/known_len), 'k_id', .~k_geno,
          quote(k_geno %in% c('1','0') & known_len > 1200))
plot.site('SITESP', "Site Specificity:\nTrue Sites / Predicted Sites (CN={0,1}, >1200 nt)\n at Fraction Overlap or Greater", 
          'Fraction of predicted DEL overlapped by known', 'Sp', 
          quote(p_geno), quote(overlap_len/pred_len), 'p_id', .~p_geno,
          quote(!is.na(p_geno) & p_geno %in% c('1','0') & pred_len > 1200))


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
ggplot(ovlp3, aes(x=start_offset, y=end_offset)) + geom_point(alpha=0.20) + xlab('Offset from True Start') + ylab('Offset from True End') + ggtitle('Boundary Offsets') + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  geom_density_2d() + xlim(-1000,1000) + ylim(-1000,1000)

