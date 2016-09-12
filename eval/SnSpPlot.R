library(dplyr)

db <- src_postgres(dbname='seq', host='localhost', port=5432, user='postgres')
print(src_tbls(db))
sn <- db %>% tbl('base_sensitivity')

ovlp <- db %>% tbl('overlap') %>% collect
site.sn <- ddply(ovlp, .(k_geno), function(df) {
  ldply(seq(0.1,1,0.1), function(thresh) {
    with(df, data.frame(coverage=thresh,
                        frac=length(unique(k_id[overlap_len/known_len >= thresh]))/(length(unique(k_id)))))
  })
})
ggplot(na.omit(site.sn), aes(x=coverage, y=frac)) + geom_point() + geom_line() + 
  facet_grid(.~k_geno) +
  ggtitle("Site Sensitivity:\nSites Found / Known Sites\n at Fraction Coverage or Greater") + 
  ylab('Sn') + xlab('Fraction of known DEL overlapped by prediction')

ovlp.srt <- ddply(arrange(subset(ovlp, !is.na(k_geno)), known_len), .(k_geno),
                  mutate,
                  tot_known_len=cumsum(known_len),
                  tot_overlap_len=cumsum(overlap_len),
                  sn=tot_overlap_len/tot_known_len)

sn.geno <- ddply(ovlp.srt, .(k_geno), summarize, sn=sum(overlap_len)/sum(known_len))




ggplot(ovlp.srt, aes(x=known_len, y=sn)) + geom_point() + 
  xlab("Length of Known DEL") + facet_grid(.~k_geno) +
  ggtitle("Base-Level Sensitivity\nup to Length") + xlim(0,5000) +
  geom_hline(aes(yintercept=sn), sn.geno) + geom_label(aes(label=sprintf("%.0f%%",sn*100), y=sn), sn.geno, x=max(ovlp.srt$known_len)/2)

ovlp2 <- db %>% tbl('overlap2') %>% filter(p_geno %in% c('1','0')) %>% collect
site.sp <- ddply(ovlp2, .(p_geno), function(df) {
  ldply(c(0.000001,seq(0.1,1,0.1)), function(thresh) {
    with(df, data.frame(coverage=thresh, 
                        frac=length(unique(p_id[overlap_len/pred_len >= thresh]))/(length(unique(p_id)))))
  })
})
ggplot(na.omit(site.sp), aes(x=coverage, y=frac)) + geom_point() + geom_line() + 
  facet_grid(.~p_geno) +
  ggtitle("Site Specificity:\nTrue Sites / Predicted Sites\nat Fraction Coverage or Greater") + 
  ylab('Sp') + xlab('Fraction of predicted DEL overlapped by known')

ovlp2.srt <- ddply(arrange(subset(ovlp2, !is.na(p_geno)), pred_len), .(p_geno),
                   mutate,
                   tot_pred_len=cumsum(pred_len),
                   tot_overlap_len=cumsum(overlap_len),
                   sp=tot_overlap_len/tot_pred_len)

sp.geno <- ddply(ovlp2.srt, .(p_geno), summarize, sp=sum(overlap_len)/sum(pred_len))
ggplot(ovlp2.srt, aes(x=pred_len, y=sp)) + geom_point() + xlab("Length of predicted DEL") + 
  ggtitle("Base-Level Specificity\nup to Length") + facet_grid(.~p_geno) +
geom_hline(aes(yintercept=sp), sp.geno) + geom_label(aes(label=sprintf("%.0f%%",sp*100), y=sp), sp.geno, x=max(ovlp2.srt$pred_len)/2)

ovlp3 <- db %>% tbl('overlap3') %>% collect
ovlp3$start.diff <- ovlp3$start_pos - ovlp3$b_start
ovlp3$end.diff <- ovlp3$end_pos - ovlp3$b_end
ovlp3b <- subset(ovlp3, abs(end.diff) < 5000 & abs(start.diff) < 5000)  # removes outliers.
qplot(ovlp3b$start.diff, ovlp3b$end.diff) + xlab('Offset from True Start') + ylab('Offset from True End') + ggtitle('Boundary Offsets') + geom_vline(xintercept=0) + geom_hline(yintercept=0)

