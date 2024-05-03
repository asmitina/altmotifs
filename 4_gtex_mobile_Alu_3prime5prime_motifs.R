setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(GenomicRanges)
library(plotrix)
library(stringr)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg38)

genome = BSgenome.Hsapiens.UCSC.hg38

revCompl = function(seq) {
  return(chartr("ATGC", "TACG", reverse(seq) ))
}

getCGContent = function(seq){
  cg.len = length(grep("C|G", strsplit(seq, "")[[1]]))
  return(round(cg.len/nchar(seq), digits = 4))
}


tr_info_new = read.delim('KGP.STR.All.tsv') # for KGP
tr_info_new = tr_info_new[tr_info_new$chr %in% paste0("chr", c(1:22)), ]

tr_info_new = read.delim('GTEx.tr.regions.samples.tsv', stringsAsFactors = F) # for GTEx
tr_info_new = tr_info_new[tr_info_new$chr %in% paste0("chr", c(1:22)), ]
tr_info_new$varid <- paste(tr_info_new$chr, tr_info_new$start, tr_info_new$end, sep="#")
tr_info_new$repeatID <- paste(tr_info_new$chr, tr_info_new$start, tr_info_new$end, tr_info_new$dominant.motif, sep="#")

repeats = read.delim('~/Downloads/UCSC_hg38_repeatMasker.tsv') # 5655664 
repeats =  repeats[repeats$genoName %in% paste0("chr", c(1:22, "M", "X", "Y")), ] # 5317291 
mobile = repeats[repeats$repFamily == 'Alu' | repeats$repClass == 'LINE',] # 2697298 

ext = 0
repClass = 'SINE'
mobile_ = mobile[mobile$repClass == repClass, ]
mobile_$newStart = ifelse(mobile_$genoStart - ext < 1, 1, mobile_$genoStart - ext)
mobile_$newEnd = mobile_$genoEnd + ext

mobile_ranges = GRanges(mobile_$genoName, IRanges(mobile_$newStart, mobile_$newEnd), mobile_$strand)
tr_info_ranges = GRanges(tr_info_new$chr, IRanges(tr_info_new$start, tr_info_new$end), "*")

overlap = findOverlaps(tr_info_ranges, mobile_ranges)
mobile_test = mobile_[overlap@to, ]# 
tr_test = tr_info_new[overlap@from,]# 

tr_test = cbind(tr_test, mobile_test[,c('genoName', 'repName', 'genoStart', 'genoEnd', 'strand')])

#### Repeats on 5 prime end of Alu elements 
#Alu dist = 0, which are on 5' from Alu - 601 repeat KGP; 940 repeat GTEx
prime5 = rbind(tr_test[tr_test$start <  mobile_test$genoStart & 
                         tr_test$end < mobile_test$genoEnd & 
                         mobile_test$strand == '+',], 
               tr_test[tr_test$start >  mobile_test$genoStart 
                       & tr_test$end > mobile_test$genoEnd & 
                         mobile_test$strand == '-',] )

#Two times repeat motif
prime5$tr_plus_ = ''
prime5$tr_minus_ = ''
for (i in 1:length(prime5$repeatID)){
  tmp.seq = getSeq(genome, prime5$chr[i], start = prime5$start[i], end = prime5$end[i])
  if(grepl(paste0(prime5$dominant.motif[i],prime5$dominant.motif[i]), as.character(reverseComplement(tmp.seq)))){
    prime5$tr_minus_[i] = '_'}
  else{
    if (grepl(paste0(prime5$dominant.motif[i],prime5$dominant.motif[i]), as.character(tmp.seq))){
      prime5$tr_plus_[i] = '+'}
  }
}

#Three times repeat motif
prime5$tr_strand = ''
i = 1
for (i in 1:length(prime5$repeatID)){
  tmp.seq = getSeq(genome, prime5$chr[i], start = prime5$start[i], end = prime5$end[i])
  if(grepl(paste0(prime5$dominant.motif[i],prime5$dominant.motif[i],prime5$dominant.motif[i]), as.character(reverseComplement(tmp.seq)))){
    prime5$tr_strand[i] = '-'}
  else{
    if (grepl(paste0(prime5$dominant.motif[i],prime5$dominant.motif[i],prime5$dominant.motif[i]), as.character(tmp.seq))){
      prime5$tr_strand[i] = '+'}
  }
}

prime5$tr_strand_ = ''
for (i in 1:length(prime5$repeatID)){
  tmp.seq = getSeq(genome, prime5$chr[i], start = prime5$start[i], end = prime5$end[i])
  if(grepl(paste0(prime5$dominant.motif[i],prime5$dominant.motif[i],prime5$dominant.motif[i]), as.character(tmp.seq))){
    prime5$tr_strand_[i] = '+'}
  else{
    if (grepl(paste0(prime5$dominant.motif[i],prime5$dominant.motif[i],prime5$dominant.motif[i]), as.character(reverseComplement(tmp.seq)))){
      prime5$tr_strand_[i] = '-'}
  }
}
table(prime5$tr_strand_)

prime5_motif = prime5[which(prime5$tr_strand == '-' | prime5$tr_strand == '+'),]
prime5_motif = prime5_motif[which(prime5_motif$tr_strand == prime5_motif$tr_strand_),]

prime5_motif$tr_motif = ''
for (i in 1:length(prime5_motif$repeatID)){
  if (prime5_motif$tr_strand[i] != ''){
    if (prime5_motif$strand[i] == '+'){
      prime5_motif$tr_motif[i] = ifelse(prime5_motif$tr_strand[i] == '+', prime5_motif$dominant.motif[i], revCompl(prime5_motif$dominant.motif[i]))}
    else{prime5_motif$tr_motif[i] = ifelse(prime5_motif$tr_strand[i] == '-', prime5_motif$dominant.motif[i], revCompl(prime5_motif$dominant.motif[i]))}
  }
}

prime5.all = as.data.frame( table(prime5_motif$dominant.motif)) 
colnames(prime5.all) = c('motif', 'freq')
prime5.all = prime5.all[order(prime5.all$freq, decreasing = T),]
prime5.all$motif = as.character(prime5.all$motif)
prime5.all$cg.perc = sapply(as.character(prime5.all$motif), getCGContent) * 100
prime5.all$col = color.gradient(prime5.all$cg.perc)

prime5.st = as.data.frame( table(prime5_motif$dominant.motif[prime5_motif$type == 'stable']))# 40; 58 GTEx
colnames(prime5.st) = c('motif', 'freq.st')
prime5.st$motif = as.character(prime5.st$motif)
prime5.st = prime5.st[order(prime5.st$freq, decreasing = T),]

for (i in 1:nrow(prime5.all)){
  if (prime5.all$motif[i] %in% prime5.st$motif){
    prime5.all$freq.st[i] = prime5.st$freq.st[prime5.st$motif == prime5.all$motif[i]]
  }else{
    prime5.all$freq.st[i] = 0
  }
}

prime5.var = as.data.frame( table(prime5_motif$dominant.motif[prime5_motif$type == 'variable']))# 44; 38 GTEx
colnames(prime5.var) = c('motif', 'freq.var')
prime5.var$motif = as.character(prime5.var$motif)
prime5.var = prime5.var[order(prime5.var$freq, decreasing = T),]

for (i in 1:nrow(prime5.all)){
  if (prime5.all$motif[i] %in% prime5.var$motif){
    prime5.all$freq.var[i] = prime5.var$freq.var[prime5.var$motif == prime5.all$motif[i]]
  }else{
    prime5.all$freq.var[i] = 0
  }
}


### PLOT - All - ST - VAR
motif.plot = prime5.all[prime5.all$freq > 1,]

pdf('output/figures/_Alu_5prime_GTEx_all_st_var_COUNT.pdf', w = 3.5, h = 5.5)
par(mar = c(4.5,4,1,1), mgp = c(3,1,0), mfrow=c(3,1)) 
my_bar = barplot(motif.plot$freq, las = 2, main = '', 
                 col = 'grey', ylim = c(0,280), ylab = 'Frequency')
text(my_bar, y = par("usr")[3]-5,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9,
     adj = 1, srt = 35)#, srt = 35
text(my_bar, motif.plot$freq + 10, paste(motif.plot$freq) ,cex=0.9) 

my_bar = barplot(motif.plot$freq.st, las = 2, main = '', # Stable STRs at 5\' of Alu elements',
                 col = alpha('springgreen2', 0.6), ylim = c(0,280), ylab = 'Frequency') 
text(my_bar, y = par("usr")[3]-5,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9,
     adj = 1, srt = 35)#, srt = 35
text(my_bar, motif.plot$freq.st + 10, paste(motif.plot$freq.st) ,cex=0.9) 

my_bar = barplot(motif.plot$freq.var, las = 2, main = '', # Variable STRs at 5\' of Alu elements',
                 col = alpha('firebrick1', 0.6), ylim = c(0,280), ylab = 'Frequency')
text(my_bar, y = par("usr")[3]-5,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9,
     adj = 1, srt = 35)#, srt = 35
text(my_bar, motif.plot$freq.var + 10, paste(motif.plot$freq.var) ,cex=0.9) 
dev.off()


pdf('output/figures/_Alu_5prime_GTEx_all_st_var_FREQ_MIR_lab_SIGN_.pdf', w = 3.5, h = 4.5)
par(mar = c(4.5,4,1,1), mgp = c(3,1,0), mfrow=c(3,1)) 
my_bar = barplot(motif.plot$freq/sum(motif.plot$freq), las = 2, main = '', # All STRs at 5\' of Alu elements',
                 col = 'grey', ylim = c(0,0.5), ylab = 'Frequency')
text(my_bar, y = par("usr")[3]-0.02,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9,
     adj = 1, srt = 35)

par(mar=c(0,5,3,3))
my_bar = barplot(motif.plot$freq.st/sum(motif.plot$freq.st), las = 2, main = '', # Stable STRs at 5\' of Alu elements',
                 col = alpha('springgreen2', 0.6), ylim = c(0,0.5), ylab = 'STABLE') 

par(mar=c(5,5,0,3))
my_bar = barplot(motif.plot$freq.var/sum(motif.plot$freq.var), las = 2, main = '', # Variable STRs at 5\' of Alu elements',
                 col = alpha('firebrick1', 0.6), ylim = c(0.5,0), ylab = 'VARIABLE')
text(my_bar-0.5, y = par("usr")[3]+0.32,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9, pos = 4,
     adj = 1, srt = 90)#, srt = 35
text(my_bar, y = par("usr")[3]+0.32,labels = ifelse(prime5.all$signif, '*', ''),xpd = NA, srt = 90,cex = 1.5,
     adj = 1, srt = 90)#, srt = 35
dev.off()

prime5.all$OR = ''
prime5.all$OR.upper = ''
prime5.all$OR.lower = ''
prime5.all$Fisher.p = ''

### Calc OR between stable and variable
for (i in 1:nrow(prime5.all)){
  test <- fisher.test(data.frame("X1" = c(prime5.all$freq.var[i], prime5.all$freq.st[i]), 
                                 "X2" = c(sum(prime5.all$freq.var) - prime5.all$freq.var[i], 
                                          sum(prime5.all$freq.st) - prime5.all$freq.st[i])))#, alternative = 'greater')
  prime5.all$OR[i] = test$estimate
  prime5.all$OR.upper[i] = test$conf.int[1]
  prime5.all$OR.lower[i] = test$conf.int[2]
  prime5.all$Fisher.p[i] = test$p.value
  # format(test$p.value, scientific = T, digits = 2)# "1.8e-208"
}

prime5.all$signif = ifelse(prime5.all$OR > 1 & prime5.all$Fisher.p < 0.05, T, F)

### BOXPLOT
prime5_motif = read.delim('_GTEx/output/alu_5prime_GTEx_motif_st_var.tsv', sep = ' ')
prime5_motif$dominant.motif = as.character(prime5_motif$dominant.motif)
prime5_motif$gc = sapply(prime5_motif$dominant.motif, getCGContent)

gc.stable = prime5_motif$gc[prime5_motif$type == 'stable']
gc.variable = prime5_motif$gc[prime5_motif$type == 'variable'] 

prime5.alt = as.data.frame(strsplit(paste(prime5_motif$motif[prime5_motif$type == 'variable'], collapse = ';'), ';')[[1]])
colnames(prime5.alt) = 'motif'
prime5.alt$gc = sapply(prime5.alt$motif, getCGContent)

test.st.dom = wilcox.test(gc.stable, gc.variable)
test.st.alt = wilcox.test(gc.stable, prime5.alt$gc)
test.dom.alt = wilcox.test(gc.variable, prime5.alt$gc)

pdf('_GTEx/output/figures/_Alu_5prime_GTEx_boxplot_st_dom_alt__PVAL_COLOR_ALL.pdf', w = 2.6, h = 4)#, h = 7
par(mar = c(4.5,4,1,1), mgp = c(3,1,0), mfrow=c(1,1))
boxplot(gc.stable, gc.var.dom, gc.var.alt, col = c(alpha(c('springgreen2', 'firebrick1', 'dodgerblue1'), 0.6)),
        xaxt = 'n', yaxt = 'n', ylab = 'GC-content')
axis(1, at = c(1,2,3), labels = c('ST', 'DOM', 'ALT'))
axis(2, at = seq(0,1,0.2), las = 2)

segments(1,0.9,2,0.9)
text(1.5,0.95, labels = paste0('p = ', format(test.st.dom$p.value, digits=2)))

segments(1,0.8,3,0.8)
text(2.5,0.85, labels = paste0('p = ', format(test.st.alt$p.value, digits=2)))

segments(2,0.7,3,0.7)
text(2.5,0.75, labels = paste0('p = ', format(test.dom.alt$p.value, digits=2, scientific = T)))

dev.off()




#### Repeats on 3 prime end of Alu elements 
#Alu dist = 0, which are on 3' from Alu - 1026 repeat KGP; 2720 GTEx
prime3 = rbind(tr_test[tr_test$start <  mobile_test$genoStart & 
                         tr_test$end < mobile_test$genoEnd & 
                         mobile_test$strand == '-',],
               tr_test[tr_test$start >  mobile_test$genoStart & 
                         tr_test$end > mobile_test$genoEnd & 
                         mobile_test$strand == '+',])

prime3$tr_strand = ''
for (i in 1:length(prime3$repeatID)){
  tmp.seq = getSeq(genome, prime3$chr[i], start = prime3$start[i], end = prime3$end[i])
  if(grepl(prime3$dominant.motif[i], as.character(tmp.seq))){
    prime3$tr_strand[i] = '+'}
  else{
    if (grepl(prime3$dominant.motif[i], as.character(reverseComplement(tmp.seq)))){
      prime3$tr_strand[i] = '-'}
  }
}


#Two times repeat motif
prime3$tr_plus_ = ''
prime3$tr_minus_ = ''
for (i in 1:length(prime3$repeatID)){
  tmp.seq = getSeq(genome, prime3$chr[i], start = prime3$start[i], end = prime3$end[i])
  if(grepl(paste0(prime3$dominant.motif[i],prime3$dominant.motif[i]), as.character(reverseComplement(tmp.seq)))){
    prime3$tr_minus_[i] = '_'}
  else{
    if (grepl(paste0(prime3$dominant.motif[i],prime3$dominant.motif[i]), as.character(tmp.seq))){
      prime3$tr_plus_[i] = '+'}
  }
}

#Three times repeat motif
prime3$tr_strand = ''
for (i in 1:length(prime3$repeatID)){
  tmp.seq = getSeq(genome, prime3$chr[i], start = prime3$start[i], end = prime3$end[i])
  if(grepl(paste0(prime3$dominant.motif[i],prime3$dominant.motif[i],prime3$dominant.motif[i]), as.character(reverseComplement(tmp.seq)))){
    prime3$tr_strand[i] = '-'}
  else{
    if (grepl(paste0(prime3$dominant.motif[i],prime3$dominant.motif[i],prime3$dominant.motif[i]), as.character(tmp.seq))){
      prime3$tr_strand[i] = '+'}
  }
}

prime3$tr_strand_ = ''
for (i in 1:length(prime3$repeatID)){
  tmp.seq = getSeq(genome, prime3$chr[i], start = prime3$start[i], end = prime3$end[i])
  if(grepl(paste0(prime3$dominant.motif[i],prime3$dominant.motif[i],prime3$dominant.motif[i]), as.character(tmp.seq))){
    prime3$tr_strand_[i] = '+'}
  else{
    if (grepl(paste0(prime3$dominant.motif[i],prime3$dominant.motif[i],prime3$dominant.motif[i]), as.character(reverseComplement(tmp.seq)))){
      prime3$tr_strand_[i] = '-'}
  }
}


prime3_motif = prime3[which(prime3$tr_strand == '-' | prime3$tr_strand == '+'),]
prime3_motif = prime3_motif[which(prime3_motif$tr_strand == prime3_motif$tr_strand_),]

prime3_motif$tr_motif = ''
for (i in 1:length(prime3_motif$repeatID)){
  if (prime3_motif$tr_strand[i] != ''){
    if (prime3_motif$strand[i] == '+'){
      prime3_motif$tr_motif[i] = ifelse(prime3_motif$tr_strand[i] == '+', prime3_motif$dominant.motif[i], revCompl(prime3_motif$dominant.motif[i]))}
    else{prime3_motif$tr_motif[i] = ifelse(prime3_motif$tr_strand[i] == '-', prime3_motif$dominant.motif[i], revCompl(prime3_motif$dominant.motif[i]))}
  }
}

prime3.all = as.data.frame( table(prime3_motif$dominant.motif)) 
colnames(prime3.all) = c('motif', 'freq')
prime3.all = prime3.all[order(prime3.all$freq, decreasing = T),]
prime3.all$motif = as.character(prime3.all$motif)
prime3.all$cg.perc = sapply(as.character(prime3.all$motif), getCGContent) * 100
prime3.all$col = color.gradient(prime3.all$cg.perc)

prime3.st = as.data.frame( table(prime3_motif$dominant.motif[prime3_motif$type == 'stable']))
colnames(prime3.st) = c('motif', 'freq.st')
prime3.st$motif = as.character(prime3.st$motif)
prime3.st = prime3.st[order(prime3.st$freq, decreasing = T),]

for (i in 1:nrow(prime3.all)){
  if (prime3.all$motif[i] %in% prime3.st$motif){
    prime3.all$freq.st[i] = prime3.st$freq.st[prime3.st$motif == prime3.all$motif[i]]
  }else{
    prime3.all$freq.st[i] = 0
  }
}

prime3.var = as.data.frame( table(prime3_motif$dominant.motif[prime3_motif$type == 'variable']))# 44;; 67 GTEx
colnames(prime3.var) = c('motif', 'freq.var')
prime3.var$motif = as.character(prime3.var$motif)
prime3.var = prime3.var[order(prime3.var$freq, decreasing = T),]

for (i in 1:nrow(prime3.all)){
  if (prime3.all$motif[i] %in% prime3.var$motif){
    prime3.all$freq.var[i] = prime3.var$freq.var[prime3.var$motif == prime3.all$motif[i]]
  }else{
    prime3.all$freq.var[i] = 0
  }
}


### PLOT - All - ST - VAR
motif.plot = prime3.all[prime3.all$freq > 3,]

pdf('output/figures/_Alu_3prime_GTEx_all_st_var_COUNT_25.pdf', w = 7, h = 6.5)#
par(mar = c(4.5,4,1,1), mgp = c(3,1,0), mfrow=c(3,1)) 
my_bar = barplot(motif.plot$freq, las = 2, main = '', # All STRs at 5\' of Alu elements',
                 col = 'grey', ylim = c(0,1250), ylab = 'Frequency')
text(my_bar, y = par("usr")[3]-50,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9,
     adj = 1, srt = 35)
text(my_bar, motif.plot$freq + 50, paste(motif.plot$freq) ,cex=0.9) 

my_bar = barplot(motif.plot$freq.st, las = 2, main = '', # Stable STRs at 5\' of Alu elements',
                 col = alpha('springgreen2', 0.6), ylim = c(0,1250), ylab = 'Frequency') 
text(my_bar, y = par("usr")[3]-50,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9,
     adj = 1, srt = 35)#, srt = 35
text(my_bar, motif.plot$freq.st + 50, paste(motif.plot$freq.st) ,cex=0.9) 

my_bar = barplot(motif.plot$freq.var, las = 2, main = '', # Variable STRs at 5\' of Alu elements',
                 col = alpha('firebrick1', 0.6), ylim = c(0,1250), ylab = 'Frequency')
text(my_bar, y = par("usr")[3]-50,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9,
     adj = 1, srt = 35)#, srt = 35
text(my_bar, motif.plot$freq.var + 50, paste(motif.plot$freq.var) ,cex=0.9) 
dev.off()



pdf('_GTEx/output/figures/_Alu_3prime_GTEx_all_st_var_FREQ_MIR_lab_SIGN_25_.pdf', w = 4, h = 5)
par(mar = c(4.5,4,1,1), mgp = c(3,1,0), mfrow=c(3,1)) 
my_bar = barplot(motif.plot$freq/sum(motif.plot$freq), las = 2, main = '', # All STRs at 5\' of Alu elements',
                 col = 'grey', ylim = c(0,0.5), ylab = 'Frequency')
text(my_bar, y = par("usr")[3]-0.02,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9,
     adj = 1, srt = 35)

par(mar=c(0,5,5,3))
my_bar = barplot(motif.plot$freq.st/sum(motif.plot$freq.st), las = 2, main = '', # Stable STRs at 5\' of Alu elements',
                 col = alpha('springgreen2', 0.6), ylim = c(0,0.5), ylab = 'STABLE') 

par(mar=c(5,5,0,3))
my_bar = barplot(motif.plot$freq.var/sum(motif.plot$freq.var), las = 2, main = '', # Variable STRs at 5\' of Alu elements',
                 col = alpha('firebrick1', 0.6), ylim = c(0.5,0), ylab = 'VARIABLE')

text(my_bar-0.5, y = par("usr")[3]+0.12,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.9, pos = 4,
     adj = 1, srt = 90)#, srt = 35

text(my_bar, y = par("usr")[3]+0.12,labels = ifelse(prime3.all$signif, '*', ''),xpd = NA, srt = 90,cex = 1.5,
     adj = 1, srt = 90)

dev.off()

prime3.all = read.delim('_GTEx/output/alu_3prime_GTEx_motif_freq_st_var.tsv', sep = ' ')
prime3.plot = prime3.all[prime3.all$motif %in% motif.plot$motif,]

prime3.alt = as.data.frame(strsplit(paste(prime3_motif$motif[prime3_motif$type == 'variable'], collapse = ';'), ';')[[1]])
colnames(prime3.alt) = 'motif'
prime3.alt$gc = sapply(prime3.alt$motif, getCGContent)

prime3.alt.freq = data.frame(table(prime3.alt$motif))
colnames(prime3.alt.freq) = c('motif', 'freq.alt')

length(prime3.all$motif %in% prime3.alt.freq$motif)

for (i in 1:nrow(prime3.all)){
  if (prime3.all$motif[i] %in% prime3.alt.freq$motif){
    prime3.all$freq.alt[i] = prime3.alt.freq$freq.alt[prime3.alt.freq$motif == prime3.all$motif[i]]
  }else{
    prime3.all$freq.alt[i] = 0
  }
}


### VAR vs ST
for (i in 1:nrow(prime3.all)){
  test <- fisher.test(data.frame("X1" = c(prime3.all$freq.var[i], prime3.all$freq.st[i]), 
                                 "X2" = c(sum(prime3.all$freq.var) - prime3.all$freq.var[i], 
                                          sum(prime3.all$freq.st) - prime3.all$freq.st[i])))#, alternative = 'greater')
  prime3.all$OR[i] = test$estimate
  prime3.all$OR.upper[i] = test$conf.int[1]
  prime3.all$OR.lower[i] = test$conf.int[2]
  prime3.all$Fisher.p[i] = test$p.value
  # format(test$p.value, scientific = T, digits = 2)# "1.8e-208"
}

### ST vs ALL
for (i in 1:nrow(prime3.all)){
  test <- fisher.test(data.frame("X1" = c(prime3.all$freq.st[i], prime3.all$freq[i]), 
                                 "X2" = c(sum(prime3.all$freq.st) - prime3.all$freq.st[i], 
                                          sum(prime3.all$freq) - prime3.all$freq[i])))#, alternative = 'greater')
  prime3.all$OR.st[i] = test$estimate
  prime3.all$OR.upper.st[i] = test$conf.int[1]
  prime3.all$OR.lower.st[i] = test$conf.int[2]
  prime3.all$Fisher.p.st[i] = test$p.value
  # format(test$p.value, scientific = T, digits = 2)# "1.8e-208"
}

### VAR vs ALL
for (i in 1:nrow(prime3.all)){
  test <- fisher.test(data.frame("X1" = c(prime3.all$freq.var[i], prime3.all$freq[i]), 
                                 "X2" = c(sum(prime3.all$freq.var) - prime3.all$freq.var[i], 
                                          sum(prime3.all$freq) - prime3.all$freq[i])))#, alternative = 'greater')
  prime3.all$OR.var[i] = test$estimate
  prime3.all$OR.upper.var[i] = test$conf.int[1]
  prime3.all$OR.lower.var[i] = test$conf.int[2]
  prime3.all$Fisher.p.var[i] = test$p.value
  # format(test$p.value, scientific = T, digits = 2)# "1.8e-208"
}

### ALT vs ALL
for (i in 1:nrow(prime3.all)){
  test <- fisher.test(data.frame("X1" = c(prime3.all$freq.alt[i], prime3.all$freq[i]), 
                                 "X2" = c(sum(prime3.all$freq.alt) - prime3.all$freq.alt[i], 
                                          sum(prime3.all$freq) - prime3.all$freq[i])))#, alternative = 'greater')
  prime3.all$OR.alt[i] = test$estimate
  prime3.all$OR.upper.alt[i] = test$conf.int[1]
  prime3.all$OR.lower.alt[i] = test$conf.int[2]
  prime3.all$Fisher.p.alt[i] = test$p.value
  # format(test$p.value, scientific = T, digits = 2)# "1.8e-208"
}

prime3.all$motif[prime3.all$OR.st > 1 & prime3.all$Fisher.p < 0.05]
prime3.all$motif[prime3.all$OR.var > 1 & prime3.all$Fisher.p < 0.05]
prime3.all$motif[prime3.all$OR.alt > 1 & prime3.all$Fisher.p < 0.05]

prime3.all$signif.var = ifelse(prime3.all$OR.var > 1 & prime3.all$Fisher.p < 0.05, T, F)
prime3.all$signif.st = ifelse(prime3.all$OR.st > 1 & prime3.all$Fisher.p < 0.05, T, F)
prime3.all$signif.alt = ifelse(prime3.all$OR.alt > 1 & prime3.all$Fisher.p.alt < 0.05, T, F)


pdf('_GTEx/output/figures/_Alu_3prime_GTEx_all_st_var_FREQ_MIR_lab_SIGN_24_ST_VAR_.pdf', w = 5, h = 4)
par(mar = c(4.5,4,1,4), mgp = c(3,1,0), mfrow=c(2,1)) 

my_bar = barplot(motif.plot$freq.st/sum(motif.plot$freq.st), las = 2, main = '', # Stable STRs at 5\' of Alu elements',
                 col = alpha('springgreen2', 0.6), ylim = c(0,0.5), ylab = 'Frequency') 
text(my_bar+1, y = par("usr")[3]-0.07,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.8, pos = 2,
     adj = 1, srt = 90)
text(my_bar, y = par("usr")[3]-0.01,labels = ifelse(prime3.all$signif.st, '*', ''),xpd = NA, srt = 90,cex = 1.5,
     adj = 1, srt = 90)

my_bar = barplot(motif.plot$freq.var/sum(motif.plot$freq.var), las = 2, main = '', # Variable STRs at 5\' of Alu elements',
                 col = alpha('firebrick1', 0.6), ylim = c(0,0.5), ylab = 'Frequency')
text(my_bar+1, y = par("usr")[3]-0.07,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.8, pos = 2,
     adj = 1, srt = 90)
text(my_bar, y = par("usr")[3]-0.01,labels = ifelse(prime3.all$signif, '*', ''),xpd = NA, srt = 90,cex = 1.5,
     adj = 1, srt = 90)
dev.off()

pdf('_GTEx/output/figures/_Alu_3prime_GTEx_all_st_var_FREQ_MIR_lab_SIGN_24_ST_VAR_ALT_versus_ALL_.pdf', w = 4, h = 4) 
par(mar = c(4.5,4,1,4), mgp = c(3,1,0), mfrow=c(3,1)) 

my_bar = barplot(motif.plot$freq.st/sum(motif.plot$freq.st), las = 2, main = '', # Stable STRs at 5\' of Alu elements',
                 col = alpha('springgreen2', 0.6), ylim = c(0,0.5), ylab = 'Frequency') 
text(my_bar+1, y = par("usr")[3]-0.07,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.8, pos = 2,
     adj = 1, srt = 90)
text(my_bar, y = par("usr")[3]-0.01,labels = ifelse(prime3.all$signif.st, '*', ''),xpd = NA, srt = 90,cex = 1.5,
     adj = 1, srt = 90)

my_bar = barplot(motif.plot$freq.var/sum(motif.plot$freq.var), las = 2, main = '', # Variable STRs at 5\' of Alu elements',
                 col = alpha('firebrick1', 0.6), ylim = c(0,0.5), ylab = 'Frequency')
text(my_bar+1, y = par("usr")[3]-0.07,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.8, pos = 2,
     adj = 1, srt = 90)
text(my_bar, y = par("usr")[3]-0.01,labels = ifelse(prime3.all$signif.var, '*', ''),xpd = NA, srt = 90,cex = 1.5,
     adj = 1, srt = 90)

my_bar = barplot(motif.plot$freq.alt/sum(motif.plot$freq.alt), las = 2, main = '', # Variable STRs at 5\' of Alu elements',
                 col = alpha('dodgerblue1', 0.6), ylim = c(0,0.5), ylab = 'Frequency')
text(my_bar+1, y = par("usr")[3]-0.07,labels = motif.plot$motif,xpd = NA, srt = 90,cex = 0.8, pos = 2,
     adj = 1, srt = 90)
text(my_bar, y = par("usr")[3]-0.01,labels = ifelse(prime3.all$signif.alt, '*', ''),xpd = NA, srt = 90,cex = 1.5,
     adj = 1, srt = 90)
dev.off()




prime3.all$OR = ''
prime3.all$OR.upper = ''
prime3.all$OR.lower = ''
prime3.all$Fisher.p = ''

### Calc OR between stable and variable
for (i in 1:nrow(prime3.all)){
  test <- fisher.test(data.frame("X1" = c(prime3.all$freq.var[i], prime3.all$freq.st[i]), 
                                 "X2" = c(sum(prime3.all$freq.var) - prime3.all$freq.var[i], 
                                          sum(prime3.all$freq.st) - prime3.all$freq.st[i])))#, alternative = 'greater')
  prime3.all$OR[i] = test$estimate
  prime3.all$OR.upper[i] = test$conf.int[1]
  prime3.all$OR.lower[i] = test$conf.int[2]
  prime3.all$Fisher.p[i] = test$p.value
  # format(test$p.value, scientific = T, digits = 2)# "1.8e-208"
}

prime3.all$motif[prime3.all$OR > 0 & prime3.all$Fisher.p < 0.05]
prime3.all$signif = ifelse(prime3.all$OR > 1 & prime3.all$Fisher.p < 0.05, T, F)
prime3.all$sign.st = ifelse(prime3.all$OR < 1 & prime3.all$Fisher.p < 0.05, T, F)

my_bar = barplot(as.numeric(prime3.all$OR[prime3.all$OR != 'Inf']), las = 2, main = '', # Variable STRs at 5\' of Alu elements',
                 col = alpha('firebrick1', 0.6), ylab = 'Freq (variable)')

prime3.all = read.delim('_GTEx/output/alu_3prime_GTEx_motif_freq_st_var.tsv', sep = ' ')# 2376


### BOXPLOT
prime3_motif = read.delim('_GTEx/output/alu_3prime_GTEx_motif_st_var.tsv', sep = ' ')
prime3_motif$dominant.motif = as.character(prime3_motif$dominant.motif)
prime3_motif$gc = sapply(prime3_motif$dominant.motif, getCGContent)

gc.stable = prime3_motif$gc[prime3_motif$type == 'stable']
gc.variable = prime3_motif$gc[prime3_motif$type == 'variable'] 

prime3.alt = as.data.frame(strsplit(paste(prime3_motif$motif[prime3_motif$type == 'variable'], collapse = ';'), ';')[[1]])
colnames(prime3.alt) = 'motif'
prime3.alt$gc = sapply(prime3.alt$motif, getCGContent)

test.st.dom = wilcox.test(gc.stable, gc.variable)
test.st.alt = wilcox.test(gc.stable, prime3.alt$gc)
test.dom.alt = wilcox.test(gc.variable, prime3.alt$gc)

pdf('_GTEx/output/figures/_Alu_3prime_GTEx_boxplot_st_dom_alt__PVAL_COLOR.pdf', w = 2.6, h = 4)#, h = 7
par(mar = c(4.5,4,1,1), mgp = c(3,1,0), mfrow=c(1,1)) 
boxplot(gc.stable, gc.var.dom, gc.var.alt, col = c(alpha(c('springgreen2', 'firebrick1', 'dodgerblue1'), 0.6)), 
        xaxt = 'n', yaxt = 'n', ylab = 'GC-content', ylim = c(0,1))
axis(1, at = c(1,2,3), labels = c('ST', 'DOM', 'ALT'))
axis(2, at = seq(0,1,0.2), las = 2)

segments(1,0.9,2,0.9)
text(1.5,0.95, labels = paste0('p = ', format(test.st.dom$p.value, digits=2)))

segments(1,0.8,3,0.8)
text(2.5,0.85, labels = paste0('p = ', format(test.st.alt$p.value, digits=2)))

segments(2,0.7,3,0.7)
text(2.5,0.75, labels = paste0('p = ', format(test.dom.alt$p.value, digits=2, scientific = T)))

dev.off()

pdf('_GTEx/output/figures/_Alu_3prime_GTEx_boxplot_st_var__PVAL.pdf', w = 2.6, h = 4)#, h = 7
par(mar = c(5,4,1,1), mgp = c(3,1,0), mfrow=c(1,1)) 
boxplot(gc.stable, gc.variable, col = c(alpha(c('springgreen2', 'firebrick1'), 0.6)), 
        xaxt = 'n', yaxt = 'n', ylab = 'GC-content', ylim = c(0,1))
axis(1, at = c(1,2), labels = c('ST', 'VAR'))
axis(2, at = seq(0,1,0.2), las = 2)

segments(1,0.9,2,0.9)
text(1.5,0.95, labels = paste0('p = ', format(test.st.dom$p.value, digits=2)))
dev.off()

write.table(prime5.all, 'output/alu_5prime_GTEx_motif_freq_st_var.tsv', quote = F, col.names = T, row.names = F)
write.table(prime5_motif, 'output/alu_5prime_GTEx_motif_st_var.tsv', quote = F, col.names = T, row.names = F)

write.table(prime3.all, 'output/alu_3prime_GTEx_motif_freq_st_var.tsv', quote = F, col.names = T, row.names = F)
write.table(prime3_motif, 'output/alu_3prime_GTEx_motif_st_var.tsv', quote = F, col.names = T, row.names = F)