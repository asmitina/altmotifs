setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(plotrix)
library(GenomicRanges)
library(scales)

getCGContent <- function(seq){
  cg.len = length(grep("C|G", strsplit(seq, "")[[1]]))
  return(round(cg.len/nchar(seq), digits = 4))
}

gtex = read.delim('../meta/GTEx_metadata.tsv', stringsAsFactors = F)
load('GO_KEGG_n_100_1000.RData')
annovar = read.delim('merged.trs.forannotation.2022-05-05.tsv_annovarIn.txt.annovar.out_rev27.6_hg38.tsv', stringsAsFactors = F)# 2447 -> 47273 Yes!
annovar$typeseq_priority = factor(annovar$typeseq_priority, 
                                  levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                             "UTR5", "UTR3", "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
annovar = annovar[order(annovar$typeseq_priority), ]
annovar = annovar[!duplicated(annovar$varid), ]# 2447 -> 45990
gtex.ann.g <- GRanges(annovar$chr, IRanges(annovar$start, annovar$end), "*")
gtex.str.g <- GRanges(gtex.str$chr, IRanges(gtex.str$start, gtex.str$end), "*")

olap <- data.frame(findOverlaps(gtex.str.g, gtex.ann.g))# 2090
annovar$typeseq_priority = as.character(annovar$typeseq_priority)
gtex.str$typeseq_priority = ''
for (i in olap$queryHits){
  gtex.str$typeseq_priority[i] = annovar$typeseq_priority[olap$subjectHits[olap$queryHits == i]]
}
dim(gtex.str[gtex.str$typeseq_priority != '',])# 18397

gtex.str = read.delim('GTEx.tr.regions.samples.tsv', stringsAsFactors = F)# 18397
gtex.str = gtex.str[gtex.str$chr %in% paste0("chr", c(1:22)), ]# 18397
gtex.str$varid <- paste(gtex.str$chr, gtex.str$start, gtex.str$end, sep="#")

gtex.exp = read.delim('merged.rare.expansions.2020-08-07.tsv', stringsAsFactors = F)# 3249
gtex.exp$varid <- paste(gtex.exp$chr, gtex.exp$start, gtex.exp$end, sep="#")
gtex.exp$gtex = ''
for (i in 1:nrow(gtex.exp)){
  samples = strsplit(gtex.exp$outliers[i], ';')[[1]]
  if (length(samples[samples %in% gtex$Sample.ID]) > 0){
    gtex.exp$gtex[i] = paste(samples[samples %in% gtex$Sample.ID], sep = ';')
  }else{next}
}
gtex.exp = gtex.exp[gtex.exp$gtex != '',]# 518 

gtex.str.g <- GRanges(gtex.str$chr, IRanges(gtex.str$start, gtex.str$end), "*")
gtex.exp.g <- GRanges(gtex.exp$chr, IRanges(gtex.exp$start, gtex.exp$end), "*")

olap <- findOverlaps(gtex.str.g, gtex.exp.g)
gtex.str$expanded <- F
gtex.str$expanded[unique(olap@from)] <- T

for (i in 1:nrow(gtex.str)){
  gtex.str$gc.cont[i] = getCGContent(gtex.str$dominant.motif[i])
  gtex.str$motif.length[i] = nchar(gtex.str$dominant.motif[i])
}

h.gtex.gc = hist(gtex.str$gc.cont, yaxt = 'n', freq = T, breaks = 20)# 24
pdf('gtex_gc-comp_unexp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.gtex.gc$counts/sum(h.gtex.gc$counts),col = alpha('thistle', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.4),yaxt = 'n', main = 'GC-content in gtex (unexp)')
axis(1, at = b[c(1,5,9,13,17,20)], labels = h.gtex.gc$breaks[c(1,5,9,13,17,21)], las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.4, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(gtex.str$gc.cont, col = alpha('thistle', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'GC-composition (%)', horizontal = T)
dev.off()

pdf('gtex_gc-comp_exp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.gtex.exp.gc$counts/sum(h.gtex.exp.gc$counts),col = alpha('orchid4', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.4),yaxt = 'n', main = 'GC-content in gtex (exp)')
axis(1, at = b[c(1,5,9,13,17)], labels = h.gtex.exp.gc$breaks[c(1,5,9,13,17)], las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.4, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(gtex.str.exp$gc.cont, col = alpha('orchid4', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'GC-composition (%)', horizontal = T)
dev.off()

pdf('gtex_motif-size_unexp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.gtex.size$counts/sum(h.gtex.size$counts),col = alpha('thistle', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n', main = 'Motif size in gtex (unexp)')
axis(1, at = b[c(1,3,5,7,11,15,18)], labels = h.gtex.size$breaks[c(1,3,5,7,11,15,19)], las = 1, padj = -1, tick = T)
# axis(1, at = c(0,h.gtex.gc$mids+0.025)[c(1,9,19)], labels = h.gtex.gc$breaks[c(1,9,19)], las = 2, padj = -1, tick = F)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(gtex.str$dominant.motif.size, col = alpha('thistle', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'Dominant motif size (bp)', horizontal = T)
dev.off()

# then take expanded and plot for expanded (I already have that info in my table!!)
gtex.str.exp = gtex.str[gtex.str$expanded == T,]# 514 -> 492 ?? first do for merged, then do for EHdn expansions
h.gtex.exp.gc = hist(gtex.str.exp$gc.cont, yaxt = 'n', freq = T, breaks = 20)

pdf('gtex_gc-comp_exp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.gtex.exp.gc$counts/sum(h.gtex.exp.gc$counts),col = alpha('orchid4', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.4),yaxt = 'n', main = 'GC-content in gtex (exp)')
# axis(1, at = b[c(3,7,11,15,18)], labels = h.gtex.exp.gc$breaks[c(3,7,11,15,19)], las = 1, padj = -1, tick = T)
# axis(1, at = b[c(4,9,14,19,20)], labels = h.gtex.exp.gc$breaks[c(4,9,14,19,21)], las = 1, padj = -1, tick = T)
axis(1, at = b[c(1,5,9,13,17)], labels = h.gtex.exp.gc$breaks[c(1,5,9,13,17)], las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.4, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(gtex.str.exp$gc.cont, col = alpha('orchid4', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'GC-composition (%)', horizontal = T)
dev.off()


h.gtex.size = hist(gtex.str.exp$dominant.motif.size, yaxt = 'n', freq = T, breaks = 20)
pdf('gtex_motif-size_exp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.gtex.size$counts/sum(h.gtex.size$counts),col = alpha('orchid4', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n', main = 'Motif size in gtex (exp)')
axis(1, at = b[c(1,3,5,7,11,15,18)], labels = h.gtex.size$breaks[c(1,3,5,7,11,15,19)], las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(gtex.str.exp$dominant.motif.size, col = alpha('orchid4', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'Dominant motif size (bp)', horizontal = T)
dev.off()

write.table(dt.out, '../_1KGP/dt.out.tsv', quote = F, sep = '\t', col.names = T, row.names = F)


# Stable vs Variable
gtex.str.var = gtex.str[gtex.str$type == 'variable',]# 2032
gtex.str.st = gtex.str[gtex.str$type == 'stable',]# 16365

h.gtex.gc.var = hist(gtex.str$gc.cont[gtex.str$type == 'variable'], yaxt = 'n', freq = T, breaks = 20)# 24
h.gtex.gc.st = hist(gtex.str$gc.cont[gtex.str$type == 'stable'], yaxt = 'n', freq = T, breaks = 20)# 24

pdf('gtex_gc-comp_unexp_ST_VAR.pdf', w = 3, h = 3.5)
par(mar = c(1,5,1,3), mgp = c(3,1,0), mfrow=c(3,1)) 
b = barplot(h.gtex.gc.st$counts/sum(h.gtex.gc.st$counts),col = alpha('springgreen2', 0.6), 
            xlab = '', ylab = 'Frequency', ylim = c(0, 0.4),yaxt = 'n',
            main = '', space = c())
axis(1, at = b[c(1,5,9,13,17,20)], labels = h.gtex.gc$breaks[c(1,5,9,13,17,21)]*100, 
     las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)

b = barplot(h.gtex.gc.var$counts/sum(h.gtex.gc.var$counts),col = alpha('firebrick1', 0.6), 
            xlab = '', ylab = 'Frequency', ylim = c(0, 0.4),yaxt = 'n',
            main = '', space = c())
axis(1, at = b[c(1,5,9,13,17,20)], labels = h.gtex.gc$breaks[c(1,5,9,13,17,21)]*100, 
     las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)

par(mar = c(3,5,3,3), mgp = c(1,1,0)) 
boxplot(gtex.str.var$gc.cont, gtex.str.st$gc.cont, col = alpha(c('firebrick1','springgreen2') , 0.6), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'GC-composition (%)', horizontal = T)
dev.off()

test.gc = wilcox.test(gtex.str$gc.cont[gtex.str$type == 'stable'], 
                      gtex.str$gc.cont[gtex.str$type == 'variable'])#, alternative = 'less')
format(test.gc$p.value, scientific=T, digits=2)# "2.5e-20"


h.gtex.size = hist(gtex.str$dominant.motif.size, yaxt = 'n', freq = T)
h.gtex.size.st = hist(gtex.str$dominant.motif.size[gtex.str$type == 'stable'], yaxt = 'n', freq = T, breaks = 20)
h.gtex.size.var = hist(gtex.str$dominant.motif.size[gtex.str$type == 'variable'], yaxt = 'n', freq = T, breaks = 20)


pdf('gtex_motif-size_unexp_ST_VAR.pdf', w = 3, h = 3.5)
par(mar = c(1,5,1,3), mgp = c(3,1,0), mfrow=c(3,1)) 
b = barplot(h.gtex.size.st$counts/sum(h.gtex.size.st$counts),col = alpha('springgreen2', 0.6), 
            xlab = '', ylab = 'Frequency', ylim = c(0, 0.4),yaxt = 'n', main = '')
axis(1, at = b[c(1,3,5,7,11,15,18)], labels = h.gtex.size.st$breaks[c(1,3,5,7,11,15,19)], las = 1, padj = -1, tick = T)
# axis(1, at = c(0,h.gtex.gc$mids+0.025)[c(1,9,19)], labels = h.gtex.gc$breaks[c(1,9,19)], las = 2, padj = -1, tick = F)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)

b = barplot(h.gtex.size.var$counts/sum(h.gtex.size.var$counts),col = alpha('firebrick1', 0.6), 
            xlab = '', ylab = 'Frequency', ylim = c(0, 0.4),yaxt = 'n', main = '')
axis(1, at = b[c(1,3,5,7,11,15,18)], labels = h.gtex.size.var$breaks[c(1,3,5,7,11,15,19)], las = 1, padj = -1, tick = T)
# axis(1, at = c(0,h.gtex.gc$mids+0.025)[c(1,9,19)], labels = h.gtex.gc$breaks[c(1,9,19)], las = 2, padj = -1, tick = F)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)

par(mar = c(3,5,3,3), mgp = c(1,1,0)) 
boxplot(gtex.str.var$dominant.motif.size, gtex.str.st$dominant.motif.size, col = alpha(c('firebrick1','springgreen2') , 0.6), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'Dominant motif size (bp)', horizontal = T)
dev.off()

test.size = wilcox.test(gtex.str$dominant.motif.size[gtex.str$type == 'stable'], 
                        gtex.str$dominant.motif.size[gtex.str$type == 'variable'])#, alternative = 'less')
format(test.size$p.value, scientific=T, digits=2)# "3.7e-103"

"#FF3030CC" "#00EE76CC"

gtex.str.exp = gtex.str[gtex.str$expanded == T,]# 514 -> 492 ?? first do for merged, then do for EHdn expansions
h.gtex.exp.gc = hist(gtex.str.exp$gc.cont, yaxt = 'n', freq = T, breaks = 20)
