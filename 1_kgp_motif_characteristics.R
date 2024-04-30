setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(plotrix)
library(GenomicRanges)
library(scales)

getCGContent <- function(seq){
  cg.len = length(grep("C|G", strsplit(seq, "")[[1]]))
  return(round(cg.len/nchar(seq), digits = 4))
}

load('GO_KEGG_n_100_1000.RData')
annovar = read.delim('../_KGP.STR.forannotation_fx.tsv.annovar.out_rev38.1m.tsv', stringsAsFactors = F)

table(annovar$typeseq_priority)
annovar$typeseq_priority = factor(annovar$typeseq_priority, 
                                  levels = c("exonic", "splicing", "exonic;splicing", "ncRNA_exonic", "ncRNA_splicing", "ncRNA_exonic;ncRNA_splicing", 
                                             "UTR5", "UTR3", "intronic", "ncRNA_intronic", "upstream", "downstream", "intergenic"))
annovar = annovar[order(annovar$typeseq_priority), ]
annovar = annovar[!duplicated(annovar$varid), ]

kgp.str = read.delim('../KGP.STR.All.info.no.genotype.Sep29.withExpansion.tsv', stringsAsFactors = F)# 11453 
kgp.str = kgp.str[kgp.str$chr %in% paste0("chr", c(1:22)), ]# 10884 - chr1-22
kgp.str = merge(kgp.str, annovar[, c("varid", "entrez_id", "typeseq_priority", 'gene_symbol')], by = "varid")
write.table(kgp.str, 'kgp.str.type.str.tsv', quote = F, row.names = F,sep = '\t')

kgp.str$gc.cont = getCGContent(kgp.str$dominant.motif)
kgp.str$motif.length = nchar(kgp.str$dominant.motif)

h.kgp.gc = hist(kgp.str$gc.cont, yaxt = 'n', freq = T, breaks = 20)
pdf('kgp_gc-comp_unexp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.kgp.gc$counts/sum(h.kgp.gc$counts),col = alpha('thistle', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n', main = 'GC-content in KGP (unexp)')
axis(1, at = b[c(3,7,11,15,18)], labels = h.kgp.gc$breaks[c(3,7,11,15,19)], las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(kgp.str$gc.cont, col = alpha('thistle', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'GC-composition (%)', horizontal = T)
dev.off()

h.kgp.size = hist(kgp.str$dominant.motif.size, yaxt = 'n', freq = T)
pdf('kgp_motif-size_unexp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.kgp.size$counts/sum(h.kgp.size$counts),col = alpha('thistle', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n', main = 'Motif size in KGP (unexp)')
axis(1, at = b[c(1,3,5,7,11,15,18)], labels = h.kgp.size$breaks[c(1,3,5,7,11,15,19)], las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(kgp.str$dominant.motif.size, col = alpha('thistle', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'Dominant motif size (bp)', horizontal = T)
dev.off()

# Expanded
kgp.str.exp = kgp.str[kgp.str$expanded == T,]# 1964
h.kgp.exp.gc = hist(kgp.str.exp$gc.cont, yaxt = 'n', freq = T, breaks = 20)

pdf('kgp_gc-comp_exp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.kgp.exp.gc$counts/sum(h.kgp.exp.gc$counts),col = alpha('orchid4', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n', main = 'GC-content in KGP (exp)')
axis(1, at = b[c(3,7,11,15,18)], labels = h.kgp.gc$breaks[c(3,7,11,15,19)], las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(kgp.str.exp$gc.cont, col = alpha('orchid4', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'GC-composition (%)', horizontal = T)
dev.off()

h.kgp.exp.size = hist(kgp.str.exp$dominant.motif.size, yaxt = 'n', freq = T, breaks = 20)
pdf('kgp_motif-size_exp.pdf', w = 4, h = 4)
par(mar = c(1,5,3,3), mgp = c(3,1,0), mfrow=c(2,1)) 
b = barplot(h.kgp.exp.size$counts/sum(h.kgp.exp.size$counts),col = alpha('orchid4', 0.8), 
            xlab = '', ylab = 'Percent of regions', ylim = c(0, 0.3),yaxt = 'n', main = 'Motif size in KGP (exp)')
axis(1, at = b[c(1,3,5,7,11,15,18)], labels = h.kgp.exp.size$breaks[c(1,3,5,7,11,15,19)], las = 1, padj = -1, tick = T)
axis(2, at = seq(0,0.3, 0.1), labels = T, las = 2)
par(mar = c(6,5,1,3), mgp = c(1,1,0)) 
boxplot(kgp.str.exp$dominant.motif.size, col = alpha('orchid4', 0.8), 
        yaxt = 'n', xaxt = 'n', ylab ='', xlab = 'Dominant motif size (bp)', horizontal = T)
dev.off()