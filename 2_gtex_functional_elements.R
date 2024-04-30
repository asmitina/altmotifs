setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(plotrix)
library(GenomicRanges)
library(scales)

### Test functional elements for all
exp = gtex.str
dt.out <- data.frame()
for(typeseq in unique(as.character(exp$typeseq_priority))){
  variable.count <- sum(exp$typeseq_priority == typeseq & exp$nmotif > 1, na.rm = T)
  stable.count <- sum(exp$typeseq_priority == typeseq & exp$nmotif == 1, na.rm = T)
  
  variable.denom <- sum(exp$nmotif > 1, na.rm = T) - variable.count
  stable.denom <- sum(exp$nmotif == 1, na.rm = T) - stable.count
  
  test <- fisher.test(data.frame("X1" = c(variable.count, stable.count),
                                 "X2" = c(variable.denom, stable.denom), stringsAsFactors = F))
  
  dt.out <- rbind(dt.out, data.frame("feature" = typeseq, "OR" = round(test$estimate, digits = 2),
                                     "OR.lower" = round(test$conf.int[1], digits = 2),
                                     "OR.upper" = round(test$conf.int[2], digits = 2),
                                     "p" = signif(test$p.value, digits = 2), stringsAsFactors = F
  ))
}

dt.out$feature <- factor(dt.out$feature, levels = dt.out$feature[order(dt.out$OR, decreasing = T)])
dt.out$fdr <- signif(p.adjust(dt.out$p, method = "BH"), digits = 2)
dt.out$fdr.range <- ifelse(dt.out$fdr < 0.2, "FDR 20%", "Not significant")
dt.out$fdr.range[dt.out$fdr < 0.20] <- "FDR 20%"
dt.out$fdr.range[dt.out$fdr < 0.15] <- "FDR 15%"
dt.out$fdr.range[dt.out$fdr < 0.10] <- "FDR 10%"
dt.out$fdr.range[dt.out$fdr < 0.05] <- "FDR 5%"
dt.out$fdr.range[dt.out$fdr < 0.01] <- "FDR 1%"
dt.out$fdr.range <- factor(dt.out$fdr.range, levels = c("FDR 1%", "FDR 5%", "FDR 10%",
                                                        "FDR 15%", "FDR 20%", "Not significant"))
# ggplot(dt.out, aes(x = feature, y = OR, color = fdr.range)) + geom_hline(yintercept = 1, lty = 2) +
#   geom_point() + geom_errorbar(aes(ymin = OR.lower, ymax = OR.upper), width =.3) + 
#   theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
#   coord_cartesian(ylim = c(0, 6.5)) + ylab("Odds ratio (variable/stable)") +
#   scale_color_manual(values = c("red", "orange", "yellow", "grey"))
#
# ggsave("output/figures/genic.element.png", width = 6, height = 4)
# write.table(dt.out, "output/genic.element.tsv", sep = "\t", row.names=F, quote=F, col.names=T)
dt.out$fdr.range = as.character(dt.out$fdr.range)
dt.out$col = unlist(sapply(dt.out$fdr.range, switch, 'Not significant'='grey', 'FDR 20%'='khaki', 'FDR 15%'='grey',
                           'FDR 10%'='yellow','FDR 5%'='orange', 'FDR 1%'='red'))
table(dt.out$fdr.range)
dt.out = dt.out[order(dt.out$OR, decreasing = T),]
genic.plot = dt.out

pdf('output/figures/genic_burden_gtex_.pdf', w = 6, h = 4)
par(mar = c(9,6,1,7), mgp = c(3,1,0), mfrow=c(1,1)) 
b = barplot(height = genic.plot$OR,
            col = alpha(genic.plot$col, alpha = 0.8),
            yaxt = 'n', ylim = c(0,7), ylab = 'Odds ratio')
arrows(b, genic.plot$OR.lower, b, genic.plot$OR.upper, col='black', lwd = 1, code = 3, angle = 90, length = 0.05)# code=1
abline(h = 1, lty = 2)
axis(2, at = c(0,1,5,7), las = 2)
axis(1, at = c(0.6:13)*1.2, las = 2, labels = F, tck=-0.05)

text(c(0.7:13)*1.2, par("usr")[3]-0.5, srt = 45, adj = 1, xpd = TRUE, cex = 0.9, 
     labels = paste(genic.plot$feature))

legend(x = "right", legend = c('FDR 1%','FDR 5%', 'FDR 20%', 'Not significant'), 
       fill = alpha(c('red', 'orange', 'khaki','grey'), alpha = 0.8), 
       inset=c(-0.4,-0.2), bty = "n", xpd = TRUE, title = 'FDR range')

dev.off()


#### COMBINED
dt.out.gtex = read.delim('_GTEx/output/genic.element.tsv', stringsAsFactors = F)
dt.out.kgp = read.delim('_1KGP/output/genic.element.tsv', stringsAsFactors = F)

dt.out.gtex = dt.out.gtex[order(match(dt.out.gtex$feature, dt.out.kgp$feature)),]
dt.out.gtex$fdr.range = as.character(dt.out.gtex$fdr.range)
dt.out.gtex$col = unlist(sapply(dt.out.gtex$fdr.range, switch, 'Not significant'='grey', 'FDR 20%'='khaki', 'FDR 15%'='grey',
                                'FDR 10%'='yellow','FDR 5%'='orange', 'FDR 1%'='red'))
table(dt.out.gtex$fdr.range)
dt.out.gtex = dt.out.gtex[order(match(dt.out.gtex$feature, dt.out.kgp$feature)),]
genic.plot.gt = dt.out.gtex

dt.out.kgp$fdr.range = as.character(dt.out.kgp$fdr.range)
dt.out.kgp$col = unlist(sapply(dt.out.kgp$fdr.range, switch, 'Not significant'='grey', 'FDR 20%'='khaki', 'FDR 15%'='grey',
                               'FDR 10%'='yellow','FDR 5%'='orange', 'FDR 1%'='red'))
table(dt.out.kgp$fdr.range)
dt.out.kgp = dt.out.kgp[order(dt.out.kgp$OR, decreasing = T),]
genic.plot.kg = dt.out.kgp

pdf('output/figures/genic_burden_GTEx_KGP_combined.pdf', w = 6.5, h = 4)
par(mar = c(9,6,1,7), mgp = c(3,1,0), mfrow=c(1,1)) 
b = barplot(c(rbind(genic.plot.gt$OR, genic.plot.kg$OR)),
            col = c(rbind(alpha(genic.plot.gt$col, alpha = 0.8), alpha(genic.plot.kg$col, alpha = 0.8))), 
            space = c(rbind(rep(1,13), rep(0,13))), yaxt = 'n', ylim = c(0,7), 
            ylab = 'Odds ratio')

arrows(b, c(rbind(genic.plot.gt$OR.lower, genic.plot.kg$OR.lower)), b, 
       c(rbind(genic.plot.gt$OR.upper, genic.plot.kg$OR.upper)), col='black', 
       lwd = 1, code = 3, angle = 90, length = 0.03)# code=1

abline(h = 1, lty = 2)
axis(2, at = c(0,1,5,7), las = 2)
axis(1, at = c(0.67:13)*3, las = 2, labels = F, tck=-0.05)
text(c(0.67:13)*3, par("usr")[3]-0.5, srt = 45, adj = 1, xpd = TRUE, cex = 0.9, 
     labels = paste(genic.plot.kg$feature))
legend(x = "right", legend = c('FDR 1%','FDR 5%', 'FDR 10%', 'FDR 20%', 'Not significant'), 
       fill = alpha(c('red', 'orange', 'yellow', 'khaki','grey'), alpha = 0.8), 
       inset=c(-0.35,-0.2), bty = "n", xpd = TRUE, title = 'FDR range')
dev.off()


pdf('output/figures/genic_burden_KGP_GTEx_combined.pdf', w = 6.5, h = 4)
par(mar = c(9,6,1,7), mgp = c(3,1,0), mfrow=c(1,1)) 
b = barplot(c(rbind(genic.plot.kg$OR, genic.plot.gt$OR)),
            col = c(rbind(alpha(genic.plot.kg$col, alpha = 0.8), alpha(genic.plot.gt$col, alpha = 0.8))), 
            space = c(rbind(rep(1,13), rep(0,13))), yaxt = 'n', ylim = c(0,7), 
            ylab = 'Odds ratio')

arrows(b, c(rbind(genic.plot.kg$OR.lower, genic.plot.gt$OR.lower)), b, 
       c(rbind(genic.plot.kg$OR.upper, genic.plot.gt$OR.upper)), col='black', 
       lwd = 1, code = 3, angle = 90, length = 0.03)

abline(h = 1, lty = 2)
axis(2, at = c(0,1,5,7), las = 2)
axis(1, at = c(0.67:13)*3, las = 2, labels = F, tck=-0.05)
text(c(0.67:13)*3, par("usr")[3]-0.5, srt = 45, adj = 1, xpd = TRUE, cex = 0.9, 
     labels = paste(genic.plot.kg$feature))
legend(x = "right", legend = c('FDR 1%','FDR 5%', 'FDR 10%', 'FDR 20%', 'Not significant'), 
       fill = alpha(c('red', 'orange', 'yellow', 'khaki','grey'), alpha = 0.8), 
       inset=c(-0.35,-0.2), bty = "n", xpd = TRUE, title = 'FDR range')
dev.off()

### Expandability
variable.count <- nrow(gtex.str[gtex.str$type == 'variable' & gtex.str$expanded == T,])# 157
stable.count <- nrow(gtex.str[gtex.str$type == 'stable' & gtex.str$expanded == T,])# 335

variable.denom <- nrow(gtex.str[gtex.str$type == 'variable' & gtex.str$expanded == F,])# 1875
stable.denom <- nrow(gtex.str[gtex.str$type == 'stable' & gtex.str$expanded == F,])# 16030

gtex.test <- fisher.test(data.frame("X1" = c(variable.count, stable.count), 
                                    "X2" = c(variable.denom, stable.denom)))#, alternative = 'greater')
format(gtex.test$p.value, scientific = T, digits = 2)# "9e-37"
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   3.273000 4.887555
# sample estimates:
#   odds ratio 
# 4.006629 