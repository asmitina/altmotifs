setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(plotrix)
library(GenomicRanges)
library(scales)
library(data.table)

#### mobile element alu sub family
####################
#### mobile elements
mle <- fread("../UCSC_hg38_repeatMasker.tsv", data.table = F)
mle <- mle[mle$repFamily == "Alu", ]

### group subfamily
mle$repName <- substr(mle$repName, 0, 4)
mle$repName <- ifelse(mle$repName %in% c("Alu", "AluY", "AluS", "AluJ"),
                      mle$repName, "Others")
###
exp = gtex.str

mle.g <- GRanges(mle$genoName, IRanges(mle$genoStart, mle$genoEnd), "*")
exp.g <- GRanges(exp$chr, IRanges(exp$start, exp$end), "*")
olap <- data.frame(findOverlaps(mle.g, exp.g))
olap$type <- exp$type[olap$subjectHits]
mle$type <- NA
mle$type[olap$queryHits] <- olap$type
dt.plot <- data.frame(t(table(mle$type, mle$repName)))
dt.plot <- dt.plot[dt.plot$Var1 %in% c("AluJ", "AluS", "AluY"), ]
names(dt.plot)[1] <- "Alu family"
dt.plot$percentage <- 0
dt.plot$percentage[dt.plot$Var2 == "stable"] <- 
  100*dt.plot$Freq[dt.plot$Var2 == "stable"]/sum(dt.plot$Freq[dt.plot$Var2 == "stable"])
dt.plot$percentage[dt.plot$Var2 == "variable"] <- 
  100*dt.plot$Freq[dt.plot$Var2 == "variable"]/sum(dt.plot$Freq[dt.plot$Var2 == "variable"])
ggplot(dt.plot, aes(x = Var2, y = percentage, fill = `Alu family`)) + xlab("") +
  geom_bar(stat = "identity", position = position_stack(), size = .2) + theme_bw() + scale_fill_brewer(palette = "Spectral")
ggsave("output/figures/alu.dis.all.png", width = 4, height = 6)

dt.out <- data.frame()
ext <- 0

for(repFam in unique(mle$repName)){
  mle.tmp <- mle[mle$repName == repFam, ]
  mle.tmp$genoStart <- ifelse(mle.tmp$genoStart - ext < 1, 1, mle.tmp$genoStart - ext)
  mle.tmp$genoEnd <- mle.tmp$genoEnd + ext
  
  mle.tmp.g <- GRanges(mle.tmp$genoName, IRanges(mle.tmp$genoStart, mle.tmp$genoEnd), "*")
  exp.g <- GRanges(exp$chr, IRanges(exp$start, exp$end), "*")
  
  olap <- findOverlaps(exp.g, mle.tmp.g)
  
  tmp.exp <- exp[unique(olap@from), ]
  
  variable.count <- sum(tmp.exp$nmotif > 1, na.rm = T)
  stable.count <- sum(tmp.exp$nmotif == 1, na.rm = T)
  
  variable.denom <- sum(exp$nmotif > 1, na.rm = T) - variable.count
  stable.denom <- sum(exp$nmotif == 1, na.rm = T) - stable.count
  
  test <- fisher.test(data.frame("X1" = c(variable.count, stable.count),
                                 "X2" = c(variable.denom, stable.denom), stringsAsFactors = F))
  
  dt.out <- rbind(dt.out, data.frame("feature" = repFam, "extension" = ext,
                                     "OR" = round(test$estimate, digits = 2),
                                     "OR.lower" = round(test$conf.int[1], digits = 2),
                                     "OR.upper" = round(test$conf.int[2], digits = 2),
                                     "p" = signif(test$p.value, digits = 2), stringsAsFactors = F
  ))
}


dt.out$extension <- paste0(dt.out$extension/1000, "KB")
dt.out$extension[dt.out$extension == "0KB"] <- "no extension"
dt.out$extension <- factor(dt.out$extension, levels = unique(dt.out$extension))
dt.out$fdr <- signif(p.adjust(dt.out$p, method = "BH"), digits = 2)
dt.out$fdr.range <- ifelse(dt.out$fdr < 0.20, "FDR 20%", "Not significant")
dt.out$fdr.range[dt.out$fdr < 0.20] <- "FDR 20%"
dt.out$fdr.range[dt.out$fdr < 0.15] <- "FDR 15%"
dt.out$fdr.range[dt.out$fdr < 0.10] <- "FDR 10%"
dt.out$fdr.range[dt.out$fdr < 0.05] <- "FDR 5%"
dt.out$fdr.range[dt.out$fdr < 0.01] <- "FDR 1%"
dt.out$fdr.range <- factor(dt.out$fdr.range, levels = c("FDR 1%", "FDR 5%", "FDR 10%",
                                                        "FDR 15%", "FDR 20%", "Not significant"))

colorscales <- list("FDR 1%" = "red3", "FDR 5%" = "red", "FDR 10%" = "darkorange",
                    "FDR 15%" = "orange", "FDR 20%" = "yellow", "Not significant" = "grey")

subfam <- unique(names(table(mle$repName))[table(mle$repName) > 1000])
# dt.out$feature <- factor(dt.out$feature, levels = dt.out$feature[order(dt.out$fdr)])
ggplot(dt.out[dt.out$feature %in% subfam, ], aes(x = feature, y = OR, color = fdr.range)) + 
  geom_hline(yintercept = 1, lty = 2) +
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = OR.lower, ymax = OR.upper), width =.3) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
  ylab("Odds ratio (variable/stable)") +
  scale_color_manual(values = colorscales)

# ggsave("../output/figures/Alu.sub.families.png", width = 11, height = 4)
ggsave("output/figures/Alu.sub.families.grouped.png", width = 5, height = 4)
write.table(dt.out, "output/Alu.sub.families.tsv", sep = "\t", row.names=F, quote=F, col.names=T)

### Alu Combined
alu.kgp = read.delim('../_1KGP/output/SINE.Alu.sub.families.tsv')
alu.kgp$cohort = 'KGP'
alu.gtex = read.delim('output/SINE.Alu.sub.families.tsv')
alu.gtex$cohort = 'GTEX'

sine = rbind(alu.kgp, alu.gtex)
sine$move = sapply(sine$cohort, switch, "GTEX" = 0.2, "KGP" = (-0.2))
sine = sine[order(sine[,1], sine[,9]),] 
sine$col = unlist(sapply(sine$fdr.range, switch, 'Not significant'='grey', 'FDR 25%'='grey', 'FDR 20%'='grey', 'FDR 15%'='khaki3',
                         'FDR 10%'='yellow','FDR 5%'='orange', 'FDR 1%'='red'))

genic.plot = sine
pdf('output/figures/Alu_families_gtex_kgp_name.pdf', w = 5, h = 4)
par(mar = c(9,6,1,8), mgp = c(3,1,0), mfrow=c(1,1)) 
plotCI((1:10)+genic.plot$move, genic.plot$OR, ylab = 'OddsRatio (variable vs stable)', xlab = '',
       main = '', ui = genic.plot$OR.upper, li = genic.plot$OR.lower, pch = NA,
       lwd = 1, xaxt="n", ylim = c(0,4), col = 'black', cex = 1, yaxt = 'n', )
points((1:10)+genic.plot$move, genic.plot$OR, pch = 21, bg = genic.plot$col, cex = 1.5, col = 'black')

axis(1, at=seq(1, 10, by=2)+0.5, labels = FALSE)
axis(2, at=seq(0, 6, by=1), labels = TRUE, las = 2)
abline(h = 1, lwd = 0.5, lty = 2)

text((1:5)*2-0.5, par("usr")[3]-0.5, srt = 90, adj = 1, xpd = TRUE, cex = 1, 
     labels = paste(unique(genic.plot$feature)))
legend('topright', pch=19, pt.cex=1.3, inset=c(-0.7,0), xpd = TRUE, bty = "n", 
       col = c('red','orange', 'yellow', 'khaki3','grey'), 
       legend=c( "FDR 1%","FDR 5%","FDR 10%","FDR 15%", "Not significant"), title = 'FDR range')
legend('topright', pch=21, pt.cex=1.3, inset=c(-0.7,0), xpd = TRUE, bty = "n", col = 'black',
       bg = c('red','orange', 'yellow', 'khaki3','grey'), 
       legend=c( "FDR 1%","FDR 5%","FDR 10%","FDR 15%", "Not significant"), title = '')
dev.off()
