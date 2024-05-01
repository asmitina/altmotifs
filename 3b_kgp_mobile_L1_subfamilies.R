setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(plotrix)
library(GenomicRanges)
library(scales)
library(data.table)

#### mobile element LINE sub family
####################
#### mobile elements
mle <- fread("../_Bank_files/UCSC_hg38_repeatMasker.tsv", data.table = F)
mle <- mle[mle$repFamily == "L1", ]

### group subfamily
mle$repName <- substr(mle$repName, 0, 4)
mle$repName <- gsub("P[0-9]", "P", mle$repName)
mle$repName <- gsub("M[0-9]", "M", mle$repName)
mle <- mle[mle$repName %in% c("L1M", "L1MA", "L1MB", "L1MC", "L1MD", "L1ME",
                              "L1P", "L1PA", "L1PB"), ]

table(mle$repName)
mle.g <- GRanges(mle$genoName, IRanges(mle$genoStart, mle$genoEnd), "*")
exp.g <- GRanges(exp$chr, IRanges(exp$start, exp$end), "*")
olap <- data.frame(findOverlaps(mle.g, exp.g))
olap$type <- exp$type[olap$subjectHits]
olap$expanded <- exp$expanded[olap$subjectHits]

mle$type <- NA
mle$type[olap$queryHits] <- olap$type

mle$expanded <- NA
mle$expanded[olap$queryHits] <- olap$expanded

dt.plot <- data.frame(t(table(mle$type, mle$repName)))
dt.plot <- dt.plot[dt.plot$Var1 %in% c("L1M", "L1MA", "L1MB", "L1MC", "L1MD", "L1ME", "L1P", "L1PA", "L1PB"), ]

names(dt.plot)[1] <- "LINE family"
dt.plot$percentage <- 0
dt.plot$percentage[dt.plot$Var2 == "stable"] <- 
  100*dt.plot$Freq[dt.plot$Var2 == "stable"]/sum(dt.plot$Freq[dt.plot$Var2 == "stable"])
dt.plot$percentage[dt.plot$Var2 == "variable"] <- 
  100*dt.plot$Freq[dt.plot$Var2 == "variable"]/sum(dt.plot$Freq[dt.plot$Var2 == "variable"])
ggplot(dt.plot, aes(x = Var2, y = percentage, fill = `LINE family`)) + xlab("") +
  geom_bar(stat = "identity", position = position_stack()) + theme_bw() + scale_fill_brewer(palette = "Spectral")
ggsave("output/figures/line.dis.top.png", width = 4, height = 5)

###

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
ggsave("output/figures/LINE.L1.sub.families.grouped.png", width = 5, height = 4)
write.table(dt.out, "output/L1.sub.families.tsv", sep = "\t", row.names=F, quote=F, col.names=T)
