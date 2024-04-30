setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(plotrix)
library(GenomicRanges)
library(scales)

### Mobile elements
####################
#### mobile elements
mle <- fread("../UCSC_hg38_repeatMasker.tsv", data.table = F)# 5,633,664
mle <- mle[mle$repFamily == "Alu" | mle$repClass == "LINE", ]# 2,870,078


exp <- kgp.str
dt.out <- data.frame()
p <- list()
set <- c()
ext = 5000
ext = 10000
ext = 20000

for(ext in c(0, 1000, 5000, 10000, 20000)){
  motif.portion <- data.frame()
  
  for(repClass in c("LINE", "SINE")){
    mle.tmp <- mle[mle$repClass == repClass, ]
    mle.tmp$genoStart <- ifelse(mle.tmp$genoStart - ext < 1, 1, mle.tmp$genoStart - ext)
    mle.tmp$genoEnd <- mle.tmp$genoEnd + ext
    
    mle.tmp.g <- GRanges(mle.tmp$genoName, IRanges(mle.tmp$genoStart, mle.tmp$genoEnd), "*")
    exp.g <- GRanges(exp$chr, IRanges(exp$start, exp$end), "*")
    
    olap <- findOverlaps(exp.g, mle.tmp.g)
    if(length(set) == 0){
      set <- olap@from
    }else{
      set <- intersect(set, olap@from)
    }
    if(ext == 0){
      
      motifs <- exp$dominant.motif[intersect(unique(olap@from), which(exp$type == "variable"))]
      # motifs <- exp$dominant.motif[intersect(unique(olap@from), which(exp$type == "variable" & exp$expanded == T))]
      motifs <- data.frame(table(motifs))
      motifs$percentage <- motifs$Freq/sum(motifs$Freq) * 100
      motifs <- motifs[order(motifs$Freq, decreasing = T), ]
      
      motif.tmp <- motifs[1:10, ]
      motif.tmp <- rbind(motif.tmp, data.frame("motifs" = "Other motif", "Freq" = NA, 
                                               "percentage" = 100 - sum(motif.tmp$percentage)))
      motif.tmp$type <- "all variable"
      motif.tmp$class <- repClass
      motif.portion <- rbind(motif.portion, motif.tmp)
      
      motifs$motif <- factor(motifs$motifs, levels = motifs$motifs)
      p[[length(p)+1]] <- ggplot(motifs[1:20, ], aes(x = motif, y = percentage))+ geom_bar(stat = "identity") +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        ggtitle(ifelse(repClass == "SINE", "Alu - variable", "LINE - variable")) + xlab("") +
        ylim(0,45)
      
      motifs <- exp$dominant.motif[intersect(unique(olap@from), which(exp$type == "variable" & exp$expanded == T))]
      motifs <- data.frame(table(motifs))
      motifs$percentage <- motifs$Freq/sum(motifs$Freq) * 100
      motifs <- motifs[order(motifs$Freq, decreasing = T), ]
      
      motif.tmp <- motifs[1:10, ]
      motif.tmp <- rbind(motif.tmp, data.frame("motifs" = "Other motif", "Freq" = NA, 
                                               "percentage" = 100 - sum(motif.tmp$percentage)))
      motif.tmp$type <- "expanded variable"
      motif.tmp$class <- repClass
      motif.portion <- rbind(motif.portion, motif.tmp)
      
      #### stable 
      
      motifs <- exp$dominant.motif[intersect(unique(olap@from), which(exp$type == "stable"))]
      
      # motifs <- exp$dominant.motif[intersect(unique(olap@from), which(exp$type == "stable" & exp$expanded == T))]
      motifs <- data.frame(table(motifs))
      motifs$percentage <- motifs$Freq/sum(motifs$Freq) * 100
      motifs <- motifs[order(motifs$Freq, decreasing = T), ]
      
      motif.tmp <- motifs[1:10, ]
      motif.tmp <- rbind(motif.tmp, data.frame("motifs" = "Other motif", "Freq" = NA, 
                                               "percentage" = 100 - sum(motif.tmp$percentage)))
      motif.tmp$type <- "all stable"
      motif.tmp$class <- repClass
      motif.portion <- rbind(motif.portion, motif.tmp)
      
      motifs$motif <- factor(motifs$motifs, levels = motifs$motifs)
      p[[length(p)+1]] <- ggplot(motifs[1:20, ], aes(x = motif, y = percentage))+ geom_bar(stat = "identity") +
        theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
        ggtitle(ifelse(repClass == "SINE", "Alu - stable", "LINE - stable")) + xlab("") +
        ylim(0,45)
      
      motifs <- exp$dominant.motif[intersect(unique(olap@from), which(exp$type == "stable" & exp$expanded == T))]
      motifs <- data.frame(table(motifs))
      motifs$percentage <- motifs$Freq/sum(motifs$Freq) * 100
      motifs <- motifs[order(motifs$Freq, decreasing = T), ]
      
      motif.tmp <- motifs[1:10, ]
      motif.tmp <- rbind(motif.tmp, data.frame("motifs" = "Other motif", "Freq" = NA, 
                                               "percentage" = 100 - sum(motif.tmp$percentage)))
      motif.tmp$type <- "expanded stable"
      motif.tmp$class <- repClass
      motif.portion <- rbind(motif.portion, motif.tmp)
    }
    tmp.exp <- exp[unique(olap@from), ]
    
    variable.count <- sum(tmp.exp$nmotif > 1, na.rm = T)
    stable.count <- sum(tmp.exp$nmotif == 1, na.rm = T)
    
    variable.denom <- sum(exp$nmotif > 1, na.rm = T) - variable.count
    stable.denom <- sum(exp$nmotif == 1, na.rm = T) - stable.count
    
    test <- fisher.test(data.frame("X1" = c(variable.count, stable.count),
                                   "X2" = c(variable.denom, stable.denom), stringsAsFactors = F))
    
    dt.out <- rbind(dt.out, data.frame("feature" = repClass, "extension" = ext,
                                       "OR" = round(test$estimate, digits = 2),
                                       "OR.lower" = round(test$conf.int[1], digits = 2),
                                       "OR.upper" = round(test$conf.int[2], digits = 2),
                                       "p" = signif(test$p.value, digits = 2), stringsAsFactors = F
    ))
  }
  
  motif.portion$class <- gsub("SINE", "Alu", motif.portion$class)
  motif.portion$group <- paste0(motif.portion$class, "-", motif.portion$type)
  if(ext == 0){
    ggplot(motif.portion, aes(x = group, y = percentage, fill = motifs)) +
      geom_bar(stat = "identity", position = "stack", color = "black") + ggtitle("LINE/Alu") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    # ggsave("../output/figures/aluline.motif.stack.png", width =12, height = 6)
    ggsave("output/figures/aluline.motif.stack.png", width =12, height = 6)
    
  }
}

cowplot::plot_grid(p[[1]], p[[3]], p[[2]], p[[4]], nrow = 2)
# ggsave("../output/figures/linealu.motif.with.expansions.png", width = 12, height = 6) ## add & expanded == T at two lines above
# ggsave("../output/figures/linealu.motif.with.png", width = 12, height = 6)
ggsave("../output/figures/linealu.motif.with.png", width = 12, height = 6)

getwd()

dt.out$extension <- paste0(dt.out$extension/1000, "KB")
# dt.out$extension[1:2] <- "no extension"
dt.out$extension <- factor(dt.out$extension, levels = unique(dt.out$extension))
dt.out$fdr <- signif(p.adjust(dt.out$p, method = "BH"), digits = 2)
dt.out$fdr.range <- ifelse(dt.out$fdr < 0.25, "FDR 25%", "Not significant")
dt.out$fdr.range[dt.out$fdr < 0.20] <- "FDR 20%"
dt.out$fdr.range[dt.out$fdr < 0.15] <- "FDR 15%"
dt.out$fdr.range[dt.out$fdr < 0.10] <- "FDR 10%"
dt.out$fdr.range[dt.out$fdr < 0.05] <- "FDR 5%"
dt.out$fdr.range[dt.out$fdr < 0.01] <- "FDR 1%"
dt.out$fdr.range <- factor(dt.out$fdr.range, levels = c("FDR 1%", "FDR 5%", "FDR 10%",
                                                        "FDR 15%", "FDR 20%", "FDR 25%", "Not significant"))
# ggplot(dt.out, aes(x = extension, y = OR, color = fdr.range, group = feature)) + 
#   geom_hline(yintercept = 1, lty = 2) +
#   geom_point(aes(shape = feature), position = position_dodge2(width = .5), size = 3) + 
#   geom_errorbar(aes(ymin = OR.lower, ymax = OR.upper), width =.3, position = position_dodge(width = .5)) + 
#   theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0)) +
#   ylab("Odds ratio (variable/stable)") +
#   scale_color_manual(values = c("red", "orange", "yellow", "grey"))
# 
# # ggsave("../output/figures/mobile.element.png", width = 6, height = 4)
# ggsave("output/figures/mobile.element.png", width = 6, height = 4)

# write.table(dt.out, "output/mobile.element.tsv", sep = "\t", row.names=F, quote=F, col.names=T)

### KGP
dt.out = read.delim('../_1KGP/output/mobile.element.tsv', stringsAsFactors = F)
dt.out$fdr.range = as.character(dt.out$fdr.range)
dt.out$col = unlist(sapply(dt.out$fdr.range, switch, 'Not significant'='grey', 'FDR 25%'='grey', 'FDR 20%'='grey', 'FDR 15%'='grey',
                           'FDR 10%'='yellow','FDR 5%'='orange', 'FDR 1%'='red'))
table(dt.out$fdr.range)
genic.plot = dt.out
genic.plot$pch = ifelse(genic.plot$feature == 'LINE', 17, 19)
genic.plot$pch = ifelse(genic.plot$feature == 'LINE', 24, 21)

genic.plot$move = sapply(genic.plot$feature, switch, "LINE" = 0.2, "SINE" = (-0.2))
genic.plot$extension[1:2] = "0   "

pdf('output/figures/mobile_elements_kgp.pdf', w = 6, h = 4)#, h = 7
par(mar = c(9,6,1,8), mgp = c(3,1,0), mfrow=c(1,1)) 
plotCI((1:10)+genic.plot$move, genic.plot$OR, ylab = 'OddsRatio (variable vs stable)', xlab = 'Extension',
       main = '', ui = genic.plot$OR.upper, li = genic.plot$OR.lower, pch = NA,
       lwd = 1, xaxt="n", ylim = c(0,3), col = 'black', cex = 1, yaxt = 'n', )
points((1:10)+genic.plot$move, genic.plot$OR, pch = genic.plot$pch, bg = genic.plot$col, cex = 1.5, col = 'black')

# plot(1:25, pch = c(1:25))
axis(1, at=seq(1, 10, by=2)+0.5, labels = FALSE)
axis(2, at=seq(0, 5, by=1), labels = TRUE, las = 2)
abline(h = 1, lwd = 0.5, lty = 2)

text((1:5)*2, par("usr")[3]-0.5, srt = 0, adj = 1, xpd = TRUE, cex = 1, 
     labels = paste(unique(genic.plot$extension)))
legend('topright', pch=19, pt.cex=1.3, inset=c(-0.45,0), xpd = TRUE, bty = "n", col = c('red','yellow','grey'), 
       legend=c( "FDR 1%","FDR 10%", "Not significant"), title = 'FDR range')
legend('topright', pch=21, pt.cex=1.3, inset=c(-0.45,0), xpd = TRUE, bty = "n", col = 'black',
       bg = c('red','yellow','grey'), legend=c( "FDR 1%","FDR 10%", "Not significant"), title = '')

legend('bottomright', pch=c(17,19), pt.cex=1.3, inset=c(-0.3,0), xpd = TRUE, bty = "n", col = 'grey', 
       legend=c("LINE","SINE"), title = 'Class')
legend('bottomright', pch=c(24,21), pt.cex=1.3, inset=c(-0.3,0), xpd = TRUE, bty = "n", col = 'black',
       bg = 'grey', legend=c( "LINE","SINE"), title = '')

dev.off()