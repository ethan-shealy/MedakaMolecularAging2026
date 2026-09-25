library(data.table)
library(tidyverse)
library(ggpubr)

setwd(".")

theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=7),
                  axis.title=element_text(size=9,face="bold"),
                  strip.text = element_text(size=9,face="bold"),
                  panel.border = element_rect(fill = NA)))

load("./Code_Data/mus_NFIXregions_RM.rda")
load("./Code_Data/medaka_NFIXregions_RM.rda")

### Process mouse data

mus.samples <- ncol(mus_allNFIX[,-1])
mus.ages <- as.numeric(unlist(mus_allNFIX[1,-ncol(mus_allNFIX)]))
mus.IDs <- colnames(mus_allNFIX[,-1])
mus.genes <- mus_allNFIX$gene[-1]
mus.windows <- nrow(mus_allNFIX[-1,])


mus_allNFIX <- mus_allNFIX[-1,-ncol(mus_allNFIX)]

#### Repeat for medaka

medaka.samples <- ncol(medaka_allNFIX[,-1])
medaka.ages <- as.numeric(unlist(medaka_allNFIX[1,-ncol(medaka_allNFIX)]))
medaka.IDs <- colnames(medaka_allNFIX[,-1])
medaka.genes <- medaka_allNFIX$gene[-1]
medaka.windows <- nrow(medaka_allNFIX[-1,])


medaka_allNFIX <- medaka_allNFIX[-1,-ncol(medaka_allNFIX)]

meta <- data.frame("Species" = rep(c("mouse", "medaka"), times = c(ncol(mus_allNFIX), ncol(medaka_allNFIX))),
                   "ID" = c(mus.IDs, medaka.IDs),
                   "Ages" = c(mus.ages, medaka.ages),
                   "Gene" = "NFIX(b)")

sampleAvgs <- list(NA)
windowAvgs <- list(NA)

mus.sampleAvgs <- colMeans(mus_allNFIX, na.rm = TRUE)
medaka.sampleAvgs <- colMeans(medaka_allNFIX, na.rm = TRUE)


meta$NFIXAvg <- c(mus.sampleAvgs, medaka.sampleAvgs)

## plot simple average across all windows
library(ggpubr)
library(tidyverse)
library(grid)
library(egg)


#### Get gene coordinates
library(psych)
library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(txdbmaker)
library(Gviz)

options(ucscChromosomeNames = FALSE)

gtrack <- GenomeAxisTrack()




## medaka

medaka.txdb <- makeTxDbFromGFF("./Code_Data/medaka.gff_nfix.gff")

txTrack <- GeneRegionTrack(
  medaka.txdb, name = "NCBI RefSeq Transcripts",
  transcriptAnnotation = "transcript", 
  chromosome = medaka.txdb$user_seqlevels,
  start = min(start(transcripts(medaka.txdb))) - 2000,
  end = max(end(transcripts(medaka.txdb))) + 2000
)

txTrack <- txTrack[which(txTrack@range@elementMetadata@listData$transcript %in% c("XM_020705512.2", 
                                                                                  "XM_020705517.2",
                                                                                  "XM_020705522.2"))]

displayPars(txTrack) <- list(fontcolor.title = "black", just.group = "above",
                             col = "black", cex.group = 0.8)

cor.medaka <- corr.test(t(medaka_allNFIX), y = meta$Ages[which(meta$Species == "medaka")], method = "spearman")

medaka.coord <- str_split_fixed(rownames(medaka_allNFIX), ":", 3)

medaka.ranges <- GRanges(medaka.coord[,1], 
                         IRanges(as.numeric(medaka.coord[,2]),
                                 as.numeric(medaka.coord[,3])))

dtrack <- DataTrack(data = as.numeric(cor.medaka$r), start = start(medaka.ranges), genome = medaka.txdb$user_seqlevels,
                    end = end(medaka.ranges), chromosome = seqnames(medaka.ranges),
                    name = "Rho")

dtrack <- dtrack[,-which(dtrack@range@ranges@start < min(start(transcripts(medaka.txdb))) - 10000 | 
                           dtrack@range@ranges@start > max(end(transcripts(medaka.txdb))) + 10000)]


displayPars(dtrack) <- list(aggregateGroups = FALSE, baseline = 0, col.baseline = "black", col.histogram = "blue", fill.histogram = "blue",
                            fontcolor.title = "black", col.axis = "black", ylim = c(-1, 1))

plotTracks(list(gtrack, dtrack, txTrack), type = "histogram", sizes = c(0.2, 0.8, 1.2), 
           collapseTranscripts = FALSE, title.width = 0.7, reverseStrand = FALSE)


map.1 <- grid.grab()

## draw it, changes optional
map.1.new <- editGrob(map.1, vp = viewport(width = unit(12, "in"), 
                                           height = unit(6, "in")))
grid.newpage()
grid.draw(map.1.new)

ggsave("./Figures/7A.1.svg",
       plot = map.1.new, device = "svg", width = 12, height = 6, units = "in")



### Plot specific promoter meth changes

unique_transcripts <- ranges(txTrack)[which(!duplicated(ranges(txTrack)$transcript))]
unique_proms <- promoters(unique_transcripts, upstream = 1000, downstream = 1000)

medaka.overlap <- findOverlaps(medaka.ranges, unique_proms)

overlaps.which <- data.frame("From" = medaka.ranges@ranges@start[medaka.overlap@from],
                             "To" = unique_proms@elementMetadata$transcript[medaka.overlap@to])

getValues <- data.frame(meta[which(meta$Species == "medaka"),], 
                        t(medaka_allNFIX[medaka.overlap@from,]))

getValues.long <- pivot_longer(getValues, cols = starts_with("NC_"), 
                               names_to = "CpG", values_to = "Methylation")

getValues.long$Start <- str_split_fixed(getValues.long$CpG, "\\.", 4)[,3]
getValues.long$Prom <- overlaps.which$To[match(getValues.long$Start, overlaps.which$From)]

#ggplot(getValues.long, aes(x = Ages, y = Methylation)) + geom_point() + geom_smooth() + 
  #facet_wrap(~Prom)

getValues.sum <- getValues.long %>% group_by(ID, Prom) %>% 
                    summarise(MethAvg = mean(Methylation, na.rm = TRUE))

getValues.sum$Ages <- getValues.long$Ages[match(getValues.sum$ID, getValues.long$ID)]

B1 <- ggplot(getValues.sum, aes(x = Ages, y = MethAvg)) + 
        geom_point() + 
        geom_smooth(method = "lm") + 
        facet_wrap(~Prom) + ylim(0, 100) + 
        ylab("Methylation Percent")

grid.newpage()
grid.draw(set_panel_size(B1, width  = unit(5, "in"), height = unit(3, "in")))


ggsave("./Figures/7B.1.svg",
       plot = B1, device = "svg", width = 5, height = 3, units = "in")


## mouse

mus.txdb <- makeTxDbFromGFF("./Code_Data/mus.gff_nfix.gff")

txTrack <- GeneRegionTrack(
  mus.txdb, name = "NCBI RefSeq Transcripts",
  transcriptAnnotation = "transcript", 
  chromosome = mus.txdb$user_seqlevels,
  start = min(start(transcripts(mus.txdb))) - 2000,
  end = max(end(transcripts(mus.txdb))) + 2000
)

txTrack <- txTrack[which(txTrack@range@elementMetadata@listData$transcript %in% c("NM_001371052.1", 
                                                                                  "NM_001371053.1"))]

displayPars(txTrack) <- list(fontcolor.title = "black", just.group = "above",
                             col = "black", cex.group = 0.8)

cor.mus <- corr.test(t(mus_allNFIX), y = meta$Ages[which(meta$Species == "mouse")], method = "spearman")

mus.coord <- str_split_fixed(rownames(mus_allNFIX), ":", 3)

mus.ranges <- GRanges(mus.coord[,1], 
                      IRanges(as.numeric(mus.coord[,2]),
                              as.numeric(mus.coord[,3])))

dtrack <- DataTrack(data = as.numeric(cor.mus$r), start = start(mus.ranges), genome = mus.txdb$user_seqlevels,
                    end = end(mus.ranges), chromosome = seqnames(mus.ranges),
                    name = "Rho")

dtrack <- dtrack[,-which(dtrack@range@ranges@start < min(start(transcripts(mus.txdb))) - 10000 | 
                           dtrack@range@ranges@start > max(end(transcripts(mus.txdb))) + 10000)]

#txTrack <- txTrack[-which(grepl("XM_", txTrack@range$symbol))]

displayPars(dtrack) <- list(aggregateGroups = FALSE, baseline = 0, col.baseline = "black", col.histogram = "blue", fill.histogram = "blue",
                            fontcolor.title = "black", col.axis = "black", ylim = c(-1, 1))

plotTracks(list(gtrack, dtrack, txTrack), type = "histogram", sizes = c(0.2, 0.8, 1.2), 
           collapseTranscripts = FALSE, title.width = 0.7, reverseStrand = TRUE)

map.2 <- grid.grab()

## draw it, changes optional
map.2.new <- editGrob(map.2, vp = viewport(width = unit(12, "in"), 
                                             height = unit(6, "in")))
grid.newpage()
grid.draw(map.2.new)

ggsave("./Figures/7A.2.svg",
       plot = map.2.new, device = "svg", width = 12, height = 6, units = "in")



### Plot specific promoter meth changes

#unique_transcripts <- ranges(txTrack)[which(!duplicated(ranges(txTrack)$transcript))]
#unique_proms <- promoters(unique_transcripts, upstream = 1000, downstream = 1000)

Inter <- GRanges("NC_000074.7", IRanges(85510000, 85515000))

#unique_proms <- c(unique_proms, Inter)

mus.overlap <- findOverlaps(mus.ranges, Inter)

overlaps.which <- data.frame("From" = mus.ranges@ranges@start[mus.overlap@from],
                             "To" = "Inter")

#overlaps.which$To[which(is.na(overlaps.which$To))] <- "Inter"

getValues <- data.frame(meta[which(meta$Species == "mouse"),], 
                        t(mus_allNFIX[mus.overlap@from,]))

getValues.long <- pivot_longer(getValues, cols = starts_with("NC_"), 
                               names_to = "CpG", values_to = "Methylation")

#getValues.long$Start <- str_split_fixed(getValues.long$CpG, "\\.", 4)[,3]
getValues.long$Prom <- "Inter"

#ggplot(getValues.long, aes(x = Ages, y = Methylation)) + geom_point() + geom_smooth() + 
#  facet_wrap(~Prom)

getValues.sum <- getValues.long %>% group_by(ID, Prom) %>% 
  summarise(MethAvg = mean(Methylation, na.rm = TRUE))

getValues.sum$Ages <- getValues.long$Ages[match(getValues.sum$ID, getValues.long$ID)]

B2 <- ggplot(getValues.sum, aes(x = Ages, y = MethAvg)) + 
        geom_point() + 
        geom_smooth(method = "lm", formula = y ~ log(x)) + 
        facet_wrap(~Prom) + #ylim(0, 100) + 
        ylab("Methylation Percent")

grid.newpage()
grid.draw(set_panel_size(B2, width  = unit(3, "in"), height = unit(3, "in")))


ggsave("./Figures/7B.2.svg",
       plot = B2, device = "svg", width = 3, height = 3, units = "in")

