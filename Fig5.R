library(tidyverse)
library(ggpubr)
library(data.table)
library(grid)
library(egg)

setwd("")


theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=6),
                  axis.title=element_text(size=9,face="bold"),
                  strip.text = element_text(size=9,face="bold"),
                  panel.border = element_rect(fill = NA)))


raw <- as.data.frame(fread("./unfiltered_genecounts.csv", check.names = FALSE))

meta <- readxl::read_xlsx("./Code_Data/2022 Aging cohort.xlsx", 
                          sheet = "Dissect", na = c("N/A", "NA", ""))[,-c(15:25)]


meta$ID <- paste0("X", meta$ID)

raw.meta <- raw[,c(1:3)]

raw.meta$V1 <- str_remove_all(raw.meta$V1, "ID")
raw.meta$V1 <- str_remove_all(raw.meta$V1, "\\.bam")

raw.meta$V1 <- paste0("X", raw.meta$V1)

raw.dat <- raw[,-c(1:3)]

rownames(raw.dat) <- raw.meta$V1

cpm.dat <- t(edgeR::cpm(t(raw.dat), log = FALSE))

cpm.filt.means <- colMeans(cpm.dat)

cpm.dat <- cpm.dat[,which(cpm.filt.means > 0.1)]

allGeneNames <- colnames(cpm.dat)


contAgeGenes <- readxl::read_excel("./Code_data/DESEQ2_all.xlsx")

library(GenomicRanges)
library(rtracklayer)

annot_full <- import.gff("./Code_Data/ASM223467v1.gff")

annot <- annot_full[which(annot_full$type == "mRNA"),]

contAgeGenes.annot <- annot[na.omit(match(contAgeGenes$Gene, annot$gene))]
contAgeGenes.annot$score <- contAgeGenes$log2FoldChange[match(contAgeGenes.annot$gene, contAgeGenes$Gene)]

contAgeGenes.bed <- as.data.frame(contAgeGenes.annot)[,c(1,2,3,8)]


### now overlap of all CpGs

CGsites <- fread("./5x_50ind_manualFilt.tab", sep = "\t", header = TRUE)

CGsites <- str_split_fixed(CGsites$coord, pattern = ":", 2)

CGsites.clean <- GRanges(seqnames = CGsites[,1], 
                         ranges = IRanges(start = as.numeric(CGsites[,2]), end = as.numeric(CGsites[,2])), strand = "*")

length(CGsites.clean) ### total measured CpGs = 2,803,340

rm(CGsites)

## now get age-related sites

ageDMCs <- read.table("./DSS/SigLogAgeCpGs.tab", sep = "\t", header = TRUE)

ageDMCs.bed <- data.frame(ageDMCs[,c(1, 2, 2)], ageDMCs$pvals*sign(ageDMCs$stat))

ageDMCs.gr <- GRanges(ageDMCs$chr, IRanges(ageDMCs$pos, ageDMCs$pos), strand = "*")

length(ageDMCs.gr) ## total ageDMCs = 23,776

### Basic test

aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.annot, maxgap = 2000L, type = "any", select = "first")

## set up data frame for plotting

chrSizes <- read.table("./Code_Data/ASM223467v1_chromSizes.txt")

chrSizes$cum <- c(0, cumsum(chrSizes$V2)[-length(chrSizes$V2)])


geneFrame <- data.frame("Name" = contAgeGenes.annot$Name,
                        "Chr" = seqnames(contAgeGenes.annot),
                        "Start" = start(contAgeGenes.annot),
                        "End" = end(contAgeGenes.annot),
                        "Overlap" = "Gene")

aDMCframe <- data.frame("Name" = paste0(seqnames(ageDMCs.gr), ":", start(ageDMCs.gr)),
                        "Chr" = seqnames(ageDMCs.gr),
                        "Start" = start(ageDMCs.gr),
                        "End" = end(ageDMCs.gr),
                        "Overlap" = ifelse(is.na(aDMCoverlap), "NotNear_aDEG", "Near_aDEG"))

table(aDMCframe$Overlap) ## 1,753 aDMCs are within 2kb of an aDEG (out of 23,776, or ~7.37%)

comboFrame <- rbind(geneFrame, aDMCframe)

comboFrame$ChrOffset <- chrSizes$cum[match(comboFrame$Chr, chrSizes$V1)] + comboFrame$Start

comboFrame$Weight <- ifelse(comboFrame$Overlap == "Gene", 1, 
                            abs(ageDMCs$stat[match(make.names(comboFrame$Name[-which(comboFrame$Overlap == "Gene")]), 
                                               ageDMCs$name)]))

library(ggridges)

comboFrame$Type <- ifelse(comboFrame$Overlap == "Gene", "Gene", "CpG")
comboFrame$ChrLength <- chrSizes$V2[match(comboFrame$Chr, chrSizes$V1)]
comboFrame$RelChromPos <- comboFrame$Start / comboFrame$ChrLength

a <- ggplot(comboFrame, aes(x = RelChromPos, y = Chr, fill = Type)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.01, alpha = 0.5) +
  scale_fill_viridis_d() +
  geom_point(data = filter(comboFrame, Overlap == "Near_aDEG"), color = "black", 
             position = position_nudge(y = 0.2), alpha = 0.2, shape = 25, show.legend = FALSE) +
  labs(x = "Chromosomal Position", y = "Chromosome", fill = "Type") +
  theme(panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_blank()) + xlim(0, 1) +
  coord_cartesian(xlim = c(0, 1))

grid.newpage()
grid.draw(set_panel_size(a, width  = unit(2.8, "in"), height = unit(7, "in")))


ggsave("./Figures/5B.svg",
       plot = a, device = "svg", width = 3.8, height = 7, units = "in")


overlapDMCs <- aDMCframe[which(aDMCframe$Overlap == "Near_aDEG"),]
overlapDMCs$Gene <- contAgeGenes.annot$gene[na.omit(aDMCoverlap)]

fwrite(overlapDMCs, "./DSS/overlapDMCs.tsv")

overlapDEGs <- table(overlapDMCs$Gene) %>% as.data.frame() %>% sort_by(~Freq, decreasing = TRUE)

fwrite(overlapDEGs, "./DSS/overlapDEGs.tsv")

#### Enrichment of overlap DEGs
library(gprofiler2)
set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e112_eg59_p19") ### Retreive archived version to get same results
get_version_info()

allgenes.ensg <- gconvert(colnames(raw.dat), organism = "olatipes")

overlapGenes.ensg <- gconvert(overlapDEGs$Var1, organism = "olatipes")

overlapGenes.go <- gost(overlapGenes.ensg$target, organism = "olatipes", 
                        custom_bg = allgenes.ensg$target, ordered_query = TRUE)

overlap.go.res <- overlapGenes.go$result

overlap.go.res$ordered <- c(1:nrow(overlap.go.res))

go.t <- ggplot(overlap.go.res, aes(x = -log10(p_value), y = reorder(term_name, ordered, decreasing = TRUE), color = source)) + 
  geom_point() + ylab("") + theme(legend.position = "top")

grid.newpage()
grid.draw(set_panel_size(go.t, width  = unit(3.2, "in"), height = unit(2.8, "in")))


ggsave("./Figures/5C.svg",
       plot = go.t, device = "svg", width = 3.6, height = 3.3, units = "in")



### Set up loop for testing by distance

distances <- c(100000L, 50000L, 20000L, 10000L, 5000L, 2000L, 1000L, 0L)

tests <- list()

for (i in 1:length(distances)){

    allCGoverlap <- findOverlaps(CGsites.clean, contAgeGenes.annot, maxgap = distances[i], type = "any", select = "first")
    
    ### now overlap age-related sites
    
    aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.annot, maxgap = distances[i], type = "any", select = "first")

    ### set up fisher test
    
    compTable <- data.frame("NonAge" = c(length(CGsites.clean) - length(na.omit(allCGoverlap)), 
                                         length(na.omit(allCGoverlap))),
                            "aDMC" = c(length(ageDMCs.gr) - length(na.omit(aDMCoverlap)),
                                       length(na.omit(aDMCoverlap))))
    rownames(compTable) <- c("NotNearaDEG", "NearaDEG")
    
    tests[[i]] <- fisher.test(compTable)        
}

tests[[6]]

ORs <- unlist(lapply(1:length(tests), function(x) tests[[x]]$estimate))
lowers <- unlist(lapply(1:length(tests), function(x) tests[[x]]$conf.int[1]))
uppers <- unlist(lapply(1:length(tests), function(x) tests[[x]]$conf.int[2]))

## 

distanceTable <- data.frame("Distance" = c("Background", "100kb", "50kb", "20kb", "10kb", "5kb", "2kb", "1kb", "Within"),
                            "OddsRatio" = c(1, ORs),
                            "Lower" = c(1, lowers),
                            "Upper" = c(1, uppers))

distanceTable$Distance <- factor(distanceTable$Distance, levels = distanceTable$Distance)
                        
b <- ggplot(distanceTable, aes(x = Distance, y = OddsRatio)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5) + 
  labs(x = "Maximum Distance from aDEG", y = "Odds ratio aDMCs vs all CpGs") + 
  scale_y_continuous(breaks = seq(0, 1.5, 0.25), limits = c(0, 1.5))

grid.newpage()
grid.draw(set_panel_size(b, width  = unit(3.2, "in"), height = unit(2.8, "in")))


ggsave("./Figures/5A.svg",
       plot = b, device = "svg", width = 3.6, height = 2, units = "in")



##### seperate into increasing and decreasing aDEGs

contAgeGenes.annot$logFC <- contAgeGenes$log2FoldChange[match(contAgeGenes.annot$gene, contAgeGenes$Gene)]
contAgeGenes.annot$Dir <- as.factor(sign(contAgeGenes.annot$logFC))

contAgeGenes.annot.dec <- contAgeGenes.annot[which(contAgeGenes.annot$logFC < 0)]
contAgeGenes.annot.inc <- contAgeGenes.annot[which(contAgeGenes.annot$logFC > 0)]

#### Now same analysis on decreasing genes only
#
#tests <- list()
#
#for (i in 1:length(distances)){
#  
#  allCGoverlap <- findOverlaps(CGsites.clean, contAgeGenes.annot.dec, maxgap = distances[i], type = "any", select = "first")
#  
#  ### now overlap age-related sites
#  
#  aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.annot.dec, maxgap = distances[i], type = "any", select = "first")
#  
#  ### set up fisher test
#  
#  compTable <- data.frame("NonAge" = c(length(CGsites.clean) - length(na.omit(allCGoverlap)), 
#                                       length(na.omit(allCGoverlap))),
#                          "aDMC" = c(length(ageDMCs.gr) - length(na.omit(aDMCoverlap)),
#                                     length(na.omit(aDMCoverlap))))
#  rownames(compTable) <- c("NotNearaDEG", "NearaDEG")
#  
#  tests[[i]] <- fisher.test(compTable)        
#}
#
#ORs.dec <- unlist(lapply(1:length(tests), function(x) tests[[x]]$estimate))
#lowers.dec <- unlist(lapply(1:length(tests), function(x) tests[[x]]$conf.int[1]))
#uppers.dec <- unlist(lapply(1:length(tests), function(x) tests[[x]]$conf.int[2]))
#
#
### repeat for increasing genes
#
#tests <- list()
#
#for (i in 1:length(distances)){
#  
#  allCGoverlap <- findOverlaps(CGsites.clean, contAgeGenes.annot.inc, maxgap = distances[i], type = "any", select = "first")
#  
#  ### now overlap age-related sites
#  
#  aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.annot.inc, maxgap = distances[i], type = "any", select = "first")
#  
#  ### set up fisher test
#  
#  compTable <- data.frame("NonAge" = c(length(CGsites.clean) - length(na.omit(allCGoverlap)), 
#                                       length(na.omit(allCGoverlap))),
#                          "aDMC" = c(length(ageDMCs.gr) - length(na.omit(aDMCoverlap)),
#                                     length(na.omit(aDMCoverlap))))
#  rownames(compTable) <- c("NotNearaDEG", "NearaDEG")
#  
#  tests[[i]] <- fisher.test(compTable)        
#}
#
#ORs.inc <- unlist(lapply(1:length(tests), function(x) tests[[x]]$estimate))
#lowers.inc <- unlist(lapply(1:length(tests), function(x) tests[[x]]$conf.int[1]))
#uppers.inc <- unlist(lapply(1:length(tests), function(x) tests[[x]]$conf.int[2]))
#
#
### 
#
#distanceTable.aDEGdir <- data.frame("Distance" = c("Background", c("20kb", "10kb", "5kb", "2kb", "1kb", "500bp", "200bp", "Within"),
#                                                   "Background", c("20kb", "10kb", "5kb", "2kb", "1kb", "500bp", "200bp", "Within")),
#                            "OddsRatio" = c(1, ORs.dec, 1, ORs.inc),
#                            "Lower" = c(1, lowers.dec, 1, lowers.inc),
#                            "Upper" = c(1, uppers.dec, 1, uppers.inc),
#                            "aDEGdir" = c(rep("Decreasing", length(distances)+1), rep("Increasing", length(distances)+1)))
#
#distanceTable.aDEGdir$Distance <- factor(distanceTable.aDEGdir$Distance, levels = distanceTable$Distance)
#
#ggplot(distanceTable.aDEGdir, aes(x = Distance, y = OddsRatio, fill = aDEGdir)) + 
#  geom_col(position = "dodge") + facet_wrap(~aDEGdir, scales = "free_x") +
#  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5, position = position_dodge(width = 1)) + 
#  labs(x = "Distance from aDEG", y = "Odds ratio aDMCs vs all CpGs")









#### Now look at likelihood ratios between increasing and decreasing aDMCs and aDEGs

ageDMCs.gr.dec <- ageDMCs.gr[which(ageDMCs$stat < 0)]
ageDMCs.gr.inc <- ageDMCs.gr[which(ageDMCs$stat > 0)]


## Decreasing aDEG, decreasing aDMC

tests <- list()

for (i in 1:length(distances)){
  
  ### now overlap age-related sites
  
  NegaDMCoverlapNegGene <- findOverlaps(ageDMCs.gr.dec, contAgeGenes.annot.dec, maxgap = distances[i], type = "any", select = "first")
  NegaDMCoverlapPosGene <- findOverlaps(ageDMCs.gr.dec, contAgeGenes.annot.inc, maxgap = distances[i], type = "any", select = "first")
  PosaDMCoverlapNegGene <- findOverlaps(ageDMCs.gr.inc, contAgeGenes.annot.dec, maxgap = distances[i], type = "any", select = "first")
  PosaDMCoverlapPosGene <- findOverlaps(ageDMCs.gr.inc, contAgeGenes.annot.inc, maxgap = distances[i], type = "any", select = "first")
  
  ### set up fisher test
  
  compTable <- data.frame("NegaDMC" = c(length(na.omit(NegaDMCoverlapNegGene)), length(na.omit(NegaDMCoverlapPosGene))),
                          "PosaDMC" = c(length(na.omit(PosaDMCoverlapNegGene)), length(na.omit(PosaDMCoverlapPosGene))))
  
  rownames(compTable) <- c("NegaDEG", "PosaDEG")
  
  tests[[i]] <- fisher.test(compTable)        
}

ORs.dir <- unlist(lapply(1:length(tests), function(x) tests[[x]]$estimate))
lowers.dir <- unlist(lapply(1:length(tests), function(x) tests[[x]]$conf.int[1]))
uppers.dir <- unlist(lapply(1:length(tests), function(x) tests[[x]]$conf.int[2]))


distanceTable.bothDir <- data.frame("Distance" = c("Background", "100kb", "50kb", "20kb", "10kb", "5kb", "2kb", "1kb", "Within"),
                                    "OddsRatio" = c(1, ORs.dir),
                                    "Lower" = c(1, lowers.dir),
                                    "Upper" = c(1, uppers.dir),
                                    "Dir" = ifelse(c(1, ORs.dir) < 0, "", ""))

distanceTable.bothDir$Distance <- factor(distanceTable.bothDir$Distance, levels = distanceTable$Distance)

ggplot(distanceTable.bothDir, aes(x = Distance, y = OddsRatio)) + 
  geom_col() + #ylim(-0.08, 0.08) + 
  geom_errorbar(aes(ymin = Lower, ymax = Upper)) + 
  labs(x = "Distance from aDEG", y = "Odds that aDMC and aDEG Direction Match") + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0))


### making sure this is correct - get data frame of overlaps with context info

aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.annot, maxgap = 2000L, type = "any", select = "first")

totalGeneTable <- data.frame("CpG" = paste0(ageDMCs$chr, ":", ageDMCs$pos),
                             "ageDir" = as.factor(sign(ageDMCs$stat)),
                             "Closest_aDEG" = contAgeGenes.annot$gene[aDMCoverlap],
                             "Closest_aDEG_Dir" = contAgeGenes.annot$Dir[aDMCoverlap])

countMatrix <- as.matrix(table(totalGeneTable$ageDir, totalGeneTable$Closest_aDEG_Dir))

colnames(countMatrix) <- c("Negative_aDEGs", "Positive_aDEGs")
rownames(countMatrix) <- c("Negative_aDMCs", "Positive_aDMCs")

library(ComplexHeatmap)



Heatmap(countMatrix, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "black"),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(countMatrix[i, j], x, y, gp = gpar(fontsize = 12))},
        name = "Count", show_heatmap_legend = FALSE, col = c("grey80", "yellow"))


c <- grid.grab()

## draw it, changes optional
c.new <- editGrob(c, vp = viewport(width = unit(3.8, "in"), 
                                             height = unit(3, "in")))
grid.newpage()
grid.draw(c.new)

ggsave("./Figures/5D.svg",
       plot = c.new, device = "svg", width = 3.8, height = 3, units = "in")




fisher.test(table(totalGeneTable$Closest_aDEG_Dir, totalGeneTable$ageDir))


linked_aDMCs <- subsetByOverlaps(ageDMCs.gr, contAgeGenes.annot, maxgap = 2000L, type = "any")

### Only count each gene once

totalGeneTable2 <- totalGeneTable[-which(duplicated(totalGeneTable$Closest_aDEG)),]

countMatrix2 <- as.matrix(table(totalGeneTable2$ageDir, totalGeneTable2$Closest_aDEG_Dir))

colnames(countMatrix2) <- c("Negative_aDEGs", "Positive_aDEGs")
rownames(countMatrix2) <- c("Negative_aDMCs", "Positive_aDMCs")


Heatmap(countMatrix2, cluster_rows = FALSE, cluster_columns = FALSE, rect_gp = gpar(col = "black"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(countMatrix2[i, j], x, y, gp = gpar(fontsize = 12))},
        name = "Count", show_heatmap_legend = FALSE, col = c("grey80", "grey60"))

fisher.test(table(totalGeneTable2$Closest_aDEG_Dir, totalGeneTable2$ageDir))



#### Within a DEGs, is there an enrichment of aDMCs in specific places?
annot_CDS <- annot_full[which(annot_full$type == "gene"),]
annot_CDS <- annot_CDS[which(annot_CDS$gene %in% colnames(cpm.dat))]


annot_exons <- annot_full[which(annot_full$type == "exon" & annot_full$gene %in% colnames(cpm.dat)),]

contAgeGenes.exons <- annot_exons[which(annot_exons$gene %in% contAgeGenes$Gene)]
contAgeGenes.exons$score <- contAgeGenes$log2FoldChange[match(contAgeGenes.exons$gene, contAgeGenes$Gene)]

#export.bed(contAgeGenes.exons, "~/exons.bed")

allGenes.introns <- setdiff(annot_CDS, annot_exons)

contAgeGenes.introns <- setdiff(contAgeGenes.annot, contAgeGenes.exons)

#export.bed(contAgeGenes.introns, "~/introns.bed")

allGenes.promoters <- promoters(annot_CDS, upstream = 2000, downstream = 2000)

contAgeGenes.promoter <- promoters(contAgeGenes.annot, upstream = 2000, downstream = 2000)

allGenes.terminators <- terminators(annot_CDS, upstream = 2000, downstream = 2000)

contAgeGenes.terminators <- terminators(contAgeGenes.annot, upstream = 2000, downstream = 2000)

#### first enrichment in exons

allCGoverlap <- findOverlaps(CGsites.clean, contAgeGenes.exons, maxgap = 0L, type = "any", select = "first")

numExonCpGs <- sum(!is.na(allCGoverlap))

### now overlap age-related sites

aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.exons, maxgap = 0L, type = "any", select = "first")

numExonDMCs <- sum(!is.na(aDMCoverlap))

### set up fisher test

compTable <- data.frame("NonAge" = c(length(CGsites.clean) - length(na.omit(allCGoverlap)), 
                                     length(na.omit(allCGoverlap))),
                        "aDMC" = c(length(ageDMCs.gr) - length(na.omit(aDMCoverlap)),
                                   length(na.omit(aDMCoverlap))))
rownames(compTable) <- c("NotWithinContext", "WithinContext")

compTable

fisher.test(compTable) ## Exon enrichment: OR = 1.27


#### enrichment in introns

allCGoverlap <- findOverlaps(CGsites.clean, contAgeGenes.introns, maxgap = 0L, type = "any", select = "first")

numIntronCpGs <- sum(!is.na(allCGoverlap))

### now overlap age-related sites

aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.introns, maxgap = 0L, type = "any", select = "first")

numIntronDMCs <- sum(!is.na(aDMCoverlap))

### set up fisher test

compTable <- data.frame("NonAge" = c(length(CGsites.clean) - length(na.omit(allCGoverlap)), 
                                     length(na.omit(allCGoverlap))),
                        "aDMC" = c(length(ageDMCs.gr) - length(na.omit(aDMCoverlap)),
                                   length(na.omit(aDMCoverlap))))
rownames(compTable) <- c("NotWithinContext", "WithinContext")

compTable

fisher.test(compTable) ## Exon enrichment: OR = 1.22


### next promoters


allCGoverlap <- findOverlaps(CGsites.clean, contAgeGenes.promoter, type = "any", select = "first")

numPromoterCpGs <- sum(!is.na(allCGoverlap))

### now overlap age-related sites

aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.promoter, type = "any", select = "first")

numPromoterDMCs <- sum(!is.na(aDMCoverlap))

### set up fisher test

compTable <- data.frame("NonAge" = c(length(CGsites.clean) - length(na.omit(allCGoverlap)), 
                                     length(na.omit(allCGoverlap))),
                        "aDMC" = c(length(ageDMCs.gr) - length(na.omit(aDMCoverlap)),
                                   length(na.omit(aDMCoverlap))))
rownames(compTable) <- c("NotWithinContext", "WithinContext")

compTable

fisher.test(compTable) ## Exon enrichment: OR = 1.28



### finally 'terminators'

allCGoverlap <- findOverlaps(CGsites.clean, contAgeGenes.terminators, maxgap = 0L, type = "any", select = "first")

numTerminatorCpGs <- sum(!is.na(allCGoverlap))

### now overlap age-related sites

aDMCoverlap <- findOverlaps(ageDMCs.gr, contAgeGenes.terminators, maxgap = 0L, type = "any", select = "first")

numTerminatorDMCs <- sum(!is.na(aDMCoverlap))

### set up fisher test

compTable <- data.frame("NonAge" = c(length(CGsites.clean) - length(na.omit(allCGoverlap)), 
                                     length(na.omit(allCGoverlap))),
                        "aDMC" = c(length(ageDMCs.gr) - length(na.omit(aDMCoverlap)),
                                   length(na.omit(aDMCoverlap))))
rownames(compTable) <- c("NotWithinContext", "WithinContext")

compTable

fisher.test(compTable) ## Exon enrichment: OR = 1.74


#### Count by context

linkedDMCs_context <- data.frame("Context" = factor(c("Promoter", "Exon", "Intron", "Terminator"), 
                                                    levels = c("Promoter", "Exon", "Intron", "Terminator")),
                                 "Count" = c(numPromoterDMCs, numExonDMCs, numIntronDMCs, numTerminatorDMCs),
                                 "TotalVolume" = c(sum(width(contAgeGenes.promoter)), sum(width(contAgeGenes.exons)),
                                                   sum(width(contAgeGenes.introns)), sum(width(contAgeGenes.terminators)))/1000)


ggplot(linkedDMCs_context, aes(x = Context, y = Count/TotalVolume, fill = Context)) + 
                geom_col() + 
                labs(x = "aDEG Genic Context", y = "Count per Kb") + 
                geom_text(aes(label = Count),
                          size = 5, fontface = "bold")

linkedDMCs_context <- data.frame("Context" = factor(c("Promoter", "Exon", "Intron", "Terminator"), 
                                                    levels = c("Promoter", "Exon", "Intron", "Terminator")),
                                 "aDMCCount" = c(numPromoterDMCs, numExonDMCs, numIntronDMCs, numTerminatorDMCs),
                                 "aDEGTotalVolume" = c(sum(width(contAgeGenes.promoter)), sum(width(contAgeGenes.exons)),
                                                   sum(width(contAgeGenes.introns)), sum(width(contAgeGenes.terminators))),
                                 "CpGCount" = c(numPromoterCpGs, numExonCpGs, numIntronCpGs, numTerminatorCpGs))

linkedDMCs_context$Ratio <- linkedDMCs_context$aDMCCount / linkedDMCs_context$CpGCount


d <- ggplot(linkedDMCs_context, aes(x = Context, y = Ratio, fill = Context)) + 
  geom_col(position = position_dodge(1)) + 
  labs(x = "aDEG Genic Context", y = "Ratio of aDMCs to total CpGs per Kb") + 
  geom_text(aes(label = aDMCCount),
            size = 8, size.unit = "pt", fontface = "bold", nudge_y = -0.0005)

d

grid.newpage()
grid.draw(set_panel_size(d, width  = unit(3.0, "in"), height = unit(2.0, "in")))


ggsave("./Figures/5E.svg",
       plot = d, device = "svg", width = 3.4, height = 2.0, units = "in")

### Fisher test for terminators specifically



term.fisher <- matrix(c(numTerminatorDMCs, numTerminatorCpGs - numTerminatorDMCs,
                        sum(linkedDMCs_context$aDMCCount)-numTerminatorDMCs, 
                        sum(linkedDMCs_context$CpGCount)-numTerminatorCpGs-numTerminatorDMCs), nrow = 2)
rownames(term.fisher) <- c("Terminator", "Other")
colnames(term.fisher) <- c("aDMC", "non-aDMC")

fisher.test(term.fisher)
