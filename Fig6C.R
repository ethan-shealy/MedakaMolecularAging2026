library(data.table)
library(ggpubr)
library(tidyverse)
library(Rsamtools)
library(GenomicAlignments)
library(txdbmaker)
library(Gviz)
library(rtracklayer)

setwd("")

options(ucscChromosomeNames=FALSE)

gtrack <- GenomeAxisTrack()

### Transcript Track

txdb = makeTxDbFromGFF("./Code_Data/NFIX_RAD23a_summarized.gff")

txTrack <- GeneRegionTrack(
  txdb, name = "NCBI RefSeq Transcripts",
  transcriptAnnotation = "transcript", 
  chromosome = "NC_019866.2",
  start = 22741733,
  end = 22952584
)


txTrack@range$gene <- ifelse(txTrack@range$gene == "rad23a", "RAD23a", "NFIXb")
txTrack@range$density

displayPars(txTrack) <- list(fontcolor.title = "black", just.group = "above",
                             col = "black", cex.group = 0.8)



### DMC track

dmcs <- import.bedGraph("./DSS/SigLogAgeCpGs.bedGraph")

dmcs.nfix <- dmcs[which(seqnames(dmcs) == "NC_019866.2" & 
                          start(dmcs) > 22747037 &
                          end(dmcs) < 22924273)]


dtrack <- DataTrack(data = dmcs.nfix$score, start = start(dmcs.nfix), genome = "ASM223467v1",
                    end = end(dmcs.nfix), chromosome = seqnames(dmcs.nfix),
                    name = "-log10(p) \naDMCS")

displayPars(dtrack) <- list(aggregateGroups = FALSE, baseline = 0, col.baseline = "black", col.histogram = "blue", fill.histogram = "blue",
                            fontcolor.title = "black", col.axis = "black")

### SNPs track

snps <- fread("./Code_Data/AgeCohort_NFIX_RAD23a.vcf")

snptrack <- AnnotationTrack(chromosome = "NC_019866.2", start = snps$POS, width = 1, 
                            name = "EA Clock SNPs")

displayPars(snptrack) <- list(fontcolor.title = "black", col = "red", fill = "red",
                              stacking = "dense", rotation.title = 0, h = 1)

### Clock sites track

clock <- read.table("./Clocks/ElasticNet_DNAm_AgeTrans_sites.bedGraph")

clock.nfix <- clock[which(clock$V1 == "NC_019866.2" & 
                          clock$V2 > 22747037 &
                          clock$V2 < 22924273),]

clocktrack <- AnnotationTrack(chromosome = "NC_019866.2", start = clock.nfix$V2, width = 1, 
                            name = "Clock CpGs")

displayPars(clocktrack) <- list(fontcolor.title = "black", col = "magenta", fill = "magenta", collapse = FALSE, 
                              stacking = "full", rotation.title = 0, cex.title = 0.6, h = 1, min.distance = 0, stackHeight = 1)


plotTracks(list(gtrack, dtrack, snptrack, clocktrack, txTrack), type = "histogram", sizes = c(0.2, 0.8, 0.25, 0.25, 1.2), 
           collapseTranscripts = FALSE, title.width = 0.7)

map.1 <- grid.grab()

## draw it, changes optional
map.1.new <- editGrob(map.1, vp = viewport(width = unit(8, "in"), 
                                           height = unit(4, "in")))
grid.newpage()
grid.draw(map.1.new)

ggsave("./Figures/6C.svg",
       plot = map.1.new, device = "svg", width = 8, height = 4, units = "in")




### Version with data track for clock sites


clock.dtrack <- DataTrack(data = clock.nfix$V4, start = clock.nfix$V2, genome = "ASM223467v1",
                    end = clock.nfix$V3, chromosome = clock.nfix$V1,
                    name = "Clock Coeff")

displayPars(clock.dtrack) <- list(aggregateGroups = FALSE, baseline = 0, col.baseline = "black", col.histogram = "magenta", fill.histogram = "magenta",
                            fontcolor.title = "black", col.axis = "black")


plotTracks(list(gtrack, dtrack, snptrack, clock.dtrack, txTrack), type = "histogram", sizes = c(0.2, 0.5, 0.15, 0.5, 1), 
           collapseTranscripts = FALSE, title.width = 0.7)

