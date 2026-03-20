library(tidyverse)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(ggpubr)
library(psych)

setwd("")

theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=6),
                  axis.title=element_text(size=9,face="bold"),
                  strip.text = element_text(size=9,face="bold"),
                  panel.border = element_rect(fill = NA)))



CpG <- fread(file = "./5x_50ind_manualFilt.tab", header = TRUE, sep = "\t", data.table = FALSE)

coords <- CpG[,1]

coords.df <- as.data.frame(str_split_fixed(coords, ":", 2))

chrSizes <- read.table("./Code_Data/ASM223467v1_chromSizes.txt")

chrSizes$cum <- c(0, cumsum(chrSizes$V2)[-length(chrSizes$V2)])

CpG <- CpG[,-1]

samples <- colnames(CpG)

## Check % of missing sites in each sample

miss <- rep(NA, length(CpG))

for (i in 1:length(CpG[1,])){
  miss[i] <- sum(is.na(CpG[,i])) / length(CpG[,i])
}

hist(miss, breaks = 20)

## Metadata matching

dat <- readxl::read_xlsx("./Code_Data/2022 Aging cohort.xlsx", na = c("", "NA", "#N/A"))

dat$K <- (100 * dat$Tot_wgt) / (dat$Tot_length^3)

dat$ID <- paste0("X", dat$ID)

dat.seq <- dat[match(samples, dat$ID),]

hist(dat.seq$Age, main = "", xlab = "Ages in months")

cpg.mat <- as.data.frame(t(CpG))


#### DMC Analysis

dmcs <- fread("./DSS/SigLogAgeCpGs.tab")

dmcs$ChrOffset <- chrSizes$cum[match(dmcs$chr, chrSizes$V1)] + dmcs$pos

manhattan <- ggplot(dmcs, aes(x = ChrOffset, y = -log10(pvals), color = chr)) + 
  geom_point(show.legend = FALSE) + labs(x = "Chromosomal Location", y = "-log10(p)") + 
  scale_color_manual(values = c("grey60", rep(c("grey20", "grey60"), 12))) + 
  scale_x_continuous(breaks = chrSizes$cum[-1], labels = chrSizes$V1[-1]) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0), axis.ticks.x = element_blank())
manhattan

library(grid)
library(egg)

grid.newpage()
grid.draw(set_panel_size(manhattan, width  = unit(7, "in"), height = unit(2, "in")))

ggsave("./Figures/4A.tif",
       plot = manhattan, device = "tiff", width = 7, height = 2, units = "in")




#rand.dmcs <- dmcs[sample(1:nrow(dmcs), 16),]
##rand.dmcs <- dmcs[order(dmcs$pvals),][1:25,]
#
#
#
#rand.mat <- data.frame(dat.seq, cpg.mat[,match(paste0(rand.dmcs$chr, ":", rand.dmcs$pos), coords)])
#colnames(rand.mat)[27:ncol(rand.mat)] <- paste0(rand.dmcs$chr, ":", rand.dmcs$pos)
#
#
#rand.mat.long <- pivot_longer(rand.mat, cols = starts_with("NC_"), 
#                              names_to = "CpG", values_to = "Methylation")
#
#
#ggplot(rand.mat.long, aes(x = Age, y = Methylation)) + geom_smooth(method = "lm", formula = y ~ log(x)) + 
#  geom_point() + facet_wrap(~CpG, scales = "free_y")

dmc.mat <- data.frame(cpg.mat[,match(paste0(dmcs$chr, ":", dmcs$pos), coords)])

missing <- rowSums(is.na(dmc.mat)); hist(missing)

library(ComplexHeatmap)


age_fun = circlize::colorRamp2(c(0, 30), c("green", "purple"))

col_ha = columnAnnotation(Age = dat.seq$Age[-which(missing > 15000)], Sex = dat.seq$Gonad_sex[-which(missing > 15000)],
                          col = list(Age = age_fun, Sex = c("F" = "salmon", "M" = "skyblue2")))

table(sign(dmcs$stat))

Heatmap(t(dmc.mat)[,-which(missing > 15000)], name = "Meth %",
        column_order = order(dat.seq$Age[-which(missing > 15000)]), 
        show_row_dend = FALSE, column_names_gp = gpar(fontsize = 8),
        row_split = as.factor(sign(dmcs$stat)), 
        row_title = c("Increasing (n = 8,348)", "Decreasing (n = 15,428)"),
        show_row_names = FALSE, bottom_annotation = col_ha, 
        raster_resize_mat = TRUE)

heat.1 <- grid.grab()

## draw it, changes optional
heat.1.new <- editGrob(heat.1, vp = viewport(width = unit(4.2, "in"), 
                   height = unit(3.6, "in")))
grid.newpage()
grid.draw(heat.1.new)

ggsave("./Figures/4B.svg",
      plot = heat.1.new, device = "svg", width = 4.2, height = 3.6, units = "in")





### Two largest DMRs

dat.seq$AgeBin <- cut(dat.seq$Age, c(0, 5, 10, 18, 30))
table(dat.seq$AgeBin)

# NFIX

cpg.nfix1 <- colMeans(CpG[which(coords.df$V1 == "NC_019866.2" & 
                                  coords.df$V2 > 22748000 &
                                  coords.df$V2 < 22749000),-1], na.rm = TRUE)

nrow(CpG[which(coords.df$V1 == "NC_019866.2" & 
                   coords.df$V2 > 22748000 &
                   coords.df$V2 < 22749000),-1])

head(cpg.nfix1)


dat.seq$NFIX1 <- cpg.nfix1[match(dat.seq$ID, names(cpg.nfix1))]

ggplot(dat.seq, aes(x = Age, y = NFIX1)) + 
  geom_point() + geom_smooth(method = "lm") + 
  labs(title = "NC_019866.2: 22,748,000 - 22,749,000") + 
  ylim(0, 100)

ggplot(dat.seq, aes(x = Age, y = NFIX1, color = Gonad_sex)) + 
  geom_point() + geom_smooth(method = "lm") + 
  labs(title = "NC_019866.2: 22,748,000 - 22,749,000")

summary(lm(NFIX1 ~ Age + Gonad_sex, data = dat.seq))

# Second promoter region


cpg.nfix2 <- colMeans(CpG[which(coords.df$V1 == "NC_019866.2" & 
                                  coords.df$V2 > 22801600 &
                                  coords.df$V2 < 22802400),-1], na.rm = TRUE)

nrow(CpG[which(coords.df$V1 == "NC_019866.2" & 
                 coords.df$V2 > 22801600 &
                 coords.df$V2 < 22802400),-1])

head(cpg.nfix2)


dat.seq$NFIX2 <- cpg.nfix2[match(dat.seq$ID, names(cpg.nfix2))]

ggplot(dat.seq, aes(x = Age, y = NFIX2)) + 
  geom_point() + geom_smooth(method = "lm") + 
  labs(title = "NC_019866.2: 22,801,600 - 22,802,400") + ylim(0,100)

ggplot(dat.seq, aes(x = Age, y = NFIX2, color = Gonad_sex)) + 
  geom_point() + geom_smooth(method = "lm") + 
  labs(title = "NC_019866.2: 22,748,500 - 22,748,800")

summary(lm(NFIX2 ~ Age + Gonad_sex, data = dat.seq))


### NFIX genotype
#
#nfixGT <- read.csv("./Code_Data/NFIX_RAD23_Genotypes.csv")
#
#dat.seq$GT <- as.factor(nfixGT$GT[match(dat.seq$ID, nfixGT$ID)])
#
#ggplot(filter(dat.seq, !is.na(GT)), aes(x = Age, y = NFIX1, color = GT)) + 
#  geom_point() + geom_smooth(method = "lm", alpha = 0.2) + 
#  coord_cartesian(ylim = c(0,105)) +
#  scale_color_manual(values = c("blue", "purple", "red"), labels = c("GG", "GA", "AA")) + 
#  ylab("NFIX Methylation")

## Gene expression


#counts <- fread("./unfiltered_genecounts.csv", sep = ",", data.table = FALSE)
#
#sampleIDs <- counts$V1 %>% str_remove_all(".bam") %>% str_replace_all("ID", "X")
#
#cpm <- as.data.frame(t(edgeR::cpm(t(counts[,-c(1:3)]))))
#
#abundances <- cpm[match(dat.seq$ID, sampleIDs),]
#
#dat.seq$NFIX_Expr <- log2(abundances$LOC101164563 + 1)
#
#ggplot(filter(dat.seq, !is.na(GT)), aes(x = Age, y = NFIX_Expr, color = GT)) + 
#  geom_point() + geom_smooth(method = "lm", alpha = 0.2) + 
#  scale_color_manual(values = c("blue", "purple", "red"), labels = c("GG", "GA", "AA")) + 
#  ylab("log2 NFIX Expression")
#
#summary(lm(NFIX_Expr ~ Age*GT, data = dat.seq))


#### PPARA / MIR
#
#cpg.ppar <- colMeans(CpG[which(coords.df$V1 == "NC_019881.2" & 
#                                  coords.df$V2 > 6761400 &
#                                  coords.df$V2 < 6762200),-1], na.rm = TRUE)
#
#
#
#nrow(CpG[which(coords.df$V1 == "NC_019881.2" & 
#                    coords.df$V2 > 6761400 &
#                    coords.df$V2 < 6762200),-1])
#
#head(cpg.ppar)
#
#
#dat.seq$PPAR <- cpg.ppar[match(dat.seq$ID, names(cpg.ppar))]
#
#ggplot(dat.seq, aes(x = Age, y = PPAR, color = Gonad_sex)) + 
#  geom_point() + geom_smooth(method = "lm", formula = y ~ log(x)) + 
#  labs(title = "NC_019881.2: 6,761,400 - 6,762,200")
#
#ggplot(dat.seq, aes(x = NFIX1, y = PPAR, color = Gonad_sex)) + 
#  geom_point() + geom_smooth(method = "lm", formula = y ~ x,
#                             aes(group = 1), color = "black")



#### CpG clocksite characteristics

clockSites <- read.csv("./Clocks/DNAm_loocvSites.csv")

clockSites$name <- str_replace_all(clockSites$sites, ":", ".")

venn <- ggVennDiagram::ggVennDiagram(list("aDMCs" = dmcs$name, "ClockSites" = clockSites$name), 
                             label = "count") + binned_scale(aesthetics = "fill",
                                                             scale_name = "stepsn", 
                                                             palette = function(x) c("grey60", "grey80", "grey80"),
                                                             breaks = c(-1, 1, 200, 60000),
                                                             guide = "colorsteps") +
  coord_flip() + theme(legend.position = "none")

grid.newpage()
grid.draw(set_panel_size(venn, width  = unit(2.8, "in"), height = unit(1.7, "in")))


ggsave("./Figures/4C.svg",
       plot = venn, device = "svg", width = 2.8, height = 1.7, units = "in")





### Promoter methylation

geneAnnot <- as.data.frame(import.gff("./Code_Data/ASM223467v1.gff"))

geneAnnot <- filter(geneAnnot, type == "gene")

proms <- fread("./DSS/Promoters/allSamples_PromGroupedCovg.tab", data.table = FALSE)

proms.meth <- data.frame(proms[,seq(3, length(proms), 2)] / proms[,seq(2, length(proms), 2)])

colnames(proms.meth) <- str_remove_all(colnames(proms[,seq(3, length(proms), 2)]), pattern = "_")
rownames(proms.meth) <- proms$Gene


#### Differentially methylated promoters analysis

DMproms <- read.table("./DSS/Promoters/SigLogAgeProms.tab",
                      sep = "\t", header = TRUE)


### Get methylation levels


mat.meta <- data.frame(dat.seq, t(proms.meth))

mat.DMproms <- data.frame(proms.meth[match(DMproms$chr, rownames(proms.meth)),])

highMissing <- which(colSums(is.na(mat.DMproms)) > 0.2*nrow(mat.DMproms))

library(ComplexHeatmap)

range(mat.DMproms, na.rm = TRUE)

col_fun = circlize::colorRamp2(c(0, 50, 100), c("blue", "white", "red"))

age_fun = circlize::colorRamp2(c(0, 30), c("green", "purple"))

col_ha = columnAnnotation(Age = dat.seq$Age[-highMissing], Sex = dat.seq$Gonad_sex[-highMissing],
                          col = list(Age = age_fun, Sex = c("F" = "salmon", "M" = "skyblue2", "Im" = "grey60")))


Heatmap(as.matrix(mat.DMproms[,-highMissing])*100, col = col_fun, na_col = "grey20",
        column_order = order(dat.seq$Age[-highMissing]), name = "Meth %",
        show_row_names = FALSE, row_split = as.factor(sign(DMproms$stat)),
        bottom_annotation = col_ha, column_names_gp = gpar(fontsize = 8),
        row_title = c("Increasing (n = 261)", "Decreasing (n = 363)"), 
        use_raster = TRUE, raster_resize_mat = TRUE)

table(as.factor(sign(DMproms$stat)))


heat.2 <- grid.grab()

## draw it, changes optional
heat.2.new <- editGrob(heat.2, vp = viewport(width = unit(4.2, "in"), 
                                             height = unit(3.6, "in")))
grid.newpage()
grid.draw(heat.2.new)

ggsave("./Figures/4G.svg",
       plot = heat.2.new, device = "svg", width = 4.2, height = 3.6, units = "in")


### Overlap with age-related genes
library(ggVennDiagram)

aDEGs <- readxl::read_excel("./Code_Data/DESEQ2_all.xlsx")

venn.2 <- ggVennDiagram(list(aDEGs$Gene, DMproms$chr), category.names = c("aDEGs", "aDMPs"), label = "count") + 
  coord_flip() + binned_scale(aesthetics = "fill", breaks = c(35), palette = function(x) c("navy", "skyblue2")) + 
  theme(legend.position = "none")

grid.newpage()
grid.draw(set_panel_size(venn.2, width  = unit(2.8, "in"), height = unit(1.7, "in")))


ggsave("./Figures/4F.svg",
       plot = venn.2, device = "svg", width = 2.8, height = 1.7, units = "in")


#### Test for enrichment (total of 16,548 genes tested for DEG)
allGenes <- fread("./genecounts_cpm.csv")[,-c(1:4)]

genes.filt <- as.data.frame(allGenes)[, colSums(allGenes > 0.5) >= 10] ### Check this number matches up with MS

allProms <- fread("./DSS/Promoters/allSamples_PromGroupedCovg.tab")

allProms.fixed <- allProms[which(make.names(allProms$Gene) %in% make.names(colnames(genes.filt))),]
allGenes.fixed <- genes.filt[,which(make.names(colnames(genes.filt)) %in% make.names(allProms$Gene))]
aDEGs.fixed <- aDEGs[which(make.names(aDEGs$Gene) %in% make.names(allProms.fixed$Gene)),]
aDMPs.fixed <- DMproms[which(make.names(DMproms$chr) %in% make.names(allProms.fixed$Gene)),]

neither <- allProms.fixed[-which(allProms.fixed$Gene %in% aDEGs$Gene | allProms.fixed$Gene %in% aDMPs.fixed$chr),]


fisher.table <- matrix(c(sum(make.names(aDEGs.fixed$Gene) %in% make.names(aDMPs.fixed$chr)),  ### aDEGs also aDMPs
                         sum(!(make.names(aDEGs.fixed$Gene) %in% make.names(aDMPs.fixed$chr))), ### aDEGs not aDMPs
                         sum(!(make.names(aDMPs.fixed$chr) %in% make.names(aDEGs.fixed$Gene))), ### aDMPs not aDEGs
                         nrow(neither) ### Total genes under consideration not aDEGs or aDMPs
                         ), nrow = 2)
colnames(fisher.table) <- c("aDEG", "non_aDEG")
rownames(fisher.table) <- c("aDMP", "non_aDMP")

fisher.test(fisher.table) ## N.S. P = 0.2486

################# Enrichments



# Age DMC enrichment ------------------------------------------------------


## First positively correlated sites 

# Import CpG island gff file generated using EMBOSS's cpgplot
islands <- import("./Code_Data/CGislands.gff", format = "gff")

# define shores, shelves, seas

shores.up <- trim(flank(islands, 2000))
shores.down <- trim(flank(islands, 2000, start = FALSE))
shores <- c(shores.up, shores.down)
shores.clean <- reduce(GenomicRanges::setdiff(shores, islands))

shelves.up <- trim(flank(shores.up, 2000))
shelves.down <- trim(flank(shores.down, 2000, start = FALSE))
shelves <- c(shelves.up, shelves.down)
shelves.no.islands <- GenomicRanges::setdiff(shelves, islands)
shelves.clean <- reduce(GenomicRanges::setdiff(shelves.no.islands, shores))


not.seas <- c(islands,shores,shelves)
seas <- gaps(not.seas)

## CGdensity enrichment 

## Find background likelihood that CpGs are found in each context

CGsites <- fread("./DSS/allSamples_Covg_5x_50L.tab", sep = "\t", header = TRUE)

CGsites <- CGsites[which(CGsites$strand == "+"),]

CGsites.clean <- GRanges(seqnames = CGsites$chr, 
                         ranges = IRanges(start = CGsites$pos, end = CGsites$pos),
                         strand = "*")

totalCG <- length(CGsites.clean)


sites.in.islands <- countOverlaps(CGsites.clean, islands, type = "within")
sum(sites.in.islands) 
cat("Proportion of CG sites that fall in islands: ", sum(sites.in.islands) / length(CGsites.clean), "\n")


sites.in.shores <- countOverlaps(CGsites.clean, shores.clean, type = "within")
sum(sites.in.shores)
cat("Proportion of CG sites that fall in shores: ", sum(sites.in.shores) / length(CGsites.clean), "\n")


sites.in.shelves <- countOverlaps(CGsites.clean, shelves.clean, type = "within")
sum(sites.in.shelves) 
cat("Proportion of CG sites that fall in shelves: ", sum(sites.in.shelves) / length(CGsites.clean), "\n")


sites.in.seas <- countOverlaps(CGsites.clean, seas, type = "within")
sum(sites.in.seas) 
cat("Proportion of CG sites that fall in seas: ", sum(sites.in.seas) / length(CGsites.clean), "\n")

# Import aDMCs w/ correlation coeffs

CPG <- read.table("./DSS/SigLogAgeCpGs.tab", 
                  sep = "\t", header = TRUE)

CPG <- GRanges(CPG$chr, IRanges(CPG$pos, CPG$pos), strand = "*", name = CPG$stat)

CPG <- CPG[which(CPG$name > 0),]

totalCpG <- length(CPG)

cpg.in.islands <- countOverlaps(CPG, islands, type = "within")
sum(cpg.in.islands) 
cat("Proportion of (+) Age-Associated CpGs that fall in islands: ", sum(cpg.in.islands) / length(CPG), "\n")


cpg.in.shores <- countOverlaps(CPG, shores.clean, type = "within")
sum(cpg.in.shores) 
cat("Proportion of (+) Age-Associated CpGs that fall in shores: ", sum(cpg.in.shores) / length(CPG), "\n")


cpg.in.shelves <- countOverlaps(CPG, shelves.clean, type = "within")
sum(cpg.in.shelves) 
cat("Proportion of (+) Age-Associated CpGs that fall in shelves: ", sum(cpg.in.shelves) / length(CPG), "\n")


cpg.in.seas <- countOverlaps(CPG, seas, type = "within")
sum(cpg.in.seas) 
cat("Proportion of (+) Age-Associated CpGs that fall in seas: ", sum(cpg.in.seas) / length(CPG), "\n")


pos.density.results <- data.frame("Context" = c("islands", "shores", "shelves", "seas","islands", "shores", "shelves", "seas"), 
                                  "Values" = c(sum(sites.in.islands), sum(sites.in.shores),
                                               sum(sites.in.shelves), sum(sites.in.seas),
                                               sum(cpg.in.islands), sum(cpg.in.shores), 
                                               sum(cpg.in.shelves), sum(cpg.in.seas)),
                                  "Proportions" = c(sum(sites.in.islands)/totalCG, sum(sites.in.shores)/totalCG,
                                                    sum(sites.in.shelves)/totalCG, sum(sites.in.seas)/totalCG,
                                                    sum(cpg.in.islands)/totalCpG, sum(cpg.in.shores)/totalCpG, 
                                                    sum(cpg.in.shelves)/totalCpG, sum(cpg.in.seas)/totalCpG),
                                  "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

pos.density.results$Context <- factor(pos.density.results$Context, levels=c("islands","shores","shelves", "seas"))

enrichment.cgdensity <- ggplot(data = pos.density.results, mapping = aes(Context, Proportions)) + 
  geom_bar(aes(fill = type), stat = "identity", position = "dodge", color = "black") + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(title = "Enrichment of CpGs with (+) Age-Association") +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black")) +
  font("xlab", size = 18, color = "black") +
  font("ylab", size = 18, color = "black") +
  font("xy.text", size = 18, color = "black") +
  font("title", size = 15, color = "black", face = "bold")

enrichment.cgdensity


## Statistical test:

#Islands
binom.test(x = sum(cpg.in.islands), n = totalCpG, p = ( sum(sites.in.islands)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p < 2.2e-16 ***

#Shores
binom.test(x = sum(cpg.in.shores), n = totalCpG, p = ( sum(sites.in.shores)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p < 2.2e-16 ***

#Shelves
binom.test(x = sum(cpg.in.shelves), n = totalCpG, p = ( sum(sites.in.shelves)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = N.S.

#Seas
binom.test(x = sum(cpg.in.seas), n = totalCpG, p = ( sum(sites.in.seas)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p < 2.2e-16 ***


## Generate genic contexts

annot <- import.gff("./Code_Data/ASM223467v1.gff")

exons <- annot[which(annot$type == "exon")]

exons.destranded <- reduce(exons, ignore.strand = TRUE)

CDS <- annot[which(annot$type == "gene")]
CDS.destranded <- reduce(CDS, ignore.strand = TRUE)
full.CDS <- GenomicRanges::union(CDS.destranded, exons.destranded)

intergenic <- gaps(CDS.destranded)
intergenic.clean <- GenomicRanges::setdiff(intergenic, exons.destranded)

exons.comp <- gaps(exons)

introns <- GenomicRanges::setdiff(CDS.destranded, exons.destranded)

promoter <- reduce(promoters(CDS), ignore.strand = TRUE)

#### Genic context Enrichment  

## Find background likelihood that CpGs are found in each context

sites.in.promo <- countOverlaps(CGsites.clean, promoter, type = "within")
sum(sites.in.promo)
cat("Proportion of CG sites that fall in promoters: ", sum(sites.in.promo) / length(CGsites.clean), "\n")


sites.in.exons <- countOverlaps(CGsites.clean, exons.destranded, type = "within")
sum(sites.in.exons) 
cat("Proportion of CG sites that fall in exons: ", sum(sites.in.exons) / length(CGsites.clean), "\n")


sites.in.introns <- countOverlaps(CGsites.clean, introns, type = "within")
sum(sites.in.introns)
cat("Proportion of CG sites that fall in introns: ", sum(sites.in.introns) / length(CGsites.clean), "\n")


sites.in.intergenic <- countOverlaps(CGsites.clean, intergenic.clean, type = "within")
sum(sites.in.intergenic) 
cat("Proportion of CG sites that fall in intergenic space: ", sum(sites.in.intergenic) / length(CGsites.clean), "\n")

# get likelihood of age DMCs in each context

cpg.in.promo <- countOverlaps(CPG, promoter, type = "within")
sum(cpg.in.promo)
cat("Proportion of (+) Age-Associated CpGs that fall in promoters: ", sum(cpg.in.promo) / length(CPG), "\n")

cpg.in.exons <- countOverlaps(CPG, exons.destranded, type = "within")
sum(cpg.in.exons)
cat("Proportion of (+) Age-Associated CpGs that fall in exons: ", sum(cpg.in.exons) / length(CPG), "\n")


cpg.in.introns <- countOverlaps(CPG, introns, type = "within")
sum(cpg.in.introns) 
cat("Proportion of (+) Age-Associated CpGs that fall in introns: ", sum(cpg.in.introns) / length(CPG), "\n")


cpg.in.intergenic <- countOverlaps(CPG, intergenic.clean, type = "within")
sum(cpg.in.intergenic) 
cat("Proportion of (+) Age-Associated CpGs that fall in intergenic space: ", sum(cpg.in.intergenic) / length(CPG), "\n")

# combine results

pos.genic.results <- data.frame("Context" = c("promoters", "exons", "introns", "intergenic", "promoters", "exons", "introns", "intergenic"), 
                                "Values" = c(sum(sites.in.promo), sum(sites.in.exons), 
                                             sum(sites.in.introns), sum(sites.in.intergenic),
                                             sum(cpg.in.promo), sum(cpg.in.exons), 
                                             sum(cpg.in.introns), sum(cpg.in.intergenic)),
                                "Proportions" = c(sum(sites.in.promo)/totalCG, sum(sites.in.exons)/totalCG, 
                                                  sum(sites.in.introns)/totalCG, sum(sites.in.intergenic)/totalCG,
                                                  sum(cpg.in.promo)/totalCpG, sum(cpg.in.exons)/totalCpG, 
                                                  sum(cpg.in.introns)/totalCpG, sum(cpg.in.intergenic)/totalCpG),
                                "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

pos.genic.results$Context <- factor(pos.genic.results$Context, levels=c("promoters", "exons", "introns", "intergenic"))


enrichment.genic <- ggplot(data = pos.genic.results, mapping = aes(Context, Proportions)) + 
  geom_bar(aes(fill = type), stat = "identity", position = "dodge", color = "black") + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(title = "Enrichment of CpGs with (+) Age-Association") +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black")) +
  font("xlab", size = 18, color = "black") +
  font("ylab", size = 18, color = "black") +
  font("xy.text", size = 18, color = "black") +
  font("title", size = 15, color = "black", face = "bold")

enrichment.genic


## Statistical test:

#Promoters
binom.test(x = sum(cpg.in.promo), n = totalCpG, p = ( sum(sites.in.promo)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.02 . (N.S after multiple testing correction)

#Exons
binom.test(x = sum(cpg.in.exons), n = totalCpG, p = ( sum(sites.in.exons)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = N.S.

#Introns
binom.test(x = sum(cpg.in.introns), n = totalCpG, p = ( sum(sites.in.introns)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 0.02 . (N.S after multiple testing correction)

#Intergenic
binom.test(x = sum(cpg.in.intergenic), n = totalCpG, p = ( sum(sites.in.intergenic)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.05 . (N.S after multiple testing correction)


## Now negatively correlated sites


CPG <- read.table("./DSS/SigLogAgeCpGs.tab", 
                  sep = "\t", header = TRUE)

CPG <- GRanges(CPG$chr, IRanges(CPG$pos, CPG$pos), strand = "*", name = CPG$stat)

CPG <- CPG[which(CPG$name < 0),]
totalCpG <- length(CPG)

cpg.in.islands <- countOverlaps(CPG, islands, type = "within")
sum(cpg.in.islands)
cat("Proportion of (-) Age-Associated CpGs that fall in islands: ", sum(cpg.in.islands) / length(CPG), "\n")


cpg.in.shores <- countOverlaps(CPG, shores.clean, type = "within")
sum(cpg.in.shores) 
cat("Proportion of (-) Age-Associated CpGs that fall in shores: ", sum(cpg.in.shores) / length(CPG), "\n")


cpg.in.shelves <- countOverlaps(CPG, shelves.clean, type = "within")
sum(cpg.in.shelves) 
cat("Proportion of (-) Age-Associated CpGs that fall in shelves: ", sum(cpg.in.shelves) / length(CPG), "\n")


cpg.in.seas <- countOverlaps(CPG, seas, type = "within")
sum(cpg.in.seas) 
cat("Proportion of (-) Age-Associated CpGs that fall in seas: ", sum(cpg.in.seas) / length(CPG), "\n")


neg.density.results <- data.frame("Context" = c("islands", "shores", "shelves", "seas","islands", "shores", "shelves", "seas"), 
                                  "Values" = c(sum(sites.in.islands), sum(sites.in.shores),
                                               sum(sites.in.shelves), sum(sites.in.seas),
                                               sum(cpg.in.islands), sum(cpg.in.shores), 
                                               sum(cpg.in.shelves), sum(cpg.in.seas)),
                                  "Proportions" = c(sum(sites.in.islands)/totalCG, sum(sites.in.shores)/totalCG,
                                                    sum(sites.in.shelves)/totalCG, sum(sites.in.seas)/totalCG,
                                                    sum(cpg.in.islands)/totalCpG, sum(cpg.in.shores)/totalCpG, 
                                                    sum(cpg.in.shelves)/totalCpG, sum(cpg.in.seas)/totalCpG),
                                  "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

neg.density.results$Context <- factor(neg.density.results$Context, levels=c("islands","shores","shelves", "seas"))


enrichment.cgdensity <- ggplot(data = neg.density.results, mapping = aes(Context, Proportions)) + 
  geom_bar(aes(fill = type), stat = "identity", position = "dodge", color = "black") + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(title = "Enrichment of CpGs with (-) Age-Association") +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black")) +
  font("xlab", size = 18, color = "black") +
  font("ylab", size = 18, color = "black") +
  font("xy.text", size = 18, color = "black") +
  font("title", size = 15, color = "black", face = "bold")
#geom_signif(stat = "identity", data = data.frame(x=c(0.65, 1.65, 3.65), 
#                                                 xend=c(1.3, 2.3, 4.3),
#                                                 y=c(0.1, 0.18, 0.73), 
#                                                 annotation=c(" *** ", "*", "   **   ")),
#            aes(x=x,xend=xend, y=y, yend=y, annotation=annotation))

enrichment.cgdensity


## Statistical test:

#Islands
binom.test(x = sum(cpg.in.islands), n = totalCpG, p = ( sum(sites.in.islands)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p < 2.2e-16 ***

#Shores
binom.test(x = sum(cpg.in.shores), n = totalCpG, p = ( sum(sites.in.shores)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p = 1.824e-10 **

#Shelves
binom.test(x = sum(cpg.in.shelves), n = totalCpG, p = ( sum(sites.in.shelves)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 0.0003584 **

#Seas
binom.test(x = sum(cpg.in.seas), n = totalCpG, p = ( sum(sites.in.seas)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p < 2.2e-16

cpg.in.promo <- countOverlaps(CPG, promoter, type = "within")
sum(cpg.in.promo) 
cat("Proportion of (+) Age-Associated CpGs that fall in promoters: ", sum(cpg.in.promo) / length(CPG), "\n")

cpg.in.exons <- countOverlaps(CPG, exons.destranded, type = "within")
sum(cpg.in.exons) 
cat("Proportion of (-) Age-Associated CpGs that fall in exons: ", sum(cpg.in.exons) / length(CPG), "\n")


cpg.in.introns <- countOverlaps(CPG, introns, type = "within")
sum(cpg.in.introns) 
cat("Proportion of (-) Age-Associated CpGs that fall in introns: ", sum(cpg.in.introns) / length(CPG), "\n")


cpg.in.intergenic <- countOverlaps(CPG, intergenic.clean, type = "within")
sum(cpg.in.intergenic) 
cat("Proportion of (-) Age-Associated CpGs that fall in intergenic space: ", sum(cpg.in.intergenic) / length(CPG), "\n")


neg.genic.results <- data.frame("Context" = c("promoters", "exons", "introns", "intergenic", "promoters", "exons", "introns", "intergenic"), 
                                "Values" = c(sum(sites.in.promo), sum(sites.in.exons), 
                                             sum(sites.in.introns), sum(sites.in.intergenic),
                                             sum(cpg.in.promo), sum(cpg.in.exons), 
                                             sum(cpg.in.introns), sum(cpg.in.intergenic)),
                                "Proportions" = c(sum(sites.in.promo)/totalCG, sum(sites.in.exons)/totalCG, 
                                                  sum(sites.in.introns)/totalCG, sum(sites.in.intergenic)/totalCG,
                                                  sum(cpg.in.promo)/totalCpG, sum(cpg.in.exons)/totalCpG, 
                                                  sum(cpg.in.introns)/totalCpG, sum(cpg.in.intergenic)/totalCpG),
                                "type" = c(rep("Background", 4), rep("Age-Associated", 4)))

neg.genic.results$Context <- factor(neg.genic.results$Context, levels=c("promoters", "exons", "introns", "intergenic"))


enrichment.genic <- ggplot(data = neg.genic.results, mapping = aes(Context, Proportions)) + 
  geom_bar(aes(fill = type), stat = "identity", position = "dodge", color = "black") + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  labs(title = "Enrichment of CpGs with (-) Age-Association") +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black")) +
  font("xlab", size = 18, color = "black") +
  font("ylab", size = 18, color = "black") +
  font("xy.text", size = 18, color = "black") +
  font("title", size = 15, color = "black", face = "bold")

enrichment.genic


## Statistical test:

# Promoters
binom.test(x = sum(cpg.in.promo), n = totalCpG, p = ( sum(sites.in.promo)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 3.32e-12 ***

#Exons
binom.test(x = sum(cpg.in.exons), n = totalCpG, p = ( sum(sites.in.exons)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 2.211e-10 **

#Introns
binom.test(x = sum(cpg.in.introns), n = totalCpG, p = ( sum(sites.in.introns)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95)  ## p < 2.2e-16 ***

#Intergenic
binom.test(x = sum(cpg.in.intergenic), n = totalCpG, p = ( sum(sites.in.intergenic)/totalCG ),
           alternative = "two.sided",
           conf.level = 0.95) ## p = 8.521e-05 *


## Now plot all of these results together


pos.density.results$sign <- paste0("+ ", pos.density.results$type)
neg.density.results$sign <- paste0("- ", pos.density.results$type)

total.density <- rbind(pos.density.results, neg.density.results)
total.density$sign[c(1,2,3,4,9,10,11,12)] <- "Background"
total.density$sign <- factor(total.density$sign, levels = c("- Age-Associated", "Background", "+ Age-Associated"))

# Plot and save fig 1F

cgdensity <- ggplot(data = total.density, aes(x = Context, y = Proportions)) + 
  geom_bar(aes(fill = sign), stat = "identity", position = "dodge", color = "black") + 
  scale_fill_manual(values=c("navy", "forestgreen", "firebrick")) +
  #geom_text(aes(group = sign, y = -0.01), position = position_dodge(0.9)) +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"))
cgdensity

grid.newpage()
grid.draw(set_panel_size(cgdensity, width  = unit(2.75, "in"), height = unit(1.75, "in")))

ggsave("./Figures/4D.svg",
       plot = cgdensity, device = "svg", width = 2.75, height = 1.75, units = "in")




#ggplot(data = total.density[-c(9:12),], aes(x = sign, y = Proportions, fill = sign, label = paste0(Context, " (n=", Values, ")"))) + 
#  geom_col(position = "stack", color = "black", ) + 
#  scale_fill_manual(values=c("skyblue4", "forestgreen", "firebrick")) +
#  ylab("Percentage of Total CpGs") +
#  geom_text(position = position_stack(vjust = 0.5), size = 4) +
#  xlab("") + 
#  theme(axis.line = element_line(colour = "black"), legend.position = "none",
#        axis.text=element_text(size=10),
#        axis.title=element_text(size=12,face="bold"))

total.density[,c(1,2,5)]


pos.genic.results$sign <- paste0("+ ", pos.genic.results$type)
neg.genic.results$sign <- paste0("- ", pos.genic.results$type)

total.genic <- rbind(pos.genic.results, neg.genic.results)
total.genic$sign[c(1,2,3,4,9,10,11,12)] <- "Background"
total.genic$sign <- factor(total.genic$sign, levels = c("- Age-Associated", "Background", "+ Age-Associated"))

# plot and save fig 1G

genic <- ggplot(data = total.genic, aes(x = Context, y = Proportions)) + 
  geom_bar(aes(fill = sign), stat = "identity", position = "dodge",   color = "black") + 
  scale_fill_manual(values=c("navy", "forestgreen", "firebrick")) +
  #geom_text(aes(group = sign, y = -0.01), position = position_dodge(0.9)) +
  ylab("Percentage of Total CpGs") +
  xlab("") +
  theme(axis.line = element_line(colour = "black"), legend.position = "none",
        axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"))

genic

grid.newpage()
grid.draw(set_panel_size(genic, width  = unit(2.75, "in"), height = unit(1.75, "in")))

ggsave("./Figures/4E.svg",
       plot = genic, device = "svg", width = 2.75, height = 1.75, units = "in")


