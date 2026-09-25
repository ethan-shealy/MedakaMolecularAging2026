library(matrixStats)
library(tidyverse)
library(data.table)
library(ggpubr)
library(grid)
library(egg)

setwd(".")

theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=7),
                  axis.title=element_text(size=9,face="bold"),
                  strip.text = element_text(size=9,face="bold"),
                  panel.border = element_rect(fill = NA)))

meta.seq <- read.csv("./Code_Data/NFIX_RAD23_Genotypes.csv")

meta.seq$AgeBin <- cut(meta.seq$Age, c(0, 10, 20, 30))

chrSizes <- read.table("./Code_Data/ASM223467v1_chromSizes.txt")

chrSizes$cum <- c(0, cumsum(chrSizes$V2)[-length(chrSizes$V2)])

CpG <- fread("./5x_50ind_manualFilt.tab", header = TRUE,
                                   colClasses = "numeric", data.table = FALSE)

coords <- CpG[,1]

rownames(CpG) <- coords
CpG <- CpG[,-1]
samples <- colnames(CpG)

CpG <- CpG[,match(meta.seq$ID, samples)]

clockSites <- read.table("./Clocks/ElasticNet_DNAm_AgeTrans_sites.bedGraph")[-1,]

clockSites$name <- paste0(clockSites$V1, ":", clockSites$V3)

clockSites.cpg <- t(CpG[match(clockSites$name, rownames(CpG)),])

coords.df <- as.data.frame(str_split_fixed(colnames(clockSites.cpg), ":", 2))

## Apply linear models to determine effect of GT*Age interaction at each site
library(broom)

results <- map(as.data.frame(clockSites.cpg), ~ lm(.x ~ Age*as.numeric(GT), data = meta.seq))

tidy_results <- map2_df(results, names(results), ~ tidy(.x) |> mutate(response = .y))

modelWide <- pivot_wider(tidy_results, names_from = "term", values_from = "estimate", id_cols = "response")

modelWide.err <- pivot_wider(tidy_results, names_from = "term", values_from = "std.error", id_cols = "response")

modelWide.sig <- pivot_wider(tidy_results, names_from = "term", values_from = "p.value", id_cols = "response")

modelWide.adj <- data.frame(modelWide.sig$response, apply(modelWide.sig[,-1], 2, p.adjust, method = "fdr"))

stopifnot(modelWide$response == colnames(clockSites.cpg))

modelWide$NFIX_Assoc <- ifelse(coords.df$V1 == "NC_019866.2" & coords.df$V2 > 22745376 & coords.df$V2 < 22926356, 
                               "Yes", "No")

modelWide$Dir <- sign(modelWide$Age) ### Split by Age direction

modelWide$GT.Error <- modelWide.err$`as.numeric(GT)`
modelWide$Int.Error <- modelWide.err$`Age:as.numeric(GT)`

#### Combine clock coeffs with Genotype effect

clockSites$ClockImp <- abs(clockSites$V4) / sum(abs(clockSites$V4))

clockSites$GenotypeEffect <- modelWide$`as.numeric(GT)`[match(clockSites$name, modelWide$response)]

clockSites$GenotypeEffectPerc <- abs(clockSites$GenotypeEffect) / sum(abs(clockSites$GenotypeEffect))

clockSites$NFIX_Assoc <- modelWide$NFIX_Assoc[match(clockSites$name, modelWide$response)]


#f1 <- ggplot(clockSites, aes(x = NFIX_Assoc, y = abs(V4), fill = NFIX_Assoc)) + 
#        stat_summary(fun = "mean", geom = "col", size = 1) +
#        geom_jitter(width = 0.2, color = "grey30") + stat_pwc(label = "p = {p}", y.position = 0.0115) +
#        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, linewidth = 2) + 
#        labs(x = "NFIX Associated?", y = "Absolute Clock Coefficient") + 
#        scale_fill_manual(values = c("blue", "red")) +  
#        theme(legend.position = "none")
#
#grid.newpage()
#grid.draw(set_panel_size(f1, width  = unit(1.5, "in"), height = unit(3, "in")))
#
#
#ggsave("./Figures/6F.1.svg",
#       plot = f1, device = "svg", width = 1.5, height = 3, units = "in")



f2 <- ggplot(clockSites, aes(x = NFIX_Assoc, y = abs(GenotypeEffect), fill = NFIX_Assoc)) + 
        stat_summary(fun = "mean", geom = "col", size = 1) +
        geom_jitter(width = 0.2, color = "grey30") + stat_pwc(label = "p = {p}") +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2, linewidth = 2) + 
        labs(x = "NFIX Associated?", y = "Absolute Genotype Effect") + 
        scale_fill_manual(values = c("blue", "red")) +  
        theme(legend.position = "none")

grid.newpage()
grid.draw(set_panel_size(f2, width  = unit(1.5, "in"), height = unit(3, "in")))


ggsave("./Figures/6F.2.svg",
       plot = f2, device = "svg", width = 1.5, height = 3, units = "in")


library(GenomicRanges)

meta <- readxl::read_xlsx("./Code_Data/2022 Aging cohort.xlsx", 
                          sheet = "Dissect", na = c("N/A", "NA", ""))[,-c(15:25)]


meta$ID <- paste0("X", meta$ID)

meta$AgeBin <- factor(cut(meta$Age, c(0, 6, 16, 30), labels = c("Young", "Mid", "Old")))

NFIXproms <- read.table("./Code_Data/NFIX_Proms.tsv", header = TRUE)

NFIXproms.gr <- GRanges(NFIXproms$Chr, IRanges(NFIXproms$Start, NFIXproms$Stop))

coords.df <- as.data.frame(str_split_fixed(coords, ":", 2))

cpg.gr <- GRanges(coords.df$V1, IRanges(as.numeric(coords.df$V2), as.numeric(coords.df$V2)))

overlaps <- findOverlaps(cpg.gr, NFIXproms.gr, type = "within")

nfixCpGs <- data.frame("Prom" = NFIXproms$Region[overlaps@to],
                       "CpG" = coords[overlaps@from], 
                       CpG[overlaps@from,])

nfixCpGs.long <- pivot_longer(nfixCpGs, cols = starts_with("X"), 
                              names_to = "ID", values_to = "Methylation") %>% na.omit()

nfixCpGs.long$AgeBin <- meta$AgeBin[match(nfixCpGs.long$ID, meta$ID)]

nfixPromMeth <- nfixCpGs.long %>% group_by(Prom, ID) %>% 
  summarise("AvgMeth" = mean(Methylation, na.rm = TRUE))

## Compare to whole gene expression

cpm <- fread("./genecounts_cpm.csv")

cpm$V1 <- str_remove_all(cpm$V1, "\\.bam") %>% str_replace_all("ID", "X")

nfixPromMeth$NFIX_Expr <- cpm$LOC101164563[match(nfixPromMeth$ID, cpm$V1)]

ggplot(nfixPromMeth, aes(x = AvgMeth, y = NFIX_Expr)) + 
  geom_point() + geom_smooth(method = "lm") + 
  labs(x = "Region Methylation", y = "NFIX CPM") +
  facet_wrap(~Prom, scales = "free_x")

nfixPromMeth$AgeBin <- meta$AgeBin[match(nfixPromMeth$ID, str_replace_all(meta$ID, "ID", "X"))]
nfixPromMeth$Age <- meta$Age[match(nfixPromMeth$ID, str_replace_all(meta$ID, "ID", "X"))]

lmProm1 <- lm(NFIX_Expr ~ Age + AvgMeth, data = filter(nfixPromMeth, Prom == "Prom1"))
lmProm2 <- lm(NFIX_Expr ~ Age + AvgMeth, data = filter(nfixPromMeth, Prom == "PromInt1"))

sjPlot::tab_model(lmProm1, lmProm2, dv.labels = c("Prom1 on NFIX", "PromInt1 on NFIX"))


library(GenomicFeatures)
library(txdbmaker)
library(GenomicAlignments)
library(Rsamtools)
library(rtracklayer)


txdb = makeTxDbFromGFF("./Code_Data/NFIX_Exons.gff")

flattenedAnnotation = exonicParts(txdb, linked.to.single.gene.only=TRUE)

names(flattenedAnnotation) =
  sprintf("%s:E%0.3d", flattenedAnnotation$gene_id, flattenedAnnotation$exonic_part)


library(DEXSeq)

### Bam files consist of just NFIXb region extracted from full STAR alignment

bam_files <- list.files("./Code_Data/splicing", pattern = ".bam$", full.names = TRUE)

sample_names <- bam_files %>% str_remove_all("ROI\\.bam") %>% str_remove_all("./Code_Data/splicing/")

## Get meta data

meta <- readxl::read_xlsx("./Code_Data/2022 Aging cohort.xlsx", 
                          sheet = "Dissect", na = c("N/A", "NA", ""))[,-c(15:25)]


meta$ID <- paste0("ID", meta$ID)

meta$AgeBin <- factor(cut(meta$Age, c(0, 6, 16, 30), labels = c("Young", "Mid", "Old")))

meta.rna <- meta[match(sample_names, meta$ID),]

# ---- Count reads per exon ----

se = summarizeOverlaps(
  flattenedAnnotation, BamFileList(bam_files), singleEnd=FALSE, fragments=TRUE, ignore.strand=TRUE)

colData(se)$condition = meta.rna$AgeBin
colData(se)$libType = factor("paired-end")
dxd = DEXSeqDataSetFromSE(se, design= ~ sample + exon + condition:exon)

countMat <- counts(dxd)[,1:85]

totCounts <- colSums(countMat)

### Compare Exon 1 Reads to total reads

raw_counts <- fread("./unfiltered_genecounts.csv")

counts.rna <- as.data.frame(raw_counts)[match(meta.rna$ID, str_remove_all(raw_counts$V1, "\\.bam")),]

exon1Expr <- data.frame("ID" = str_replace_all(meta.rna$ID, "ID", "X"),
                        "Exon1Abs" = countMat[1,],
                        "Exon1CPM" = countMat[1,] / rowSums(counts.rna[,-c(1:4)]),
                        "Exon1Perc" = 100 * countMat[1,] / totCounts)


nfixPromMeth$Exon1CPM <- exon1Expr$Exon1CPM[match(nfixPromMeth$ID, exon1Expr$ID)]
nfixPromMeth$Exon1Perc <- exon1Expr$Exon1Perc[match(nfixPromMeth$ID, exon1Expr$ID)]


ggplot(nfixPromMeth, aes(x = AvgMeth, y = Exon1Perc)) + 
  geom_point() + geom_smooth(method = "lm") + 
  labs(x = "Region Methylation", y = "Exon 1 Usage \n(as % of all NFIX reads)") +
  facet_wrap(~Prom, scales = "free_x")


lmProm1 <- lm(Exon1Perc ~ Age*AvgMeth, data = filter(nfixPromMeth, Prom == "Prom1"))
lmProm2 <- lm(Exon1Perc ~ Age*AvgMeth, data = filter(nfixPromMeth, Prom == "PromInt1"))

sjPlot::tab_model(lmProm1, lmProm2, dv.labels = c("Prom1 on Exon1", "PromInt1 on Exon1"))

#### Test relationship between age and Exon usage

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )

dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition")

dxr1 = DEXSeqResults( dxd )

res <- as.data.frame(dxr1)

sjPlot::tab_df(res[,1:12])

#res$transcripts[which(res$padj < 0.1)]

#### Test relationship between GVAN and Exon usage

meta.seq <- meta[match(sample_names, meta$ID),]

nfix_gt <- read.csv("./Code_Data/NFIX_RAD23_Genotypes.csv")

meta.seq$GT <- factor(nfix_gt$GT[match(meta.seq$ID, str_replace_all(nfix_gt$ID, "X", "ID"))])

nfixPromMeth$GT <- meta.seq$GT[match(nfixPromMeth$ID, str_replace_all(meta.seq$ID, "ID", "X"))]

I <- ggplot(filter(nfixPromMeth, !is.na(GT)), aes(x = Age, y = Exon1Perc, color = GT)) + 
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Age", y = "Exon 1 Usage \n(as % of all NFIX reads)", color = "GT")  +
  scale_color_manual(values = c("blue", "purple", "red"), labels = c("GG", "GA", "AA"))


grid.newpage()
grid.draw(set_panel_size(I, width  = unit(3.5, "in"), height = unit(3, "in")))


ggsave("./Figures/6I.svg",
       plot = I, device = "svg", width = 3.5, height = 3, units = "in")


lmGT <- lm(Exon1Perc ~ Age + as.numeric(GT), data = filter(nfixPromMeth, Prom == "Prom1"))

summary(lmGT)

sjPlot::tab_model(lmGT)


bam_files <- bam_files[-which(is.na(meta.seq$GT))]
meta.seq <- meta.seq[-which(is.na(meta.seq$GT)),]

se = summarizeOverlaps(flattenedAnnotation, BamFileList(bam_files), singleEnd=FALSE, fragments=TRUE, ignore.strand=TRUE)

colData(se)$condition = meta.seq$GT
colData(se)$libType = factor("paired-end")
dxd = DEXSeqDataSetFromSE(se, design= ~ sample + exon + condition:exon)

countMat <- counts(dxd)[,1:66]

totCounts <- colSums(countMat)

#keepWidth <- width(dxd@rowRanges) > 30
#keepDepth <- rowMeans(countMat) > 1

#dxd <- dxd[keepWidth & keepDepth,]

dxd = estimateSizeFactors( dxd )
dxd = estimateDispersions( dxd )
dxd = testForDEU( dxd )

dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition")

dxr1 = DEXSeqResults( dxd )

res <- as.data.frame(dxr1)


#### GVAN effects on promoter methylation

H <- ggplot(filter(nfixPromMeth, !is.na(GT), 
              Prom %in% c("Prom1", "Prom2")), aes(x = Age, y = AvgMeth, color = GT)) + 
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Age", y = "Promoter Methylation", color = "GT")  +
  facet_wrap(~Prom, labeller = labeller("Prom" = c("Prom1" = "DMR 1", "Prom2" = "DMR 2"))) + 
  coord_cartesian(ylim = c(0, 100)) + theme(legend.position = "top") +
  scale_color_manual(values = c("blue", "purple", "red"), labels = c("GG", "GA", "AA"))


grid.newpage()
grid.draw(set_panel_size(H, width  = unit(4, "in"), height = unit(3, "in")))


ggsave("./Figures/6H.svg",
       plot = H, device = "svg", width = 4, height = 3, units = "in")


lmGT.prom2 <- lm(AvgMeth ~ Age + as.numeric(GT), data = filter(nfixPromMeth, Prom == "Prom2"))
summary(lmGT.prom2)

lmGT.prom1 <- lm(AvgMeth ~ Age + as.numeric(GT), data = filter(nfixPromMeth, Prom == "Prom1"))
summary(lmGT.prom1)

### NFIX GT effect on NFIX / RAD23a expression

#NFIXb
ggplot(filter(nfixPromMeth, !is.na(GT), 
                   Prom == "Prom1"), aes(x = Age, y = NFIX_Expr, color = GT)) + 
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Age", y = "NFIX Expression", color = "GT")  +
  coord_cartesian(ylim = c(0, 100)) + theme(legend.position = "top") +
  scale_color_manual(values = c("blue", "purple", "red"), labels = c("GG", "GA", "AA"))


lmGT.nfixExpr <- lm(NFIX_Expr ~ Age + as.numeric(GT), data = filter(nfixPromMeth, Prom == "Prom1"))
summary(lmGT.nfixExpr)

nfixPromMeth$RAD23a_Expr <- cpm$rad23a[match(nfixPromMeth$ID, cpm$V1)]

#Rad23a
ggplot(filter(nfixPromMeth, !is.na(GT), 
              Prom == "Prom1"), aes(x = Age, y = RAD23a_Expr, color = GT)) + 
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Age", y = "RAD23a Expression", color = "GT")  +
  coord_cartesian(ylim = c(0, 100)) + theme(legend.position = "top") +
  scale_color_manual(values = c("blue", "purple", "red"), labels = c("GG", "GA", "AA"))


lmGT.rad23Expr <- lm(RAD23a_Expr ~ Age + as.numeric(GT), data = filter(nfixPromMeth, Prom == "Prom1"))
summary(lmGT.rad23Expr)



#### NFIX effect on other genes

cor.nfix <- psych::corr.test(as.matrix(cpm[,-c(1:4)]), y = cpm$LOC101164563, method = "spearman",
                              ci = FALSE, use = "pairwise")
cor.nfix.res <- data.frame(cor.nfix[1:6])

cor.nfix.sig <- cor.nfix.res[which(cor.nfix.res$p.adj < 0.05),]

cor.nfix.res[which(rownames(cor.nfix.res) == "mtor"),]

#### GVAN effect on mtor 

nfixPromMeth$mtor <- cpm$mtor[match(nfixPromMeth$ID, cpm$V1)]

ggplot(filter(nfixPromMeth, !is.na(GT), 
              Prom %in% c("Prom1")), aes(x = NFIX_Expr, y = mtor)) + 
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "NFIXb CPM", y = "mTOR CPM")

J <- ggplot(filter(nfixPromMeth, !is.na(GT), 
              Prom %in% c("Prom1")), aes(x = NFIX_Expr, y = mtor, color = GT)) + 
        geom_point() + geom_smooth(method = "lm") +
        labs(x = "NFIXb CPM", y = "mTOR CPM", color = "GT")  +
        scale_color_manual(values = c("blue", "purple", "red"), labels = c("GG", "GA", "AA"))


grid.newpage()
grid.draw(set_panel_size(J, width  = unit(3.5, "in"), height = unit(3, "in")))


ggsave("./Figures/6J.svg",
       plot = J, device = "svg", width = 3.5, height = 3, units = "in")


nfixPromMeth$Gonad_sex <- meta$Gonad_sex[match(nfixPromMeth$ID, str_replace_all(meta$ID, "ID", "X"))]

mtor.lm <- lm(mtor ~ NFIX_Expr*as.numeric(GT) + Gonad_sex + Age, data = filter(nfixPromMeth, !is.na(GT), Prom %in% c("Prom1")))
summary(mtor.lm)

