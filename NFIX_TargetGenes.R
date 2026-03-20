library(tidyverse)
library(ggpubr)
library(data.table)
library(grid)
library(egg)

setwd("")


theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=7),
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

meta.seq <- meta[match(rownames(raw.dat), meta$ID),]



cpm.dat <- t(log2(edgeR::cpm(t(raw.dat)) + 1))


allGeneNames <- colnames(cpm.dat)


contAgeGenes <- readxl::read_excel("./Code_Data/DESEQ2_all.xlsx")

library(GenomicRanges)
library(rtracklayer)
library(BSgenome)
library(memes)


oryLat <- readDNAStringSet("./Code_Data/ASM223467v1.fa", format="fasta", use.names=TRUE)
oryLat@ranges@NAMES <- str_remove_all(oryLat@ranges@NAMES, " Oryzias.*$")


annot_full <- import.gff("./Code_Data/ASM223467v1.gff")

geneDat <- annot_full[which(annot_full$type == "gene"),]

#genes_measured <- geneDat[which(geneDat$gene %in% allGeneNames)]

geneDat.aDEGs <- geneDat[which(geneDat$gene %in% contAgeGenes$Gene)]

aDEG.proms <- promoters(geneDat.aDEGs, upstream = 5000, downstream = 2000)

AgeProms <- get_sequence(aDEG.proms, oryLat)
AgeProms
write_fasta(AgeProms, path = "./DSS/aDEG_seq7000.fa")

## MEME Result:

## MA1528.1	NFIX	p = 3.08e-1	q = 3.08e-1	e = 3.08e-1	TP = 591 / 1385 (42.7%)	FP = 577 / 1385 (41.7%)	Enrichment = 1.02	Thresh = 8.68



geneDat.randGenes <- geneDat[sample(1:length(geneDat), length(AgeProms))]

rand.proms <- promoters(geneDat.randGenes, upstream = 5000, downstream = 2000)

RandProms <- get_sequence(rand.proms, oryLat)
RandProms
write_fasta(RandProms, path = "~/Parrott_Lab/medaka_EMseq/DSS/RandomGene_seq7000.fa")


### Same test for aDMCs

ageSites <- fread("./DSS/SigLogAgeCpGs.tab")

ageSites.gr <- GRanges(ageSites$chr, IRanges(ageSites$pos - 250, ageSites$pos + 250), strand = "*")


AgeSeqs <- get_sequence(ageSites.gr, oryLat)
AgeSeqs
write_fasta(AgeSeqs, path = "./DSS/aDMC_seq250.fa")

## MEME result:

## MA1528.1	NFIX	p = 2.81e-7	q = 2.81e-7	e = 2.81e-7	TP = 564 / 21399 (2.6%)	FP = 409 / 21399 (1.9%)	Enrichment = 1.38	Thresh = 10.93




### random CpGs


#randSites <- fread("E:/medaka/MedClock2/5x_60ind_manualFilt.tab")[,1]
#
#randSites.coord <- as.data.frame(str_split_fixed(randSites$coord, ":", 2))[sample(1:nrow(randSites), nrow(ageSites)),]
#
#randSites.gr <- GRanges(randSites.coord$V1, 
#                        IRanges(as.numeric(randSites.coord$V2) - 250, as.numeric(randSites.coord$V2) + 250), 
#                        strand = "*")
#
#
#RandSeqs <- get_sequence(randSites.gr, oryLat)
#RandSeqs
#write_fasta(RandSeqs, path = "./DSS/Random_seq250.fa")



### repeat again with overlap DMCs

overlapSites <- fread("./DSS/overlapDMCs.tsv")

overlapSites.gr <- GRanges(overlapSites$Chr, IRanges(overlapSites$Start - 250, overlapSites$Start + 250), strand = "*")


overlapSeqs <- get_sequence(overlapSites.gr, oryLat)
overlapSeqs
write_fasta(overlapSeqs, path = "./DSS/overlapDMCs_seq250.fa")

## MEME result:

## MA1528.1	NFIX	p = 3.08e-1	e = 3.08e-1	q = 3.08e-1	TP = 591 / 1385 (42.7%)	FP = 577 / 1385 (41.7%)	Enrichment = 1.02	Thresh = 8.68

### overlapDEGs

overlapgenes <- fread("./DSS/overlapDEGs.tsv")

overlapgenes <- geneDat.aDEGs[which(geneDat.aDEGs$gene %in% overlapgenes$Var1)]

overlapgenesProms <- promoters(overlapgenes, upstream = 5000, downstream = 2000)

overlapGeneSeqs <- get_sequence(overlapgenesProms, oryLat)
overlapGeneSeqs
write_fasta(overlapGeneSeqs, path = "./DSS/overlapDEGs_seq7000.fa")

## MEME result: 

## MA1528.1	NFIX	p = 9.30e-2	q = 9.30e-2	e = 9.30e-2	TP = 234 / 526 (44.5%)	FP = 568 / 1385 (41.0%)	Ernichment = 1.09	Thresh = 8.71


### Split aDMCs into increasing and decreasing

ageSites.gr.dec <- ageSites.gr[which(ageSites$stat < 0)]

AgeSeqs.dec <- get_sequence(ageSites.gr.dec, oryLat)
AgeSeqs.dec
write_fasta(AgeSeqs.dec, path = "./DSS/aDMC_decreasing_seq250.fa")

## MEME result: MA1528.1	NFIX	p = 6.87e-17	q = 6.87e-17	e = 6.87e-17	TP = 501 / 13886 (3.6%)	FP = 454 / 21399 (2.1%)	Enrichment = 1.70	Threshold = 10.63


ageSites.gr.inc <- ageSites.gr[which(ageSites$stat > 0)]

AgeSeqs.inc <- get_sequence(ageSites.gr.inc, oryLat)
AgeSeqs.inc
write_fasta(AgeSeqs.inc, path = "./DSS/aDMC_increasing_seq250.fa")

## MEME result: MA1528.1	NFIX	p = 5.38e-6	q = 5.38e-6	e = 5.38e-6	TP = 1973 / 7514 (26.3%)	FP = 5072 / 21399 (23.7%)	Enrichment = 1.11	Thresh = 0.11



#### Read in all sites with NFIX motifs

sites <- fread("./Code_data/aDMC_NFIsites.tsv")
seqs <- fread("./Code_data/aDMC_NFIsequences.tsv")

seqs$Stat <- ageSites$stat[match(seqs$seq_ID, names(AgeSeqs))]
seqs$sig <- ageSites$pvals[match(seqs$seq_ID, names(AgeSeqs))]


table(sign(seqs$Stat))
plot(-log10(seqs$sig), seqs$seq_Score)

seqs.coord <- str_split_fixed(seqs$seq_ID, "[:-]", 3)

seqs <- cbind(seqs, seqs.coord)
seqs$V2 <- as.numeric(seqs$V2)
seqs$V3 <- as.numeric(seqs$V3)


chrSizes <- read.table("./Code_Data/ASM223467v1_chromSizes.txt")

chrSizes$cum <- c(0, cumsum(chrSizes$V2)[-length(chrSizes$V2)])


seqs$ChrOffset <- chrSizes$cum[match(seqs$V1, chrSizes$V1)] + seqs$V2

seqs.bed <- data.frame(seqs$V1, seqs$V2, seqs$V3)

#fwrite(seqs.bed, file = "~/NFIX_Binding.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

ggplot(seqs, aes(x = ChrOffset, y = seq_Score, color = V1)) + 
  geom_point(show.legend = FALSE) + labs(x = "Chromosomal Location", y = "Motif Score") + 
  scale_color_manual(values = c("grey60", rep(c("grey20", "grey60"), 12))) + 
  scale_x_continuous(breaks = chrSizes$cum[-1], labels = chrSizes$V1[-1]) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0), axis.ticks.x = element_blank())


ggplot(seqs, aes(x = V2)) + 
  geom_density(show.legend = FALSE) + labs(x = "Chromosomal Location", y = "Motif Score") + 
  scale_color_manual(values = c("grey60", rep(c("grey20", "grey60"), 12))) + 
  facet_wrap(~V1, nrow = 24, scales = "free_x", strip.position = "left")

library(ggridges)

ggplot(seqs, aes(x = V2, y = V1, fill = V1)) + 
  geom_density_ridges(rel_min_height = 0.005) + 
  scale_fill_manual(values = c("grey60", rep(c("grey20", "grey60"), 12))) + 
  labs(x = "Chromosomal Location", y = "Motif Density")


ggplot(seqs, aes(x = V2, y = V1, fill = V1, height = after_stat(density))) + 
  geom_density_ridges(rel_min_height = 0.005, stat = "density") + 
  scale_fill_manual(values = c("grey60", rep(c("grey20", "grey60"), 12))) + 
  labs(x = "Chromosomal Location", y = "Motif Density")








#### Does NFIX genotype influence methylation at NFIX-targeted aDMCs? Compared to non-NFIX targeted aDMCs?

nfix_aDMCs <- na.omit(ageSites[match(seqs$seq_ID, names(AgeSeqs))])

non_nfix_aDMCs <- ageSites[-which(ageSites$name %in% nfix_aDMCs$name),]

meth <- fread("./5x_50ind_manualFilt.tab", sep = "\t", header = TRUE)
meth$coord <- str_replace_all(meth$coord, ":", ".")


meth.DMCs <- as.data.frame(t(meth[match(ageSites$name, meth$coord),-1]))
colnames(meth.DMCs) <- ageSites$name

meta.gt <- read.csv("./Code_Data/NFIX_RAD23_Genotypes.csv")
meta.gt <- meta.gt[-which(is.na(meta.gt$GT)),]

meth.DMCs <- meth.DMCs[match(meta.gt$ID, rownames(meth.DMCs)),]

meta.DMCs <- data.frame(meta.gt, meth.DMCs) %>%
  pivot_longer(cols = starts_with("NC_"), names_to = "CpG", values_to = "Methylation")

## Apply linear models to determine effect of GT*Age interaction at each site
library(broom)

results <- map(meth.DMCs, ~ lm(.x ~ Age*as.numeric(GT), data = meta.gt))

tidy_results <- map2_df(results, names(results), ~ tidy(.x) |> mutate(response = .y))

modelWide <- pivot_wider(tidy_results, names_from = "term", values_from = "estimate", id_cols = "response")

modelWide.err <- pivot_wider(tidy_results, names_from = "term", values_from = "std.error", id_cols = "response")

modelWide.sig <- pivot_wider(tidy_results, names_from = "term", values_from = "p.value", id_cols = "response")


stopifnot(modelWide$response == colnames(meth.DMCs))

modelWide$NFIX_Target <- ifelse(modelWide$response %in% nfix_aDMCs$name, "Yes", "No")

modelWide$Dir <- sign(modelWide$Age) ### Split by Age direction

modelWide$GT.Error <- modelWide.err$`as.numeric(GT)`
modelWide$Int.Error <- modelWide.err$`Age:as.numeric(GT)`
modelWide$GT.p <- modelWide.sig$`as.numeric(GT)`
modelWide$Int.p <- modelWide.sig$`Age:as.numeric(GT)`

modelWide <- modelWide[-which(modelWide$Age == 0),]

### Look for differences in Age effect across target status

ggplot(modelWide, aes(x = NFIX_Target, y = Age, fill = NFIX_Target))  +
  geom_jitter(width = 0.2, color = "grey30", alpha = 0.1) + geom_boxplot() + 
  stat_summary(fun = "mean", geom = "point", size = 3, color = "black") +
  stat_pwc(label = "p = {p}") + labs(x = "NFIX Motif?", y = "Effect of Age") + 
  scale_fill_manual(values = c("blue", "red")) + 
  facet_wrap(~Dir, nrow = 1, scales = "free_x",
             labeller = labeller("Dir" = c("-1" = "Decreasing with Age", "1" = "Increasing with Age")))

### Look for differences in Genotype effect across target status

ggplot(modelWide, aes(x = NFIX_Target, y = `as.numeric(GT)`, fill = NFIX_Target))  +
  geom_jitter(width = 0.2, color = "grey30", alpha = 0.1) + geom_boxplot() + 
  stat_summary(fun = "mean", geom = "point", size = 3, color = "black") +
  stat_pwc(label = "p = {p}") + labs(x = "NFIX Motif?", y = "Genotype Effect") + 
  scale_fill_manual(values = c("blue", "red")) 


ggplot(modelWide, aes(x = NFIX_Target, y = `as.numeric(GT)`, fill = NFIX_Target))  +
  facet_wrap(~Dir, nrow = 1, scales = "free_y",
             labeller = labeller("Dir" = c("-1" = "Decreasing with Age", "1" = "Increasing with Age"))) + 
  geom_jitter(width = 0.2, color = "grey30", alpha = 0.1) + 
  geom_boxplot(outliers = FALSE) + 
  stat_summary(fun = "mean", geom = "point", size = 3, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_pwc(label = "p = {p}") + labs(x = "NFIX Motif?", y = "Genotype Effect") + 
  scale_fill_manual(values = c("blue", "red")) + theme(legend.position = "none")

modelWide$sigGT <- ifelse(modelWide$GT.p < 0.01, "Y", "N")

table(modelWide$sigGT, modelWide$NFIX_Target)
table(modelWide$sigGT, modelWide$NFIX_Target, modelWide$Dir)


### Repeat for interaction effect


ggplot(modelWide, aes(x = NFIX_Target, y = Age, fill = NFIX_Target))  +
  facet_wrap(~Dir, nrow = 1, scales = "free_y",
             labeller = labeller("Dir" = c("-1" = "Decreasing with Age", "1" = "Increasing with Age"))) + 
  #geom_jitter(width = 0.2, color = "grey30", alpha = 0.1) + 
  geom_boxplot(outliers = FALSE) + 
  stat_summary(fun = "mean", geom = "point", size = 3, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_pwc(label = "p = {p}", y.position = 0.8, bracket.shorten = 0.5, vjust = 2) + labs(x = "NFIX Motif?", y = "Rate of Change in Ref. Hom.") + 
  scale_fill_manual(values = c("blue", "red")) + theme(legend.position = "none")


ggplot(modelWide, aes(x = NFIX_Target, y = Age + 2*`Age:as.numeric(GT)`, fill = NFIX_Target))  +
  facet_wrap(~Dir, nrow = 1, scales = "free_y",
             labeller = labeller("Dir" = c("-1" = "Decreasing with Age", "1" = "Increasing with Age"))) + 
  #geom_jitter(width = 0.2, color = "grey30", alpha = 0.1) + 
  geom_boxplot(outliers = FALSE) + 
  stat_summary(fun = "mean", geom = "point", size = 3, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_pwc(label = "p = {p}", y.position = 0.8, bracket.shorten = 0.5, vjust = 2) + labs(x = "NFIX Motif?", y = "Rate of Change in Alt. Hom.") + 
  scale_fill_manual(values = c("blue", "red")) + theme(legend.position = "none")


modelWide$sigInt <- ifelse(modelWide$Int.p < 0.01, "Y", "N")

table(modelWide$sigInt, modelWide$NFIX_Target); fisher.test(table(modelWide$sigInt, modelWide$NFIX_Target))
table(modelWide$sigInt, modelWide$NFIX_Target, modelWide$Dir)

#### Compare hom. ref. and alt. ref. in same plot

modelWide$RefRate <- modelWide$Age
modelWide$AltRate <- modelWide$Age + 2*modelWide$`Age:as.numeric(GT)`

modelWide.long <- pivot_longer(modelWide, cols = ends_with("Rate"), 
                               names_to = "GT", values_to = "Rate")

modelWide.long$GT <- factor(modelWide.long$GT, levels = c("RefRate", "AltRate"))

ggplot(modelWide.long, aes(x = GT, y = Rate, fill = NFIX_Target))  +
  facet_wrap(~Dir, nrow = 1, scales = "free_y",
             labeller = labeller("Dir" = c("-1" = "Decreasing with Age", "1" = "Increasing with Age"))) + 
  #geom_jitter(width = 0.2, color = "grey30", alpha = 0.1) + 
  scale_x_discrete(labels = c("Hom. Ref.", "Hom. Alt.")) +
  geom_boxplot(outliers = FALSE, color = "black") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_pwc(hide.ns = TRUE, label = "p = {p}{p.signif}", y.position = 0.9, bracket.shorten = 0.7, tip.length = 0.01) + 
  labs(x = "NFIX Motif?", y = "Rate of Change at aDMCs \n(%methylation / month)") + 
  scale_fill_viridis_d(name = "NFIX Target?") + theme(legend.position = "top")
