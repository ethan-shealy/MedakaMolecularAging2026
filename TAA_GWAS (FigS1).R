library(matrixStats)
library(tidyverse)
library(data.table)
library(ggpubr)

setwd("")

theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=7),
                  axis.title=element_text(size=9,face="bold"),
                  strip.text = element_text(size=9,face="bold"),
                  panel.border = element_rect(fill = NA)))

### Get metadata

meta <- readxl::read_xlsx("./Code_Data/2022 Aging cohort.xlsx", 
                          sheet = "Dissect", na = c("N/A", "NA", ""))[,-c(15:25)]

meta$ID <- paste0("X", meta$ID)

expr_res <- readxl::read_excel("./Clocks/LOOCV_Performance.xlsx")
expr_res$ID <- str_remove_all(expr_res$ID, "\\.bam") %>% str_replace_all("ID", "X")

expr_res.error <- lm(Predicted_age ~ age, data = expr_res)

meta$TAA <- residuals(expr_res.error)[match(meta$ID, expr_res$ID)]
meta$TA_Est <- expr_res$Predicted_age[match(meta$ID, expr_res$ID)]

plot(meta$Age, meta$TAA)

meta <- meta[-which(meta$TAA > 40),]
meta <- meta[-which(meta$Gonad_sex == "Im"),]

ggplot(meta, aes(x = Age, y = TA_Est)) + geom_point() + geom_smooth(method = "lm")

jitterer <- position_jitter(width = .25, seed = 1999)

ggplot(meta, aes(x = Age, y = TA_Est)) + 
  geom_smooth(method = "lm") + 
  geom_segment(aes(yend = TA_Est - TAA, xend = Age), color = "red", 
               linewidth = 2, alpha = 0.25) + geom_point() + 
  labs(x = "Actual Age (months)", y = "DNAm LOOCV Estimated Age")

resids <- ggplot(meta, aes(x = Age, y = TA_Est)) + 
  geom_smooth(method = "lm", color = "black", se = FALSE) + 
  geom_segment(aes(yend = TA_Est - TAA, xend = Age, color = TAA), 
               linewidth = 1.5, alpha = 0.75, position = jitterer) + 
  geom_point(position = jitterer) + 
  labs(x = "Actual Age (months)", y = "DNAm LOOCV Estimated Age") + 
  scale_color_gradient2(low = "blue", mid = "grey", high = "red") + 
  theme_bw()

library(grid)
library(egg)

grid.newpage()
grid.draw(set_panel_size(resids, width  = unit(2.5, "in"), height = unit(2, "in")))

ggsave("./Figures/TAA_LOOCV_Resids.svg",
       plot = resids, device = "svg", width = 3.5, height = 2.2, units = "in")

hist(meta$TAA, breaks = 12)

### Get genotype matrix
library(vcfR)
library(adegenet)

var <- read.vcfR("./AgeCohort_SNPs_variant_final.vcf.gz")

gl <- vcfR2genlight(var)

gl@ind.names <- str_remove_all(gl@ind.names, "\\.bam") %>% str_replace_all("ID", "X")
gl@ind.names

meta.seq <- meta[which(meta$ID %in% gl@ind.names),]
meta.seq <- meta.seq[-which(is.na(meta.seq$TAA)),]

gl <- gl[match(meta.seq$ID, gl@ind.names),]

codedDF <- data.frame("SNPID" = paste0(gl@chromosome, "_", gl@position),
                      t(tab(gl, NA.method = "zero")))

fwrite(codedDF, "./Code_Data/genotype_matrix.txt",
       sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

SNPcoords <- data.frame("SNPID" = paste0(gl@chromosome, "_", gl@position),
                        "chr" = gl@chromosome, "pos" = gl@position)

fwrite(SNPcoords, "./Code_Data/snp_locations.txt",
       sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# File paths 
genotype_file <- "./Code_Data/genotype_matrix.txt"     # SNPs x samples (0, 1, 2)
snp_loc_file <- "./Code_Data/snp_locations.txt"        # SNPID, chr, pos


# Load genotype data
geno_dt <- fread(genotype_file)
head(geno_dt)

geno_mat <- t(geno_dt[,-1])
rownames(geno_mat) <- colnames(geno_dt)[-1]
colnames(geno_mat) <- geno_dt$SNPID

snp_loc <- fread(snp_loc_file)
head(snp_loc)

### Proceed with analysis
library(rrBLUP)
library(qqman)

K <- A.mat(geno_mat)  # Additive genetic relationship matrix
str(K)

heatmap(K)

stopifnot(rownames(K) == rownames(geno_mat))

geno_df <- data.frame(snp_loc, as.matrix(geno_dt[,-1]) - 1) ### Convert 0, 1, 2 genotype notation to -1, 0, 1 by subtracting 1

gwas_results <- GWAS(
  pheno = as.data.frame(meta.seq[,c(1, 4, 6, 15)]),
  geno = geno_df,
  fixed = c("Age", "Gonad_sex"),
  K = K,  
  n.PC = 0,
  min.MAF = 0.1,
  plot = TRUE
)

gwas_bed <- data.frame("chr" = gwas_results$chr,
                       "start" = gwas_results$pos - 1,
                       "end" = gwas_results$pos,
                       "score" = gwas_results$TAA)

fwrite(gwas_bed, "./TAA_GWAS.bedGraph", sep = "\t",
       col.names = FALSE, quote = FALSE)

chrSizes <- read.table("./Code_Data/ASM223467v1_chromSizes.txt")

chrSizes$cum <- c(0, cumsum(chrSizes$V2)[-length(chrSizes$V2)])

gwas_bed$ChrOffset <- chrSizes$cum[match(gwas_bed$chr, chrSizes$V1)] + gwas_bed$start

manhat <- ggplot(filter(gwas_bed, chr != "NC_004387.1"), aes(x = ChrOffset, y = score, color = chr)) + 
  geom_point(show.legend = FALSE) + 
  scale_x_continuous(breaks = chrSizes$cum[-1], labels = seq(1, length(chrSizes$V1[-1]), 1)) + 
  scale_color_manual(values = rep(c("grey70", "grey20"), 12)) + 
  theme(panel.grid.minor = element_blank()) + 
  labs(x = "Chromosomal Location", y = "-log10(p) for TAA Association")


grid.newpage()
grid.draw(set_panel_size(manhat, width  = unit(7, "in"), height = unit(2, "in")))

ggsave("./Figures/TAA_manhattan.tiff",
       plot = manhat, device = "tiff", width = 7.25, height = 2.3, units = "in")


gwas_sorted <- gwas_results[order(gwas_results$TAA, decreasing = TRUE),]

head(gwas_sorted)

gwas_filt <- gwas_sorted[which(gwas_sorted$TAA > 0),]
gwas_filt$TAA_p <- 10^-gwas_filt$TAA

## Genomic inflation test
hist(gwas_filt$TAA_p, breaks = 99)

# Convert p-values to chi-squared statistics
gwas_chisq <- qchisq(1 - gwas_filt$TAA_p, df = 1)

# Calculate lambda (genomic inflation factor)
lambda_gc <- median(gwas_chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
print(lambda_gc)

qq <- ggplot(gwas_filt, aes(x = -log10(ppoints(length(TAA_p))), y = -log10(TAA_p))) + 
  geom_point(size = 0.7) + geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + 
  annotate(x = 3.8, y = 1, geom = "text", label = paste0("Lambda \nGC = ", round(lambda_gc, 3)), size = 2.5) + 
  labs(x = "Expected \n-log10(p)", y = "Observed \n-log10(p)") + theme(axis.text = element_text(size = 8))

grid.newpage()
grid.draw(set_panel_size(qq, width  = unit(2.5, "in"), height = unit(2, "in")))

ggsave("./Figures/TAA_qq.svg",
       plot = qq, device = "svg", width = 3, height = 2.2, units = "in")

ggsave("./Figures/TAA_qq_panel.tiff",
       plot = qq, device = "tiff", width = 1.5, height = 1.2, units = "in")

## Get GTs with missing values present

gl <- vcfR2genlight(var)

gl@ind.names <- str_remove_all(gl@ind.names, "\\.bam") %>% str_replace_all("ID", "X")
gl@ind.names

gl <- gl[match(meta.seq$ID, gl@ind.names),]

codedDF <- data.frame("SNPID" = paste0(gl@chromosome, "_", gl@position),
                      t(tab(gl, NA.method = "asis")))

stopifnot(colnames(codedDF)[-1] == meta.seq$ID)


SNP_to_plot <- "NC_019874.2_13037322"
var@fix[which(paste0(var@fix[,1], "_", var@fix[,2]) == SNP_to_plot),]

meta.seq$GT <- as.factor(unlist(codedDF[which(colnames(geno_mat) == SNP_to_plot),-1]))

GT_TAA <- ggplot(filter(meta.seq, !is.na(GT)), aes(x = GT, y = TAA, color = GT)) + 
  geom_boxplot(outliers = FALSE) + geom_jitter(width = 0.1) + 
  scale_color_manual(values = c("blue", "red"), labels = c("CC", "CA"))


grid.newpage()
grid.draw(set_panel_size(GT_TAA, width  = unit(2, "in"), height = unit(3, "in")))

ggsave("./Figures/GT_TAA.svg",
       plot = GT_TAA, device = "svg", width = 2.5, height = 3, units = "in")


GT_EA <- ggplot(filter(meta.seq, !is.na(GT)), aes(x = Age, y = TA_Est, color = GT)) + 
  geom_point() + geom_smooth(method = "lm", alpha = 0.2) + 
  scale_color_manual(values = c("blue", "red"), labels = c("CC", "CA"))


grid.newpage()
grid.draw(set_panel_size(GT_EA, width  = unit(4, "in"), height = unit(3, "in")))

ggsave("./Figures/GT_TA_Est.svg",
       plot = GT_EA, device = "svg", width = 5, height = 3, units = "in")


lm.test <- lm(TA_Est ~ Age + Gonad_sex + as.numeric(as.character(GT)), data = meta.seq)
summary(lm.test)


### What genes DO show association in expression with the identified genotype?

######### All genes' expression levels - which ones do the SNP affect?

raw <- as.data.frame(fread("./unfiltered_genecounts.csv", check.names = FALSE))

raw$V1 <- str_remove_all(raw$V1, "ID") %>% 
  str_remove_all("\\.bam") %>%
  paste0("X", .)

rownames(raw) <- raw$V1

raw.dat <- na.omit(raw[match(meta$ID, raw$V1),-c(1,2,3)])

dds <- DESeqDataSetFromMatrix(countData = t(raw.dat),
                              colData = meta.seq,
                              design = ~ 1)

# Variance stabilizing transformation
library(DESeq2)

vsd <- vst(dds, blind = TRUE)
expr_matrix <- data.frame("gene" = rownames(vsd), assay(vsd))

age_genes <- readxl::read_excel("./Code_Data/DESEQ2_all.xlsx")

expr_matrix <- t(expr_matrix[match(age_genes$Gene, expr_matrix$gene),-1])
expr_matrix[1:10,1:10]
colnames(expr_matrix) <- age_genes$Gene

## Linear model approach
meta.modeling <- meta.seq %>% mutate(GT = as.numeric(as.character(GT))) %>% filter(!is.na(meta.seq$GT))

missingGenes <- unique(which(is.na(expr_matrix), arr.ind = TRUE)[,2])

expr_matrix <- expr_matrix[match(meta.modeling$ID, rownames(expr_matrix)),-missingGenes]

modelResults <- list()

for (i in 1:ncol(expr_matrix)){
  if(i %% 1000 == 0) {print(i)}
  modelResults[[i]] <- lm(expr_matrix[,i] ~ GT*Age + Gonad_sex, data = meta.modeling)
}

extract_summ <- function(model, extractVariable) {
  summ <- summary(model)
  coefs <- summ$coefficients
  if (extractVariable %in% rownames(coefs)) {
    return(data.frame(beta = coefs[extractVariable, "Estimate"], 
                      p_value = coefs[extractVariable, "Pr(>|t|)"]))
  } else {
    return(data.frame(beta = NA, 
                      p_value = NA))
  }
}

GTeffect_expr <- lapply(modelResults, extract_summ, extractVariable = "GT:Age")

GTeffect_expr.df <- data.frame("Gene" = age_genes$Gene[-missingGenes],
                               do.call(rbind, GTeffect_expr))

GTeffect_expr.df$p_adjusted <- p.adjust(GTeffect_expr.df$p_value, method = "fdr")

qq(GTeffect_expr.df$p_value)

GTeffect_mat.df <- data.frame(meta.modeling, 
                              expr_matrix[,which(GTeffect_expr.df$p_adjusted < 0.05)])

ggplot(GTeffect_mat.df, aes(x = as.factor(GT), y = expr, color = as.factor(GT))) + geom_boxplot() + geom_point(color = "black") + 
  labs(x = "Genotype", y = "IGF2 log2(CPM)", color = "GT") +
  scale_color_manual(values = c("blue", "purple", "red"), labels = c("CC", "CA", "AA"))

ggplot(GTeffect_mat.df, aes(x = Age, y = expr, color = as.factor(GT))) + 
  geom_point(size = 2) + geom_smooth(method = "lm") + 
  labs(x = "Age", y = "IGF2 log2(CPM)", color = "GT") + 
  scale_color_manual(values = c("blue", "purple", "red"), labels = c("CC", "CA", "AA"))



######## All aDMC sites - are any affected by identified genotype?

methData <- fread("./5x_50ind_manualFilt.tab", sep = "\t", 
                  header = TRUE, data.table = FALSE)

methData$coord <- make.names(methData$coord)

methData <- methData[,c(1, which(colnames(methData) %in% meta.seq$ID))]

## now get age-related sites

ageDMCs <- read.table("./DSS/SigLogAgeCpGs.tab", sep = "\t", header = TRUE)

meth_dt <- methData[which(methData$coord %in% ageDMCs$name),]
colnames(meth_dt)[1] <- "CpG"

meth_matrix <- t(meth_dt[,-1]); meth_matrix[1:10,1:10]
colnames(meth_matrix) <- meth_dt$CpG

meth_matrix <- meth_matrix[match(meta.modeling$ID, rownames(meth_matrix)),]

modelResults <- list()

for (i in 1:ncol(meth_matrix)){
  if(i %% 1000 == 0) {print(i)}
  modelResults[[i]] <- lm(meth_matrix[,i] ~ GT*Age + Gonad_sex, data = meta.modeling)
}

extract_summ <- function(model, extractVariable) {
  summ <- summary(model)
  coefs <- summ$coefficients
  if (extractVariable %in% rownames(coefs)) {
    return(data.frame(beta = coefs[extractVariable, "Estimate"], 
                      p_value = coefs[extractVariable, "Pr(>|t|)"]))
  } else {
    return(data.frame(beta = NA, 
                      p_value = NA))
  }
}

GTeffect <- lapply(modelResults, extract_summ, extractVariable = "GT:Age")

GTeffect.df <- data.frame("CpG" = meth_dt$CpG,
                          do.call(rbind, GTeffect))

GTeffect.df$p_adjusted <- p.adjust(GTeffect.df$p_value, method = "fdr")

qq(GTeffect.df$p_value)

GTeffect_mat.df <- data.frame(meta.modeling, 
                              meth_matrix[,which(GTeffect.df$p_adjusted < 0.05)]) %>%
  pivot_longer(cols = starts_with("NC_"), names_to = "CpG", values_to = "Methylation")

ggplot(GTeffect_mat.df, aes(x = as.factor(GT), y = Methylation, color = as.factor(GT))) + geom_boxplot() + geom_point(color = "black") + 
  labs(x = "Genotype", y = "Methylation %", color = "GT") +
  facet_wrap(~CpG, scales = "free_y") +
  scale_color_manual(values = c("blue", "red"), labels = c("CC", "CA"))

TAA_SNP_effect_Meth <- ggplot(GTeffect_mat.df, aes(x = Age, y = Methylation, color = as.factor(GT))) + 
                          geom_point(size = 2, alpha = 0.6, position = position_jitter(height = 2)) + geom_smooth(method = "lm") + 
                          labs(x = "Age", y = "Methylation %", color = "GT") + 
                          facet_wrap(~CpG, scales = "free_y") +
                          scale_color_manual(values = c("blue", "red"), labels = c("CC", "CA")) + 
                          theme(legend.position = "inside", legend.position.inside = c(0.8, 0.15))

grid.newpage()
grid.draw(set_panel_size(TAA_SNP_effect_Meth, width  = unit(1, "in"), height = unit(0.8, "in")))

ggsave("./Figures/GT_TA_MethEffect.svg",
       plot = TAA_SNP_effect_Meth, device = "svg", width = 3.5, height = 3, units = "in")
