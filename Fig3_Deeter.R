#TASK 3: SPEARMAN CORRELATION MODEL
#Load packages
library(psych)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(clipr)
library(edgeR)
library(DESeq2)
library(ggrepel)
library(psych)
library(tidyr)
library(openxlsx) 
library(statmod) #for robust dispersion
library(chromoMap)


#Load Transcriptome_KnownAge_CountData_Dec2024.RData from file path
load("./Code_Data/Transcriptome_KnownAge_CountData_Dec2024.RData")

################################################################################
#TASK 4: DESEQ2 MODEL - CONTINUOUS MODEL 2

#Order of operations:
#1. Run GEX ~ age + sex + interaction
#2. Find genes with significant interaction
#3. Remove significant interaction terms from df
#4. Re-run model without interaction

library(DESeq2)
###
#1. Run GEX ~ age + sex + interaction
###
#4.1: Prepare data
View(counts_raw) #use the raw counts, NOT the CPM
samples <- rownames(genecounts_raw) 
View(samples)

#4.2: Create colData (df containing sample names, age_scaled and sex)
colData <- data.frame(sample = samples,
                      age_log = genecounts_raw$age_log, 
                      sex= genecounts_raw$sex) 

colData$sex <- factor(colData$sex, levels = c("M","F"))
colData$age_log <- as.numeric(colData$age_log)
View(colData)

#4.3 Create DESeqDataSet object
dds_age_log <-  DESeqDataSetFromMatrix(countData = counts_raw,
                                       colData = colData, 
                                       design = ~age_log + sex + sex*age_log) 


#4.4 Run DESeq2
dds_age_log <- DESeq(dds_age_log) #performs normalization, estimation of dispersion, and differential expression analysis
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

?DESeq()

#4.5 Results
resultsNames(dds_age_log) #View coefficient names so that they can be specified properly during interpretation

#4.5.1: AGE_LOG
results_age_log_model2 <- results(dds_age_log, name = "age_log")
results_age_log_model2 <- as.data.frame(results_age_log_model2)
results_age_log_sig_model2 <- results_age_log_model2[which(results_age_log_model2$padj < 0.05), ] #325 genes
results_age_log_sig_model2$Gene <- rownames(results_age_log_sig_model2)
View(results_age_log_sig_model2) #Confirmed

#Filter genes increasing, decreasing with age
results_age_log_model2_increasing <- results_age_log_sig_model2[which(results_age_log_sig_model2$log2FoldChange > 0), ]
View(results_age_log_model2_increasing) #192 genes increase with age, confirmed
results_age_log_model2_decreasing <- results_age_log_sig_model2[which(results_age_log_sig_model2$log2FoldChange < 0), ]
View(results_age_log_model2_decreasing) #133 genes decrease with age, confirmed

#4.5.2: SEX (sex_F_vs_M)
results_sex_model2 <- results(dds_age_log, name="sex_F_vs_M")
results_sex_model2 <- as.data.frame(results_sex_model2)
results_sex_model2_significant <- results_sex_model2[which(results_sex_model2$padj < 0.05), ] #90 genes

#Filter genes higher in females, higher in males
results_sex_model2_significant_females_greater <- results_sex_model2_significant[which(results_sex_model2_significant$log2FoldChange > 0), ]
View(results_sex_model2_significant_females_greater) #73 genes higher in females
results_sex_model2_significant_males_greater <- results_sex_model2_significant[which(results_sex_model2_significant$log2FoldChange < 0), ]
View(results_sex_model2_significant_males_greater) #17 genes higher in males

###
#2. Find genes with significant interaction
###

#4.5.3: AGE_LOG*SEX INTERACTION
results_interaction_model2 <- results(dds_age_log, name="age_log.sexF")
results_interaction_model2 <- as.data.frame(results_interaction_model2)
results_interaction_model2_significant <- results_interaction_model2[which(results_interaction_model2$padj < 0.05), ]
View(results_interaction_model2_significant) #5 genes have a significant age*sex interaction

#The 5 genes are: (1) LOC101160321 (2) LOC101165969 (3) LOC101173743 (4) mb (5) slc25a48

###
#3. Remove the following genes from counts_raw
counts_raw2 <- trans_readcounts_pass_raw #
counts_raw2 <- lapply(counts_raw, as.numeric)
counts_raw2 <- as.data.frame(counts_raw)
rownames(counts_raw2) <- rownames(trans_readcounts_pass_raw)
colnames(counts_raw2) <- NULL
View(counts_raw2)

genecounts_raw <- as.data.frame(trans_readcounts_pass_raw)
genecounts_raw <- t(genecounts_raw) #rows= samples, columns = genes
genecounts_raw <- as.data.frame(genecounts_raw)
View(genecounts_raw)
#Remove columns of the following 5 genes
genecounts_raw2 <- genecounts_raw[ ,!colnames(genecounts_raw) %in% c("LOC101160321", "LOC101165969", "LOC101173743", "mb", "slc25a48")]

genecounts_raw2 <- data.frame(sex=sex_column, genecounts_raw2)
genecounts_raw2 <- data.frame(age_log = age_log_column, genecounts_raw2)
genecounts_raw2 <- data.frame(age = age_months_column, genecounts_raw2)


#Remove rows titled the above 5 genes:
# Remove rows where row names are in the specified list
counts_raw2 <- counts_raw[!rownames(counts_raw) %in% c("LOC101160321", "LOC101165969", "LOC101173743", "mb", "slc25a48"), ] #Confirmed, 16,543 observations


####
#4. Re-run model without interaction

#4.1: Prepare data
View(counts_raw2) #use the raw counts, NOT the CPM
samples2 <- rownames(genecounts_raw2) 
View(samples2)


#4.2: Create colData (df containing sample names, age_scaled and sex)
colData2 <- data.frame(sample = samples2,
                       age_log = genecounts_raw2$age_log, 
                       sex= genecounts_raw2$sex) 

colData2$sex <- factor(colData$sex, levels = c("M","F"))
colData2$age_log <- as.numeric(colData$age_log)
View(colData2)

#4.3 Create DESeqDataSet object
dds_age_log <-  DESeqDataSetFromMatrix(countData = counts_raw2,
                                       colData = colData2, 
                                       design = ~age_log + sex) #!!!!!!!!!! NO INTERACTION TERM !!!!!!!!!!!!


#4.4 Run DESeq2
dds_age_log <- DESeq(dds_age_log) #performs normalization, estimation of dispersion, and differential expression analysis
#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing

#4.5 Results
resultsNames(dds_age_log) #View coefficient names so that they can be specified properly during interpretation

#4.5.1: AGE_LOG
results_age_log_model2 <- results(dds_age_log, name = "age_log")
results_age_log_model2 <- as.data.frame(results_age_log_model2)
results_age_log_sig_model2 <- results_age_log_model2[which(results_age_log_model2$padj < 0.05), ] #1,547 genes
results_age_log_sig_model2$Gene <- rownames(results_age_log_sig_model2)
View(results_age_log_sig_model2)
#> resultsNames(dds_age_log) #View coefficient names so that they can be specified properly during interpretation
#[1] "Intercept"  "age_log"    "sex_F_vs_M"

write_xlsx(results_age_log_model2, "./Code_Data/DESEQ2_all.xlsx")



library(ComplexHeatmap)
library(ggpubr)
library(grid)
library(egg)

load("./Code_Data/heatmap_object_aDEGs.RData")

age_fun <- circlize::colorRamp2(c(0,30), c("green","purple"))
col_ha <- columnAnnotation(Age=genecounts_cpm$age, Sex= genecounts_cpm$sex,
                           col= list(Age=age_fun, Sex=c("F"="salmon","M"="skyblue2")))

figure1c1 <- Heatmap(t(genecounts_heatmap_scaled),
                     name= "ScaledCPM",
                     row_split = continuous_aDEGs$group,
                     column_order = order(genecounts_cpm$age),
                     show_row_dend = FALSE,
                     column_names_gp = gpar(fontsize=8),
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     bottom_annotation = col_ha,
                     use_raster = TRUE)

figure1c1


heat.1 <- grid.grab()

## draw it, changes optional
heat.1.new <- editGrob(heat.1, vp = viewport(width = unit(5.2, "in"), 
                                             height = unit(3.6, "in")))
grid.newpage()
grid.draw(heat.1.new)

ggsave("./Figures/3A.svg",
       plot = heat.1.new, device = "svg", width = 6.2, height = 4.2, units = "in")



