library(data.table)


#### SET PARAMETERS

setwd("")

methBase <- "./methylBase_allSites_5x.txt"

covgThreshold <- 5

sampleThreshold <- 50


#### Run Script

methDat <- fread(methBase, skip = 14)

head(methDat)

methMeta <- fread(methBase, nrows = 13, sep = "\t")

sampleIDs <- unlist(stringr::str_split(methMeta[3,], ";"))

sampleIDs[1] <- stringr::str_remove(sampleIDs[1], "#SI:")

#sampleIDs <- paste0("S", sampleIDs)

coords <- methDat[,c(1:4)]

dat <- methDat[,-c(1:4)]

rm(methDat)

cvgCols <- seq(1, length(dat), 3)

covg <- dat[, ..cvgCols]

range(covg, na.rm = TRUE)

rm(dat)

gc()

head(covg)

AvgCovg <- rowMeans(covg, na.rm = TRUE)

## Look at missingness before changing threshold

missing <- rowSums(is.na(covg))

hist(missing, breaks = 96, xlab = "Sites Missing in X Samples")

thresh <- which(covg < covgThreshold, arr.ind = TRUE)

covg[thresh] <- NA

## Look at new missingness

missing <- rowSums(is.na(covg))

hist(missing, breaks = 36)

badSites <-  which(missing > length(sampleIDs)-sampleThreshold)

gc()

if (any(missing > length(sampleIDs)-sampleThreshold)){
  covg_passed <- covg[-badSites,]
} else {
  covg_passed <- covg
}

rm(covg)

gc()

hist(rowMeans(covg_passed, na.rm = TRUE), breaks = 100)
hist(log(rowMeans(covg_passed, na.rm = TRUE)), breaks = 100)

range(covg_passed, na.rm = TRUE)


methDat <- fread(methBase, skip = 14)

tCols <- seq(7, length(methDat), 3)

nC_cov <- methDat[, -..tCols]

nC_passed <- nC_cov[-badSites,]
rm(nC_cov)

colnames(nC_passed)[1:4] <- c("chr", "pos", "pos2", "strand")
colnames(nC_passed)[seq(5, length(nC_passed), 2)] <- paste0("N_", sampleIDs)
colnames(nC_passed)[seq(6, length(nC_passed), 2)] <- paste0("X_", sampleIDs)

head(nC_passed)
gc()

coords.short <- paste0(nC_passed$chr, ":", nC_passed$pos)

nC_passed <- as.data.frame(nC_passed)

fwrite(nC_passed, file = "./DSS/allSamples_Covg_5x_50L.tab",
       col.names = TRUE, sep = "\t", quote = FALSE)

for (i in 1:length(sampleIDs)) {
  print(sampleIDs[i])
  df <- na.omit(data.frame("chr" = nC_passed$chr, 
                   "pos" = nC_passed$pos,
                   "N" = nC_passed[seq(5, length(nC_passed), 2)][,i],
                   "X" = nC_passed[seq(6, length(nC_passed), 2)][,i]))
  
  fwrite(df, file = paste0("./DSS/dataFiles/", sampleIDs[i], ".dss"),
         col.names = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)
}

rm(list=ls())


### Now run DSS on these
library(tidyverse)
library(data.table)
library(DSS)

meta <- readxl::read_xlsx("./Code_Data/2022 Aging cohort.xlsx", 
                          sheet = "Dissect", na = c("N/A", "NA", ""))[,-c(15:25)]

methMeta <- fread(methBase, nrows = 13, sep = "\t")

sampleIDs <- unlist(stringr::str_split(methMeta[3,], ";"))

sampleIDs[1] <- stringr::str_remove(sampleIDs[1], "#SI:")

meta <- meta[which(meta$ID %in% sampleIDs),]

ImInd <- which(meta$Gonad_sex == "Im")

meta <- meta[-ImInd,]
sampleIDs <- sampleIDs[-ImInd]

X <- model.matrix(~ Age + Gonad_sex + Age:Gonad_sex, meta)

X

design <- data.frame("ID" = sampleIDs,
                     "Age" = meta$Age,
                     "Sex" = meta$Gonad_sex)

datFiles <- as.list(paste0("./DSS/dataFiles/", sampleIDs, ".dss"))

## Get values for methylation sites at age-related sites

methData <- fread("./DSS/allSamples_Covg_5x_50L.tab",
                  header = TRUE, sep = "\t", data.table = FALSE)


N <- as.matrix(methData[,seq(5, length(methData), 2)])
C <- as.matrix(methData[,seq(6, length(methData), 2)])

mat <- t(C*100 / N)

rownames(mat) <- stringr::str_remove_all(rownames(mat), "X_")

colnames(mat) <- paste0(methData$chr, ".", methData$pos)


for (i in 1:length(sampleIDs)) {
  print(sampleIDs[i])
  assign(paste0("X", sampleIDs[i]), fread(input = file.path(datFiles[i]), header=TRUE, sep = "\t"))
}



BSobj <- makeBSseqData(list(X1,X2,X3,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,X22,X24,X25,X26,X27,
                            X28,X29,X30,X31,X32,X33,X34,X35,X36,X37,X38,X39,X40,X41,X42,X43,X44,X46,X47,X49,X50,
                            X51,X52,X53,X54,X56,X59,X60,X61,X65,X68,X70,X71,X72,X73,X74,X75,X76,X77,X78,X81,X82,X83,
                            X85,X87,X88,X90,X91,X94,X95,X98,X99,X100,X101,X102,X103,X108,X110,X111,X114,X115,X117,
                            X119,X120,X123,X125,X130,X132,X133,X134,X137,X140), 
                       sampleNames = as.list(paste0("X", sampleIDs)))

rm(X1,X2,X3,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,X22,X24,X25,X26,X27,
   X28,X29,X30,X31,X32,X33,X34,X35,X36,X37,X38,X39,X40,X41,X42,X43,X44,X46,X47,X49,X50,
   X51,X52,X53,X54,X56,X59,X60,X61,X65,X68,X70,X71,X72,X73,X74,X75,X76,X77,X78,X81,X82,X83,
   X85,X87,X88,X90,X91,X94,X95,X98,X99,X100,X101,X102,X103,X108,X110,X111,X114,X115,X117,
   X119,X120,X123,X125,X130,X132,X133,X134,X137,X140)

#### first run model with interaction term

design$Sex <- factor(design$Sex, levels = c("F", "M"))


DMLfit = DMLfit.multiFactor(BSobj = BSobj, design=design, formula=~ log(Age) + Sex + log(Age):Sex,
                            smoothing = TRUE, smoothing.span = 100)


### Now perform tests

DMLtest.int <- DMLtest.multiFactor(DMLfit, coef=4)
sig.int.sites <- DMLtest.int[which(DMLtest.int$fdrs < 0.05),]
sig.int.sites$name <- paste0(sig.int.sites$chr, ".", sig.int.sites$pos)

write.table(sig.int.sites, file = "./DSS/intSites_smoothed100.tab",
            sep = "\t", col.names = TRUE)

intSites <- read.table("./DSS/intSites_smoothed100.tab", 
                       sep = "\t", header = TRUE)

BSobj.coords <- paste0(decode(BSobj@rowRanges@seqnames), ".", start(BSobj@rowRanges@ranges))

BSobj <- BSobj[-which(BSobj.coords %in% intSites$name),]

#### Now run the model with no interaction term

DMLfit = DMLfit.multiFactor(BSobj = BSobj, design=design, formula=~ log(Age) + Sex,
                            smoothing = TRUE, smoothing.span = 100)

## Now perform tests

head(DMLfit$X)

DMLtest.age <- DMLtest.multiFactor(DMLfit, coef=2)
sig.age.sites <- DMLtest.age[which(DMLtest.age$fdrs < 0.05),]
sig.age.sites$name <- paste0(sig.age.sites$chr, ".", sig.age.sites$pos)


DMLtest.sex <- DMLtest.multiFactor(DMLfit, coef = 3)
sig.sex.sites <- DMLtest.sex[which(DMLtest.sex$fdrs < 0.05),]
sig.sex.sites$name <- paste0(sig.sex.sites$chr, ".", sig.sex.sites$pos)

ggVennDiagram::ggVennDiagram(list("AgeSites" = sig.age.sites$name, "SexSites" = sig.sex.sites$name)) + ggplot2::coord_flip()

fwrite(sig.age.sites, file = "./DSS/SigLogAgeCpGs.tab", sep = "\t")


