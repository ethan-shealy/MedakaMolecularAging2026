library(data.table)


#### SET PARAMETERS

setwd("")

methBase <- "./methylBase_allSites_5x.txt"

covgThreshold <- 5

sampleThreshold <- 30


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
colnames(nC_passed)[seq(5, 196, 2)] <- paste0("N_", sampleIDs)
colnames(nC_passed)[seq(6, 196, 2)] <- paste0("X_", sampleIDs)

head(nC_passed)
gc()

coords.short <- paste0(nC_passed$chr, ":", nC_passed$pos)

### Import and overlap
library(rtracklayer)

genes <- import.gff3("./Code_Data/ASM223467v1.gff")
head(genes)

genes <- genes[which(genes$type == "gene")]

library(dplyr)

promoters <- promoters(genes, upstream = 1000, downstream = 1000)

meth.gr <- GRanges(seqnames = nC_passed$chr, IRanges(start = nC_passed$pos, end = nC_passed$pos2), strand = "*")

head(meth.gr)

meth.genes <- findOverlaps(meth.gr, promoters, type = "within")

head(meth.genes)

meth.overlap <- data.frame("CpG" = coords.short[meth.genes@from],
                           "Gene" = genes$Name[meth.genes@to],
                           nC_passed[meth.genes@from,-c(1:4)])

meth.grouped <- meth.overlap %>%
  group_by(Gene) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE) %>%
  as.data.frame()

fwrite(meth.grouped, file = "./DSS/Promoters/allSamples_PromGroupedCovg.tab",
       col.names = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)

for (i in 1:length(sampleIDs)) {
  print(sampleIDs[i])
  df <- na.omit(data.frame("chr" = meth.grouped$Gene, 
                   "pos" = 1,
                   "N" = meth.grouped[seq(2, 193, 2)][,i],
                   "X" = meth.grouped[seq(3, 193, 2)][,i]))
  
  fwrite(df, file = paste0("./DSS/Promoters/", sampleIDs[i], ".dss"),
         col.names = TRUE, quote = FALSE, sep = "\t", row.names = FALSE)
}

rm(list=ls())

#### Now run DSS on these
library(data.table)
library(DSS)

methBase <- "./methylBase_allSites_5x.txt"


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
                     "Gonad_sex" = meta$Gonad_sex)

datFiles <- as.list(paste0("./DSS/Promoters/", sampleIDs, ".dss"))

for (i in 1:length(sampleIDs)) {
  assign(paste0("X", sampleIDs[i]), fread(input = file.path(datFiles[i]), header=TRUE, sep = "\t"))
}


BSobj <- makeBSseqData(list(X1,X2,X3,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,X22,X24,X25,X26,X27,
                       X28,X29,X30,X31,X32,X33,X34,X35,X36,X37,X38,X39,X40,X41,X42,X43,X44,X46,X47,X49,X50,
                       X51,X52,X53,X54,X56,X59,X60,X61,X65,X68,X70,X71,X72,X73,X74,X75,X76,X77,X78,X81,X82,X83,
                       X85,X87,X88,X90,X91,X94,X95,X98,X99,X100,X101,X102,X103,X108,X110,X111,X114,X115,X117,
                       X119,X120,X123,X125,X130,X132,X133,X134,X137,X140), 
                       sampleNames = as.list(paste0("X", sampleIDs)))


DMLfit = DMLfit.multiFactor(BSobj = BSobj, design=design, formula=~ log(Age) + Gonad_sex + log(Age):Gonad_sex)


### Now perform tests

DMLtest.int <- DMLtest.multiFactor(DMLfit, coef=4)
sig.int.sites <- DMLtest.int[which(DMLtest.int$fdrs < 0.1),]
write.table(sig.int.sites, file = "./DSS/Promoters/SigLogIntProms.tab",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)


BSobj <- BSobj[-which(DMLtest.int$fdrs < 0.1),]

DMLfit = DMLfit.multiFactor(BSobj = BSobj, design=design, formula=~ log(Age) + Gonad_sex)


DMLtest.age <- DMLtest.multiFactor(DMLfit, coef=2)
sig.age.sites <- DMLtest.age[which(DMLtest.age$fdrs < 0.1),]
write.table(sig.age.sites, file = "./DSS/Promoters/SigLogAgeProms.tab",
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
