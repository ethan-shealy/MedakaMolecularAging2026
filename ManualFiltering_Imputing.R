library(data.table)


#### SET PARAMETERS

setwd("")

methBase <- "./methylBase_allSites_5x.txt"

covgThreshold <- 5

sampleThreshold <- 60


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


#ggplot(data.frame("AvgCovg" = AvgCovg), aes(x = AvgCovg)) + geom_histogram(color = "black") + 
#  scale_x_log10(n.breaks = 6) + xlab("Average Coverage at CpGs (among samples present)")

#hist(log(AvgCovg), xlab = "Log-Average Coverage at CpGs (among samples present)")

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

cCols <- seq(6, length(methDat), 3)
tCols <- seq(7, length(methDat), 3)

nC <- methDat[, ..cCols]
nT <- methDat[, ..tCols]

rm(methDat)

percMeth <- nC*100 / (nC+nT)

AvgMeth <- rowMeans(percMeth, na.rm = TRUE)

rm(nC, nT)

covStats <- data.frame("chr" = coords$V1,
                       "pos" = coords$V2,
                       "AvgMeth" = AvgMeth,
                       "AvgCov" = AvgCovg)


#ggplot(covStats, aes(x = AvgMeth, y = AvgCov)) + 
 #         geom_point() + geom_density2d() + geom_smooth(method = "lm")

#covStats.long <- tidyr::pivot_longer(covStats, cols = c(3,4), names_to = "Measure", values_to = "Value")

#ggplot(covStats, aes(x = pos, y = AvgCov)) + 
 #       scattermore::geom_scattermore() + facet_wrap(~chr, scales = "free_x")

covBed <- data.frame("chr" = coords$V1,
                     "start" = coords$V2,
                     "end" = coords$V2,
                     "AvgCov" = AvgCovg)

#fwrite(covBed, file = paste0("./", covgThreshold, "x_", sampleThreshold, "ind_coverageMap.bedGraph"),
  #     sep = "\t", col.names = FALSE, quote = FALSE)

methBed <- data.frame("chr" = coords$V1,
                     "start" = coords$V2,
                     "end" = coords$V2,
                     "AvgMeth" = AvgMeth)

#fwrite(methBed, file = paste0("./", covgThreshold, "x_", sampleThreshold, "ind_methMap.bedGraph"),
 #    sep = "\t", col.names = FALSE, quote = FALSE)

colnames(percMeth) <- sampleIDs

hist(rowMeans(percMeth, na.rm = TRUE), breaks = 99, 
     main = "Covg Cutoff = 5x, 0%", xlab = "Avg Methylation")

final_passed <- data.frame("coord" = paste0(coords$V1[-badSites], ":", coords$V2[-badSites]), percMeth[-badSites,])

hist(rowMeans(percMeth[-badSites,], na.rm = TRUE), breaks = 99, 
     main = paste0("Covg Cutoff = ", covgThreshold, "x, ", sampleThreshold, "L"), 
     xlab = "Avg Methylation")

fwrite(final_passed, file = paste0("./", covgThreshold, "x_", sampleThreshold, "ind_manualFilt.tab"), sep = "\t",
       col.names = TRUE, quote = FALSE)

### Imputation

missingBySample <- colSums(is.na(final_passed)) / nrow(final_passed)

hist(missingBySample, breaks = 12)

highMissingSamples <- which(missingBySample > 0.4)

final_passed <- final_passed[,-highMissingSamples]

missingBySample2 <- colSums(is.na(final_passed)) / nrow(final_passed)

hist(missingBySample2, breaks = 12, xlim = c(0, 1))

## Total missingness
sum(is.na(final_passed)) ## 3,915,594

# % missing 
sum(is.na(final_passed)) / (nrow(final_passed)*(ncol(final_passed)-1)) ## 13.5%

## mean imputation

meanimpute <- function(x) ifelse(is.na(x), mean(x, na.rm=T), x)
methImp <- apply(final_passed[,-1], 1, meanimpute)

dat.imp <- data.frame(final_passed$coord, t(methImp))

colnames(dat.imp) <- colnames(final_passed)

fwrite(dat.imp, file = paste0("./", covgThreshold, "x_", sampleThreshold, "ind_meanImputed.tab"), sep = "\t",
       col.names = TRUE, quote = FALSE)
