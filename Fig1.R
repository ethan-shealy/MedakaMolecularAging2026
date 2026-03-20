library(tidyverse)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(ggpubr)
library(psych)
library(grid)
library(egg)

setwd("")

theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=6),
                  axis.title=element_text(size=9,face="bold"),
                  strip.text = element_text(size=9,face="bold"),
                  panel.border = element_rect(fill = NA)))

### Just methylation clock

CpG.imputed <- as.data.frame(fread("./5x_60ind_meanImputed.tab", header = TRUE,
                                   colClasses = "numeric"))

coords <- CpG.imputed[,1]

coords.df <- as.data.frame(str_split_fixed(coords, ":", 2))

chrSizes <- read.table("./Code_Data/ASM223467v1_chromSizes.txt")

chrSizes$cum <- c(0, cumsum(chrSizes$V2)[-length(chrSizes$V2)])

CpG <- CpG.imputed[,-1]
rownames(CpG) <- CpG.imputed$coord
samples <- colnames(CpG)

## Metadata matching

dat <- readxl::read_xlsx("./Code_Data/2022 Aging cohort.xlsx", na = c("", "NA", "#N/A"))[,1:14]

dat$K <- (100 * dat$Tot_wgt) / (dat$Tot_length^3)

dat$ID <- paste0("X", dat$ID)

### Define age transformation after Horvath (2013)

adult.thresh <- 10

age_ranges <- c(1:60)

train.age.fun <- ifelse(age_ranges <= adult.thresh, 
                        log(age_ranges)-log(adult.thresh+0.01),
                        (age_ranges - adult.thresh) / (adult.thresh+0.01))

plot(age_ranges, train.age.fun)

inv_ages = approxfun(x=train.age.fun, y=age_ranges)

dat$AgeTrans <- train.age.fun[match(dat$Age, age_ranges)]

dat.seq <- dat[match(samples, dat$ID),]

cpg.mat <- data.frame(t(CpG))

set.seed(1998)
gc()

library(glmnet)

### Train / test

#index <- sample(1:length(medaka_age), 16) # Random sampling
index <- order(dat.seq$Age)[seq(1, length(dat.seq$Age), 5)] # Somewhat random sampling along range of ages

#dividing into training and test sets
test.CpG <- t(CpG)[index,]
test.age <- dat.seq$AgeTrans[index]
hist(test.age)

train.CpG <- t(CpG)[-index,] ## take all that are not in the test set
train.age <- dat.seq$AgeTrans[-index]
hist(train.age)

gc()

#using 10 fold cross validation to estimate the lambda parameter 
age.cv <- cv.glmnet(train.CpG, 
                    as.matrix(train.age),
                    alpha=0.5, nfolds= 10, family="gaussian", 
                    trace.it = TRUE)

plot(age.cv)

#define the lambda parameter
best_lambda <- age.cv$lambda.min

#coefficients for fitted model
age_coeff <- coef(age.cv, s=best_lambda)

#fitted values
fitted <- inv_ages(predict(age.cv, train.CpG, s=best_lambda))
actual.fitted <- data.frame(dat.seq[-index,], "s1" = fitted)

#calculate R2 for fitted data
R2.train <- cor(fitted, train.age)^2

##get the error associated with each predicted age
CpG.train_error <- mean(abs(actual.fitted$s1 - actual.fitted$Age))

#predicted values for test set
predicted <- inv_ages(predict(age.cv, test.CpG, s=best_lambda))

#actual by predicted values for training data
actual.predicted <- data.frame(dat.seq[index,], "s1" = predicted)

##get the error associated with each predicted age
CpG.predicted_age <- actual.predicted$s1
CpG.testerror <- (actual.predicted$s1 - actual.predicted$Age)
CpG.abs_error <- abs(CpG.testerror)

CpG.AAE <- mean(CpG.abs_error)
CpG.AAE

##PULLING OUT CLOCK SITES
coeffs <- coef(age.cv, s = best_lambda)
coeffs.df <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x)

coeffs.format <- data.frame(stringr::str_split_fixed(coeffs.df$name, ":", 2))[-1,]

options(scipen = 99)

coeffs.bed <- data.frame(coeffs.format$X1, as.numeric(coeffs.format$X2)-1, coeffs.format$X2, coeffs.df$coefficient[-1])
coeffs.bed <- rbind(c("Intercept", NA, NA, coeffs.df$coefficient[1]), coeffs.bed)


write.table(coeffs.bed, file = "./Clocks/ElasticNet_DNAm_AgeTrans_sites.bedGraph",
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

options(scipen = 3)

#### PLOT CLOCKs

##linear model for R2 and pvalue

trainlm <- lm(actual.fitted$s1 ~ actual.fitted$Age, data = actual.fitted)
summary(trainlm)

testlm <- lm(actual.predicted$s1 ~ actual.predicted$Age, data = actual.predicted)
summary(testlm)

##CpG Test

plot.CpG.test <- ggplot(data=actual.predicted, aes(x=Age, y=s1)) +
  geom_point() + 
  #geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="black") + 
  labs(x = "Actual Age", y = "Predicted Epigenetic Age") +
  ylim(2, 29) + 
  annotate("text", x = 8, y = 25, label = paste0("MAE: ", round(CpG.AAE, 2), " months"),
           size = 7, size.unit = "pt")

plot.CpG.test


grid.newpage()
grid.draw(set_panel_size(plot.CpG.test, width  = unit(1.5, "in"), height = unit(1.5, "in")))

ggsave("./Figures/1D.svg",
       plot = plot.CpG.test, device = "svg", width = 1.5, height = 1.5, units = "in")

#### Representative clock blow-up image

plot.CpG.test <- ggplot(data=actual.predicted, aes(x=Age, y=s1)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() + 
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="black") + 
  labs(x = "Actual Age", y = "Predicted Epigenetic Age") +
  ylim(2, 29)

plot.CpG.test


grid.newpage()
grid.draw(set_panel_size(plot.CpG.test, width  = unit(3.5, "in"), height = unit(3, "in")))

ggsave("./Figures/DNAmTest_BlowUp.svg",
       plot = plot.CpG.test, device = "svg", width = 4, height = 3.2, units = "in")


## Write output

actual.total <- rbind(actual.fitted, actual.predicted)
actual.total$Set <- c(rep("Train", nrow(actual.fitted)), rep("Test", nrow(actual.predicted)))

ggplot(data=actual.total, aes(x=Age, y=s1, color = Set)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(size = 2) + 
  geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, alpha = 0.2) + 
  labs(x = "Actual Age", y = "Predicted Epigenetic Age") +
  ylim(2, 29)


write.csv(actual.total, "./Clocks/AgeTrans_RepClock_res.csv")


############# DNAm LOOCV
rownames(CpG) <- coords

sites <- c()
Predicted_age <- rep(NA, nrow(dat.seq))

### Set age transformation

### Use 10-fold cross validation to estimate lambda

set.seed(1998) #setting seed
age.cv <- cv.glmnet(t(CpG), 
                    dat.seq$AgeTrans,
                    alpha=0.5, nfolds=10, family="gaussian",
                    trace.it = TRUE)

plot(age.cv)
#define the lambda parameter
best_lambda <- age.cv$lambda.min

##elastic net for loop -- modified from SG script
for (i in 1:nrow(dat.seq)){
  print(paste0("Running model ", i, " of ", nrow(dat.seq)))
  train.CpG <- t(CpG)[-i,] ##train is all except the 1 test sample per iteration
  test.CpG <- t(CpG)[i,] ##test is left out sample
  test.age <- dat.seq$AgeTrans[c(i)]
  train.age <- dat.seq$AgeTrans[c(-i)]
  
  ##run elastic net
  model <- glmnet(as.matrix(train.CpG), 
                  as.matrix(train.age),
                  alpha=0.5, family="gaussian", nlambda=100)
  #coefficients for fitted model
  age_coeff <- coef(model, s=best_lambda)
  #fitted values
  fitted <- inv_ages(predict(model, as.matrix(train.CpG), s=best_lambda))
  actual.fitted <- data.frame(fitted, train.age)
  #predicted values for test set
  predicted <- inv_ages(predict(model, t(test.CpG), s=best_lambda))
  #actual by predicted values for training data
  actual.predicted <- data.frame(predicted, test.age)
  ##get the error associated with each predicted age
  Predicted_age[i] <- actual.predicted$predicted
  age.coeffs <- coef(model, s =best_lambda)
  sites <- append(sites,(name = age.coeffs@Dimnames[[1]][age.coeffs@i + 1]))
}

##merge output into one dataset
cpgClock.test <- data.frame(dat.seq, "DNAmAge" = Predicted_age)

cpgClock.test$Test_Error <- cpgClock.test$DNAmAge - cpgClock.test$Age

write.csv(cpgClock.test, "./Clocks/AgeTrans_LOOCV_res.csv")

mean(abs(cpgClock.test$Test_Error)) # 3.17

c <- ggplot(cpgClock.test, aes(x = Age, y = DNAmAge)) + 
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() + labs(x = "Age (months)", y = "Predicted Age (CpG LOOCV)") +
  #geom_smooth(color = "black", method = "lm", alpha = 0.2) + 
  ylim(range(dat.seq$Age)) + 
  annotate(geom = "text", x = 9, y = 26, label = "MAE: 3.17mo", size = 9, size.unit = "pt")
c

grid.newpage()
grid.draw(set_panel_size(c, width  = unit(1.5, "in"), height = unit(1.5, "in")))

ggplot(cpgClock.test, aes(x = Age, y = DNAmAge, color = Gonad_sex)) + 
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point(size = 2) + labs(x = "Age (months)", y = "Predicted Age (CpG LOOCV)") +
  #geom_smooth(color = "black", method = "lm", alpha = 0.2) + 
  ylim(range(dat.seq$Age)) 

cpgClock.test$Resids <- residuals(lm(cpgClock.test$DNAmAge ~ cpgClock.test$Age))

ggplot(cpgClock.test, aes(x = Gonad_sex, y = Resids, fill = Gonad_sex)) + 
  geom_boxplot() + stat_pwc(label = "Wilcoxon, p = {p}") + labs(x = "Gonad Sex", y = "LOOCV Residuals")

ggsave("./Figures/1C.svg",
       plot = c, device = "svg", width = 1.5, height = 1.5, units = "in")


# Get all sites from LOOCV

sites.tab <- as.data.frame(table(sites))[-1,]

write.csv(sites.tab, "./Clocks/DNAm_loocvSites.csv", row.names = FALSE)



##### Gene-expression based clock

gexClockEstimates <- readxl::read_excel("./Clocks/test_age_estimates_8020_Clock.xlsx")

gexClockLOOCV <- readxl::read_excel("./Clocks/LOOCV_Performance.xlsx")

mean(abs(gexClockEstimates$predicted_age - gexClockEstimates$age), na.rm = TRUE) ## MAE: 3.54

gex.test <- ggplot(gexClockEstimates, aes(x = age, y = predicted_age)) + 
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() + labs(x = "Age (months)", y = "Predicted Age") +
  annotate(geom = "text", x = 9, y = 25, label = "MAE: 3.54mo", size = 9, size.unit = "pt") + 
  ylim(1, 29) + xlim(1, 29)

mean(abs(gexClockLOOCV$Predicted_age - gexClockLOOCV$age), na.rm = TRUE) ## MAE: 4.36

grid.newpage()
grid.draw(set_panel_size(gex.test, width  = unit(1.5, "in"), height = unit(1.5, "in")))

ggsave("./Figures/1B.svg",
       plot = gex.test, device = "svg", width = 1.5, height = 1.5, units = "in")


ggplot(gexClockLOOCV, aes(x = age, y = Predicted_age)) + 
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() + labs(x = "Age (months)", y = "Predicted Age") + #ylim(0, 32)
  annotate(geom = "text", x = 9, y = 25, label = "MAE: 4.36mo", size = 9, size.unit = "pt")

gex.loocv <- ggplot(gexClockLOOCV, aes(x = age, y = Predicted_age)) + 
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() + labs(x = "Age (months)", y = "Predicted Age") + ylim(0, 41) +
  annotate(geom = "text", x = 9, y = 35, label = "MAE: 4.36mo", size = 9, size.unit = "pt")

ggsave("./Figures/1A.svg",
       plot = gex.loocv, device = "svg", width = 1.5, height = 1.5, units = "in")


### Now a clock combining both predictors

cpm <- as.data.frame(fread("./genecounts_cpm.csv", check.names = FALSE))

cpm.meta <- cpm[,c(1:4)]

cpm.meta$V1 <- str_remove_all(cpm.meta$V1, "ID")
cpm.meta$V1 <- str_remove_all(cpm.meta$V1, "\\.bam")

cpm.meta$V1 <- paste0("X", cpm.meta$V1)

cpm.dat <- cpm[,-c(1:4, length(cpm))]

rownames(cpm.dat) <- cpm.meta$V1

rownames(cpg.mat)
rownames(cpm.dat)

lowVar <- which(apply(cpg.mat, 2, sd) == 0) #### Remove invariable sites so scaling can function

combined.mat <- merge(scale(cpg.mat[,-lowVar]), scale(cpm.dat), by = 0)

rownames(combined.mat) <- combined.mat$Row.names
combined.mat <- combined.mat[,-1]

combined.meta <- dat[match(rownames(combined.mat), dat$ID),]

# LOOCV
combo.cv <- cv.glmnet(as.matrix(combined.mat), 
                      combined.meta$AgeTrans,
                      alpha=0.5, nfolds= nrow(combined.meta), family="gaussian", 
                      trace.it = TRUE, keep = TRUE)

plot(combo.cv)

loocv.combo.res <- matrix(inv_ages(combo.cv$fit.preval), nrow = nrow(combined.mat))

combo.lambda.min <- which(combo.cv$lambda == combo.cv$lambda.min)


comboClock.test <- data.frame(combined.meta, "comboAge" = loocv.combo.res[,combo.lambda.min])

i <- ggplot(comboClock.test, aes(x = Age, y = comboAge)) + 
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() + labs(x = "Age (months)", y = "Predicted Age (Combined LOOCV)") + ylim(0, 30) + 
  annotate(geom = "text", x = 9, y = 29, label = "MAE: 3.45mo", size = 9, size.unit = "pt")

grid.newpage()
grid.draw(set_panel_size(i, width  = unit(1.5, "in"), height = unit(1.5, "in")))


ggsave("./Figures/1I.svg",
       plot = i, device = "svg", width = 1.5, height = 1.5, units = "in")


mean(abs(comboClock.test$comboAge - comboClock.test$Age), na.rm = TRUE) ## MAE: 3.45

write.csv(comboClock.test, "./Clocks/ComboClock_LOOCV_res.csv")


##PULLING OUT CLOCK SITES
coeffs <- coef(combo.cv, s = combo.cv$lambda.min)
coeffs.df <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x)[-1,]

coeffs.df$Type <- ifelse(grepl("NC_", coeffs.df$name), "CpG", "Gene")

j <- ggplot(coeffs.df, aes(x = reorder(name, coefficient), color = Type)) + 
  theme(axis.text.x = element_blank(), 
        legend.position = "none", panel.grid = element_blank()) + 
  labs(y = "Scaled Exponentiated Betas", x = "Predictor") + 
  scale_color_manual(values = c("green", "purple")) + #xlim(0.89, 1.11) + 
  geom_segment(aes(y = coefficient, yend = 0, x = reorder(name, coefficient)), linewidth = 2,) +
  geom_text(aes(y = ifelse(coefficient < 0, 0.01, -0.035), 
                label = ifelse(!grepl("NC_", name), name, "")), 
            angle = 90, color = "black", hjust = 0.2, size = 7, size.unit = "pt") 

grid.newpage()
grid.draw(set_panel_size(j, width  = unit(4.8, "in"), height = unit(1.5, "in")))


ggsave("./Figures/1J.svg",
       plot = j, device = "svg", width = 4.8, height = 1.5, units = "in")







######################### Promoter clock

proms <- fread("./DSS/Promoters/allSamples_PromGroupedCovg.tab", 
               data.table = FALSE)

proms.mat <- 100*t(data.frame(proms[,seq(3, length(proms), 2)] / proms[,seq(2, length(proms), 2)]))

samples <- rownames(proms.mat) %>% str_remove_all("_")

# filter by missingness

highMissing <- which(colSums(is.na(proms.mat)) > 0.2*nrow(proms.mat))

# impute by median

medianImpute <- function(x) {ifelse(is.na(x), median(x, na.rm = TRUE), x)}

proms.imp <- t(apply(proms.mat[,-highMissing], 2, medianImpute))


colnames(proms.imp) <- samples
rownames(proms.imp) <- proms$Gene[-highMissing]

head(proms.imp)

dat.clean <- dat[match(colnames(proms.imp), dat$ID),]

stopifnot(dat.clean$ID == colnames(proms.imp))

##break into age and methylation values


set.seed(1998)

index <- order(dat.clean$Age)[seq(1, length(dat.clean$Age), 5)] # Somewhat random sampling along range of ages

#dividing into training and test sets
test.prom <- t(proms.imp)[index,]
test.age <- dat.clean$AgeTrans[index]
hist(test.age)

train.prom <- t(proms.imp)[-index,] ## take all that are not in the test set
train.age <- dat.clean$AgeTrans[-index]
hist(train.age)

gc()

#using 10 fold cross validation to estimate the lambda parameter 
age.cv <- cv.glmnet(train.prom, 
                    as.matrix(train.age),
                    alpha=0.5, nfolds= 10, family="gaussian", 
                    trace.it = TRUE)

plot(age.cv)

#define the lambda parameter
best_lambda <- age.cv$lambda.min

#coefficients for fitted model
age_coeff <- coef(age.cv, s=best_lambda)

#fitted values
fitted <- inv_ages(predict(age.cv, train.prom, s=best_lambda))
actual.fitted <- data.frame(dat.clean[-index,], fitted)

#calculate R2 for fitted data
R2.train <- cor(fitted, train.age)^2

##get the error associated with each predicted age
actual.fitted$trainerrors <- (actual.fitted$fitted - actual.fitted$Age)
actual.fitted$abserror <- abs(actual.fitted$trainerrors)
prom.train_error <- mean(actual.fitted$abserror)

#predicted values for test set
predicted <- inv_ages(predict(age.cv, test.prom, s=best_lambda))

#actual by predicted values for training data
actual.predicted <- data.frame(dat.clean[index,], predicted)

##get the error associated with each predicted age
prom.predicted_age <- actual.predicted$predicted
prom.testerror <- (actual.predicted$predicted - actual.predicted$Age)
prom.abs_error <- abs(prom.testerror)


prom.AAE <- mean(prom.abs_error)
prom.AAE

##PULLING OUT CLOCK SITES
coeffs <- coef(age.cv, s = best_lambda)
coeffs.df <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x)

coeffs.format <- data.frame(stringr::str_split_fixed(coeffs.df$name, ":", 2))[-1,]

options(scipen = 99)

coeffs.bed <- data.frame(coeffs.format$X1, as.numeric(coeffs.format$X2)-1, coeffs.format$X2, coeffs.df$coefficient[-1])
coeffs.bed <- rbind(c("Intercept", NA, NA, coeffs.df$coefficient[1]), coeffs.bed)

#### GO analysis
#library(gprofiler2)
#
#geneBackground <- gconvert(rownames(proms.imp), organism = "olatipes", target = "ENSG")
#
#clockgenes.ENSG <- gconvert(coeffs.df$name[-1], organism = "olatipes",
#                            target = "ENSG")
#
#clockgenes.go <- gost(clockgenes.ENSG$target, organism = "olatipes", 
#                      custom_bg = geneBackground$target)$result
#
#clockgenes.go
#

##linear model for R2 and pvalue

trainlm <- lm(actual.fitted$fitted ~ actual.fitted$Age, data = actual.fitted)
summary(trainlm)

testlm <- lm(actual.predicted$predicted ~ actual.predicted$Age, data = actual.predicted)
summary(testlm)

##prom Test

plot.prom.test <- ggplot(data=actual.predicted, aes(x=Age, y=predicted)) +
  #geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="blue") + 
  geom_point() + ylim(1, 29) +
  theme(panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor = element_line(colour = "grey90"), axis.line = element_line(colour = "black")) + 
  labs(x = "Actual Age", y = "Predicted Promoter Age") +  
  annotate("text", x = 8, y = 25, label = paste0("MAE: ", round(prom.AAE, 2), "mo"),
           size = 8, size.unit = "pt")

plot.prom.test

grid.newpage()
grid.draw(set_panel_size(plot.prom.test, width  = unit(1.5, "in"), height = unit(1.5, "in")))

ggsave("./Figures/1F.svg",
       plot = plot.prom.test, device = "svg", width = 1.5, height = 1.5, units = "in")


### Plot out coeffs

ggplot(coeffs.df[-1,], aes(x = reorder(name, coefficient), y = coefficient)) + 
  geom_col(aes(fill = as.factor(sign(coefficient))), color = "black") + 
  scale_fill_manual(values = c("navy", "firebrick")) + 
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4)) + 
  #scale_y_continuous(limits = c(-0.0045, 0.012), breaks = seq(-0.004, 0.012, 0.001)) +
  labs(x = "Predictors (n = 56)", y = "Beta coefficient")


################################################################ LOOCV


dataset<-as.matrix(1:length(dat.clean$ID))
##create lists to input data values into
R2trains <- c()
Train_MAE <- c()
Test_Error <- c()
sites <- c()
Predicted_age <- c()
sample_names <- colnames(proms.imp)


### Use 10-fold cross validation to estimate lambda

set.seed(1998) #setting seed
age.cv <- cv.glmnet(t(proms.imp), 
                    dat.clean$AgeTrans,
                    alpha=0.5, nfolds=10, family="gaussian",
                    trace.it = TRUE)

plot(age.cv)
#define the lambda parameter
best_lambda <- age.cv$lambda.min

##elastic net for loop -- modified from SG script
for (i in 1:nrow(dataset)){
  print(paste0("Running model ", i, " of ", nrow(dataset)))
  train.prom <- t(proms.imp)[-i,] ##train is all except the 1 test sample per iteration
  test.prom <- t(proms.imp)[i,] ##test is left out sample
  test.age <- dat.clean$AgeTrans[c(i)]
  train.age <- dat.clean$AgeTrans[c(-i)]
  
  ##run elastic net
  model <- glmnet(as.matrix(train.prom), 
                  as.matrix(train.age),
                  alpha=0.5, family="gaussian", nlambda=100)
  #coefficients for fitted model
  age_coeff <- coef(model, s=best_lambda)
  #fitted values
  fitted <- inv_ages(predict(model, as.matrix(train.prom), s=best_lambda))
  actual.fitted <- data.frame(fitted, train.age)
  #calculate R2 for fitted data
  R2.train <- cor(fitted, train.age)^2
  R2trains[i] <- R2.train ##write R2 from each iteration to list
  #predicted values for test set
  predicted <- inv_ages(predict(model, t(test.prom), s=best_lambda))
  #actual by predicted values for training data
  actual.predicted <- data.frame(predicted, test.age)
  ##get the error associated with each predicted age
  Predicted_age[i] <- actual.predicted$predicted
  age.coeffs <- coef(model, s =best_lambda)
  sites <- append(sites,(name = age.coeffs@Dimnames[[1]][age.coeffs@i + 1]))
}

##merge output into one dataset
merged_info <- data.frame(dat.clean, R2trains, Predicted_age)

merged_info$Test_Error <- merged_info$Predicted_age - merged_info$Age
mean(abs(merged_info$Test_Error))

write.csv(merged_info, "./Clocks/PromClock_LOOCV_res.csv")

loocv.test.lm <- lm(Predicted_age ~ Age + Tot_wgt + Gonad_sex, data = merged_info)
summary(loocv.test.lm)

## Plot Test Data
f <- ggplot(merged_info, aes(x = Age, y = Predicted_age)) + 
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() + labs(x = "Age (months)", y = "Predicted Age") +
  #geom_smooth(color = "black", method = "lm", alpha = 0.2) + 
  ylim(1, 29) + 
  annotate(geom = "text", x = 9, y = 26, label = "MAE: 4.32mo", size = 9, size.unit = "pt")
f

grid.newpage()
grid.draw(set_panel_size(f, width  = unit(1.5, "in"), height = unit(1.5, "in")))


ggsave("./Figures/1E.svg",
       plot = f, device = "svg", width = 1.5, height = 1.5, units = "in")


################### Entropy-based clocks


entropy <- function(x, delta = 0.0001){
  z <- x / 100
  (-z * log2(z+delta)) - ( (1 - z)*log2(1-z + delta) )
}


cpg.entropy <- entropy(CpG.imputed[,-1])

set.seed(1998)

index <- order(dat.seq$Age)[seq(1, length(dat.seq$Age), 5)] # Somewhat random sampling along range of ages

#dividing into training and test sets
test.entropy <- t(cpg.entropy)[index,]
test.age <- dat.seq$AgeTrans[index]
hist(test.age)

train.entropy <- t(cpg.entropy)[-index,] ## take all that are not in the test set
train.age <- dat.seq$AgeTrans[-index]
hist(train.age)

gc()

#using 10 fold cross validation to estimate the lambda parameter 
age.cv <- cv.glmnet(train.entropy, 
                    as.matrix(train.age),
                    alpha=0.5, nfolds= 10, family="gaussian", 
                    trace.it = TRUE)

plot(age.cv)

#define the lambda parameter
best_lambda <- age.cv$lambda.min

#coefficients for fitted model
age_coeff <- coef(age.cv, s=best_lambda)

#fitted values
fitted <- inv_ages(predict(age.cv, train.entropy, s=best_lambda))
actual.fitted <- data.frame(dat.seq[-index,], fitted)

#calculate R2 for fitted data
R2.train <- cor(fitted, train.age)^2

##get the error associated with each predicted age
actual.fitted$trainerrors <- (actual.fitted$fitted - actual.fitted$Age)
actual.fitted$abserror <- abs(actual.fitted$trainerrors)
entropy.train_error <- mean(actual.fitted$abserror)

#predicted values for test set
predicted <- inv_ages(predict(age.cv, test.entropy, s=best_lambda))

#actual by predicted values for training data
actual.predicted <- data.frame(dat.seq[index,], predicted)

##get the error associated with each predicted age
entropy.predicted_age <- actual.predicted$predicted
entropy.testerror <- (actual.predicted$predicted - actual.predicted$Age)
entropy.abs_error <- abs(entropy.testerror)


entropy.AAE <- mean(entropy.abs_error)
entropy.AAE

##PULLING OUT CLOCK SITES
coeffs <- coef(age.cv, s = best_lambda)
coeffs.df <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x)

coeffs.format <- data.frame(stringr::str_split_fixed(coeffs.df$name, ":", 2))[-1,]

options(scipen = 99)

coeffs.bed <- data.frame(coeffs.format$X1, as.numeric(coeffs.format$X2)-1, coeffs.format$X2, coeffs.df$coefficient[-1])
coeffs.bed <- rbind(c("Intercept", NA, NA, coeffs.df$coefficient[1]), coeffs.bed)


options(scipen = 3)

#### PLOT CLOCKs
library(ggpubr)

testlm <- lm(actual.predicted$predicted ~ actual.predicted$Age, data = actual.predicted)
summary(testlm)

##entropy Test

plot.entropy.test <- ggplot(data=actual.predicted, aes(x=Age, y=predicted)) +
  geom_point() + 
  #geom_smooth(method=lm, na.rm = TRUE, fullrange= TRUE, aes(group=1),colour="black") + 
  labs(x = "Actual Age", y = "Predicted Epigenetic Age") +
  ylim(2, 29) + 
  annotate("text", x = 8, y = 25, label = paste0("MAE: ", round(entropy.AAE, 2), " months"),
           size = 7, size.unit = "pt")

plot.entropy.test


grid.newpage()
grid.draw(set_panel_size(plot.entropy.test, width  = unit(1.5, "in"), height = unit(1.5, "in")))

ggsave("./Figures/1H.svg",
       plot = plot.entropy.test, device = "svg", width = 1.5, height = 1.5, units = "in")



############# Entropy LOOCV
age.cv <- cv.glmnet(t(cpg.entropy), 
                    dat.seq$AgeTrans,
                    alpha=0.5, nfolds= nrow(dat.seq), family="gaussian", 
                    trace.it = TRUE, keep = TRUE)

plot(age.cv)


loocv.res <- matrix(inv_ages(age.cv$fit.preval), nrow = ncol(cpg.entropy))

entropy.lambda.min <- which(age.cv$lambda == age.cv$lambda.min)


entropyClock.test <- data.frame(dat.seq, "EntropyAge" = loocv.res[,entropy.lambda.min])

write.csv(entropyClock.test, "./Clocks/EntropyClock_LOOCV_res.csv")

g <- ggplot(entropyClock.test, aes(x = Age, y = EntropyAge)) + 
  #geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_point() + labs(x = "Age (months)", y = "Predicted Age (entropy LOOCV)") +
  #geom_smooth(color = "black", method = "lm", alpha = 0.2) + 
  ylim(range(dat.seq$Age)) + 
  annotate(geom = "text", x = 9, y = 26, label = "MAE: 4.82mo", size = 9, size.unit = "pt")
g

grid.newpage()
grid.draw(set_panel_size(g, width  = unit(1.5, "in"), height = unit(1.5, "in")))


ggsave("./Figures/1G.svg",
       plot = g, device = "svg", width = 1.5, height = 1.5, units = "in")

# Mean absolute error

mean(abs(entropyClock.test$EntropyAge - entropyClock.test$Age), na.rm = TRUE) ## 4.82
