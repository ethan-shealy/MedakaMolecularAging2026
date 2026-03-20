###################
#LOAD PACKAGES
###################

library(rtracklayer)
library(ggpubr)
library(data.table)
library(glmnet)
library(ggfortify)
library(openxlsx)

###################
#SET WD
setwd("")

###################

theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=12),
                  axis.title=element_text(size=16,face="bold"),
                  strip.text = element_text(size=16,face="bold"),
                  panel.border = element_rect(fill = NA)))


###################
#PART ONE: TRANSCRIPTOME CLOCKS
###################


###################
#CLOCK 1: 80/20 TRAIN/TEST METHOD
###################


#### USE CPM, normalized by library size
load("./Code_Data/Transcriptome_KnownAge_CountData_Dec2024.RData")
genecounts_clock <- genecounts_cpm
genecounts_clock$id <- rownames(genecounts_clock) 
genecounts_clock <- genecounts_clock[, c("id", setdiff(names(genecounts_clock), "id"))]
clock_meta <- genecounts_clock[,1:4]
genecounts_clock_genes <- genecounts_clock[, -c(1:4)]
genecounts_clock_genes <- t(genecounts_clock_genes)
colnames(genecounts_clock_genes) <- NULL
colnames(genecounts_clock_genes) <- clock_meta$id


##create index for randomly selecting ages
set.seed(123) #generate random numbers, but each time it will be the same, FOR REPRODUCIBILITY
medaka_age1 <- genecounts_clock$age
#index1 <- sample(1:length(medaka_age1), 17)  # Random sampling of 17 elements
index1 <- order(medaka_age1)[seq(1, length(medaka_age1), 5)] # Somewhat random sampling along range of ages

print(index1) #CONFIRMED: 31 79 51 14 67 42 50 43 81 25 69 57  9 26  7 83 72

#divide genecounts_clock into training and test sets
test.genes1 <- t(genecounts_clock_genes)[index1,] #matrix of 17 test individuals

test.age <- medaka_age1[index1]
print(test.age) #CONFIRMED: 15 16  3  4 16 19  4 15  9  5 10  9 14  4 16  7  5
hist(test.age)

train.genes1 <- t(genecounts_clock_genes)[-index1,]

train.age <- medaka_age1[-index1] #matrix of 67 training individuals
print(train.age) #[1] 29  5  4  2 28 19 15 29  8  8 25  3 21 28 17 16 14 10 10  8 27  3  2 28 19 14 21 10  7  7  5  4  3 21 19 16 10  9 28  7  5 29 17 17 14 14  9  6  8
#[50]  5  4  2 27 17 17 14  9  8 25 23  2 27 19 17 14  9  7
hist(train.age)

#use 10 fold cross validation to estimate the lambda parameter
age.cv <- cv.glmnet(train.genes1,
                    log(as.matrix(train.age)),
                    alpha=0.5, nfolds=10, family="gaussian",
                    trace.it = TRUE)
plot(age.cv)



#define the lambda parameter- regularization value that yields best model performance
best_lambda <- age.cv$lambda.min #CONFIRMED: 0.0168

#coefficients for fitted model
age_coeff <- coef(age.cv, s=best_lambda) 

#fitted values
fitted <- exp(predict(age.cv, train.genes1, s=best_lambda))
actual.fitted <- data.frame(clock_meta[-index1,], fitted)
plot <- ggplot(data=actual.fitted, aes(x=age, y=fitted))+
  geom_point()
plot


#calculate R2 for fitted data
R2.train <- cor(fitted, train.age)^2 #0.998

##get the error associated with each predicted age
actual.fitted$trainerrors <- (actual.fitted$s1 - actual.fitted$age)
actual.fitted$abserror <- abs(actual.fitted$trainerrors)
gene.train_error <- mean(actual.fitted$abserror)

#predicted values for test set
predicted <- exp(predict(age.cv, (test.genes1), s=best_lambda))

#actual predicted values for training data
actual.predicted <- data.frame(clock_meta[index1,], predicted, test.age)


actual.predicted_train <- data.frame(clock_meta[-index1,], fitted, train.age)


#For Ethan (4/3/25)
test.age.estimates.df.80_20_Model <- actual.predicted
test.age.estimates.df.80_20_Model <- test.age.estimates.df.80_20_Model[, -6] #remove extra column of actual ages
colnames(test.age.estimates.df.80_20_Model)[5] <- "predicted_age"

write.csv(test.age.estimates.df.80_20_Model, "./Clocks/MED_test_age_estimates.csv")

##get the error associated with each predicted age
gene.predicted_age <- actual.predicted$s1
gene.testerror <- (actual.predicted$s1 - actual.predicted$test.age)
gene.testabs_error <- abs(gene.testerror)
gene.MedianAE <- median(gene.testabs_error)
gene.MedianAE #CONFIRMED: Median Absolute Error: 3.39
gene.MeanAbsoluteError <- mean(gene.testabs_error) #CONFIRMED: Mean Absolute Error: 3.35
gene.MeanAbsoluteError


##PULLING OUT THE CLOCK GENES !!!!!!!!!!!!!!
coeffs <- coef(age.cv, s = best_lambda)
coeffs.df <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x)


#make clock genes dataframe
clock_genes_model1 <- coeffs.df


#### PLOT CLOCKS
library(ggpubr)

#Gene Training
plot.gene.train <- ggplot(data=actual.fitted, aes(x=train.age, y=s1)) +
  geom_point(size = 4, alpha=0.8)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange = TRUE,
              aes(group=1),colour="black")+ 
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Actual Age", y = "Predicted Age", title = "Train Set Predictions")

plot.gene.train

#Linear model for R2 and pvalue
trainlm <- lm(actual.fitted$s1 ~ actual.fitted$age, data=actual.fitted)
summary(trainlm)

#Gene Test
plot.gene.test <- ggplot(data=actual.predicted, aes(x=age, y=s1))+
  geom_point(size = 4, alpha=0.8)+
  geom_smooth(method=lm, na.rm = TRUE, fullrange = TRUE,
              aes(group=1),colour="black")+ 
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  labs(x = "Actual Age", y = "Predicted Epigenetic Age", title = "Test Set Predictions (20-80 Test-Train)") +
  annotate("text", x = 5, y = 25, label = "Median Absolute Error - 3.39 mo", 
           size = 5, color = "black", hjust = 0)

plot.gene.test

#Linear model for R2 and pvalue
testlm <- lm(actual.predicted$s1 ~ actual.predicted$age, data=actual.predicted)
summary(testlm) 
#Adjusted R2=0.435, p-value:0.002343, F-stat: 13.36 on 1 and 15 df



###################
#CLOCK 2: Leave One Out METHOD
####################

#Load library
library(dplyr)
library(ggplot2)
library(glmnet)

#Prepare count data
#### USE CPM, normalized by library size
load("./Code_Data/Transcriptome_KnownAge_CountData_Dec2024.RData")
genecounts_clock <- genecounts_cpm
genecounts_clock$id <- rownames(genecounts_clock) 
genecounts_clock <- genecounts_clock[, c("id", setdiff(names(genecounts_clock), "id"))]
genecounts_clock_genes <- genecounts_clock[, -c(1:4)]
genecounts_clock_genes <- t(genecounts_clock_genes)
colnames(genecounts_clock_genes) <- NULL
colnames(genecounts_clock_genes) <- paste0("X", 1:84)

dataset <- as.matrix(1:length(genecounts_clock$id))
##create lists to input data values into
R2trains <- rep(NA, nrow(dataset))
Train_MAE <- rep(NA, nrow(dataset))
Test_Error <- rep(NA, nrow(dataset))
genes <- c()
Predicted_age <- rep(NA, nrow(dataset))
sample_names <- colnames(genecounts_clock$id)

### Use 10fold cross validation to estimate lambda
set.seed(123)
age.cv <- cv.glmnet(t(genecounts_clock_genes),
                    log(as.matrix(genecounts_clock$age)),
                    alpha=0.5, nfolds=10, family="gaussian",
                    trace.it = TRUE)

plot(age.cv)
#define the lambda parameter
best_lambda <- age.cv$lambda.min #0.01

##elastic net for loop
for (i in 1:nrow(dataset)){
  print(paste0("Running model ", i, " of ", nrow(dataset)))
  train.gene <- t(genecounts_clock_genes)[-i,] ##training is all EXCEPT the 1 test sample per iteration
  test.gene <- t(genecounts_clock_genes)[i,] ##test is the one sample left out
  test.age <- genecounts_clock$age[c(i)]
  train.age <- genecounts_clock$age[c(-i)]
  
  ##run elastic net
  model <- glmnet(as.matrix(train.gene),
                  log(as.matrix(train.age)),
                  alpha=0.5, family="gaussian",nlambda = 100)
  #coefficients for fitted model
  age_coeff <- coef(model, s=best_lambda)
  #fitted values
  fitted <- exp(predict(model, as.matrix(train.gene), s=best_lambda))
  actual.fitted <- data.frame(fitted, train.age)
  #calculate R2 for fitted data
  R2.train <- cor(fitted, train.age)^2
  R2trains[i] <- R2.train ##write R2 from each iteration to list
  ##get the error associated with each predicted age
  actual.fitted$trainerrors <- (actual.fitted$s1 - actual.fitted$train.age)
  actual.fitted$abserror <- abs(actual.fitted$trainerrors)
  train_error <- mean(actual.fitted$abserror)
  Train_MAE[i] <- train_error
  #predicted values for test set
  predicted <- exp(predict(model, t(test.gene), s=best_lambda))
  #actual by predicted values for training data
  actual.predicted <- data.frame(predicted, test.age)
  ##get the error associated with each predicted age
  Predicted_age[i] <- actual.predicted$s1
  testerror <- (actual.predicted$s1 - actual.predicted$s1)
  Test_Error[i] <- testerror
  age.coeffs <- coef(model, s=best_lambda)
  genes <- append(genes,(name=age.coeffs@Dimnames[[1]][age.coeffs@i + 1]))
}

##merge output into one dataset
merged_info <- data.frame(genecounts_clock$id, genecounts_clock$age, genecounts_clock$sex, R2trains, Train_MAE, Predicted_age, Test_Error)
View(merged_info)
colnames(merged_info)[colnames(merged_info) == "genecounts_clock.id"] <- "ID"
colnames(merged_info)[colnames(merged_info) == "genecounts_clock.age"] <- "age"
colnames(merged_info)[colnames(merged_info) == "genecounts_clock.sex"] <- "sex"

gene_table <- as.data.frame(table(genes))

#LINEAR MODEL - LOOCV
loocv.test.lm <- lm(age ~ log(Predicted_age), data=merged_info)
summary(loocv.test.lm) 
#Adjusted R-squared: 0.6321
#F-stat: 143.6 on 1 and 82 df
#p-value: 2.2e-16

loocv_transcriptome_residuals <- residuals(loocv.test.lm)
loocv_transcriptome_residuals <- as.data.frame(loocv_transcriptome_residuals)
loocv_transcriptome_residuals$age <- merged_info$age
loocv_transcriptome_residuals$id <- merged_info$id
loocv_transcriptome_residuals$sex <- merged_info$sex

#Calculate predicted age - age
View(merged_info)
merged_info$discordance = merged_info$Predicted_age - merged_info$age
merged_info$absdisc = abs(merged_info$discordance)
MAE <- median(merged_info$absdisc)
MAE

#Save gene table and merged_info
write.xlsx(gene_table, file="./Clocks/ClockGenes_LOOCV_Model.xlsx")
write.xlsx(merged_info, file="./Clocks/LOOCV_Performance.xlsx")
