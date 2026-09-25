library(tidyverse)
library(ggpubr)
library(data.table)
library(splines)
library(grid)
library(egg)

setwd(".")


theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=6),
                  axis.title=element_text(size=9,face="bold"),
                  strip.text = element_text(size=9,face="bold"),
                  panel.border = element_rect(fill = NA)))



DNAm_res <- read.csv("./Clocks/AgeTrans_LOOCV_res.csv", 
                      header = TRUE)


DNAm.test.lm <- lm(DNAmAge ~ Age + Gonad_sex, data = DNAm_res)
summary(DNAm.test.lm)

DNAm_res$resids <- DNAm.test.lm$residuals



entropy_res <- read.csv("./Clocks/EntropyClock_LOOCV_res.csv", 
                      header = TRUE)


entropy.test.lm <- lm(EntropyAge ~ Age + Gonad_sex, data = entropy_res)
summary(entropy.test.lm)

entropy_res$resids <- entropy.test.lm$residuals


promoter_res <- read.csv("./Clocks/PromClock_LOOCV_res.csv", 
                      header = TRUE)[-96,]

prom.test.lm <- lm(Predicted_age ~ Age + Gonad_sex, data = promoter_res)
summary(prom.test.lm)

promoter_res$resids <- prom.test.lm$residuals

promoter_res <- promoter_res[match(DNAm_res$ID, promoter_res$ID),]



gex_res <- readxl::read_excel("./Clocks/LOOCV_Performance.xlsx")

gex.test.lm <- lm(Predicted_age ~ age + sex, data = gex_res)
summary(gex.test.lm)

gex_res$resids <- gex.test.lm$residuals

gex_res$ID <- str_remove_all(gex_res$ID, "\\.bam") %>% str_replace_all("ID", "X")

gex_res <- gex_res[match(DNAm_res$ID, gex_res$ID),]


total <- data.frame(DNAm_res[,-c(16:21)], 
                    "CpGAge" = DNAm_res$DNAmAge, "CpGResids" = DNAm_res$resids,
                    "EntropyAge" = entropy_res$EntropyAge, "EntropyResids" = entropy_res$resids,
                    "PromoterAge" = promoter_res$Predicted_age, "PromoterResids" = promoter_res$resids,
                    "GEXAge" = gex_res$Predicted_age, "GEXResids" = gex_res$resids)

gex_dnam_resids <- ggplot(total, aes(x = CpGResids, y = GEXResids, color = Age)) + 
  geom_smooth(aes(group = 1), method = "lm", alpha = 0.2, color = "black") + 
  geom_point() + scale_color_viridis_b() + ylim(-15, 15)

grid.newpage()
grid.draw(set_panel_size(gex_dnam_resids, 
                         width  = unit(3, "in"), 
                         height = unit(3, "in")))


ggsave("./Figures/2C.svg",
       plot = gex_dnam_resids, device = "svg", width = 3, height = 2.5, units = "in")


prom_dnam_resids <- ggplot(total, aes(x = CpGResids, y = PromoterResids, color = Age)) + 
  geom_smooth(aes(group = 1), method = "lm", alpha = 0.2, color = "black") + 
  geom_point() + scale_color_viridis_b()#+ ylim(-15, 15)

grid.newpage()
grid.draw(set_panel_size(prom_dnam_resids, 
                         width  = unit(3, "in"), 
                         height = unit(3, "in")))


ggsave("./Figures/2B.svg",
       plot = prom_dnam_resids, device = "svg", width = 3, height = 2.5, units = "in")

entropy_dnam_resids <- ggplot(total, aes(x = CpGResids, y = EntropyResids)) + 
  geom_point() + geom_smooth(method = "lm") + ylim(-7, 7)

grid.newpage()
grid.draw(set_panel_size(entropy_dnam_resids, 
                         width  = unit(3, "in"), 
                         height = unit(3, "in")))


library(ComplexHeatmap)

#col_fun = circlize::colorRamp2(c(0, 0.4, 1), c("red", "white", "black"))
#
#Heatmap(cor(total[,c(4, 17, 19, 21, 23)], method = "spearman", use = "pairwise"), 
#        name = "Spearman Cor \nSignificance", col = col_fun, rect_gp = gpar(col = "black"),
#        cell_fun = function(j, i, x, y, width, height, fill) {
#          grid.text(matrix(round(psych::corr.test(total[,c(4, 17, 19, 21, 23)], 
#                                 method = "spearman", use = "pairwise", adjust = "none")$p, 3),
#                           nrow = 5)[i, j], x, y, gp = gpar(fontsize = 10))})

col_fun = circlize::colorRamp2(c(-.4, 0, .4), c("blue", "white", "red"))

cor <- psych::corr.test(total[,c(5, 17, 19, 21, 23)], 
                        method = "spearman", use = "pairwise")

cor_mat <- paste(round(cor$r, 2), round(cor$p, 2), sep = "\n")
cor_mat

Heatmap(cor(total[,c(5, 17, 19, 21, 23)], method = "spearman", use = "pairwise"), 
        name = "Spearman \nCor Coeff", col = col_fun, rect_gp = gpar(col = "black"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(matrix(cor_mat,
                           nrow = 5)[i, j], x, y, gp = gpar(fontsize = 10))})


heat.1 <- grid.grab()

## draw it, changes optional
heat.1.new <- editGrob(heat.1, vp = viewport(width = unit(6, "in"), 
                                             height = unit(4.5, "in")))
grid.newpage()
grid.draw(heat.1.new)

ggsave("./Figures/2A.svg",
       plot = heat.1.new, device = "svg", width = 6, height = 4.5, units = "in")




## Calculate CI for correlation relationship using Fisher transformation

rho <- cor.test(total$GEXResids, total$PromoterResids, method = "spearman")$estimate
n <- length(total$GEXResids)
delta <- 1.96 / sqrt(n - 3)

#Lower
tanh(atanh(rho) - delta) # -0.01
tanh(atanh(rho) + delta) # 0.41

## More conservative estimate

spearman_CI <- function(x, y, alpha = 0.05){
  rs <- cor(x, y, method = "spearman", use = "complete.obs")
  n <- sum(complete.cases(x, y))
  sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
}


spearman_CI(total$GEXResids, total$PromoterResids)
