library(matrixStats)
library(tidyverse)
library(data.table)
library(ggpubr)
library(grid)
library(egg)


setwd("")

theme_set(theme_bw() + 
            theme(legend.title = element_text(size=16,face="bold"),
                  axis.text=element_text(size=12),
                  axis.title=element_text(size=16,face="bold"),
                  strip.text = element_text(size=16,face="bold"),
                  panel.border = element_rect(fill = NA)))

testMeta.seq <- read.csv("./Code_Data/Validation_GTs.csv")

testMeta.seq$GT <- factor(testMeta.seq$GT)

lm.0 <- lm(EA_Est ~ Treat + Age + Final_sex + GT_Num, data = testMeta.seq)
summary(lm.0)

## Get Age residuals

testMeta.resids <- testMeta.seq[-which(is.na(testMeta.seq$EA_Est)),]

testMeta.resids$EA_Resids <- residuals(lm(EA_Est ~ Age, data = testMeta.resids))

e <- ggplot(testMeta.resids, aes(x = GT, y = EA_Resids/30.5, color = GT)) + 
        geom_boxplot() + geom_jitter(width = 0.1) + 
        scale_color_manual(values = c("blue", "purple", "red")) +
        scale_x_discrete(labels = c("TT", "TC", "CC")) +
        labs(x = "GT", y = "EAA") + 
        theme_bw() + 
        theme(legend.title = element_text(size=16,face="bold"),
              axis.text=element_text(size=7),
              axis.title=element_text(size=9,face="bold"),
              strip.text = element_text(size=9,face="bold"),
              panel.border = element_rect(fill = NA), 
              legend.position = "none")


grid.newpage()
grid.draw(set_panel_size(e, width  = unit(1.5, "in"), height = unit(3, "in")))


ggsave("./Figures/6E.svg",
       plot = e, device = "svg", width = 2, height = 3, units = "in")



ggplot(testMeta.resids, aes(x = Age, y = EA_Est, color = GT)) + 
  geom_point() + geom_smooth(method = "lm", alpha = 0.2) +
  scale_color_manual(values = c("blue", "purple", "red"), labels = c("TT", "TC", "CC")) +
  labs(x = "Age (Days)", y = "EA Estimate (Days)", color = "GT")

sjPlot::tab_model(lm(EA_Est ~ Age + Treat + Final_sex + I(GT_Num*30.5), data = testMeta.seq), 
                  pred.labels = c("Intercept", "Age", "TreatHigh", "TreatLow", "TreatMed", "Sex [F]", "Genotype"), 
                  dv.labels = "Epigenetic Age Estimate")
