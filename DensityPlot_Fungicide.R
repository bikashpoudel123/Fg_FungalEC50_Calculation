library(dplyr)
##SplitPlotHistogram
getwd()
setwd("/Users/bikashpoudel/Library/CloudStorage/OneDrive-NorthDakotaUniversitySystem/Dr. Shaobin Zhong/FHB/Fungicide Sensitivity Manuscript/")

####Data Importation#####
#PROTHIOCONAZOLE
PROTHIOCONAZOLE <- read.csv(file = "Fg_Fungicide_Sensitivity_PRO.csv", sep = ",", header = TRUE)
str(PROTHIOCONAZOLE)
PROTHIOCONAZOLE$Isolate <- as.factor(PROTHIOCONAZOLE$Isolate)
PROTHIOCONAZOLE$Rep <- as.factor(PROTHIOCONAZOLE$Rep)
PROTHIOCONAZOLE$FUNGICIDE <- as.factor(PROTHIOCONAZOLE$FUNGICIDE)
PROTHIOCONAZOLE$TRIAL <- as.factor(PROTHIOCONAZOLE$TRIAL)
PROTHIOCONAZOLE$EC50 <- as.numeric(PROTHIOCONAZOLE$EC50)

#TEBUCONAZOLE
TEBUCONAZOLE <- read.csv(file = "TEB_TRIAL1_drc.csv", sep = ",", header = TRUE)
str(TEBUCONAZOLE)
TEBUCONAZOLE$Trial <- as.factor(TEBUCONAZOLE$Trial)
TEBUCONAZOLE$Isolate <- as.factor(TEBUCONAZOLE$Isolate)
TEBUCONAZOLE$Year <- as.factor(TEBUCONAZOLE$Year)
TEBUCONAZOLE$Fungicide <- as.factor(TEBUCONAZOLE$Fungicide)
TEBUCONAZOLE$Model <- as.factor(TEBUCONAZOLE$Model)
TEBUCONAZOLE$Std.error <- as.numeric(TEBUCONAZOLE$Std.error)
TEBUCONAZOLE$EC50.estimate <- as.numeric(TEBUCONAZOLE$EC50.estimate)
TEBUCONAZOLE$AIC <- as.numeric(TEBUCONAZOLE$AIC)
TEBUCONAZOLE$Lower <- as.numeric(TEBUCONAZOLE$Lower)
TEBUCONAZOLE$Upper <- as.numeric(TEBUCONAZOLE$Upper)

#COMBINED
COMBINED <- read.csv(file = "EC50_drc_TEBandPRO_output.csv", sep = ",", header = TRUE)
str(COMBINED)
COMBINED$Trial <- as.factor(COMBINED$Trial)
COMBINED$Isolate <- as.factor(COMBINED$Isolate)
COMBINED$Year <- as.factor(COMBINED$Year)
COMBINED$Fungicide <- as.factor(COMBINED$Fungicide)
COMBINED$Model <- as.factor(COMBINED$Model)
COMBINED$AIC <- as.numeric(COMBINED$AIC)
COMBINED$EC50 <- as.numeric(COMBINED$EC50)

Teb <- data.frame(subset(COMBINED, (Fungicide=="Tebuconazole")))
str(Teb)

Pro <- data.frame(subset(COMBINED, (Fungicide=="Prothioconazole")))
str(Pro)
##Loading Libraries###
install.packages("rcompanion")
library(rcompanion)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(viridis)
library(RColorBrewer)
library(car)

#DensityPlots
tiff(filename = "DensityPlot_TEB_combinedtrials_drc.tiff", width = 6.4, height = 3, units = "in", res = 300)
P <- ggplot(COMBINED, aes(x=EC50, fill=Fungicide)) + 
  geom_histogram(aes(y=..density..), alpha=0.2, position = "identity", bins = 45)+
  geom_density(alpha=0.25)
#with title
#P+theme_classic()+theme(legend.position = "top", plot.title = element_text(hjust = 0.5))+scale_x_continuous(name = "EC50 values", limits = c(0,30), breaks = seq(0,30,5))+
 # scale_y_continuous(name = "Density", limits = c(0,0.5))+labs(title = "Distribution of EC50 among Fg Isolates: Tebuconazole vs Prothioconazole")
#without title
shapiro.test(Teb$EC50)
levene
P+theme_classic()+theme(legend.position = "top", plot.title = element_text(hjust = 0.5))+scale_x_continuous(name = "EC50 values", limits = c(0,20), breaks = seq(0,25,2))+
  scale_y_continuous(name = "Density", limits = c(0,0.7))
dev.off()
p<-ggplot(df, aes(x=weight, fill=sex)) +
  geom_density(alpha=0.4)

ggplot(df, aes(x=weight, color=sex, fill=sex)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2) 
#overlaid histograms
ggplot(PROTHIOCONAZOLE, aes(x=EC50, color=TRIAL)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")

ggplot(PROTHIOCONAZOLE, aes(x=EC50, color=TRIAL)) +
  geom_histogram(fill="white")

#Regression
ggp <- ggplot(COMBINED, aes(Fungicide, EC50)) +           
  geom_point()
ggp
ggp +                                     
  stat_smooth(method = "lm",
              formula = EC50 ~ Fungicide,
              geom = "smooth")

#Histogram with normal curve
tiff("EC50 values distribution among Fg Isolates_Tebuconazole.tiff", height = 5, width = 5, units = "in", res = 300)
plotNormalHistogram(Teb$EC50, col = "yellow", linecol = "green",xlab="EC50 values of Isolates: Tebuconazole", ylab="Number of Isolates", xlim=c(0,12))
legend(5, 1000, legend = c("Shapiro.wilk test", "p-value<2.2e-16"), cex = 1)
dev.off()

tiff("EC50 values distribution among Fg Isolates_Prothioconazole.tiff", height = 5, width = 5, units = "in", res = 300)
plotNormalHistogram(Pro$EC50, col = "yellow", linecol = "green",xlab="EC50 values of Isolates: Prothioconazole", ylab="Number of Isolates", xlim=c(0,30))
legend(5, 1000, legend = c("Shapiro.wilk test", "p-value<2.2e-16"), cex = 1)
dev.off()

#Levene's test for Homogeneity of variances
leveneTest(EC50~ Trial, Pro)

#Facet plot with facet_wrap in ggplot2
tiff("EC50 values distribution among Fg Isolates.tiff", height = 4.06, width = 7, units = "in", res = 300)
COMBINED %>%
  ggplot(aes(x=EC50, fill=Fungicide))+
  geom_histogram(aes(y=..density..), position = "identity", bins = 30)+
  geom_density(alpha=0.5)+
  facet_wrap(~Fungicide, ncol = 1)+
  facet_grid(~factor(Fungicide, levels = c("Tebuconazole","Prothioconazole"))) +
  theme_prism()+scale_x_continuous(name = "EC50 values", limits = c(0,12.5), breaks = seq(0,12.5,2.5))+
  theme(legend.position = "None")+
  labs(x="EC50 values")+
  geom_vline(data = COMBINED, aes(xintercept=mean(EC50), color=Fungicide))
dev.off()
