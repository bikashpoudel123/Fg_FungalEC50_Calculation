#http://www.darrenkoppel.com/2020/09/04/dose-response-modelling-and-model-selection-in-r/
#EC50 calculation for 2008 isolates
#**************************Load library and load input dataset ***************************#
library(tidyverse)
library(drc)
TEB_EC50 <- read.csv(file = "FgEC50eval_TEB_drc_R.csv", sep = ",", header = TRUE)
str(TEB_EC50)
TEB_EC50$Trial <- as.factor(TEB_EC50$Trial)
TEB_EC50$Isolate <- as.factor(TEB_EC50$Isolate)
TEB_EC50$Rep <- as.factor(TEB_EC50$Rep)
TEB_EC50$Fungicide <- as.factor(TEB_EC50$Fungicide)
#*******************Step1: Data input*****************************************************#
#start by converting all the responses into a percent of the control response
toxdata <- toxdata %>% 
  mutate(percent_response = rootl/(mean(toxdata$rootl[toxdata$conc==0]))*100)
TEB_EC50_Percentresponse <- TEB_EC50 %>%
  mutate(percent_response = Growth/(mean(TEB_EC50$Growth[TEB_EC50$Conc==0]))*100)
PH1OCT14.percent <- subset(TEB_EC50_Percentresponse, (Isolate=="PH-1_OCT14"))
F8_59.percent <- subset(TEB_EC50_Percentresponse, (Isolate=="F8_59"))

#check that we now have a response out of 100
head(R1261.percent)

#********************************* Step2: Choose a model*********************************#
#Choosing a model
#The drc package is kind enouh to have a function made to compare models, 
#‘mselect’ which is a part of the drc package. The key inputs include a model 
#(so run any one first), a list of models you want to compare the initial model to “fctList”, 
#and whether you want to chuck in a linear regression “linreg” (as well as cubic and quadratic). 
#Because we’re using percent response as our response variable, I will fix three parameter 
#models to start at 100.
modelF8_59.1.LL3<- drm(percent_response~Conc, data=F8_59.percent, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(modelF8_59.1.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 

modelPH1OCT14.1.LL3<- drm(percent_response~Conc, data=PH1OCT14.percent, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(modelPH1OCT14.1.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 

#*******Step3: picking up best model and plotting dose response curve******************#

modelPH1.1.W13 <-  drm(percent_response~Conc, data=PH1.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF8_1.1.W13 <-  drm(percent_response~Conc, data=F8_1.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF8_2.1.W13 <-  drm(percent_response~Conc, data=F8_2.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF8_3.1.W14 <-  drm(percent_response~Conc, data=F8_3.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF8_4.1.W14 <-  drm(percent_response~Conc, data=F8_4.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope","Lower Limit","Upper Limit",  "ED50")))
modelF8_5.1.W14 <-  drm(percent_response~Conc, data=F8_5.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope", "Lower Limit","Upper Limit",  "ED50")))
modelF8_7.1.W14 <-  drm(percent_response~Conc, data=F8_7.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF8_9.1.W13 <-  drm(percent_response~Conc, data=F8_9.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF8_10.1.W13 <-  drm(percent_response~Conc, data=F8_10.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF8_11.1.W14 <-  drm(percent_response~Conc, data=F8_11.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF8_12.1.W13 <-  drm(percent_response~Conc, data=F8_12.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelPH1OCT14.1.LL3 <-  drm(percent_response~Conc, data=PH1OCT14.percent, fct=LL.3(fixed=c(NA,100, NA), names = c("Slope","Upper Limit",  "ED50")))
modelF8_13.1.W14 <-  drm(percent_response~Conc, data=F8_13.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF8_20.1.W14 <-  drm(percent_response~Conc, data=F8_20.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope","Lower Limit","Upper Limit",  "ED50")))
modelF8_21.1.W14 <-  drm(percent_response~Conc, data=F8_21.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF8_28.1.W14 <-  drm(percent_response~Conc, data=F8_28.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF8_29.1.W14 <-  drm(percent_response~Conc, data=F8_29.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF8_30.1.W14 <-  drm(percent_response~Conc, data=F8_30.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF8_31.1.W14 <-  drm(percent_response~Conc, data=F8_31.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope","Lower Limit","Upper Limit",  "ED50")))
modelF8_32.1.W14 <-  drm(percent_response~Conc, data=F8_32.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope","Lower Limit","Upper Limit",  "ED50")))
modelF8_26.1.W14 <-  drm(percent_response~Conc, data=F8_26.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope","Lower Limit","Upper Limit",  "ED50")))
modelF8_33.1.W13 <-  drm(percent_response~Conc, data=F8_33.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF8_34.1.W14 <-  drm(percent_response~Conc, data=F8_34.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF8_35.1.W14 <-  drm(percent_response~Conc, data=F8_35.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope","Lower Limit","Upper Limit", "ED50")))
modelF8_36.1.W13 <-  drm(percent_response~Conc, data=F8_36.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF8_37.1.W14 <-  drm(percent_response~Conc, data=F8_37.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF8_38.1.W14 <-  drm(percent_response~Conc, data=F8_38.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF8_39.1.W14 <-  drm(percent_response~Conc, data=F8_39.percent, fct=W1.4(fixed=c(NA,0,100, NA), names = c("Slope", "Lower Limit","Upper Limit",  "ED50")))
modelF8_44.1.W14 <-  drm(percent_response~Conc, data=F8_44.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF8_47.1.W14 <-  drm(percent_response~Conc, data=F8_47.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF8_59.1.W14 <-  drm(percent_response~Conc, data=F8_59.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))

tiff("DRC_2008Isolates_TEB.tiff", units = "in", height = 9, width = 7, res = 300)
plot(modelPH1.1.W13, main="Dose Response Curve: Tebuconazole: 2008-Isolates",xlab="Fungicide Concentration (mg/L)", ylab="Percent Growth Inhibition",  type='all', pch =1,lty=1, lwd=2)
plot(modelPH1SEP27.1.W14, main="Dose Response Curve: Tebuconazole: 2008-Isolates",xlab="Fungicide Concentration (mg/L)", ylab="Percent Growth Inhibition",  ylim = c(0,100), type='all', pch = 16,lty=1, lwd=2)
plot(modelF8_1.1.W13, add = TRUE, col = "orange", pch =2, lty=2, lwd=2)
plot(modelF8_2.1.W13, add = TRUE, col = "blue", pch =3, lty=3, lwd=2)
plot(modelF8_3.1.W14, add = TRUE, col = "purple", pch =4, lty=4, lwd=2)
plot(modelF8_4.1.W14, add = TRUE, col = "magenta",pch = 5, lty = 5, lwd=2)
plot(modelF8_5.1.W14, add = TRUE, col = "aquamarine",pch = 6, lty = 6, lwd=2)
plot(modelF8_7.1.W14, add = TRUE, col = "azure3", pch =7, lty=7, lwd=2)
plot(modelF8_9.1.W13, add = TRUE, col = "brown", pch =8, lty=8, lwd=2)
plot(modelF8_10.1.W13, add = TRUE, col = "burlywood",pch = 9, lty = 9, lwd=2)
plot(modelF8_11.1.W14, add = TRUE, col = "chartreuse", pch =10, lty=10, lwd=2)
plot(modelF8_12.1.W13, add = TRUE, col = "cadetblue", pch =11, lty=11, lwd=2)
plot(modelF8_13.1.W14, add = TRUE, col = "coral", pch =12, lty=12, lwd=2)
plot(modelF8_20.1.W14, add = TRUE, col = "chocolate",pch = 13, lty = 13, lwd=2)
plot(modelF8_21.1.W14, add = TRUE, col = "antiquewhite", pch =14, lty=14, lwd=2)
plot(modelF8_28.1.W14, add = TRUE, col = "bisque", pch =15, lty=15, lwd=2)
plot(modelF8_29.1.W14, add = TRUE, col = "turquoise", pch =16, lty=16, lwd=2)
plot(modelF8_30.1.W14, add = TRUE, col = "wheat4", pch =17, lty=17, lwd=2)
plot(modelF8_31.1.W14, add = TRUE, col = "peachpuff3",pch = 18, lty = 18, lwd=2)
plot(modelF8_32.1.W14, add = TRUE, col = "darkgoldenrod1", pch =19, lty=19, lwd=2)
plot(modelF8_26.1.W14, add = TRUE, col = "mediumorchid2",pch = 20, lty = 20, lwd=2)
plot(modelF8_33.1.W13, add = TRUE, col = "cornflowerblue", pch =21, lty=21, lwd=2)
plot(modelF8_34.1.W14, add = TRUE, col = "slategray2",pch = 22, lty = 22, lwd=2)
plot(modelF8_35.1.W14, add = TRUE, col = "yellow", pch =23, lty=23, lwd=2)
plot(modelF8_36.1.W13, add = TRUE, col = "green", pch =24, lty=24, lwd=2)
plot(modelF8_37.1.W14, add = TRUE, col = "violet", pch =25, lty=25, lwd=2)
plot(modelF8_38.1.W14, add = TRUE, col = "red",pch = 1, lty = 1, lwd=2)
plot(modelF8_39.1.W14, add = TRUE, col = "cyan", pch =2, lty=2, lwd=2)
plot(modelF8_44.1.W14, add = TRUE, col = "aliceblue", pch =3, lty=3, lwd=2)
plot(modelF8_47.1.W14, add = TRUE, col = "yellowgreen", pch =4, lty=4, lwd=2)
plot(modelF8_59.1.W14, add = TRUE, col = "pink",pch = 5, lty = 5, lwd=2)

legend(20,100, legend = c("PH-1","F8-1","F8-2","F8-3","F8-4","F8-5","F8-7","F8-9","F8-10","F8-11","F8-12","F8-13","F8-20","F8-21","F8-28","F8-29","F8-30","F8-31","F8-32","F8-26","F8-33","F8-34","F8-35","F8-36","F8-37","F8-38","F8-39","F8-44","F8-47","F8-59"), col = c("Black","Orange","Blue","Purple","Magenta","Aquamarine","Azure3","Brown","Burlywood","Chartreuse","Cadetblue","Coral","Chocolate","antiquewhite","Bisque","turquoise","Wheat4","Peachpuff3","Darkgoldenrod1","Mediumorchid2","Cornflowerblue","Slategray2","Yellow","Green","Violet","Red","cyan","Aliceblue","Yellowgreen","Pink"), lty = 1:30, cex = 0.65)
dev.off()

ED(modelF8_59.1.W14, 50, interval = "delta")
ED(modelPH1OCT14.1.W14, 50, interval = "delta")

