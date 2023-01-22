#http://www.darrenkoppel.com/2020/09/04/dose-response-modelling-and-model-selection-in-r/
#EC50 calculation for 2010 isolates
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
F10_4H2_1.percent <- subset(TEB_EC50_Percentresponse, (Isolate=="4H2_1"))
PH1JAN10.percent <- subset(TEB_EC50_Percentresponse, (Isolate=="PH-1_JAN10"))

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
modelF10_4H2_1.1.LL3<- drm(percent_response~Conc, data=F10_4H2_1.percent, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(modelF10_4H2_1.1.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 

modelPH1JAN10.1.LL3<- drm(percent_response~Conc, data=PH1JAN10.percent, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(modelPH1JAN10.1.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 

#*******Step3: picking up best model and plotting dose response curve******************#
modelALLPH1_SEP27.1.W14 <-  drm(percent_response~Conc, data=AllPH1_SEP27.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope", "Lower Limit","Upper Limit",  "ED50")))
modelF10_4H2_1.1.W14 <-  drm(percent_response~Conc, data=F10_4H2_1.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope","Lower Limit","Upper Limit",  "ED50")))
modelF10_55_1.1.W14 <-  drm(percent_response~Conc, data=F10_55_1.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF10_79_1.1.W14 <-  drm(percent_response~Conc, data=F10_79_1.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit","Upper Limit",  "ED50")))
modelF10_120_2.1.W13 <-  drm(percent_response~Conc, data=F10_120_2.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF10_121_03.1.W14 <-  drm(percent_response~Conc, data=F10_121_03.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF10_124_1.1.W14 <-  drm(percent_response~Conc, data=F10_124_1.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF10_127_3.1.W14 <-  drm(percent_response~Conc, data=F10_127_3.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope","Lower Limit", "Upper Limit",  "ED50")))
modelF10_135_2.1.LL3 <-  drm(percent_response~Conc, data=F10_135_2.percent, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelF10_136_04.1.W14 <-  drm(percent_response~Conc, data=F10_136_04.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF10_137_5.1.W14 <-  drm(percent_response~Conc, data=F10_137_5.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF10_140H1_3.1.W14 <-  drm(percent_response~Conc, data=F10_140H1_3.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF10_148H1_3.1.W14 <-  drm(percent_response~Conc, data=F10_148H1_3.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelF10_156H1_1.1.LL4 <-  drm(percent_response~Conc, data=F10_156H1_1.percent, fct=LL.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF10_160H2_2.1.W14 <-  drm(percent_response~Conc, data=F10_160H2_2.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF10_541_1.1.LL4 <-  drm(percent_response~Conc, data=F10_541_1.percent, fct=LL.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelPH1JAN10.1.W13 <-  drm(percent_response~Conc, data=PH1JAN10.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope","Upper Limit",  "ED50")))
modelF10_6H1_2.1.W14 <-  drm(percent_response~Conc, data=F10_6H1_2.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelF10_166H1_3.1.W14 <-  drm(percent_response~Conc, data=F10_166H1_3.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))

tiff("DRC_2010Isolates_TEB.tiff", units = "in", height = 9, width = 7, res = 300)
plot(modelPH1_SEP27.1.W14, main="Dose Response Curve: Tebuconazole: 2010-Isolates",xlab="Fungicide Concentration (mg/L)", ylab="Percent Growth Inhibition",  type='all', pch =1,lty=1, lwd=2)
plot(modelF10_4H2_1.1.W14, add = TRUE, col = "orange",pch = 2, lty = 2, lwd=2)
plot(modelF10_55_1.1.W14, add = TRUE, col = "blue",pch = 3, lty = 3, lwd=2)
plot(modelF10_79_1.1.W14, add = TRUE, col = "purple", pch =4, lty=4, lwd=2)
plot(modelF10_120_2.1.W13, add = TRUE, col = "magenta", pch =5, lty=5, lwd=2)
plot(modelF10_121_03.1.W14, add = TRUE, col = "aquamarine", pch =6, lty=6, lwd=2)
plot(modelF10_124_1.1.W14, add = TRUE, col = "azure3", pch =7, lty=7, lwd=2)
plot(modelF10_127_3.1.W14, add = TRUE, col = "brown", pch =8, lty=8, lwd=2)
#plot(modelF10_135_2.1.W14, add = TRUE, col = "burlywood", pch =9, lty=9, lwd=2)
plot(modelF10_136_04.1.W14, add = TRUE, col = "chartreuse", pch =10, lty=10, lwd=2)
plot(modelF10_137_5.1.W14, add = TRUE, col = "cadetblue", pch =11, lty=11, lwd=2)
plot(modelF10_140H1_3.1.W14, add = TRUE, col = "coral", pch =12, lty=12, lwd=2)
plot(modelF10_148H1_3.1.W14, add = TRUE, col = "chocolate", pch =13, lty=13, lwd=2)
plot(modelF10_156H1_1.1.LL4, add = TRUE, col = "antiquewhite", pch =14, lty=14, lwd=2)
plot(modelF10_160H2_2.1.W14, add = TRUE, col = "bisque", pch =15, lty=15, lwd=2)
#plot(modelF10_541_1.1.W14, add = TRUE, col = "turquoise", pch =16, lty=16, lwd=2)
plot(modelF10_6H1_2.1.W14, add = TRUE, col = "wheat4", pch =17, lty=17, lwd=2)
plot(modelF10_166H1_3.1.W14, add = TRUE, col = "peachpuff3", pch =18, lty=18, lwd=2)

legend(20,100, legend = c("PH-1","F10_4H2_1","F10_55_1","F10_79_1","F10_120_2","F10_121_03","F10_124_1","F10_127_3","Fg10_136_04","F10_137_5","F10_140H1_3","F10_148H1_3","F10_156H1_1","Fg10_160H2_2","F10_6H1_2","F10_166H1_3"), col = c("Black","Orange","Blue","Purple","Magenta","Aquamarine","Azure3","Brown","Chartreuse","Cadetblue","Coral","Chocolate","Antiquewhite","Bisque","Wheat4","Peachpuff3"), lty = 1:18, cex = 0.65)
dev.off()

ED(modelF10_4H2_1.1.W14, 50, interval = "delta")
ED(modelPH1JAN10.1.W13, 50, interval = "delta")
