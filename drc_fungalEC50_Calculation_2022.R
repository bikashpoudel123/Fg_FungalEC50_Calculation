#http://www.darrenkoppel.com/2020/09/04/dose-response-modelling-and-model-selection-in-r/
#Risolates

library(tidyverse)
library(drc)
TEB_EC50 <- read.csv(file = "FgEC50eval_TEB_drc_R.csv", sep = ",", header = TRUE)
str(TEB_EC50)
TEB_EC50$Trial <- as.factor(TEB_EC50$Trial)
TEB_EC50$Isolate <- as.factor(TEB_EC50$Isolate)
TEB_EC50$Rep <- as.factor(TEB_EC50$Rep)
TEB_EC50$Fungicide <- as.factor(TEB_EC50$Fungicide)
#**********************************************************************************************#
#General format for EC50 evaluation: model picks the best fit
model<- drm(rootl~conc, data=ryegrass, fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
#you don't need the 'names = ' argument but it's useful to label the b, c, d, and e parameters until you're familiar with
plot(model, type="all")
#model f
summary(modelPH1)
#Calculate EC50
ED(model, c(10,20,50), interval="delta")
#**********************************************************************************************#

#start by converting all the responses into a percent of the control response
toxdata <- toxdata %>% 
  mutate(percent_response = rootl/(mean(toxdata$rootl[toxdata$conc==0]))*100)
TEB_EC50_Percentresponse <- TEB_EC50 %>%
  mutate(percent_response = Growth/(mean(TEB_EC50$Growth[TEB_EC50$Conc==0]))*100)
R1707.percent <- subset(TEB_EC50_Percentresponse, (Isolate=="R1707"))
PH1SEP27.percent <- subset(TEB_EC50_Percentresponse, (Isolate=="PH-1_SEP27"))


#R801_1.percent <- subset(TEB_EC50_Percentresponse, (Isolate=="R801-1" & Trial=="3"))
#PH1SEP26.percent <- subset(TEB_EC50_Percentresponse, (Isolate=="PH-1_SEP26" & Trial=="3"))

#check that we now have a response out of 100
head(R1261.percent)

#*****************************************************************************************#
#when fixing model parameters, we use "NA" to indicate that it needs to be calculated 
model_fixed<- drm(percent_response~conc, data=toxdata, 
                  fct=LL.4(fixed=c(NA, 0, 100, NA),
                           names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
#fit multiple charts in one page
par(mfrow = c(4, 2)) 
plot(model_fixed, main="LL.4(fixed=c(NA, 0, 100, NA))")
#******************************************************************************************#

#Choosing a model
#The drc package is kind enouh to have a function made to compare models, 
#‘mselect’ which is a part of the drc package. The key inputs include a model 
#(so run any one first), a list of models you want to compare the initial model to “fctList”, 
#and whether you want to chuck in a linear regression “linreg” (as well as cubic and quadratic). 
#Because we’re using percent response as our response variable, I will fix three parameter 
#models to start at 100.
#example dataset
model.LL3<- drm(percent_response~conc, data=toxdata, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(model.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 

#Dataset
modelR1707.1.LL3<- drm(percent_response~Conc, data=R1707.percent, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(modelR1707.1.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 

modelPH1SEP27.1.LL3<- drm(percent_response~Conc, data=PH1SEP27.percent, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(modelPH1SEP27.1.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 

#************************************************************************************#
#for example data, W2.3 was the best fit model
model.W23 <-  drm(percent_response~conc, data=toxdata, fct=W2.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
plot(model.W23, add=TRUE,col="orange",lty=1, lwd=2)
#EC50 calculation with 95% CI
ED(model.W23, 50, interval = "delta")
#***********************************************************************************#

#for PH1 data, W1.3 was the best model  (best fit)
modelPH1SEP27.1.W14 <-  drm(percent_response~Conc, data=PH1_SEP27.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope","Lower Limit", "Upper Limit",  "ED50")))
modelR804_1.1.W14 <-  drm(percent_response~Conc, data=R804_1.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1171.1.W14 <-  drm(percent_response~Conc, data=R1171.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1213.1.W14 <-  drm(percent_response~Conc, data=R1214.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1214.1.W14 <-  drm(percent_response~Conc, data=R1214.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope","Lower Limit", "Upper Limit",  "ED50")))
modelR1217_1.1.W14 <-  drm(percent_response~Conc, data=R1217_1.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1218_1.1.W14 <-  drm(percent_response~Conc, data=R1218_1.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope","Lower Limit", "Upper Limit",  "ED50")))
modelR1226.1.W14 <-  drm(percent_response~Conc, data=R1226.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope","Lower Limit", "Upper Limit",  "ED50")))
modelR1231.1.W14 <-  drm(percent_response~Conc, data=R1231.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope","Lower Limit", "Upper Limit",  "ED50")))
modelR1236.1.W14 <-  drm(percent_response~Conc, data=R1236.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1237.1.W14 <-  drm(percent_response~Conc, data=R1237.percent, fct=W1.4(fixed=c(NA,0,100, NA), names = c("Slope", "Lower Limit","Upper Limit",  "ED50")))
modelR1238.1.W14 <-  drm(percent_response~Conc, data=R1238.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1240_1.1.W14 <-  drm(percent_response~Conc, data=R1240_1.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1246.1.W14 <-  drm(percent_response~Conc, data=R1246.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1247.1.W13 <-  drm(percent_response~Conc, data=R1247.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelR1250.1.W14<-  drm(percent_response~Conc, data=R1250.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1261.1.W13 <-  drm(percent_response~Conc, data=R1261.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelR1262_1.1.LL4 <-  drm(percent_response~Conc, data=R1262_1.percent, fct=LL.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR1265_1.1.W14 <-  drm(percent_response~Conc, data=R1265_1.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelR1266_1.1.W14 <-  drm(percent_response~Conc, data=R1266_1.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelPH1SEP26.1.W14 <-  drm(percent_response~Conc, data=PH1SEP26.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope","Lower Limit","Upper Limit",  "ED50")))
modelR366_1.1.W14 <-  drm(percent_response~Conc, data=R366_1.percent, fct=W1.4(fixed=c(NA,0, 100, NA), names = c("Slope","Lower Limit","Upper Limit",  "ED50")))
modelR370.1.W13 <-  drm(percent_response~Conc, data=R370.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelR372.1.W13 <-  drm(percent_response~Conc, data=R372.percent, fct=W1.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
modelR397.1.W14 <-  drm(percent_response~Conc, data=R397.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR407_1.1.W14 <-  drm(percent_response~Conc, data=R407_1.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR409_1.1.W14 <-  drm(percent_response~Conc, data=R409_1.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR458.1.W14 <-  drm(percent_response~Conc, data=R458.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR684.1.W14 <-  drm(percent_response~Conc, data=R684.percent, fct=W1.4(fixed=c(NA, 0, 100, NA), names = c("Slope", "Lower Limit", "Upper Limit",  "ED50")))
modelR801_1.1.W14 <-  drm(percent_response~Conc, data=R801_1.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelPH1OCT11.1.W14 <-  drm(percent_response~Conc, data=PH1OCT11.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))
modelR1707.1.W14 <-  drm(percent_response~Conc, data=R1707.percent, fct=W1.4(fixed=c(NA, 0,100, NA), names = c("Slope", "Lower Limit" ,"Upper Limit",  "ED50")))

#**********************DOSE_RESPONSE_CURVE_FOR_R_ISOLATES**********************************#


#All R isolates into one dose-response curve
tiff("DRC_RIsolates_TEB.tiff", units = "in", height = 9, width = 7, res = 300)
plot(modelPH1SEP27.1.W14, main="Dose Response Curve: Tebuconazole: R-Isolates",xlab="Fungicide Concentration (mg/L)", ylab="Percent Growth Inhibition",  ylim = c(0,100), type='all', pch = 16,lty=1, lwd=2) +scale_x_continuous(levels=c("0","0.5","1","2.5","5","10","25","50","100"))
plot(modelR804_1.1.W14, add = TRUE, col = "orange",pch = 2, lty=2, lwd=2)
plot(modelR1171.1.W14, add = TRUE, col = "blue",pch = 3, lty=3, lwd=2)
plot(modelR1213.1.W14, add = TRUE, col = "purple",pch = 4, lty = 4, lwd=2)
plot(modelR1214.1.W14, add = TRUE, col = "magenta",pch = 5, lty = 5, lwd=2)
plot(modelR1217_1.1.W14, add = TRUE, col = "aquamarine",pch = 6, lty = 6, lwd=2)
plot(modelR1218_1.1.W14, add = TRUE, col = "azure3",pch = 7, lty = 7, lwd=2)
plot(modelR1226.1.W14, add = TRUE, col = "brown",pch = 8, lty = 8, lwd=2)
plot(modelR1231.1.W14, add = TRUE, col = "burlywood",pch = 9, lty = 9, lwd=2)
plot(modelR1236.1.W14, add = TRUE, col = "chartreuse",pch = 10, lty = 10, lwd=2)
plot(modelR1237.1.W14, add = TRUE, col = "cadetblue",pch = 11, lty = 11, lwd=2)
plot(modelR1238.1.W14, add = TRUE, col = "coral",pch = 12, lty = 12, lwd=2)
plot(modelR1240.1.W14, add = TRUE, col = "chocolate",pch = 13, lty = 13, lwd=2)
plot(modelR1246.1.W14, add = TRUE, col = "antiquewhite",pch = 14, lty = 14, lwd=2)
plot(modelR1247.1.W13, add = TRUE, col = "bisque",pch = 15, lty = 15, lwd=2)
plot(modelR1250.1.W14, add = TRUE, col = "turquoise",pch = 16, lty = 16, lwd=2)
plot(modelR1261.1.W13, add = TRUE, col = "wheat4",pch = 17, lty = 17, lwd=2)
plot(modelR1262_1.1.LL4, add = TRUE, col = "peachpuff3",pch = 18, lty = 18, lwd=2)
plot(modelR1265_1.1.W14, add = TRUE, col = "darkgoldenrod1",pch = 19, lty = 19, lwd=2)
plot(modelR1266_1.1.W14, add = TRUE, col = "mediumorchid2",pch = 20, lty = 20, lwd=2)
plot(modelR366_1.1.W14, add = TRUE, col = "cornflowerblue",pch = 21, lty = 21, lwd=2)
plot(modelR370.1.W13, add = TRUE, col = "slategray2",pch = 22, lty = 22, lwd=2)
plot(modelR372.1.W13, add = TRUE, col = "yellow",pch = 23, lty = 23, lwd=2)
plot(modelR397.1.W14, add = TRUE, col = "green",pch = 24, lty = 24, lwd=2)
plot(modelR407_1.1.W14, add = TRUE, col = "violet",pch = 25, lty = 25, lwd=2)
plot(modelR409_1.1.W14, add = TRUE, col = "red",pch = 1, lty = 1, lwd=2)
plot(modelR458.1.W14, add = TRUE, col = "cyan",pch = 2, lty = 2, lwd=2)
plot(modelR684.1.W14, add = TRUE, col = "aliceblue",pch = 3, lty = 3, lwd=2)
plot(modelR801_1.1.W14, add = TRUE, col = "yellowgreen",pch = 4, lty = 4, lwd=2)
plot(modelR1707.1.W14, add = TRUE, col = "pink",pch = 5, lty = 5, lwd=2)
legend(20,100, legend = c("PH-1","R804-1","R1171","R1213","R1214","R1217-1","R1218-1","R1226","R1231","R1236","R1237","R1238","R1240-1","R1246","R1247","R1250","R1261","R1262-1","R1265-1","R1266-1","R366-1","R370","R372","R397","R407-1","R409","R458","R684","R801","R1707"), col = c("Black","Orange","Blue","Purple","Magenta","Aquamarine","Azure3","Brown","Burlywood","Chartreuse","Cadetblue","Coral","Chocolate","Antiquewhite","Bisque","Turquoise","wheat4","Peachpuff3","Darkgoldenrod1","Mediumorchid2","Cornflowerblue","Slategray2","yellow","Green","Violet","Red","cyan","Aliceblue","Yellowgreen","Pink"), lty = 1:30, cex = 0.65)
dev.off()

#Calculating ED50
ED(modelR1707.1.W14, 50, interval = "delta")
ED(modelPH1SEP27.1.W14, 50, interval = "delta")
