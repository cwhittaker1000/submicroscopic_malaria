#Access relevant packages
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

#Assign whichever dataset being used to data_frame
data_frame <- as.data.frame(Data.4th.February)

#Manipulate data_frame in prep for analysis
data_frame$Old_or_New <- as.factor(data_frame$Old_or_New)
data_frame$Global_Region <- as.factor(data_frame$Global_Region)

#Retaining Zeroes and Changing, per Lucy's suggestions:
#Changes PCR 0s
for (i in 1:341) {
  if (data_frame[i,8] == 0) {
    data_frame[i, 9] = (data_frame[i, 8]+0.1)/(data_frame[i, 7]) 
  }
}

#Changes microscopy 0s
for (i in 1:341) {
  if (data_frame[i,11] == 0) {
    data_frame[i, 12] = (data_frame[i, 11]+0.1)/(data_frame[i, 10]) 
  }
}

#Subsetting for low PCR places
subset_low_PCR <- data_frame[data_frame$PCR_Prev < 0.1,]

#Subsetting for low micro places
subset_low_micro <- data_frame[data_frame$Micro_Prev < 0.1,]

###INVESTIGATING SUBSET LOW PCR###

#Create individual datasets for each global region
South_America <- subset_low_PCR[subset_low_PCR$Global_Region == "South America",]
West_Africa <- subset_low_PCR[subset_low_PCR$Global_Region == "West Africa",]
East_Africa <- subset_low_PCR[subset_low_PCR$Global_Region == "East Africa",]
Asia <- subset_low_PCR[subset_low_PCR$Global_Region == "Asia",]

#Microscopy Onto PCR Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(West_Africa$PCR_Prev ~ 0 + West_Africa$Micro_Prev))
East_Africa_lm_zeroint <- summary(lm(East_Africa$PCR_Prev ~ 0 + East_Africa$Micro_Prev))
Asia_lm_zeroint <- summary(lm(Asia$PCR_Prev ~ 0 + Asia$Micro_Prev))
South_America_lm_zeroint <- summary(lm(South_America$PCR_Prev ~ 0 + South_America$Micro_Prev))

LM <- seq(0,0.2,0.001)

plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 0.1), ylim = c(0,0.12), pch = 20, col = "blue")
lines(South_America_lm_zeroint$coefficients[1]*LM, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 0.1), ylim = c(0,0.12), pch = 20, col = "red")
lines(West_Africa_lm_zeroint$coefficients[1]*LM, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,0.1), ylim = c(0,0.12), pch = 20, col = "green")
lines(East_Africa_lm_zeroint$coefficients[1]*LM, LM, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,0.2), ylim = c(0,0.12), pch = 20, col = "black")
lines(Asia_lm_zeroint$coefficients[1]*LM, LM, col = "black", lwd = 2)

legend("topleft", legend = c("Asia", "West Africa", "East Africa", "South America"), 
       col = c("black", "red", "green", "blue"), lwd = c(2,2,2,2))#Create individual datasets for each global region


#PCR Onto Microscopy Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(West_Africa$Micro_Prev ~ 0 + West_Africa$PCR_Prev))
East_Africa_lm_zeroint <- summary(lm(East_Africa$Micro_Prev ~ 0 + East_Africa$PCR_Prev ))
Asia_lm_zeroint <- summary(lm(Asia$Micro_Prev ~ 0 + Asia$PCR_Prev))
South_America_lm_zeroint <- summary(lm(South_America$Micro_Prev ~ 0 + South_America$PCR_Prev))

PCR <- seq(0,0.2,0.001)

plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 0.1), ylim = c(0,0.12), pch = 20, col = "blue")
lines(PCR, South_America_lm_zeroint$coefficients[1]*PCR, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 0.1), ylim = c(0,0.12), pch = 20, col = "red")
lines(PCR, West_Africa_lm_zeroint$coefficients[1]*PCR, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,0.1), ylim = c(0,0.12), pch = 20, col = "green")
lines(PCR, East_Africa_lm_zeroint$coefficients[1]*PCR, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,0.1), ylim = c(0,0.12), pch = 20, col = "black")
lines(PCR, Asia_lm_zeroint$coefficients[1]*PCR, col = "black", lwd = 2)

legend("topleft", legend = c("Asia", "West Africa", "East Africa", "South America"), 
       col = c("black", "red", "green", "blue"), lwd = c(2,2,2,2))#Create individual datasets for each global region


###INVESTIGATING SUBSET LOW MICROSCOPY###

#Create individual datasets for each global region
South_America <- subset_low_micro[subset_low_micro$Global_Region == "South America",]
West_Africa <- subset_low_micro[subset_low_micro$Global_Region == "West Africa",]
East_Africa <- subset_low_micro[subset_low_micro$Global_Region == "East Africa",]
Asia <- subset_low_micro[subset_low_micro$Global_Region == "Asia",]

#Microscopy Onto PCR Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(West_Africa$PCR_Prev ~ 0 + West_Africa$Micro_Prev))
East_Africa_lm_zeroint <- summary(lm(East_Africa$PCR_Prev ~ 0 + East_Africa$Micro_Prev))
Asia_lm_zeroint <- summary(lm(Asia$PCR_Prev ~ 0 + Asia$Micro_Prev))
South_America_lm_zeroint <- summary(lm(South_America$PCR_Prev ~ 0 + South_America$Micro_Prev))

LM <- seq(0,0.2,0.001)

plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 0.1), ylim = c(0,0.1), pch = 20, col = "blue")
lines(South_America_lm_zeroint$coefficients[1]*LM, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 0.1), ylim = c(0,0.1), pch = 20, col = "red")
lines(West_Africa_lm_zeroint$coefficients[1]*LM, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,0.1), ylim = c(0,0.1), pch = 20, col = "green")
lines(East_Africa_lm_zeroint$coefficients[1]*LM, LM, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,0.2), ylim = c(0,0.1), pch = 20, col = "black")
lines(Asia_lm_zeroint$coefficients[1]*LM, LM, col = "black", lwd = 2)

legend("topleft", legend = c("Asia", "West Africa", "East Africa", "South America"), 
       col = c("black", "red", "green", "blue"), lwd = c(2,2,2,2))#Create individual datasets for each global region


#PCR Onto Microscopy Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(West_Africa$Micro_Prev ~ 0 + West_Africa$PCR_Prev))
East_Africa_lm_zeroint <- summary(lm(East_Africa$Micro_Prev ~ 0 + East_Africa$PCR_Prev ))
Asia_lm_zeroint <- summary(lm(Asia$Micro_Prev ~ 0 + Asia$PCR_Prev))
South_America_lm_zeroint <- summary(lm(South_America$Micro_Prev ~ 0 + South_America$PCR_Prev))

PCR <- seq(0,0.2,0.001)

plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "blue")
lines(PCR, South_America_lm_zeroint$coefficients[1]*PCR, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 0.2), ylim = c(0,0.2), pch = 20, col = "red")
lines(PCR, West_Africa_lm_zeroint$coefficients[1]*PCR, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "green")
lines(PCR, East_Africa_lm_zeroint$coefficients[1]*PCR, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,0.2), ylim = c(0,0.2), pch = 20, col = "black")
lines(PCR, Asia_lm_zeroint$coefficients[1]*PCR, col = "black", lwd = 2)

legend("topleft", legend = c("Asia", "West Africa", "East Africa", "South America"), 
       col = c("black", "red", "green", "blue"), lwd = c(2,2,2,2))#Create individual datasets for each global region

