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
data_frame <- data_frame[order(data_frame$PCR_Percentage),]

#Create individual datasets for each global region
South_America <- data_frame[data_frame$Global_Region == "South America",]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa",]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa",]
Asia <- data_frame[data_frame$Global_Region == "Asia",]

#for trying with all the zeros removed
South_America <- data_frame[data_frame$Global_Region == "South America" & data_frame$Was_Initially_Zero. == "N",]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa" & data_frame$Was_Initially_Zero. == "N",]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa" & data_frame$Was_Initially_Zero. == "N",]
Asia <- data_frame[data_frame$Global_Region == "Asia" & data_frame$Was_Initially_Zero. == "N",]

#Defining LM variable for plotting
LM <- seq(0,1,0.001)

#Fitting Linear Models To Data For All Data Data Points,
#NO TRANSFORMATION
#Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(West_Africa$PCR_Prev ~ 0 + West_Africa$Micro_Prev))
East_Africa_lm_zeroint <- summary(lm(East_Africa$PCR_Prev ~ 0 + East_Africa$Micro_Prev))
Asia_lm_zeroint <- summary(lm(Asia$PCR_Prev ~ 0 + Asia$Micro_Prev))
South_America_lm_zeroint <- summary(lm(South_America$PCR_Prev ~ 0 + South_America$Micro_Prev))

plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(South_America_lm_zeroint$coefficients[1]*LM, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_zeroint$coefficients[1]*LM, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_zeroint$coefficients[1]*LM, LM, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_zeroint$coefficients[1]*LM, LM, col = "black", lwd = 2)

legend("topleft", legend = c("Asia", "West Africa", "East Africa", "South America"), 
       col = c("black", "red", "green", "blue"), lwd = c(2,2,2,2))

#Intercept Not Set to 0
West_Africa_lm_nonzeroint <- summary(lm(West_Africa$PCR_Prev ~ West_Africa$Micro_Prev))
East_Africa_lm_nonzeroint <- summary(lm(East_Africa$PCR_Prev ~ East_Africa$Micro_Prev))
Asia_lm_nonzeroint <- summary(lm(Asia$PCR_Prev ~ Asia$Micro_Prev))
South_America_lm_nonzeroint <- summary(lm(South_America$PCR_Prev ~ South_America$Micro_Prev))

plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(South_America_lm_nonzeroint$coefficients[1,1]+ South_America_lm_nonzeroint$coefficients[2,1]*LM, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_nonzeroint$coefficients[1,1]+ West_Africa_lm_nonzeroint$coefficients[2,1]*LM, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_nonzeroint$coefficients[1,1] + East_Africa_lm_nonzeroint$coefficients[2,1]*LM, LM, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_nonzeroint$coefficients[1,1]+ Asia_lm_nonzeroint$coefficients[2,1]*LM, LM, col = "black", lwd = 2)

#Fitting Linear Models To Data For All Data Data Points
#SQRT TRANSFORMATION
#Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(sqrt(West_Africa$PCR_Prev) ~ 0 + sqrt(West_Africa$Micro_Prev)))
East_Africa_lm_zeroint <- summary(lm(sqrt(East_Africa$PCR_Prev) ~ 0 + sqrt(East_Africa$Micro_Prev)))
Asia_lm_zeroint <- summary(lm(sqrt(Asia$PCR_Prev) ~ 0 + sqrt(Asia$Micro_Prev)))
South_America_lm_zeroint <- summary(lm(sqrt(South_America$PCR_Prev) ~ 0 + sqrt(South_America$Micro_Prev)))

#sqrt scale plotting
plot(sqrt(South_America$PCR_Prev), sqrt(South_America$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(South_America_lm_zeroint$coefficients[1]*sqrt(LM), sqrt(LM), col = "blue", lwd = 2)
points(sqrt(West_Africa$PCR_Prev), sqrt(West_Africa$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_zeroint$coefficients[1]*sqrt(LM), sqrt(LM), col = "red", lwd = 2)
points(sqrt(East_Africa$PCR_Prev), sqrt(East_Africa$Micro_Prev), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_zeroint$coefficients[1]*sqrt(LM), sqrt(LM), col = "green", lwd = 2)
points(sqrt(Asia$PCR_Prev), sqrt(Asia$Micro_Prev), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_zeroint$coefficients[1]*sqrt(LM), sqrt(LM), col = "black", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines((South_America_lm_zeroint$coefficients[1]*sqrt(LM))^2, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines((West_Africa_lm_zeroint$coefficients[1]*sqrt(LM))^2, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines((East_Africa_lm_zeroint$coefficients[1]*sqrt(LM))^2, LM, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines((Asia_lm_zeroint$coefficients[1]*sqrt(LM))^2, LM, col = "black", lwd = 2)

#Intercept Not Set to 0
West_Africa_lm_nonzeroint <- summary(lm(sqrt(West_Africa$PCR_Prev) ~ sqrt(West_Africa$Micro_Prev)))
East_Africa_lm_nonzeroint <- summary(lm(sqrt(East_Africa$PCR_Prev) ~ sqrt(East_Africa$Micro_Prev)))
Asia_lm_nonzeroint <- summary(lm(sqrt(Asia$PCR_Prev) ~ sqrt(Asia$Micro_Prev)))
South_America_lm_nonzeroint <- summary(lm(sqrt(South_America$PCR_Prev) ~ sqrt(South_America$Micro_Prev)))

#sqrt scale plotting
plot(sqrt(South_America$PCR_Prev), sqrt(South_America$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(South_America_lm_nonzeroint$coefficients[1,1] + South_America_lm_nonzeroint$coefficients[2,1]*sqrt(LM), sqrt(LM), col = "blue", lwd = 2)
points(sqrt(West_Africa$PCR_Prev), sqrt(West_Africa$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_nonzeroint$coefficients[1,1] + West_Africa_lm_nonzeroint$coefficients[2,1]*sqrt(LM), sqrt(LM), col = "red", lwd = 2)
points(sqrt(East_Africa$PCR_Prev), sqrt(East_Africa$Micro_Prev), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_nonzeroint$coefficients[1,1] + East_Africa_lm_nonzeroint$coefficients[2,1]*sqrt(LM), sqrt(LM), col = "green", lwd = 2)
points(sqrt(Asia$PCR_Prev), sqrt(Asia$Micro_Prev), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_nonzeroint$coefficients[1,1] + Asia_lm_nonzeroint$coefficients[2,1]*sqrt(LM), sqrt(LM), col = "black", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines((South_America_lm_nonzeroint$coefficients[1,1] + South_America_lm_nonzeroint$coefficients[2,1]*sqrt(LM))^2, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines((West_Africa_lm_nonzeroint$coefficients[1,1] + West_Africa_lm_nonzeroint$coefficients[2,1]*sqrt(LM))^2, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines((East_Africa_lm_nonzeroint$coefficients[1,1] + East_Africa_lm_nonzeroint$coefficients[2,1]*sqrt(LM))^2, LM, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines((Asia_lm_nonzeroint$coefficients[1,1] + Asia_lm_nonzeroint$coefficients[2,1]*sqrt(LM))^2, LM, col = "black", lwd = 2)


#Fitting Linear Models To Data For All Data Data Points
#LOGIT TRANSFORMATION
data_frame <- Data.4th.February

#Removing Zeroes
data_frame <- data_frame[!data_frame$Was_Initially_Zero == "Y",]

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

#Manipulate data_frame in prep for analysis
data_frame$Old_or_New <- as.factor(data_frame$Old_or_New)
data_frame$Global_Region <- as.factor(data_frame$Global_Region)

#Create individual datasets for each global region
South_America <- data_frame[data_frame$Global_Region == "South America",]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa",]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa",]
Asia <- data_frame[data_frame$Global_Region == "Asia",]

#Defining LM variable for plotting
LM <- seq(0,1,0.001)

#Intercept Set to 0
West_Africa_Micro_zeroint <- summary(lm(logit(West_Africa$PCR_Prev) ~ 0 + logit(West_Africa$Micro_Prev)))
East_Africa_Micro_zeroint <- summary(lm(logit(East_Africa$PCR_Prev) ~ 0 + logit(East_Africa$Micro_Prev)))
Asia_Micro_zeroint <- summary(lm(logit(Asia$PCR_Prev) ~ 0 + logit(Asia$Micro_Prev)))
South_America_Micro_zeroint <- summary(lm(logit(South_America$PCR_Prev) ~ 0 + logit(South_America$Micro_Prev)))

#logit scale plotting
plot(logit(South_America$PCR_Prev), logit(South_America$Micro_Prev), xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(South_America_Micro_zeroint$coefficients[1]*logit(LM), logit(LM), col = "blue", lwd = 2)
points(logit(West_Africa$PCR_Prev), logit(West_Africa$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_Micro_zeroint$coefficients[1]*logit(LM), logit(LM), col = "red", lwd = 2)
points(logit(East_Africa$PCR_Prev), logit(East_Africa$Micro_Prev), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_Micro_zeroint$coefficients[1]*logit(LM), logit(LM), col = "green", lwd = 2)
points(logit(Asia$PCR_Prev), logit(Asia$Micro_Prev), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_Micro_zeroint$coefficients[1]*logit(LM), logit(LM), col = "black", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(expit(South_America_Micro_zeroint$coefficients[1]*logit(LM)), expit(logit(LM)), col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(expit(West_Africa_Micro_zeroint$coefficients[1]*logit(LM)), expit(logit(LM)), col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(expit(East_Africa_Micro_zeroint$coefficients[1]*logit(LM)), LM, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(Asia_Micro_zeroint$coefficients[1]*logit(LM)), LM, col = "black", lwd = 2)

#Intercept Not Set to 0
West_Africa_Micro_nonzeroint <- summary(lm(logit(West_Africa$PCR_Prev) ~ logit(West_Africa$Micro_Prev)))
East_Africa_Micro_nonzeroint <- summary(lm(logit(East_Africa$PCR_Prev) ~ logit(East_Africa$Micro_Prev)))
Asia_Micro_nonzeroint <- summary(lm(logit(Asia$PCR_Prev) ~ logit(Asia$Micro_Prev)))
South_America_Micro_nonzeroint <- summary(lm(logit(South_America$PCR_Prev) ~ logit(South_America$Micro_Prev)))

#logit scale plotting
plot(logit(South_America$PCR_Prev), logit(South_America$Micro_Prev), xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(South_America_Micro_nonzeroint$coefficients[1,1] + South_America_Micro_nonzeroint$coefficients[2,1]*logit(LM), logit(LM), col = "blue", lwd = 2)
points(logit(West_Africa$PCR_Prev), logit(West_Africa$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_Micro_nonzeroint$coefficients[1,1] + West_Africa_Micro_nonzeroint$coefficients[2,1]*logit(LM), logit(LM), col = "red", lwd = 2)
points(logit(East_Africa$PCR_Prev), logit(East_Africa$Micro_Prev), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_Micro_nonzeroint$coefficients[1,1] + East_Africa_Micro_nonzeroint$coefficients[2,1]*logit(LM), logit(LM), col = "green", lwd = 2)
points(logit(Asia$PCR_Prev), logit(Asia$Micro_Prev), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_Micro_nonzeroint$coefficients[1,1] + Asia_Micro_nonzeroint$coefficients[2,1]*logit(LM), logit(LM), col = "black", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(expit(South_America_Micro_nonzeroint$coefficients[1,1] + South_America_Micro_nonzeroint$coefficients[2,1]*logit(LM)), LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(expit(West_Africa_Micro_nonzeroint$coefficients[1,1] + West_Africa_Micro_nonzeroint$coefficients[2,1]*logit(LM)), LM, col = "red", lwd = 2)
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(expit(East_Africa_Micro_nonzeroint$coefficients[1,1] + East_Africa_Micro_nonzeroint$coefficients[2,1]*logit(LM)), LM, col = "green", lwd = 2)
points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(Asia_Micro_nonzeroint$coefficients[1,1] + Asia_Micro_nonzeroint$coefficients[2,1]*logit(LM)), LM, col = "black", lwd = 2)



###REVERSING THE REGRESION

#Intercept Not Set to 0
West_Africa_Micro_nonzeroint <- summary(lm(logit(West_Africa$Micro_Prev) ~ logit(West_Africa$PCR_Prev)))
East_Africa_Micro_nonzeroint <- summary(lm(logit(East_Africa$Micro_Prev) ~ logit(East_Africa$PCR_Prev)))
Asia_Micro_nonzeroint <- summary(lm(logit(Asia$Micro_Prev) ~ logit(Asia$PCR_Prev)))
South_America_Micro_nonzeroint <- summary(lm(logit(South_America$Micro_Prev) ~ logit(South_America$PCR_Prev)))

logit_PCR <- logit(seq(0,1,0.001))

#logit scale plotting
plot(logit(South_America$PCR_Prev), logit(South_America$Micro_Prev)) #, xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(logit_PCR, South_America_Micro_nonzeroint$coefficients[1,1] + South_America_Micro_nonzeroint$coefficients[2,1]*logit_PCR, col = "blue", lwd = 2)

plot(logit(West_Africa$PCR_Prev), logit(West_Africa$Micro_Prev)) #, xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(logit_PCR, West_Africa_Micro_nonzeroint$coefficients[1,1] + West_Africa_Micro_nonzeroint$coefficients[2,1]*logit_PCR, col = "blue", lwd = 2)

plot(logit(East_Africa$PCR_Prev), logit(East_Africa$Micro_Prev)) #, xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(logit_PCR, East_Africa_Micro_nonzeroint$coefficients[1,1] + East_Africa_Micro_nonzeroint$coefficients[2,1]*logit_PCR, col = "blue", lwd = 2)

plot(logit(Asia$PCR_Prev), logit(Asia$Micro_Prev)) #, xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(logit_PCR, Asia_Micro_nonzeroint$coefficients[1,1] + Asia_Micro_nonzeroint$coefficients[2,1]*logit_PCR, col = "blue", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0,1), ylim = c(0,1)) #, xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(expit(logit_PCR), expit(South_America_Micro_nonzeroint$coefficients[1,1] + South_America_Micro_nonzeroint$coefficients[2,1]*logit_PCR), col = "blue", lwd = 2)

points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(expit(logit_PCR), expit(West_Africa_Micro_nonzeroint$coefficients[1,1] + West_Africa_Micro_nonzeroint$coefficients[2,1]*logit_PCR), col = "red", lwd = 2)

points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(expit(logit_PCR), expit(East_Africa_Micro_nonzeroint$coefficients[1,1] + East_Africa_Micro_nonzeroint$coefficients[2,1]*logit_PCR), col = "green", lwd = 2)

points(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(logit_PCR), expit(Asia_Micro_nonzeroint$coefficients[1,1] + Asia_Micro_nonzeroint$coefficients[2,1]*logit_PCR), col = "black", lwd = 2)

      