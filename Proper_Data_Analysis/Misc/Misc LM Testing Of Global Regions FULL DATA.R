#Access relevant packages
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

#Assign whichever dataset being used to data_frame
data_frame <- as.data.frame(Rough_Preliminary_Data_Exploration_For_R)

#Manipulate data_frame in prep for analysis
data_frame$Old_Or_New <- as.factor(data_frame$Old_Or_New)
data_frame$Global_Region <- as.factor(data_frame$Global_Region)
data_frame <- data_frame[order(data_frame$PCR_Percentage),]

#Create individual datasets for each global region
South_America <- data_frame[data_frame$Global_Region == "South America",]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa",]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa",]
Asia <- data_frame[data_frame$Global_Region == "Asia",]

#Defining LM variable for plotting
LM <- seq(0,1,0.001)

#Fitting Linear Models To Data For All Data Data Points,
#NO TRANSFORMATION
#Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(West_Africa$PCR_Percentage ~ 0  West_Africa$LM_Percentage))
East_Africa_lm_zeroint <- summary(lm(East_Africa$PCR_Percentage ~ 0 + East_Africa$LM_Percentage))
Asia_lm_zeroint <- summary(lm(Asia$PCR_Percentage ~ 0 + Asia$LM_Percentage))
South_America_lm_zeroint <- summary(lm(South_America$PCR_Percentage ~ 0 + South_America$LM_Percentage))

plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(South_America_lm_zeroint$coefficients[1]*LM, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_zeroint$coefficients[1]*LM, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_zeroint$coefficients[1]*LM, LM, col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_zeroint$coefficients[1]*LM, LM, col = "black", lwd = 2)

#Intercept Not Set to 0
West_Africa_lm_nonzeroint <- summary(lm(West_Africa$PCR_Percentage ~ West_Africa$LM_Percentage))
East_Africa_lm_nonzeroint <- summary(lm(East_Africa$PCR_Percentage ~ East_Africa$LM_Percentage))
Asia_lm_nonzeroint <- summary(lm(Asia$PCR_Percentage ~ Asia$LM_Percentage))
South_America_lm_nonzeroint <- summary(lm(South_America$PCR_Percentage ~ South_America$LM_Percentage))

plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(South_America_lm_nonzeroint$coefficients[1,1]+ South_America_lm_nonzeroint$coefficients[2,1]*LM, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_nonzeroint$coefficients[1,1]+ West_Africa_lm_nonzeroint$coefficients[2,1]*LM, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_nonzeroint$coefficients[1,1] + East_Africa_lm_nonzeroint$coefficients[2,1]*LM, LM, col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_nonzeroint$coefficients[1,1]+ Asia_lm_nonzeroint$coefficients[2,1]*LM, LM, col = "black", lwd = 2)

#Fitting Linear Models To Data For All Data Data Points
#SQRT TRANSFORMATION
#Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(sqrt(West_Africa$PCR_Percentage) ~ 0 + sqrt(West_Africa$LM_Percentage)))
East_Africa_lm_zeroint <- summary(lm(sqrt(East_Africa$PCR_Percentage) ~ 0 + sqrt(East_Africa$LM_Percentage)))
Asia_lm_zeroint <- summary(lm(sqrt(Asia$PCR_Percentage) ~ 0 + sqrt(Asia$LM_Percentage)))
South_America_lm_zeroint <- summary(lm(sqrt(South_America$PCR_Percentage) ~ 0 + sqrt(South_America$LM_Percentage)))

#sqrt scale plotting
plot(sqrt(South_America$PCR_Percentage), sqrt(South_America$LM_Percentage), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(South_America_lm_zeroint$coefficients[1]*sqrt(LM), sqrt(LM), col = "blue", lwd = 2)
points(sqrt(West_Africa$PCR_Percentage), sqrt(West_Africa$LM_Percentage), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_zeroint$coefficients[1]*sqrt(LM), sqrt(LM), col = "red", lwd = 2)
points(sqrt(East_Africa$PCR_Percentage), sqrt(East_Africa$LM_Percentage), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_zeroint$coefficients[1]*sqrt(LM), sqrt(LM), col = "green", lwd = 2)
points(sqrt(Asia$PCR_Percentage), sqrt(Asia$LM_Percentage), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_zeroint$coefficients[1]*sqrt(LM), sqrt(LM), col = "black", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines((South_America_lm_zeroint$coefficients[1]*sqrt(LM))^2, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines((West_Africa_lm_zeroint$coefficients[1]*sqrt(LM))^2, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines((East_Africa_lm_zeroint$coefficients[1]*sqrt(LM))^2, LM, col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines((Asia_lm_zeroint$coefficients[1]*sqrt(LM))^2, LM, col = "black", lwd = 2)

#Intercept Not Set to 0
West_Africa_lm_nonzeroint <- summary(lm(sqrt(West_Africa$PCR_Percentage) ~ sqrt(West_Africa$LM_Percentage)))
East_Africa_lm_nonzeroint <- summary(lm(sqrt(East_Africa$PCR_Percentage) ~ sqrt(East_Africa$LM_Percentage)))
Asia_lm_nonzeroint <- summary(lm(sqrt(Asia$PCR_Percentage) ~ sqrt(Asia$LM_Percentage)))
South_America_lm_nonzeroint <- summary(lm(sqrt(South_America$PCR_Percentage) ~ sqrt(South_America$LM_Percentage)))

#sqrt scale plotting
plot(sqrt(South_America$PCR_Percentage), sqrt(South_America$LM_Percentage), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(South_America_lm_nonzeroint$coefficients[1,1] + South_America_lm_nonzeroint$coefficients[2,1]*sqrt(LM), sqrt(LM), col = "blue", lwd = 2)
points(sqrt(West_Africa$PCR_Percentage), sqrt(West_Africa$LM_Percentage), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_nonzeroint$coefficients[1,1] + West_Africa_lm_nonzeroint$coefficients[2,1]*sqrt(LM), sqrt(LM), col = "red", lwd = 2)
points(sqrt(East_Africa$PCR_Percentage), sqrt(East_Africa$LM_Percentage), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_nonzeroint$coefficients[1,1] + East_Africa_lm_nonzeroint$coefficients[2,1]*sqrt(LM), sqrt(LM), col = "green", lwd = 2)
points(sqrt(Asia$PCR_Percentage), sqrt(Asia$LM_Percentage), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_nonzeroint$coefficients[1,1] + Asia_lm_nonzeroint$coefficients[2,1]*sqrt(LM), sqrt(LM), col = "black", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 1), ylim = c(0,.1), pch = 20, col = "blue")
lines((South_America_lm_nonzeroint$coefficients[1,1] + South_America_lm_nonzeroint$coefficients[2,1]*sqrt(LM))^2, LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines((West_Africa_lm_nonzeroint$coefficients[1,1] + West_Africa_lm_nonzeroint$coefficients[2,1]*sqrt(LM))^2, LM, col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines((East_Africa_lm_nonzeroint$coefficients[1,1] + East_Africa_lm_nonzeroint$coefficients[2,1]*sqrt(LM))^2, LM, col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines((Asia_lm_nonzeroint$coefficients[1,1] + Asia_lm_nonzeroint$coefficients[2,1]*sqrt(LM))^2, LM, col = "black", lwd = 2)


#Fitting Linear Models To Data For All Data Data Points
#LOGIT TRANSFORMATION

#Removing Zeroes
data_frame <- Rough_Preliminary_Data_Exploration_For_R[!(Rough_Preliminary_Data_Exploration_For_R$PCR_Percentage == 0 | 
                                                           Rough_Preliminary_Data_Exploration_For_R$LM_Percentage == 0),]


#Manipulate data_frame in prep for analysis
data_frame$Old_Or_New <- as.factor(data_frame$Old_Or_New)
data_frame$Global_Region <- as.factor(data_frame$Global_Region)
data_frame <- data_frame[order(data_frame$PCR_Percentage),]

#Create individual datasets for each global region
South_America <- data_frame[data_frame$Global_Region == "South America",]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa",]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa",]
Asia <- data_frame[data_frame$Global_Region == "Asia",]

#Defining LM variable for plotting
LM <- seq(0,1,0.001)

#Intercept Set to 0
West_Africa_lm_zeroint <- summary(lm(logit(West_Africa$PCR_Percentage) ~ 0 + logit(West_Africa$LM_Percentage)))
East_Africa_lm_zeroint <- summary(lm(logit(East_Africa$PCR_Percentage) ~ 0 + logit(East_Africa$LM_Percentage)))
Asia_lm_zeroint <- summary(lm(logit(Asia$PCR_Percentage) ~ 0 + logit(Asia$LM_Percentage)))
South_America_lm_zeroint <- summary(lm(logit(South_America$PCR_Percentage) ~ 0 + logit(South_America$LM_Percentage)))

#logit scale plotting
plot(logit(South_America$PCR_Percentage), logit(South_America$LM_Percentage), xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(South_America_lm_zeroint$coefficients[1]*logit(LM), logit(LM), col = "blue", lwd = 2)
points(logit(West_Africa$PCR_Percentage), logit(West_Africa$LM_Percentage), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_zeroint$coefficients[1]*logit(LM), logit(LM), col = "red", lwd = 2)
points(logit(East_Africa$PCR_Percentage), logit(East_Africa$LM_Percentage), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_zeroint$coefficients[1]*logit(LM), logit(LM), col = "green", lwd = 2)
points(logit(Asia$PCR_Percentage), logit(Asia$LM_Percentage), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_zeroint$coefficients[1]*logit(LM), logit(LM), col = "black", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(expit(South_America_lm_zeroint$coefficients[1]*logit(LM)), expit(logit(LM)), col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(expit(West_Africa_lm_zeroint$coefficients[1]*logit(LM)), expit(logit(LM)), col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(expit(East_Africa_lm_zeroint$coefficients[1]*logit(LM)), LM, col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(Asia_lm_zeroint$coefficients[1]*logit(LM)), LM, col = "black", lwd = 2)

#Intercept Not Set to 0
West_Africa_lm_nonzeroint <- summary(lm(logit(West_Africa$PCR_Percentage) ~ logit(West_Africa$LM_Percentage)))
East_Africa_lm_nonzeroint <- summary(lm(logit(East_Africa$PCR_Percentage) ~ logit(East_Africa$LM_Percentage)))
Asia_lm_nonzeroint <- summary(lm(logit(Asia$PCR_Percentage) ~ logit(Asia$LM_Percentage)))
South_America_lm_nonzeroint <- summary(lm(logit(South_America$PCR_Percentage) ~ logit(South_America$LM_Percentage)))

#logit scale plotting
plot(logit(South_America$PCR_Percentage), logit(South_America$LM_Percentage), xlim = c(-10, 10), ylim = c(-10,10), pch = 20, col = "blue")
lines(South_America_lm_nonzeroint$coefficients[1,1] + South_America_lm_nonzeroint$coefficients[2,1]*logit(LM), logit(LM), col = "blue", lwd = 2)
points(logit(West_Africa$PCR_Percentage), logit(West_Africa$LM_Percentage), xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(West_Africa_lm_nonzeroint$coefficients[1,1] + West_Africa_lm_nonzeroint$coefficients[2,1]*logit(LM), logit(LM), col = "red", lwd = 2)
points(logit(East_Africa$PCR_Percentage), logit(East_Africa$LM_Percentage), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(East_Africa_lm_nonzeroint$coefficients[1,1] + East_Africa_lm_nonzeroint$coefficients[2,1]*logit(LM), logit(LM), col = "green", lwd = 2)
points(logit(Asia$PCR_Percentage), logit(Asia$LM_Percentage), xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(Asia_lm_nonzeroint$coefficients[1,1] + Asia_lm_nonzeroint$coefficients[2,1]*logit(LM), logit(LM), col = "black", lwd = 2)

#natural scale plotting
plot(South_America$PCR_Percentage, South_America$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(expit(South_America_lm_nonzeroint$coefficients[1,1] + South_America_lm_nonzeroint$coefficients[2,1]*logit(LM)), LM, col = "blue", lwd = 2)
points(West_Africa$PCR_Percentage, West_Africa$LM_Percentage, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "red")
lines(expit(West_Africa_lm_nonzeroint$coefficients[1,1] + West_Africa_lm_nonzeroint$coefficients[2,1]*logit(LM)), LM, col = "red", lwd = 2)
points(East_Africa$PCR_Percentage, East_Africa$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "green")
lines(expit(East_Africa_lm_nonzeroint$coefficients[1,1] + East_Africa_lm_nonzeroint$coefficients[2,1]*logit(LM)), LM, col = "green", lwd = 2)
points(Asia$PCR_Percentage, Asia$LM_Percentage, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(Asia_lm_nonzeroint$coefficients[1,1] + Asia_lm_nonzeroint$coefficients[2,1]*logit(LM)), LM, col = "black", lwd = 2)