#Access relevant packages
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

#Assign whichever dataset being used to data_frame
data_frame <- Data.4th.February.Dec.Rounded

hist(logit(data_frame$PCR_Prev))
hist(logit(data_frame$Micro_Prev))

####SQUARE ROOT TRANSFORMATION STUFF####

#Subset the data into old and new
old_data <- data_frame[data_frame$Old_or_New == "Old",]
new_data <- data_frame[data_frame$Old_or_New == "New",]

###Basic linear model for old data with sqrt transformation
##OLD DATA
#With intercept set to 0
old_data_basic_model_intercept_zero <- summary(lm(sqrt(old_data$PCR_Prev) ~ 0 + sqrt(old_data$Micro_Prev)))
LM_full <- seq(0,1,0.001)
plot(sqrt(old_data$PCR_Prev), sqrt(old_data$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20)
lines(old_data_basic_model_intercept_zero$coefficients[1]*sqrt(LM_full), sqrt(LM_full), col = "black", lwd = 2)

plot(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines((old_data_basic_model_intercept_zero$coefficients[1]*sqrt(LM_full))^2, sqrt(LM_full)^2, col = "blue", lwd = 2)

#Without intercept set to 0
old_data_basic_model_intercept_any <- summary(lm(sqrt(old_data$PCR_Prev) ~ sqrt(old_data$Micro_Prev)))
LM_full <- seq(0,1,0.001)
plot(sqrt(old_data$PCR_Prev), sqrt(old_data$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20)
lines(old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full), sqrt(LM_full), col = "black", lwd = 2)

plot(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines((old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full))^2, sqrt(LM_full)^2, col = "blue", lwd = 2)

##NEW DATA
#With intercept set to 0
new_data_basic_model_intercept_zero <- summary(lm(sqrt(new_data$PCR_Prev) ~ 0 + sqrt(new_data$Micro_Prev)))
LM_full <- seq(0,1,0.001)
plot(sqrt(new_data$PCR_Prev), sqrt(new_data$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20)
lines(new_data_basic_model_intercept_zero$coefficients[1]*sqrt(LM_full), sqrt(LM_full), col = "black", lwd = 2)

plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines((new_data_basic_model_intercept_zero$coefficients[1]*sqrt(LM_full))^2, sqrt(LM_full)^2, col = "blue", lwd = 2)

#Without intercept set to 0
new_data_basic_model_intercept_any <- summary(lm(sqrt(new_data$PCR_Prev) ~ sqrt(new_data$Micro_Prev)))
LM_full <- seq(0,1,0.001)
plot(sqrt(new_data$PCR_Prev), sqrt(new_data$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20)
lines(new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full), sqrt(LM_full), col = "black", lwd = 2)

plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines((new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full))^2, sqrt(LM_full)^2, col = "blue", lwd = 2)

##ALL DATA
#With intercept set to 0
data_frame_basic_model_intercept_zero <- summary(lm(sqrt(data_frame$PCR_Prev) ~ 0 + sqrt(data_frame$Micro_Prev)))
LM_full <- seq(0,1,0.001)
plot(sqrt(data_frame$PCR_Prev), sqrt(data_frame$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20)
lines(data_frame_basic_model_intercept_zero$coefficients[1]*sqrt(LM_full), sqrt(LM_full), col = "black", lwd = 2)

plot(data_frame$PCR_Prev, data_frame$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines((data_frame_basic_model_intercept_zero$coefficients[1]*sqrt(LM_full))^2, sqrt(LM_full)^2, col = "blue", lwd = 2)

#Without intercept set to 0
data_frame_basic_model_intercept_any <- summary(lm(sqrt(data_frame$PCR_Prev) ~ sqrt(data_frame$Micro_Prev)))
LM_full <- seq(0,1,0.001)
plot(sqrt(data_frame$PCR_Prev), sqrt(data_frame$Micro_Prev), xlim = c(0, 1), ylim = c(0,1), pch = 20)
lines(data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full), sqrt(LM_full), col = "black", lwd = 2)

plot(data_frame$PCR_Prev, data_frame$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines((data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full))^2, sqrt(LM_full)^2, col = "blue", lwd = 2)


###summary of sqrt transformation stuff
plot(data_frame$PCR_Prev, data_frame$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines((data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full))^2,
      sqrt(LM_full)^2, col = "black", lwd = 2)
points(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines((new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full))^2,
      sqrt(LM_full)^2, col = "blue", lwd = 2)
points(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "orange")
lines((old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*sqrt(LM_full))^2,
      sqrt(LM_full)^2, col = "orange", lwd = 2)

legend("topleft", legend = c("Old", "New", "Combined"), col = c("orange", "blue", "black"),
       lwd = c(2,2,2))

####LOGIT TRANSFORMATION STUFF####
#NOTE REQUIRES REMOVAL OF ZEROES ATM- NEED TO ADDRESS THIS

#Assign whichever dataset being used to data_frame
data_frame <- Data.4th.February

#Remove 0s from the data_frame so that logit works
data_frame <- data_frame[!data_frame$Was_Initially_Zero == "Y",]

#Changes the 0s to Lucy's suggestion
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


#Subset the data into old and new
old_data <- data_frame[data_frame$Old_or_New == "Old",]
new_data <- data_frame[data_frame$Old_or_New == "New",]

##OLD DATA
#Basic linear model with logit transformation
#With intercept set to 0
old_data_basic_model_intercept_zero <- summary(lm(logit(old_data$PCR_Prev) ~ 0 + logit(old_data$Micro_Prev)))
LM_full <- seq(0,1,0.00001)
plot(logit(old_data$PCR_Prev), logit(old_data$Micro_Prev), pch = 20) # xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(old_data_basic_model_intercept_zero$coefficients[1]*logit(LM_full), logit(LM_full), col = "black", lwd = 2)

plot(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(old_data_basic_model_intercept_zero$coefficients[1]*logit(LM_full)), expit(logit(LM_full)), col = "blue", lwd = 2)

#Without intercept set to 0
old_data_basic_model_intercept_any <- summary(lm(logit(old_data$PCR_Prev) ~ logit(old_data$Micro_Prev)))
LM_full <- seq(0,1,0.00001)
plot(logit(old_data$PCR_Prev), logit(old_data$Micro_Prev), pch = 20) #, xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*logit(LM_full), logit(LM_full), col = "black", lwd = 2)

plot(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*logit(LM_full)), expit(logit(LM_full)), col = "blue", lwd = 2)

##NEW DATA
#With intercept set to 0
new_data_basic_model_intercept_zero <- summary(lm(logit(new_data$PCR_Prev) ~ 0 + logit(new_data$Micro_Prev)))
LM_full <- seq(0,1,0.00001)
plot(logit(new_data$PCR_Prev), logit(new_data$Micro_Prev), pch = 20) #, xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(new_data_basic_model_intercept_zero$coefficients[1]*logit(LM_full), logit(LM_full), col = "black", lwd = 2)

plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(new_data_basic_model_intercept_zero$coefficients[1]*logit(LM_full)), expit(logit(LM_full)), col = "blue", lwd = 2)

#Without intercept set to 0
new_data_basic_model_intercept_any <- summary(lm(logit(new_data$PCR_Prev) ~ logit(new_data$Micro_Prev)))
LM_full <- seq(0,1,0.000001)
plot(logit(new_data$PCR_Prev), logit(new_data$Micro_Prev), pch = 20, xlim = c(-10, 3), ylim = c(-12, 3)) # xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*logit(LM_full), logit(LM_full), col = "black", lwd = 2)

plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*logit(LM_full)), expit(logit(LM_full)), col = "blue", lwd = 2)

##ALL DATA
#With intercept set to 0
data_frame_basic_model_intercept_zero <- summary(lm(logit(data_frame$PCR_Prev) ~ 0 + logit(data_frame$Micro_Prev)))
LM_full <- seq(0,1,0.001)
plot(logit(data_frame$PCR_Prev), logit(data_frame$Micro_Prev), xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(data_frame_basic_model_intercept_zero$coefficients[1]*logit(LM_full), logit(LM_full), col = "black", lwd = 2)

plot(data_frame$PCR_Prev, data_frame$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(data_frame_basic_model_intercept_zero$coefficients[1]*logit(LM_full)), expit(logit(LM_full)), col = "blue", lwd = 2)

#Without intercept set to 0
data_frame_basic_model_intercept_any <- summary(lm(logit(data_frame$PCR_Prev) ~ logit(data_frame$Micro_Prev)))
LM_full <- seq(0,1,0.0001)
plot(logit(data_frame$PCR_Prev), logit(data_frame$Micro_Prev), xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*logit(LM_full), logit(LM_full), col = "black", lwd = 2)

plot(data_frame$PCR_Prev, data_frame$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*logit(LM_full)), expit(logit(LM_full)), col = "blue", lwd = 2)


###summary of logit transformation stuff
plot(data_frame$PCR_Prev, data_frame$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*logit(LM_full)),
      expit(logit(LM_full)), col = "black", lwd = 2)
points(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(expit(new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*logit(LM_full)),
      expit(logit(LM_full)), col = "blue", lwd = 2)
points(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "orange")
lines(expit(old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*logit(LM_full)),
      expit(logit(LM_full)), col = "orange", lwd = 2)

x <- seq(0,1,0.001)
y <- seq(0,1,0.001)
lines(x,y)

### PCR Regressed Onto Microscopy

#Old Data Without intercept set to 0
old_data_basic_model_intercept_any <- summary(lm(logit(old_data$Micro_Prev) ~ logit(old_data$PCR_Prev)))
PCR_full <- seq(0,1,0.0001)
plot(logit(old_data$PCR_Prev), logit(old_data$Micro_Prev), pch = 20) #, xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(logit(PCR_full), old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full), col = "black", lwd = 2)

plot(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(PCR_full, expit(old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full)), col = "blue", lwd = 2)

#New Data without intercept set to 0
new_data_basic_model_intercept_any <- summary(lm(logit(new_data$Micro_Prev) ~ logit(new_data$PCR_Prev)))
PCR_full <- seq(0,1,0.0001)
plot(logit(new_data$PCR_Prev), logit(new_data$Micro_Prev), pch = 20, xlim = c(-10, 3), ylim = c(-12, 3)) # xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(logit(PCR_full), new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full), col = "orange", lwd = 2)

plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "orange")
lines(PCR_full, expit(new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full)), col = "orange", lwd = 2)

#Full Data
#Without intercept set to 0
data_frame_basic_model_intercept_any <- summary(lm(logit(data_frame$Micro_Prev) ~ logit(data_frame$PCR_Prev)))
PCR_full <- seq(0,1,0.0001)
plot(logit(data_frame$PCR_Prev), logit(data_frame$Micro_Prev), xlim = c(-5, 5), ylim = c(-7,7), pch = 20)
lines(logit(PCR_full), data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full), col = "black", lwd = 2)

plot(data_frame$PCR_Prev, data_frame$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(PCR_full, expit(data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full)), col = "black", lwd = 2)

#Altogether now
plot(data_frame$PCR_Prev, data_frame$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "black")
lines(PCR_full, expit(data_frame_basic_model_intercept_any$coefficients[1,1]+data_frame_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full)), col = "black", lwd = 2)
points(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "blue")
lines(PCR_full, expit(old_data_basic_model_intercept_any$coefficients[1,1]+old_data_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full)), col = "blue", lwd = 2)
points(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0,1), pch = 20, col = "orange")
lines(PCR_full, expit(new_data_basic_model_intercept_any$coefficients[1,1]+new_data_basic_model_intercept_any$coefficients[2,1]*logit(PCR_full)), col = "orange", lwd = 2)

legend("topleft", legend = c("Old", "New", "Full"), lwd = c(2,2,2), col = c("blue","orange","black"))

