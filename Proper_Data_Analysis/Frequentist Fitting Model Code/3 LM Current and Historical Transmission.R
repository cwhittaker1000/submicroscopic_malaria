#Access relevant packages
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

#Assign whichever dataset being used to data_frame
data_frame <- Data.4th.February

#Subset and manipulate the data 
data_frame$Old_or_New <- as.factor(data_frame$Old_or_New)
data_frame$Trans_Hist <- data_frame$Transmission_Setting_History
data_frame$Trans_Hist <- as.factor(data_frame$Trans_Hist)
data_frame$Trans_Now <- data_frame$Transmission_Setting_Current
data_frame$Trans_Now <- as.factor(data_frame$Trans_Now)

#Add columns for composite components
data_frame$high_low <-
data_frame$high_high
data_frame$low_low
data_frame$


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

#Changing NAs to 0 for places where I don't have data
for (i in 1:341) {
  if (is.na(data_frame[i,15]) == T) {
    data_frame[i, 15] <- 0 
  }
}

for (i in 1:341) {
  if (is.na(data_frame[i,16]) == T) {
    data_frame[i, 16] <- 0 
  }
}

#Removing Zeroes
data_frame <- data_frame[!data_frame$Was_Initially_Zero == "Y",]

#Transmission History Assessment
#Note: ask them about this, but I struggled to subset using !NA- kept the NAs in there for some reason. Had to adapt
# and do the below instead!
#if 0s not removed
data_frame <- data_frame[order(data_frame$Trans_Hist),]
data_frame <- data_frame[1:234,]

#if 0s removed
data_frame <- data_frame[order(data_frame$Trans_Hist),]
data_frame <- data_frame[1:229,]

#Transmission Setting History, LOGIT TRANSFORMATION
#NOTE REQUIRES REMOVAL OF ZEROES ATM- NEED TO ADDRESS THIS

#HISTORICALLY LOW
subset_hist_low <- data_frame[data_frame$Trans_Hist == 0,]
subset_hist_low <- subset_hist_low[order(subset_hist_low$Micro_Prev),]

#Logit transformed data, no 0 intercept
#NOTE REGRESSING MICROSCOPY ONTO PCR, NOT OTHER WAY AROUND
low_nonzeroint <- summary(lm(logit(subset_hist_low$Micro_Prev) ~ logit(subset_hist_low$PCR_Prev)))
logit_PCR <- logit(seq(0,1,0.001))
plot(logit(subset_hist_low$Micro_Prev), logit(subset_hist_low$PCR_Prev))#, xlim =c(-10,10), ylim = c(-10,10))
lines(low_nonzeroint$coefficients[1,1] + low_nonzeroint$coefficients[2,1]*logit_PCR, logit_PCR)

plot(subset_hist_low$PCR_Prev, subset_hist_low$Micro_Prev, pch = 20, col ="black")
lines(expit(logit_PCR), expit(low_nonzeroint$coefficients[1,1] + low_nonzeroint$coefficients[2,1]*logit_PCR))

#NOTE THAT REVERSING THE REGRESSION IE PCR ONTO MICROSCOPY, CHANGES THE RESULTS QUITE DRAMATICALLY
# low_nonzeroint <- summary(lm(logit(subset_hist_low$PCR_Prev)~ logit(subset_hist_low$Micro_Prev)))
# logit_LM <- logit(seq(0,1,0.001))
# plot(logit(subset_hist_low$Micro_Prev), logit(subset_hist_low$PCR_Prev),  xlim =c(-10,10), ylim = c(-10,10))
# lines(logit_LM, low_nonzeroint$coefficients[1,1] + low_nonzeroint$coefficients[2,1]*logit_LM)
# 
# plot(subset_hist_low$PCR_Prev, subset_hist_low$Micro_Prev, pch = 20, col ="black")
# lines(expit(low_nonzeroint$coefficients[1,1] + low_nonzeroint$coefficients[2,1]*logit_LM), expit(logit_LM))

#HISTORICALLY HIGH
subset_hist_high <- data_frame[data_frame$Trans_Hist == 1,]
subset_hist_high <- subset_hist_high[order(subset_hist_high$Micro_Prev),]

#Logit transformed data, no 0 intercept
high_nonzeroint <- summary(lm(logit(subset_hist_high$Micro_Prev) ~ logit(subset_hist_high$PCR_Prev)))
logit_PCR <- logit(seq(0,1,0.001))
plot(logit(subset_hist_high$Micro_Prev), logit(subset_hist_high$PCR_Prev)) #, xlim =c(-10,10), ylim = c(-10,10))
lines(high_nonzeroint$coefficients[1,1] + high_nonzeroint$coefficients[2,1]*logit_PCR, logit_PCR)

plot(subset_hist_high$PCR_Prev, subset_hist_high$Micro_Prev)
lines(expit(logit_PCR), expit(high_nonzeroint$coefficients[1,1] + high_nonzeroint$coefficients[2,1]*logit_PCR))

plot(subset_hist_low$PCR_Prev, subset_hist_low$Micro_Prev, xlim = c(0,1), ylim = c(0,1), col = "limegreen", pch =20)
lines(expit(logit_PCR), expit(low_nonzeroint$coefficients[1,1] + low_nonzeroint$coefficients[2,1]*logit_PCR), col = "limegreen", lwd = 2)
points(subset_hist_high$PCR_Prev, subset_hist_high$Micro_Prev, pch = 20, col = "lightslateblue")
lines(expit(logit_PCR), expit(high_nonzeroint$coefficients[1,1] + high_nonzeroint$coefficients[2,1]*logit_PCR), lwd =2, col ="lightslateblue")

legend("topleft", legend = c("High Historically", "Low Historically"), 
       col = c("lightslateblue", "limegreen"), lwd = c(2,2))

#Transmission Setting Curent, LOGIT TRANSFORMATION
#NOTE REQUIRES REMOVAL OF ZEROES ATM- NEED TO ADDRESS THIS

#CURRENTLY LOW
subset_curr_low <- data_frame[data_frame$Trans_Now == 0,]
subset_curr_low <- subset_curr_low[order(subset_curr_low$Micro_Prev),]

#Logit transformed data, no 0 intercept
#NOTE REGRESSING MICROSCOPY ONTO PCR, NOT OTHER WAY AROUND
low_nonzeroint <- summary(lm(logit(subset_curr_low$Micro_Prev) ~ logit(subset_curr_low$PCR_Prev)))
logit_PCR <- logit(seq(0,1,0.001))
plot(logit(subset_curr_low$Micro_Prev), logit(subset_curr_low$PCR_Prev))#, xlim =c(-10,10), ylim = c(-10,10))
lines(low_nonzeroint$coefficients[1,1] + low_nonzeroint$coefficients[2,1]*logit_PCR, logit_PCR)

plot(subset_curr_low$PCR_Prev, subset_curr_low$Micro_Prev, xlim = c(0,1), ylim = c(0,1))
lines(expit(logit_PCR), expit(low_nonzeroint$coefficients[1,1] + low_nonzeroint$coefficients[2,1]*logit_PCR))

#CURRENTLY HIGH
subset_curr_high <- data_frame[data_frame$Trans_Now == 1,]
subset_curr_high <- subset_curr_high[order(subset_curr_high$Micro_Prev),]

#Logit transformed data, no 0 intercept
high_nonzeroint <- summary(lm(logit(subset_curr_high$Micro_Prev) ~ logit(subset_curr_high$PCR_Prev)))
logit_PCR <- logit(seq(0,1,0.001))
plot(logit(subset_curr_high$Micro_Prev), logit(subset_curr_high$PCR_Prev), xlim =c(-10,10), ylim = c(-10,10))
lines(high_nonzeroint$coefficients[1,1] + high_nonzeroint$coefficients[2,1]*logit_PCR, logit_PCR)

plot(subset_curr_high$PCR_Prev, subset_curr_high$Micro_Prev)
lines(expit(logit_PCR), expit(high_nonzeroint$coefficients[1,1] + high_nonzeroint$coefficients[2,1]*logit_PCR))

plot(subset_curr_low$PCR_Prev, subset_curr_low$Micro_Prev, xlim = c(0,1), ylim = c(0,1), col = "purple", pch =20)
lines(expit(logit_PCR), expit(low_nonzeroint$coefficients[1,1] + low_nonzeroint$coefficients[2,1]*logit_PCR), col = "purple", lwd = 2)
points(subset_curr_high$PCR_Prev, subset_curr_high$Micro_Prev, pch = 20, col = "tomato")
lines(expit(logit_PCR), expit(high_nonzeroint$coefficients[1,1] + high_nonzeroint$coefficients[2,1]*logit_PCR), lwd =2, col ="tomato")


legend("topleft", legend = c("High Current", "Low Current"), 
       col = c("tomato", "purple"), lwd = c(2,2))

# summary(lm(subset_hist_high$PCR_Prev ~ subset_hist_high$Micro_Prev))
# MICRO_1 <- seq(0,1,0.001)
# PCR_1 <- 0.10336 + 1.16735*MICRO_1
# plot(subset_hist_high$PCR_Prev, subset_hist_high$Micro_Prev)
# lines(PCR_1, MICRO_1)
# 
# summary(lm(subset_hist_high$Micro_Prev ~ subset_hist_high$PCR_Prev))
# PCR_2 <- seq(0,1,0.001)
# MICRO_2 <- -0.007481 + 0.593509*PCR_2
# plot(subset_hist_high$PCR_Prev, subset_hist_high$Micro_Prev)
# lines(PCR_2, MICRO_2)
