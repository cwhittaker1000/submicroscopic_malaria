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

###OPTIONAL
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

###OPTIONAL
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

#DO THIS IF NAS HAVEN'T BEEN CHANGED
#if 0s not removed
data_frame <- data_frame[order(data_frame$Trans_Hist),]
data_frame <- data_frame[1:234,]

# checking whether including Lucy's data changes the conclusions
data_frame <- data_frame[data_frame$Old_or_New == "New",]


#if 0s removed
data_frame <- data_frame[order(data_frame$Trans_Hist),]
data_frame <- data_frame[1:229,]

#Add columns for composite components
for (i in 1:length(data_frame[,1])) {
  if (data_frame$Trans_Hist[i] == 1 & data_frame$Trans_Now[i] == 1) {
    data_frame$high_high[i] <- 1 }
    else {
      data_frame$high_high[i] <- 0
    }
  if (data_frame$Trans_Hist[i] == 1 & data_frame$Trans_Now[i] == 0) {
    data_frame$high_low[i] <- 1 }
    else {
      data_frame$high_low[i] <- 0 
      }
  if (data_frame$Trans_Hist[i] == 0 & data_frame$Trans_Now[i] == 1) {
    data_frame$low_high[i] <- 1 }
    else {
      data_frame$low_high[i] <- 0
    }
  if (data_frame$Trans_Hist[i] == 0 & data_frame$Trans_Now[i] == 0) {
    data_frame$low_low[i] <- 1 }
    else {
      data_frame$low_low[i] <- 0 
      }
}


#Creating Data Subsets
subset_high_high <- data_frame[data_frame$high_high == 1,]
subset_high_low <- data_frame[data_frame$high_low == 1,]
subset_low_high <- data_frame[data_frame$low_high == 1,]
subset_low_low <- data_frame[data_frame$low_low == 1,]

#Logit transformed data, no 0 intercept
#NOTE REGRESSING MICROSCOPY ONTO PCR, NOT OTHER WAY AROUND
high_high_lm <- summary(lm(logit(subset_high_high$Micro_Prev) ~ logit(subset_high_high$PCR_Prev)))
high_low_lm <- summary(lm(logit(subset_high_low$Micro_Prev) ~ logit(subset_high_low$PCR_Prev)))
low_high_lm <- summary(lm(logit(subset_low_high$Micro_Prev) ~ logit(subset_low_high$PCR_Prev)))
low_low_lm <- summary(lm(logit(subset_low_low$Micro_Prev) ~ logit(subset_low_low$PCR_Prev)))

logit_PCR <- logit(seq(0,1,0.001))

###PLOTTING ON THE LOGIT SCALE###

#High high
plot(logit(subset_high_high$Micro_Prev), logit(subset_high_high$PCR_Prev))#, xlim =c(-10,10), ylim = c(-10,10))
lines(high_high_lm$coefficients[1,1] + high_high_lm$coefficients[2,1]*logit_PCR, logit_PCR)

#High low
plot(logit(subset_high_low$Micro_Prev), logit(subset_high_low$PCR_Prev))#, xlim =c(-10,10), ylim = c(-10,10))
lines(high_low_lm$coefficients[1,1] + high_low_lm$coefficients[2,1]*logit_PCR, logit_PCR)

#Low high
plot(logit(subset_low_high$Micro_Prev), logit(subset_low_high$PCR_Prev))#, xlim =c(-10,10), ylim = c(-10,10))
lines(low_high_lm$coefficients[1,1] + low_high_lm$coefficients[2,1]*logit_PCR, logit_PCR)

#Low low
plot(logit(subset_low_low$Micro_Prev), logit(subset_low_low$PCR_Prev))#, xlim =c(-10,10), ylim = c(-10,10))
lines(low_low_lm$coefficients[1,1] + low_low_lm$coefficients[2,1]*logit_PCR, logit_PCR)


###PLOTTING ON THE NATURAL SCALE###

#High high
plot(subset_high_high$PCR_Prev, subset_high_high$Micro_Prev)
lines(expit(logit_PCR), expit(high_high_lm$coefficients[1,1] + high_high_lm$coefficients[2,1]*logit_PCR))

#High low
plot(subset_high_low$PCR_Prev, subset_high_low$Micro_Prev)
lines(expit(logit_PCR), expit(high_low_lm$coefficients[1,1] + high_low_lm$coefficients[2,1]*logit_PCR))

#Low high
plot(subset_low_high$PCR_Prev, subset_low_high$Micro_Prev)
lines(expit(logit_PCR), expit(low_high_lm$coefficients[1,1] + low_high_lm$coefficients[2,1]*logit_PCR))

#Low low
plot(subset_low_low$PCR_Prev, subset_low_low$Micro_Prev)
lines(expit(logit_PCR), expit(low_low_lm$coefficients[1,1] + low_low_lm$coefficients[2,1]*logit_PCR))

#Altogether now
plot(subset_high_high$PCR_Prev, subset_high_high$Micro_Prev, pch = 20, col = "slateblue4")
lines(expit(logit_PCR), expit(high_high_lm$coefficients[1,1] + high_high_lm$coefficients[2,1]*logit_PCR), lwd = 2, col = "slateblue4")
points(subset_high_low$PCR_Prev, subset_high_low$Micro_Prev, pch = 20, col = "tan1")
lines(expit(logit_PCR), expit(high_low_lm$coefficients[1,1] + high_low_lm$coefficients[2,1]*logit_PCR), lwd = 2, col = "tan1")
points(subset_low_high$PCR_Prev, subset_low_high$Micro_Prev, pch = 20, col = "maroon1")
lines(expit(logit_PCR), expit(low_high_lm$coefficients[1,1] + low_high_lm$coefficients[2,1]*logit_PCR), lwd = 2, col = "maroon1")
points(subset_low_low$PCR_Prev, subset_low_low$Micro_Prev, pch = 20, col = "green3")
lines(expit(logit_PCR), expit(low_low_lm$coefficients[1,1] + low_low_lm$coefficients[2,1]*logit_PCR), lwd = 2, col = "green3")

legend("topleft", legend = c("High High", "High Low", "Low High", "Low Low"), 
       col = c("slateblue4", "tan1", "maroon1", "green3"), lwd = c(2,2,2,2))

#Investigating a linear fit
plot(subset_low_low$PCR_Prev, subset_low_low$Micro_Prev, pch = 20, col = "green3")
low_low_linear_lm <- summary(lm(subset_low_low$Micro_Prev ~ 0 + subset_low_low$PCR_Prev))
lines(expit(logit_PCR), low_low_linear_lm$coefficients[1,1]*expit(logit_PCR), lwd = 2, col = "green3")




  