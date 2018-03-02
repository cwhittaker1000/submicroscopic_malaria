#Access relevant packages
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

#Assign whichever dataset being used to data_frame
data_frame <- Data.4th.February

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

#Subsetting by PCR Prevalence
subset_low <- data_frame[data_frame$PCR_Prev < 0.2,]
subset_mid <- data_frame[data_frame$PCR_Prev < 0.5 & data_frame$PCR_Prev >= 0.2,]
subset_high <- data_frame[data_frame$PCR_Prev > 0.5,]

#For remvoing zeroes
subset_low <- subset_low[!(subset_low$PCR_Prev == 0 | subset_low$Micro_Prev == 0),]
subset_mid <- subset_mid[!(subset_mid$PCR_Prev == 0 | subset_mid$Micro_Prev == 0),]
subset_high <- subset_high[!(subset_high$PCR_Prev == 0 | subset_high$Micro_Prev == 0),]


#Creating the Linear Models, Microscopy Regressed Onto PCR
low_PCR <- summary(lm(logit(subset_low$PCR_Prev) ~ logit(subset_low$Micro_Prev)))
mid_PCR <- summary(lm(logit(subset_mid$PCR_Prev) ~ logit(subset_mid$Micro_Prev)))
high_PCR <- summary(lm(logit(subset_high$PCR_Prev) ~ logit(subset_high$Micro_Prev)))

#Defining logit of microscopy
logit_LM <- logit(seq(0,1,0.001))
LM <- expit(logit_LM)

#Plotting on the logit scale
plot(logit(subset_low$Micro_Prev), logit(subset_low$PCR_Prev),  xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "black")
lines(logit_LM, low_PCR$coefficients[1,1] + low_PCR$coefficients[2,1]*logit_LM)
points(logit(subset_mid$Micro_Prev), logit(subset_mid$PCR_Prev),  xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(logit_LM, mid_PCR$coefficients[1,1] + mid_PCR$coefficients[2,1]*logit_LM, col = "blue")
points(logit(subset_high$Micro_Prev), logit(subset_high$PCR_Prev),  xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(logit_LM, high_PCR$coefficients[1,1] + high_PCR$coefficients[2,1]*logit_LM, col = "green")

#Plotting on the natural scale
plot(subset_low$PCR_Prev, subset_low$Micro_Prev, xlim =c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(low_PCR$coefficients[1,1] + low_PCR$coefficients[2,1]*logit_LM), LM)
points(subset_mid$PCR_Prev, subset_mid$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(expit(mid_PCR$coefficients[1,1] + mid_PCR$coefficients[2,1]*logit_LM), LM, col = "blue")
points(subset_high$PCR_Prev, subset_high$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(expit(high_PCR$coefficients[1,1] + high_PCR$coefficients[2,1]*logit_LM), LM, col = "green")

#Creating the Linear Models, PCR Regressed Onto Microscopy
low_Micro <- summary(lm(logit(subset_low$Micro_Prev) ~ logit(subset_low$PCR_Prev)))
mid_Micro <- summary(lm(logit(subset_mid$Micro_Prev) ~ logit(subset_mid$PCR_Prev)))
high_Micro <- summary(lm(logit(subset_high$Micro_Prev) ~ logit(subset_high$PCR_Prev)))

#Defining logit of microscopy
logit_PCR <- logit(seq(0,1,0.001))
PCR <- expit(logit_PCR)

#Plotting on the logit scale
plot(logit(subset_low$PCR_Prev), logit(subset_low$Micro_Prev), xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "black")
lines(logit_PCR, low_Micro$coefficients[1,1] + low_Micro$coefficients[2,1]*logit_PCR)
points(logit(subset_mid$PCR_Prev), logit(subset_mid$Micro_Prev),   xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(logit_PCR, mid_Micro$coefficients[1,1] + mid_Micro$coefficients[2,1]*logit_PCR, col = "blue")
points(logit(subset_high$PCR_Prev), logit(subset_high$Micro_Prev),  xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(logit_PCR, high_Micro$coefficients[1,1] + high_PCR$coefficients[2,1]*logit_PCR, col = "green")

#Plotting on the natural scale, full lines
plot(subset_low$PCR_Prev, subset_low$Micro_Prev, xlim =c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(PCR, expit(low_Micro$coefficients[1,1] + low_Micro$coefficients[2,1]*logit_PCR), lwd = 2)
points(subset_mid$PCR_Prev, subset_mid$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(PCR, expit(mid_Micro$coefficients[1,1] + mid_Micro$coefficients[2,1]*logit_PCR), col = "blue", lwd = 2)
points(subset_high$PCR_Prev, subset_high$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(PCR, expit(high_Micro$coefficients[1,1] + high_Micro$coefficients[2,1]*logit_PCR), col = "green", lwd = 2)

#Plotting on the natural scale, lines for each section's bit only
PCR_low <- logit(seq(0,0.2,0.001))
PCR_mid <- logit(seq(0.2,0.5,0.001))
PCR_high <- logit(seq(0.5,1,0.001))

plot(subset_low$PCR_Prev, subset_low$Micro_Prev, xlim =c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(PCR_low), expit(low_Micro$coefficients[1,1] + low_Micro$coefficients[2,1]*PCR_low), lwd = 2)
points(subset_mid$PCR_Prev, subset_mid$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(expit(PCR_mid), expit(mid_Micro$coefficients[1,1] + mid_Micro$coefficients[2,1]*PCR_mid), col = "blue", lwd = 2)
points(subset_high$PCR_Prev, subset_high$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(expit(PCR_high), expit(high_Micro$coefficients[1,1] + high_Micro$coefficients[2,1]*PCR_high), col = "green", lwd = 2)

legend("topleft", legend = c("Low", "Mid", "High"), 
       col = c("black", "blue", "green"), lwd = c(2,2,2))

#Assign whichever dataset being used to data_frame
data_frame <- Collated_Data


########Subsetting by Microscopy Prevalence
subset_low <- data_frame[data_frame$Micro_Prev < 0.2,]
subset_mid <- data_frame[data_frame$Micro_Prev < 0.5 & data_frame$Micro_Prev >= 0.2,]
subset_high <- data_frame[data_frame$Micro_Prev > 0.5,]

#Removing zeroes if needs be
subset_low <- subset_low[!(subset_low$PCR_Prev == 0 | subset_low$Micro_Prev == 0),]
subset_mid <- subset_mid[!(subset_mid$PCR_Prev == 0 | subset_mid$Micro_Prev == 0),]
subset_high <- subset_high[!(subset_high$PCR_Prev == 0 | subset_high$Micro_Prev == 0),]

#Creating the Linear Models, Microscopy Regressed Onto PCR
low_Micro <- summary(lm(logit(subset_low$Micro_Prev) ~ logit(subset_low$PCR_Prev)))
mid_Micro <- summary(lm(logit(subset_mid$Micro_Prev) ~ logit(subset_mid$PCR_Prev)))
high_Micro <- summary(lm(logit(subset_high$Micro_Prev) ~ logit(subset_high$PCR_Prev)))

#Defining logit of PCR
logit_PCR <- logit(seq(0,1,0.001))
PCR <- expit(logit_PCR)

#Plotting on the logit scale
plot(logit(subset_low$Micro_Prev), logit(subset_low$PCR_Prev)) #,  xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "black")
lines(low_Micro$coefficients[1,1] + low_Micro$coefficients[2,1]*logit_PCR, logit_PCR)
plot(logit(subset_mid$Micro_Prev), logit(subset_mid$PCR_Prev)) #,  xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(mid_Micro$coefficients[1,1] + mid_Micro$coefficients[2,1]*logit_PCR, logit_PCR)
plot(logit(subset_high$Micro_Prev), logit(subset_high$PCR_Prev)) #,  xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(high_Micro$coefficients[1,1] + high_Micro$coefficients[2,1]*logit_PCR, logit_PCR)

#Plotting on the natural scale
plot(subset_low$PCR_Prev, subset_low$Micro_Prev, xlim =c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(PCR, expit(low_Micro$coefficients[1,1] + low_Micro$coefficients[2,1]*logit_PCR))
points(subset_mid$PCR_Prev, subset_mid$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(PCR, expit(mid_Micro$coefficients[1,1] + mid_Micro$coefficients[2,1]*logit_PCR), col = "blue")
points(subset_high$PCR_Prev, subset_high$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(PCR, expit(high_Micro$coefficients[1,1] + high_Micro$coefficients[2,1]*logit_PCR), col = "green")

#Creating the Linear Models, PCR Regressed Onto Microscopy
low_Micro <- summary(lm(logit(subset_low$Micro_Prev) ~ logit(subset_low$PCR_Prev)))
mid_Micro <- summary(lm(logit(subset_mid$Micro_Prev) ~ logit(subset_mid$PCR_Prev)))
high_Micro <- summary(lm(logit(subset_high$Micro_Prev) ~ logit(subset_high$PCR_Prev)))

#Defining logit of microscopy
logit_PCR <- logit(seq(0,1,0.001))
PCR <- expit(logit_PCR)

#Plotting on the logit scale
plot(logit(subset_low$PCR_Prev), logit(subset_low$Micro_Prev), xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "black")
lines(logit_PCR, low_Micro$coefficients[1,1] + low_Micro$coefficients[2,1]*logit_PCR)
points(logit(subset_mid$PCR_Prev), logit(subset_mid$Micro_Prev),   xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(logit_PCR, mid_Micro$coefficients[1,1] + mid_Micro$coefficients[2,1]*logit_PCR, col = "blue")
points(logit(subset_high$PCR_Prev), logit(subset_high$Micro_Prev),  xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(logit_PCR, high_Micro$coefficients[1,1] + high_PCR$coefficients[2,1]*logit_PCR, col = "green")

#Plotting on the natural scale
plot(subset_low$PCR_Prev, subset_low$Micro_Prev, xlim =c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(PCR, expit(low_Micro$coefficients[1,1] + low_Micro$coefficients[2,1]*logit_PCR))
points(subset_mid$PCR_Prev, subset_mid$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(PCR, expit(mid_Micro$coefficients[1,1] + mid_Micro$coefficients[2,1]*logit_LM), col = "blue")
points(subset_high$PCR_Prev, subset_high$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(PCR, expit(high_Micro$coefficients[1,1] + high_Micro$coefficients[2,1]*logit_PCR), col = "green")

#Plotting on the natural scale, lines for each section's bit only
PCR_low <- logit(seq(0,0.6,0.001))
PCR_mid <- logit(seq(0.2,0.85,0.001))
PCR_high <- logit(seq(0.62,1,0.001))

plot(subset_low$PCR_Prev, subset_low$Micro_Prev, xlim =c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(PCR_low), expit(low_Micro$coefficients[1,1] + low_Micro$coefficients[2,1]*PCR_low), lwd = 2)
points(subset_mid$PCR_Prev, subset_mid$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "blue")
lines(expit(PCR_mid), expit(mid_Micro$coefficients[1,1] + mid_Micro$coefficients[2,1]*PCR_mid), col = "blue", lwd = 2)
points(subset_high$PCR_Prev, subset_high$Micro_Prev, xlim =c(-10,10), ylim = c(-10,10), pch = 20, col = "green")
lines(expit(PCR_high), expit(high_Micro$coefficients[1,1] + high_Micro$coefficients[2,1]*PCR_high), col = "green", lwd = 2)

