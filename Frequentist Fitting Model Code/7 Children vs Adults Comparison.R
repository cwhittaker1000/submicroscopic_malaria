# Access relevant packages
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

# Assign whichever dataset being used to data_frame
data_frame <- as.data.frame(Data.4th.February)

# Manipulate data_frame in prep for analysis
data_frame$Age_Group <- as.factor(data_frame$Age_Group)

# Create individual datasets for each global region
children <- data_frame[data_frame$Age_Group == "0-15", ]
adults <- data_frame[data_frame$Age_Group == "15+", ]

# Without the zeroes
children <- data_frame[data_frame$Age_Group == "0-15" & data_frame$Was_Initially_Zero. == "N", ]
adults <- data_frame[data_frame$Age_Group == "15+" & data_frame$Was_Initially_Zero. == "N", ]

# Running the linear regression, with PCR regressed onto Microscopy. 
children_lm <- summary(lm(logit(children$Micro_Prev) ~ logit(children$PCR_Prev)))

#Logit
adults_lm <- summary(lm(logit(adults$Micro_Prev) ~ logit(adults$PCR_Prev)))

#Linear
adults_lm_linear <- summary(lm(adults$Micro_Prev ~ 0 + adults$PCR_Prev))

# Defining logitPCR for plotting purposes
logit_PCR <- logit(seq(0,1,0.001))

# Logit scale plotting
plot(logit(children$PCR_Prev), logit(children$Micro_Prev))
lines(logit_PCR, children_lm$coefficients[1,1] + children_lm$coefficients[2,1]*logit_PCR, lwd = 2)

plot(logit(adults$PCR_Prev), logit(adults$Micro_Prev))
lines(logit_PCR, adults_lm$coefficients[1,1] + adults_lm$coefficients[2,1]*logit_PCR, lwd = 2)

# Natural scale plotting
plot(children$PCR_Prev, children$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(logit_PCR), expit(children_lm$coefficients[1,1] + children_lm$coefficients[2,1]*logit_PCR), col = "black", lwd = 2)

points(adults$PCR_Prev, adults$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "orange")
lines(expit(logit_PCR), expit(adults_lm$coefficients[1,1] + adults_lm$coefficients[2,1]*logit_PCR), col = "orange", lwd = 2)
lines(expit(logit_PCR), adults_lm_linear$coefficients[1,1]*expit(logit_PCR), col = "orange", lwd = 2)

