# Access relevant packages
library(rjags)
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)
library(binom)

# Load in the dataset
data_frame <- Data.4th.February
data_frame <- Data.4th.February.Dec.Rounded

# Subset the data into old and new
old_data <- data_frame[data_frame$Old_or_New == "Old",]
new_data <- data_frame[data_frame$Old_or_New == "New",]
full_data <- data_frame

# When removing the studies with zero values
old_data <- data_frame[data_frame$Old_or_New == "Old" & data_frame$Was_Initially_Zero. == "N", ]
new_data <- data_frame[data_frame$Old_or_New == "New" & data_frame$Was_Initially_Zero. == "N", ]
full_data <- data_frame[data_frame$Was_Initially_Zero. == "N", ]

## Specifying the data in the format required for input into RJAGS

#### NOTE THAT THERE ARE CURRENTLY DECIMALS IN THE DATA WHERE I USED LUCY'S TRICK- UNSURE HOW TO INCORPORATE
#### THESE CURRENTLY, MUST CHECK WITH HER. 

# Data Lucy provided me with
old_microscopy_PCR_comparison <- list(prev_pcr = old_data$PCR_N_Positive, ## number positive by PCR,
                                  prev_microscopy = old_data$Microscopy_N_Positive, ## number positive by microscopy,
                                  total_pcr = old_data$PCR_N_Tested, ## number tested by PCR,
                                  total_microscopy = old_data$Microscopy_N_Tested, ## number tested by microscopy,
                                  N = 98) #106) zeroes included                # Total sample data overall

# Data and results from the new systematic review
new_microscopy_PCR_comparison <- list(prev_pcr = new_data$PCR_N_Positive, ## number positive by PCR,
                                        prev_microscopy = new_data$Microscopy_N_Positive,## number positive by microscopy,
                                        total_pcr = new_data$PCR_N_Tested, ## number tested by PCR,
                                        total_microscopy = new_data$Microscopy_N_Tested, ## number tested by microscopy,
                                      N = 217) #235) zeroes included

# Combined old and new data together
full_microscopy_PCR_comparison <- list(prev_pcr = full_data$PCR_N_Positive, ## number positive by PCR,
                                        prev_microscopy = full_data$Microscopy_N_Positive, ## number positive by microscopy,
                                        total_pcr = full_data$PCR_N_Tested, ## number tested by PCR,
                                        total_microscopy = full_data$Microscopy_N_Tested, ## number tested by microscopy,
                                        N = 315) #344) zeroes included    

# Specifying the parameters of interest that RJAGS will track and output, and the 
# initial values to start the chain with
params <- c("delt", "beta", "taud")
inits <- list(beta=0.0001,delt=0.0001,taud=0.0001)

# Specifying and initialising the RJAGS model- creates a JAGS model object
old_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                  data = old_microscopy_PCR_comparison, n.chains = 1, inits = inits)

new_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # NEW DATA
                                  data = new_microscopy_PCR_comparison, n.chains = 1, inits = inits)

full_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # ALL DATA
                                  data = full_microscopy_PCR_comparison, n.chains = 1, inits = inits)

# Running, updating and iterating the RJAGS model

# OLD DATA
update(old_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
old_micr_PCR_comp_model <- coda.samples(old_micr_PCR_comp_model, params, 
                                    n.iter = 35000, thin = 70) # Model updating
summary(old_micr_PCR_comp_model)
plot(old_micr_PCR_comp_model)
autocorr.plot(old_micr_PCR_comp_model)

old_beta_mean <- mean(as.array(old_micr_PCR_comp_model[, 1]))
old_delt_mean <- mean(as.array(old_micr_PCR_comp_model[, 2]))

# NEW DATA
update(new_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
new_micr_PCR_comp_model <- coda.samples(new_micr_PCR_comp_model, params, 
                                        n.iter = 35000, thin = 70) # Model updating
summary(new_micr_PCR_comp_model)
plot(new_micr_PCR_comp_model)
autocorr.plot(new_micr_PCR_comp_model)

new_beta_mean <- mean(as.array(new_micr_PCR_comp_model[, 1]))
new_delt_mean <- mean(as.array(new_micr_PCR_comp_model[, 2]))

# FULL DATA
update(full_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
full_micr_PCR_comp_model <- coda.samples(full_micr_PCR_comp_model, params, 
                                        n.iter = 35000, thin = 70) # Model updating
summary(full_micr_PCR_comp_model)
plot(full_micr_PCR_comp_model)
autocorr.plot(full_micr_PCR_comp_model)

full_beta_mean <- mean(as.array(full_micr_PCR_comp_model[, 1]))
full_delt_mean <- mean(as.array(full_micr_PCR_comp_model[, 2]))

# Specifying PCR prevalence to plot against
PCR_prevalence <- seq(0.001,1,0.001)
logit_PCR_prevalence <- logit(PCR_prevalence)

# Calculating the fitted values specifying the best fit line on the logit scale
# DON'T FORGET IT'S delta' + (1 + BETA) * logit_PCR_prevalence
old_fitted_logit_microscopy <- old_delt_mean + (1 + old_beta_mean)  * logit_PCR_prevalence
new_fitted_logit_microscopy <- new_delt_mean + (1 + new_beta_mean) * logit_PCR_prevalence
full_fitted_logit_microscopy <- full_delt_mean + (1 + full_beta_mean)  * logit_PCR_prevalence

# Converting them to the natural scale
old_fitted_microscopy <- expit(old_fitted_logit_microscopy)
new_fitted_microscopy <- expit(new_fitted_logit_microscopy)
full_fitted_microscopy <- expit(full_fitted_logit_microscopy)

# Plotting the results
plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "orange",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "blue")
lines(seq(0,1,0.01), seq(0,1,0.01))
lines(PCR_prevalence, new_fitted_microscopy, col = "orange", lwd = 3)
lines(PCR_prevalence, old_fitted_microscopy, col = "blue", lwd = 3)
lines(PCR_prevalence, full_fitted_microscopy, col = "black", lwd = 3)

legend("topleft", 
       legend = c("Old Data", "New Data",  "Old and New Data"), 
       col = c("blue", "orange", "black"),
       lty = 1,
       lwd = 2,
       pt.cex = 2, 
       cex = 0.8)

# Plotting the results separately

### OLD DATA

# Plotting the data for microscopy and PCR prevalence and the fitted line
plot(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, old_fitted_microscopy, col = "black", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 1)

# Plotting the confidence intervals for the data
old_data_Lower_Micro_Confint <- binom.confint(old_data$Microscopy_N_Positive, old_data$Microscopy_N_Tested, 
                                          conf.level = 0.95, methods = "exact")[, 5]
old_data_Upper_Micro_Confint <- binom.confint(old_data$Microscopy_N_Positive, old_data$Microscopy_N_Tested,
                                          conf.level = 0.95, methods = "exact")[, 6]
old_data_Lower_PCR_Confint <- binom.confint(old_data$PCR_N_Positive, old_data$PCR_N_Tested, 
                                        conf.level = 0.95, methods = "exact")[, 5]
old_data_Upper_PCR_Confint <- binom.confint(old_data$PCR_N_Positive, old_data$PCR_N_Tested, 
                                        conf.level = 0.95, methods = "exact")[, 6]

# Plotting the confidence intervals
arrows(old_data_Lower_PCR_Confint, old_data$Micro_Prev, old_data_Upper_PCR_Confint, 
       old_data$Micro_Prev, length=0, angle=90, code=1)
arrows(old_data$PCR_Prev, old_data_Lower_Micro_Confint, old_data$PCR_Prev, 
       old_data_Upper_Micro_Confint, length=0, angle=90, code=1)

# Plotting the sensitivity
plot(old_data$PCR_Prev, old_data$Micro_Prev/old_data$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, old_fitted_microscopy/PCR_prevalence, col = "black", lwd = 3)


### NEW DATA

# Plotting the data for microscopy and PCR prevalence and the fitted line
plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, new_fitted_microscopy, col = "black", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 1)

# Plotting the confidence intervals for the data
new_data_Lower_Micro_Confint <- binom.confint(new_data$Microscopy_N_Positive, new_data$Microscopy_N_Tested, 
                                              conf.level = 0.95, methods = "exact")[, 5]
new_data_Upper_Micro_Confint <- binom.confint(new_data$Microscopy_N_Positive, new_data$Microscopy_N_Tested,
                                              conf.level = 0.95, methods = "exact")[, 6]
new_data_Lower_PCR_Confint <- binom.confint(new_data$PCR_N_Positive, new_data$PCR_N_Tested, 
                                            conf.level = 0.95, methods = "exact")[, 5]
new_data_Upper_PCR_Confint <- binom.confint(new_data$PCR_N_Positive, new_data$PCR_N_Tested, 
                                            conf.level = 0.95, methods = "exact")[, 6]

# Plotting the confidence intervals
arrows(new_data_Lower_PCR_Confint, new_data$Micro_Prev, new_data_Upper_PCR_Confint, 
       new_data$Micro_Prev, length=0, angle=90, code=1)
arrows(new_data$PCR_Prev, new_data_Lower_Micro_Confint, new_data$PCR_Prev, 
       new_data_Upper_Micro_Confint, length=0, angle=90, code=1)

# Plotting the sensitivity
plot(new_data$PCR_Prev, new_data$Micro_Prev/new_data$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, new_fitted_microscopy/PCR_prevalence, col = "black", lwd = 3)

### FULL DATA

# Plotting the data for microscopy and PCR prevalence and the fitted line
plot(full_data$PCR_Prev, full_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, full_fitted_microscopy, col = "black", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 1)

# Plotting the confidence intervals for the data
full_data_Lower_Micro_Confint <- binom.confint(full_data$Microscopy_N_Positive, full_data$Microscopy_N_Tested, 
                                              conf.level = 0.95, methods = "exact")[, 5]
full_data_Upper_Micro_Confint <- binom.confint(full_data$Microscopy_N_Positive, full_data$Microscopy_N_Tested,
                                              conf.level = 0.95, methods = "exact")[, 6]
full_data_Lower_PCR_Confint <- binom.confint(full_data$PCR_N_Positive, full_data$PCR_N_Tested, 
                                            conf.level = 0.95, methods = "exact")[, 5]
full_data_Upper_PCR_Confint <- binom.confint(full_data$PCR_N_Positive, full_data$PCR_N_Tested, 
                                            conf.level = 0.95, methods = "exact")[, 6]

# Plotting the confidence intervals
arrows(full_data_Lower_PCR_Confint, full_data$Micro_Prev, full_data_Upper_PCR_Confint, 
       full_data$Micro_Prev, length=0, angle=90, code=1)
arrows(full_data$PCR_Prev, full_data_Lower_Micro_Confint, full_data$PCR_Prev, 
       full_data_Upper_Micro_Confint, length=0, angle=90, code=1)

# Plotting the sensitivity
plot(full_data$PCR_Prev, full_data$Micro_Prev/full_data$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, full_fitted_microscopy/PCR_prevalence, col = "black", lwd = 3)
