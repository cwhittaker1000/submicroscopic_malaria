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

#Manipulate data_frame in prep for analysis
data_frame$Old_or_New <- as.factor(data_frame$Old_or_New)
data_frame$Global_Region <- as.factor(data_frame$Global_Region)

#Create individual datasets for each global region
Asia <- data_frame[data_frame$Global_Region == "Asia",]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa",]
South_America <- data_frame[data_frame$Global_Region == "South America",]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa",]

# When removing the studies with zero values
Asia <- data_frame[data_frame$Global_Region == "Asia" & data_frame$Was_Initially_Zero. == "N", ]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa" & data_frame$Was_Initially_Zero. == "N", ]
South_America <- data_frame[data_frame$Global_Region == "South America" & data_frame$Was_Initially_Zero. == "N", ]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa" & data_frame$Was_Initially_Zero. == "N", ]

## Specifying the data in the format required for input into RJAGS

#### NOTE THAT THERE ARE CURRENTLY DECIMALS IN THE DATA WHERE I USED LUCY'S TRICK- UNSURE HOW TO INCORPORATE
#### THESE CURRENTLY, MUST CHECK WITH HER. 

# Asian data
Asia_microscopy_PCR_comparison <- list(prev_pcr = Asia$PCR_N_Positive, ## number positive by PCR,
                                      prev_microscopy = Asia$Microscopy_N_Positive, ## number positive by microscopy,
                                      total_pcr = Asia$PCR_N_Tested, ## number tested by PCR,
                                      total_microscopy = Asia$Microscopy_N_Tested, ## number tested by microscopy,
                                      N = 34) #39) zeroes included                # Total sample data overall

# East African data
East_Africa_microscopy_PCR_comparison <- list(prev_pcr = East_Africa$PCR_N_Positive, ## number positive by PCR,
                                      prev_microscopy = East_Africa$Microscopy_N_Positive,## number positive by microscopy,
                                      total_pcr = East_Africa$PCR_N_Tested, ## number tested by PCR,
                                      total_microscopy = East_Africa$Microscopy_N_Tested, ## number tested by microscopy,
                                      N = 156) #162) zeroes included

# South American data
South_America_microscopy_PCR_comparison <- list(prev_pcr = South_America$PCR_N_Positive, ## number positive by PCR,
                                       prev_microscopy = South_America$Microscopy_N_Positive, ## number positive by microscopy,
                                       total_pcr = South_America$PCR_N_Tested, ## number tested by PCR,
                                       total_microscopy = South_America$Microscopy_N_Tested, ## number tested by microscopy,
                                       N = 18) #30) zeroes included

# West African data
West_Africa_microscopy_PCR_comparison <- list(prev_pcr = West_Africa$PCR_N_Positive, ## number positive by PCR,
                                                prev_microscopy = West_Africa$Microscopy_N_Positive, ## number positive by microscopy,
                                                total_pcr = West_Africa$PCR_N_Tested, ## number tested by PCR,
                                                total_microscopy = West_Africa$Microscopy_N_Tested, ## number tested by microscopy,
                                                N = 77) #79) zeroes included

# Specifying the parameters of interest that RJAGS will track and output, and the 
# initial values to start the chain with
params <- c("delt", "beta", "taud")
inits <- list(beta=0.0001,delt=0.0001,taud=0.0001)

# Specifying and initialising the RJAGS model- creates a JAGS model object
Asia_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                      data = Asia_microscopy_PCR_comparison, n.chains = 1, inits = inits)

East_Africa_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # NEW DATA
                                      data = East_Africa_microscopy_PCR_comparison, n.chains = 1, inits = inits)

South_America_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # ALL DATA
                                       data = South_America_microscopy_PCR_comparison, n.chains = 1, inits = inits)

West_Africa_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # ALL DATA
                                                data = West_Africa_microscopy_PCR_comparison, n.chains = 1, inits = inits)


# Running, updating and iterating the RJAGS model

# ASIA DATA
update(Asia_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
Asia_micr_PCR_comp_model <- coda.samples(Asia_micr_PCR_comp_model, params, 
                                        n.iter = 35000, thin = 70) # Model updating
summary(Asia_micr_PCR_comp_model)
plot(Asia_micr_PCR_comp_model)
autocorr.plot(Asia_micr_PCR_comp_model)

Asia_beta_mean <- mean(as.array(Asia_micr_PCR_comp_model[, 1]))
Asia_delt_mean <- mean(as.array(Asia_micr_PCR_comp_model[, 2]))

# East Africa DATA
update(East_Africa_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
East_Africa_micr_PCR_comp_model <- coda.samples(East_Africa_micr_PCR_comp_model, params, 
                                        n.iter = 35000, thin = 70) # Model updating
summary(East_Africa_micr_PCR_comp_model)
plot(East_Africa_micr_PCR_comp_model)
autocorr.plot(East_Africa_micr_PCR_comp_model)

East_Africa_beta_mean <- mean(as.array(East_Africa_micr_PCR_comp_model[, 1]))
East_Africa_delt_mean <- mean(as.array(East_Africa_micr_PCR_comp_model[, 2]))

# South America DATA
update(South_America_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
South_America_micr_PCR_comp_model <- coda.samples(South_America_micr_PCR_comp_model, params, 
                                         n.iter = 35000, thin = 70) # Model updating
summary(South_America_micr_PCR_comp_model)
plot(South_America_micr_PCR_comp_model)
autocorr.plot(South_America_micr_PCR_comp_model)

South_America_beta_mean <- mean(as.array(South_America_micr_PCR_comp_model[, 1]))
South_America_delt_mean <- mean(as.array(South_America_micr_PCR_comp_model[, 2]))

# West Africa DATA
update(West_Africa_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
West_Africa_micr_PCR_comp_model <- coda.samples(West_Africa_micr_PCR_comp_model, params, 
                                                  n.iter = 35000, thin = 70) # Model updating
summary(West_Africa_micr_PCR_comp_model)
plot(West_Africa_micr_PCR_comp_model)
autocorr.plot(West_Africa_micr_PCR_comp_model)

West_Africa_beta_mean <- mean(as.array(West_Africa_micr_PCR_comp_model[, 1]))
West_Africa_delt_mean <- mean(as.array(West_Africa_micr_PCR_comp_model[, 2]))

# Specifying PCR prevalence to plot against
PCR_prevalence <- seq(0.001,1,0.001)
logit_PCR_prevalence <- logit(PCR_prevalence)

# Calculating the fitted values specifying the best fit line on the logit scale
# DON'T FORGET IT'S delta' + (1 + BETA) * logit_PCR_prevalence
Asia_fitted_logit_microscopy <- Asia_delt_mean + (1 + Asia_beta_mean)  * logit_PCR_prevalence
East_Africa_fitted_logit_microscopy <- East_Africa_delt_mean + (1 + East_Africa_beta_mean) * logit_PCR_prevalence
South_America_fitted_logit_microscopy <- South_America_delt_mean + (1 + South_America_beta_mean)  * logit_PCR_prevalence
West_Africa_fitted_logit_microscopy <- West_Africa_delt_mean + (1 + West_Africa_beta_mean) * logit_PCR_prevalence

# Converting them to the natural scale
Asia_fitted_microscopy <- expit(Asia_fitted_logit_microscopy)
East_Africa_fitted_microscopy <- expit(East_Africa_fitted_logit_microscopy)
South_America_fitted_microscopy <- expit(South_America_fitted_logit_microscopy)
West_Africa_fitted_microscopy <- expit(West_Africa_fitted_logit_microscopy)

# Plotting the results altoether
plot(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
points(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "green")
points(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "blue")
points(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "red")

lines(seq(0,1,0.01), seq(0,1,0.01))
lines(PCR_prevalence, Asia_fitted_microscopy, col = "black", lwd = 3)
lines(PCR_prevalence, East_Africa_fitted_microscopy, col = "green", lwd = 3)
lines(PCR_prevalence, South_America_fitted_microscopy, col = "blue", lwd = 3)
lines(PCR_prevalence, West_Africa_fitted_microscopy, col = "red", lwd = 3)

legend("topleft", 
       legend = c("Asia", "East Africa", "South America", "West Africa"), 
       col = c("black", "green", "blue", "red"),
       lty = 1,
       lwd = 2,
       pt.cex = 2, 
       cex = 0.8)

# Plotting the results separately

### ASIA

# Plotting the data for microscopy and PCR prevalence and the fitted line
plot(Asia$PCR_Prev, Asia$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, Asia_fitted_microscopy, col = "black", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 1)

# Plotting the confidence intervals for the data
Asia_Lower_Micro_Confint <- binom.confint(Asia$Microscopy_N_Positive, Asia$Microscopy_N_Tested, 
                                          conf.level = 0.95, methods = "exact")[, 5]
Asia_Upper_Micro_Confint <- binom.confint(Asia$Microscopy_N_Positive, Asia$Microscopy_N_Tested,
                                          conf.level = 0.95, methods = "exact")[, 6]
Asia_Lower_PCR_Confint <- binom.confint(Asia$PCR_N_Positive, Asia$PCR_N_Tested, 
                                        conf.level = 0.95, methods = "exact")[, 5]
Asia_Upper_PCR_Confint <- binom.confint(Asia$PCR_N_Positive, Asia$PCR_N_Tested, 
                                        conf.level = 0.95, methods = "exact")[, 6]

# Plotting the confidence intervals
arrows(Asia_Lower_PCR_Confint, Asia$Micro_Prev, Asia_Upper_PCR_Confint, 
       Asia$Micro_Prev, length=0, angle=90, code=1)
arrows(Asia$PCR_Prev, Asia_Lower_Micro_Confint, Asia$PCR_Prev, 
       Asia_Upper_Micro_Confint, length=0, angle=90, code=1)

# Plotting the sensitivity
plot(Asia$PCR_Prev, Asia$Micro_Prev/Asia$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, Asia_fitted_microscopy/PCR_prevalence, col = "black", lwd = 3)

### East Africa

# Plotting the data for microscopy and PCR prevalence and the fitted line
plot(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "green",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, East_Africa_fitted_microscopy, col = "green", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 1)

# Plotting the confidence intervals for the data
East_Africa_Lower_Micro_Confint <- binom.confint(East_Africa$Microscopy_N_Positive, East_Africa$Microscopy_N_Tested, 
                                          conf.level = 0.95, methods = "exact")[, 5]
East_Africa_Upper_Micro_Confint <- binom.confint(East_Africa$Microscopy_N_Positive, East_Africa$Microscopy_N_Tested,
                                          conf.level = 0.95, methods = "exact")[, 6]
East_Africa_Lower_PCR_Confint <- binom.confint(East_Africa$PCR_N_Positive, East_Africa$PCR_N_Tested, 
                                        conf.level = 0.95, methods = "exact")[, 5]
East_Africa_Upper_PCR_Confint <- binom.confint(East_Africa$PCR_N_Positive, East_Africa$PCR_N_Tested, 
                                        conf.level = 0.95, methods = "exact")[, 6]

# Plotting the confidence intervals
arrows(East_Africa_Lower_PCR_Confint, East_Africa$Micro_Prev, East_Africa_Upper_PCR_Confint, 
       East_Africa$Micro_Prev, length=0, angle=90, code=1)
arrows(East_Africa$PCR_Prev, East_Africa_Lower_Micro_Confint, East_Africa$PCR_Prev, 
       East_Africa_Upper_Micro_Confint, length=0, angle=90, code=1)

# Plotting the sensitivity
plot(East_Africa$PCR_Prev, East_Africa$Micro_Prev/East_Africa$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "green",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, East_Africa_fitted_microscopy/PCR_prevalence, col = "green", lwd = 3)


### South America

# Plotting the data for microscopy and PCR prevalence and the fitted line
plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "blue",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, South_America_fitted_microscopy, col = "blue", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 1)

# Plotting the confidence intervals for the data
South_America_Lower_Micro_Confint <- binom.confint(South_America$Microscopy_N_Positive, South_America$Microscopy_N_Tested, 
                                                 conf.level = 0.95, methods = "exact")[, 5]
South_America_Upper_Micro_Confint <- binom.confint(South_America$Microscopy_N_Positive, South_America$Microscopy_N_Tested,
                                                 conf.level = 0.95, methods = "exact")[, 6]
South_America_Lower_PCR_Confint <- binom.confint(South_America$PCR_N_Positive, South_America$PCR_N_Tested, 
                                               conf.level = 0.95, methods = "exact")[, 5]
South_America_Upper_PCR_Confint <- binom.confint(South_America$PCR_N_Positive, South_America$PCR_N_Tested, 
                                               conf.level = 0.95, methods = "exact")[, 6]

# Plotting the confidence intervals
arrows(South_America_Lower_PCR_Confint, South_America$Micro_Prev, South_America_Upper_PCR_Confint, 
       South_America$Micro_Prev, length=0, angle=90, code=1)
arrows(South_America$PCR_Prev, South_America_Lower_Micro_Confint, South_America$PCR_Prev, 
       South_America_Upper_Micro_Confint, length=0, angle=90, code=1)

# Plotting the sensitivity
plot(South_America$PCR_Prev, South_America$Micro_Prev/South_America$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "blue",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, South_America_fitted_microscopy/PCR_prevalence, col = "blue", lwd = 3)

### West Africa

# Plotting the data for microscopy and PCR prevalence and the fitted line
plot(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "red",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, West_Africa_fitted_microscopy, col = "red", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 1)

# Plotting the confidence intervals for the data
West_Africa_Lower_Micro_Confint <- binom.confint(West_Africa$Microscopy_N_Positive, West_Africa$Microscopy_N_Tested, 
                                                   conf.level = 0.95, methods = "exact")[, 5]
West_Africa_Upper_Micro_Confint <- binom.confint(West_Africa$Microscopy_N_Positive, West_Africa$Microscopy_N_Tested,
                                                   conf.level = 0.95, methods = "exact")[, 6]
West_Africa_Lower_PCR_Confint <- binom.confint(West_Africa$PCR_N_Positive, West_Africa$PCR_N_Tested, 
                                                 conf.level = 0.95, methods = "exact")[, 5]
West_Africa_Upper_PCR_Confint <- binom.confint(West_Africa$PCR_N_Positive, West_Africa$PCR_N_Tested, 
                                                 conf.level = 0.95, methods = "exact")[, 6]

# Plotting the confidence intervals
arrows(West_Africa_Lower_PCR_Confint, West_Africa$Micro_Prev, West_Africa_Upper_PCR_Confint, 
       West_Africa$Micro_Prev, length=0, angle=90, code=1)
arrows(West_Africa$PCR_Prev, West_Africa_Lower_Micro_Confint, West_Africa$PCR_Prev, 
       West_Africa_Upper_Micro_Confint, length=0, angle=90, code=1)

# Plotting the sensitivity
plot(West_Africa$PCR_Prev, West_Africa$Micro_Prev/West_Africa$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "red",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence, West_Africa_fitted_microscopy/PCR_prevalence, col = "red", lwd = 3)
