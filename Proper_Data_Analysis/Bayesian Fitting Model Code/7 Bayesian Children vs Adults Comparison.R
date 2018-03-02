# Access relevant packages
library(rjags)
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

# Load in the dataset
data_frame <- Data.4th.February
data_frame <- Data.4th.February.Dec.Rounded

# Manipulate data_frame in prep for analysis
data_frame$Age_Group <- as.factor(data_frame$Age_Group)

# Create individual datasets for each global region
children <- data_frame[data_frame$Age_Group == "0-15", ]
adults <- data_frame[data_frame$Age_Group == "15+", ]

# Without the zeroes
children <- data_frame[data_frame$Age_Group == "0-15" & data_frame$Was_Initially_Zero. == "N", ]
adults <- data_frame[data_frame$Age_Group == "15+" & data_frame$Was_Initially_Zero. == "N", ]

## Specifying the data in the format required for input into RJAGS

#### NOTE THAT THERE ARE CURRENTLY DECIMALS IN THE DATA WHERE I USED LUCY'S TRICK- UNSURE HOW TO INCORPORATE
#### THESE CURRENTLY, MUST CHECK WITH HER. 

# Data Lucy provided me with
child_microscopy_PCR_comparison <- list(prev_pcr = children$PCR_N_Positive, ## number positive by PCR,
                                      prev_microscopy = children$Microscopy_N_Positive, ## number positive by microscopy,
                                      total_pcr = children$PCR_N_Tested, ## number tested by PCR,
                                      total_microscopy = children$Microscopy_N_Tested, ## number tested by microscopy,
                                      N = 143) #144) zeroes included                # Total sample data overall

# Data and results from the new systematic review
adult_microscopy_PCR_comparison <- list(prev_pcr = adults$PCR_N_Positive, ## number positive by PCR,
                                      prev_microscopy = adults$Microscopy_N_Positive,## number positive by microscopy,
                                      total_pcr = adults$PCR_N_Tested, ## number tested by PCR,
                                      total_microscopy = adults$Microscopy_N_Tested, ## number tested by microscopy,
                                      N = 15) #18) zeroes included



# Specifying the parameters of interest that RJAGS will track and output, and the 
# initial values to start the chain with
params <- c("delt", "beta", "taud")
inits <- list(beta=0.0001,delt=0.0001,taud=0.0001)

# Specifying and initialising the RJAGS model- creates a JAGS model object
children_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                      data = child_microscopy_PCR_comparison, n.chains = 1, inits = inits)

adult_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # NEW DATA
                                      data = adult_microscopy_PCR_comparison, n.chains = 1, inits = inits)


# Running, updating and iterating the RJAGS model

# CHILD DATA
update(children_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
children_micr_PCR_comp_model <- coda.samples(children_micr_PCR_comp_model, params, 
                                        n.iter = 35000, thin = 70) # Model updating
summary(children_micr_PCR_comp_model)
plot(children_micr_PCR_comp_model)
autocorr.plot(children_micr_PCR_comp_model)

child_beta_mean <- mean(as.array(children_micr_PCR_comp_model[, 1]))
child_delt_mean <- mean(as.array(children_micr_PCR_comp_model[, 2]))

# ADULT DATA
update(adult_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
adult_micr_PCR_comp_model <- coda.samples(adult_micr_PCR_comp_model, params, 
                                        n.iter = 35000, thin = 70) # Model updating
summary(adult_micr_PCR_comp_model)
plot(adult_micr_PCR_comp_model)
autocorr.plot(adult_micr_PCR_comp_model)

adult_beta_mean <- mean(as.array(adult_micr_PCR_comp_model[, 1]))
adult_delt_mean <- mean(as.array(adult_micr_PCR_comp_model[, 2]))

# Specifying PCR prevalence to plot against
PCR_prevalence <- seq(0.001,1,0.001)
logit_PCR_prevalence <- logit(PCR_prevalence)

# Calculating the fitted values specifying the best fit line on the logit scale
# DON'T FORGET IT'S delta' + (1 + BETA) * logit_PCR_prevalence
child_fitted_logit_microscopy <- child_delt_mean + (1 + child_beta_mean)  * logit_PCR_prevalence
adult_fitted_logit_microscopy <- adult_delt_mean + (1 + adult_beta_mean) * logit_PCR_prevalence

# Converting them to the natural scale
child_fitted_microscopy <- expit(child_fitted_logit_microscopy)
adult_fitted_microscopy <- expit(adult_fitted_logit_microscopy)

# Plotting the results
plot(children$PCR_Prev, children$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "orange",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(adults$PCR_Prev, adults$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "red")
lines(seq(0,1,0.01), seq(0,1,0.01))
lines(PCR_prevalence, child_fitted_microscopy, col = "orange", lwd = 3)
lines(PCR_prevalence, adult_fitted_microscopy, col = "red", lwd = 3)

legend("topleft", 
       legend = c("Children", "Adults"), 
       col = c("orange", "red"),
       lty = 1,
       lwd = 2,
       pt.cex = 2, 
       cex = 0.8)
