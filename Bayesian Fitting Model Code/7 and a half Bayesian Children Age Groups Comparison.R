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
data_frame$Age_Group_2 <- as.factor(data_frame$Age_Group_2)

# Create individual datasets for each global region
young_children <- data_frame[data_frame$Age_Group_2 == "0-5", ]
old_children <- data_frame[data_frame$Age_Group_2 == "5-15", ]
adults <- data_frame[data_frame$Age_Group == "15+",]

# Without the zeroes
young_children <- data_frame[data_frame$Age_Group_2 == "0-5" & data_frame$Was_Initially_Zero. == "N", ]
old_children <- data_frame[data_frame$Age_Group_2 == "5-15" & data_frame$Was_Initially_Zero. == "N", ]
adults <- data_frame[data_frame$Age_Group == "15+" & data_frame$Was_Initially_Zero. == "N",]

## Specifying the data in the format required for input into RJAGS

#### NOTE THAT THERE ARE CURRENTLY DECIMALS IN THE DATA WHERE I USED LUCY'S TRICK- UNSURE HOW TO INCORPORATE
#### THESE CURRENTLY, MUST CHECK WITH HER. 

young_child_microscopy_PCR_comparison <- list(prev_pcr = young_children$PCR_N_Positive, ## number positive by PCR,
                                        prev_microscopy = young_children$Microscopy_N_Positive, ## number positive by microscopy,
                                        total_pcr = young_children$PCR_N_Tested, ## number tested by PCR,
                                        total_microscopy = young_children$Microscopy_N_Tested, ## number tested by microscopy,
                                        N = 24) #24) zeroes included                # Total sample data overall

old_child_microscopy_PCR_comparison <- list(prev_pcr = old_children$PCR_N_Positive, ## number positive by PCR,
                                              prev_microscopy = old_children$Microscopy_N_Positive, ## number positive by microscopy,
                                              total_pcr = old_children$PCR_N_Tested, ## number tested by PCR,
                                              total_microscopy = old_children$Microscopy_N_Tested, ## number tested by microscopy,
                                              N = 31) #31) zeroes included                # Total sample data overall

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
young_children_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                           data = young_child_microscopy_PCR_comparison, n.chains = 1, inits = inits)

old_children_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                                 data = old_child_microscopy_PCR_comparison, n.chains = 1, inits = inits)

adult_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # NEW DATA
                                        data = adult_microscopy_PCR_comparison, n.chains = 1, inits = inits)


# Running, updating and iterating the RJAGS model

# YOUNG CHILD DATA
update(young_children_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
young_children_micr_PCR_comp_model <- coda.samples(young_children_micr_PCR_comp_model, params, 
                                             n.iter = 35000, thin = 70) # Model updating
summary(young_children_micr_PCR_comp_model)
plot(young_children_micr_PCR_comp_model)
autocorr.plot(young_children_micr_PCR_comp_model)

young_child_beta_mean <- mean(as.array(young_children_micr_PCR_comp_model[, 1]))
young_child_delt_mean <- mean(as.array(young_children_micr_PCR_comp_model[, 2]))

# OLD CHILD DATA
update(old_children_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
old_children_micr_PCR_comp_model <- coda.samples(old_children_micr_PCR_comp_model, params, 
                                                   n.iter = 35000, thin = 70) # Model updating
summary(old_children_micr_PCR_comp_model)
plot(old_children_micr_PCR_comp_model)
autocorr.plot(old_children_micr_PCR_comp_model)

old_child_beta_mean <- mean(as.array(old_children_micr_PCR_comp_model[, 1]))
old_child_delt_mean <- mean(as.array(old_children_micr_PCR_comp_model[, 2]))

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
young_child_fitted_logit_microscopy <- young_child_delt_mean + (1 + young_child_beta_mean)  * logit_PCR_prevalence
old_child_fitted_logit_microscopy <- old_child_delt_mean + (1 + old_child_beta_mean)  * logit_PCR_prevalence
adult_fitted_logit_microscopy <- adult_delt_mean + (1 + adult_beta_mean) * logit_PCR_prevalence

# Converting them to the natural scale
young_child_fitted_microscopy <- expit(young_child_fitted_logit_microscopy)
old_child_fitted_microscopy <- expit(old_child_fitted_logit_microscopy)
adult_fitted_microscopy <- expit(adult_fitted_logit_microscopy)

# Plotting the results
plot(young_children$PCR_Prev, young_children$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#EAE600",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(old_children$PCR_Prev, old_children$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "orange")
points(adults$PCR_Prev, adults$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "red")
lines(seq(0,1,0.01), seq(0,1,0.01))
lines(PCR_prevalence, young_child_fitted_microscopy, col = "#EAE600", lwd = 3)
lines(PCR_prevalence, old_child_fitted_microscopy, col = "orange", lwd = 3)
lines(PCR_prevalence, adult_fitted_microscopy, col = "red", lwd = 3)

legend("topleft", 
       legend = c("Young Children", "Older Children", "Adults"), 
       col = c("#EAE600", "orange", "red"),
       lty = 1,
       lwd = 2,
       pt.cex = 2, 
       cex = 0.8)

##plotting using other colours
plot(young_children$PCR_Prev, young_children$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#F75C03",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(old_children$PCR_Prev, old_children$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#D80230")
points(adults$PCR_Prev, adults$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#820263")
lines(seq(0,1,0.01), seq(0,1,0.01))
lines(PCR_prevalence, young_child_fitted_microscopy, col = "#F75C03", lwd = 3)
lines(PCR_prevalence, old_child_fitted_microscopy, col = "#D80230", lwd = 3)
lines(PCR_prevalence, adult_fitted_microscopy, col = "#820263", lwd = 3)

legend("topleft", 
       legend = c("Young Children", "Older Children", "Adults"), 
       col = c("#F75C03", "#D80230", "#820263"),
       lty = 1,
       lwd = 2,
       pt.cex = 2, 
       cex = 0.8)

