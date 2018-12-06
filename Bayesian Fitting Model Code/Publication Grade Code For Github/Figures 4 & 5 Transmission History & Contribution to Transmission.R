# Whittaker et al., 2019 : Variation in the Prevalence of Submicroscopic Malaria Infections: Historical Transmission 
#                          Intensity and Age as Key Determinants

# This paper represents an update on the systematic reviews of submicroscopic malaria infections  published
# by Okell et al in 2009 (https://academic.oup.com/jid/article/200/10/1509/879741) and by Okell et al in 2012
# (https://www.nature.com/articles/ncomms2241). In it, statistical analyses of both data from these previous 
# reviews and new data collected during the updating process are undertaken and results displayed graphically.

# The code below is responsible for the analyses and plotting that produced Figure 1 of the paper. Any questions, 
# queries, comments, or mistakes, please feel free to get in touch at charles.whittaker16@imperial.ac.uk :) 

# Analyses Responsible for Producing Figure 4 In the Paper

# Access and load various required packages for the analyses
library(rjags); library(moments); library(gtools); library(nortest); library(ggplot2); library(LMest); library(ssa); library(binom);
library(tidyverse)

# Load in the dataset and subset the data by the survey region's transmission history (African surveys only):
data_frame <- Whittaker.et.al.R.Import
data_frame$Trans_Hist <- data_frame$Transmission_Setting_History
data_frame$Trans_Hist <- as.factor(data_frame$Trans_Hist)
data_frame$Trans_Now <- data_frame$Transmission_Setting_Current
data_frame$Trans_Now <- as.factor(data_frame$Trans_Now)
data_frame <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ] #Remove age disaggregated data

data_frame <- data_frame[order(data_frame$Trans_Hist), ] # orders by transmission history status (puts NAs at end)
data_frame <- data_frame[1:164, ] # removes NAs 

# Add columns assigning survey to High High, High Low, or Low Low
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

high_high_subset <- data_frame[data_frame$high_high == 1, ]
high_low_subset <- data_frame[data_frame$high_low == 1, ]
low_low_subset <- data_frame[data_frame$low_low == 1, ]


# Specifying the Data in the Format Required for Input into RJAGS
# High High Data
high_high_microscopy_PCR_comparison <- list(prev_pcr = high_high_subset$PCR_N_Positive, ## number positive by PCR,
                                            prev_microscopy = high_high_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                            total_pcr = high_high_subset$PCR_N_Tested, ## number tested by PCR,
                                            total_microscopy = high_high_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                            N = 71) 

# High Low Data
high_low_microscopy_PCR_comparison <- list(prev_pcr = high_low_subset$PCR_N_Positive, ## number positive by PCR,
                                           prev_microscopy = high_low_subset$Microscopy_N_Positive,## number positive by microscopy,
                                           total_pcr = high_low_subset$PCR_N_Tested, ## number tested by PCR,
                                           total_microscopy = high_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                           N = 65) 
# Low Low Data
low_low_microscopy_PCR_comparison <- list(prev_pcr = low_low_subset$PCR_N_Positive, ## number positive by PCR,
                                          prev_microscopy = low_low_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                          total_pcr = low_low_subset$PCR_N_Tested, ## number tested by PCR,
                                          total_microscopy = low_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                          N = 28) 


# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud")
inits <- list(beta = 0.0001, delt = 0.0001, taud= 0.0001)

# Specifying and initialising the RJAGS model- creates a JAGS model object
High_high_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                            data = high_high_microscopy_PCR_comparison, n.chains = 1, inits = inits)

High_low_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # NEW DATA
                                           data = high_low_microscopy_PCR_comparison, n.chains = 1, inits = inits)

Low_low_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # ALL DATA
                                          data = low_low_microscopy_PCR_comparison, n.chains = 1, inits = inits)

# Running, updating and iterating the RJAGS model

# High High DATA
update(High_high_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
High_high_micr_PCR_comp_model <- coda.samples(High_high_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating

# Exploring MCMC output
summary(High_high_micr_PCR_comp_model)
plot(High_high_micr_PCR_comp_model)
autocorr.plot(High_high_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
High_high_beta_mean <- mean(as.array(High_high_micr_PCR_comp_model[, 1]))
High_high_delt_mean <- mean(as.array(High_high_micr_PCR_comp_model[, 2]))

# High low DATA
update(High_low_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
High_low_micr_PCR_comp_model <- coda.samples(High_low_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating

# Exploring MCMC output
summary(High_low_micr_PCR_comp_model)
plot(High_low_micr_PCR_comp_model)
autocorr.plot(High_low_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
High_low_beta_mean <- mean(as.array(High_low_micr_PCR_comp_model[, 1]))
High_low_delt_mean <- mean(as.array(High_low_micr_PCR_comp_model[, 2]))

# Low low DATA
update(Low_low_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
Low_low_micr_PCR_comp_model <- coda.samples(Low_low_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating

# Exploring MCMC output
summary(Low_low_micr_PCR_comp_model)
plot(Low_low_micr_PCR_comp_model)
autocorr.plot(Low_low_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
Low_low_beta_mean <- mean(as.array(Low_low_micr_PCR_comp_model[, 1]))
Low_low_delt_mean <- mean(as.array(Low_low_micr_PCR_comp_model[, 2]))

# Calculating the Predicted Values of Microscopy Prevalence Given a PCR Prevalence
  # Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_high_high <- seq(0.03,0.95,0.001)
logit_PCR_prevalence_high_high <- logit(PCR_prevalence_high_high)
High_high_fitted_logit_microscopy <- High_high_delt_mean + (1 + High_high_beta_mean)  * logit_PCR_prevalence_high_high
High_high_fitted_microscopy <- expit(High_high_fitted_logit_microscopy)

PCR_prevalence_high_low <- seq(0.004,0.9,0.001)
logit_PCR_prevalence_high_low <- logit(PCR_prevalence_high_low)
High_low_fitted_logit_microscopy <- High_low_delt_mean + (1 + High_low_beta_mean) * logit_PCR_prevalence_high_low
High_low_fitted_microscopy <- expit(High_low_fitted_logit_microscopy)

PCR_prevalence_low_low <- seq(0.01,0.55,0.001)
logit_PCR_prevalence_low_low <- logit(PCR_prevalence_low_low)
Low_low_fitted_logit_microscopy <- Low_low_delt_mean + (1 + Low_low_beta_mean)  * logit_PCR_prevalence_low_low
Low_low_fitted_microscopy <- expit(Low_low_fitted_logit_microscopy)

# 95% Credible Interval Calculation for Prevalence and Sensitivity
  # High_High
High_high_chains <- High_high_micr_PCR_comp_model[[1]]
High_high_pred_mean_dist <- matrix(NA, nrow = nrow(High_high_chains), ncol = length(PCR_prevalence_high_high))
for (i in 1:nrow(High_high_pred_mean_dist)){
  High_high_value_logit_scale <- High_high_chains[i, "delt"] + 
    (1 + High_high_chains[i, "beta"]) * logit(PCR_prevalence_high_high)
  High_high_pred_mean_dist[i, ] <- expit(High_high_value_logit_scale)
}
High_high_credible_lower <- apply(High_high_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
High_high_credible_upper <- apply(High_high_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
high_high_sensitivity_upper <- High_high_credible_upper / PCR_prevalence_high_high
high_high_sensitivity_lower <- High_high_credible_lower / PCR_prevalence_high_high

  # High_Low
High_low_chains <- High_low_micr_PCR_comp_model[[1]]
High_low_pred_mean_dist <- matrix(NA, nrow = nrow(High_low_chains), ncol = length(PCR_prevalence_high_low))
for (i in 1:nrow(High_low_pred_mean_dist)){
  High_low_value_logit_scale <- High_low_chains[i, "delt"] + 
    (1 + High_low_chains[i, "beta"]) * logit(PCR_prevalence_high_low)
  High_low_pred_mean_dist[i, ] <- expit(High_low_value_logit_scale)
}
High_low_credible_lower <- apply(High_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
High_low_credible_upper <- apply(High_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
high_low_sensitivity_upper <- High_low_credible_upper / PCR_prevalence_high_low
high_low_sensitivity_lower <- High_low_credible_lower / PCR_prevalence_high_low

  # Low_Low
Low_low_chains <- Low_low_micr_PCR_comp_model[[1]]
Low_low_pred_mean_dist <- matrix(NA, nrow = nrow(Low_low_chains), ncol = length(PCR_prevalence_low_low))
for (i in 1:nrow(Low_low_pred_mean_dist)){
  Low_low_value_logit_scale <- Low_low_chains[i, "delt"] + 
    (1 + Low_low_chains[i, "beta"]) * logit(PCR_prevalence_low_low)
  Low_low_pred_mean_dist[i, ] <- expit(Low_low_value_logit_scale)
}
Low_low_credible_lower <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
Low_low_credible_upper <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
low_low_sensitivity_upper <- Low_low_credible_upper / PCR_prevalence_low_low
low_low_sensitivity_lower <- Low_low_credible_lower / PCR_prevalence_low_low

# NOTE FIGURE 4A WAS MADE IN ADOBE ILLUSTRATOR AND SO NO CODE FO IT IS FEATURED HERE. 

# Figure 4B Plotting - Sensitivity Against PCR Prevalence for All 3 Transmission Settings - Data & Modelled Relationship
plot(high_high_subset$PCR_Prev, high_high_subset$Micro_Prev/high_high_subset$PCR_Prev, 
     xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#00A600FF", xlab = "PCR Prevalence", ylab = "Sensitivity")
lines(PCR_prevalence_high_high, High_high_fitted_microscopy/PCR_prevalence_high_high, col = "#00A600FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_high, rev(PCR_prevalence_high_high)), 
        y = c(high_high_sensitivity_upper, rev(high_high_sensitivity_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)

plot(high_low_subset$PCR_Prev, high_low_subset$Micro_Prev/high_low_subset$PCR_Prev, 
     xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#ECB176FF", xlab = "PCR Prevalence", ylab = "Sensitivity")
lines(PCR_prevalence_high_low, High_low_fitted_microscopy/PCR_prevalence_high_low, col = "#ECB176FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_low, rev(PCR_prevalence_high_low)), 
        y = c(high_low_sensitivity_upper, rev(high_low_sensitivity_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)

plot(low_low_subset$PCR_Prev, low_low_subset$Micro_Prev/low_low_subset$PCR_Prev, 
     xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "darkgrey", xlab = "PCR Prevalence", ylab = "Sensitivity")
lines(PCR_prevalence_low_low, Low_low_fitted_microscopy/PCR_prevalence_low_low, col = "darkgrey", lwd = 3)
polygon(x = c(PCR_prevalence_low_low, rev(PCR_prevalence_low_low)), 
        y = c(low_low_sensitivity_upper, rev(low_low_sensitivity_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)

# Figure 4C Plotting - Microscopy Prevalence Against PCR Prevalence for All 3 Transmission Settings - Data & Modelled Relationship
plot(high_high_subset$PCR_Prev, high_high_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#00A600FF",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(high_low_subset$PCR_Prev, high_low_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#ECB176FF")
points(low_low_subset$PCR_Prev, low_low_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "darkgrey")
lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)
lines(PCR_prevalence_high_high, High_high_fitted_microscopy, col = "#00A600FF", lwd = 3)
lines(PCR_prevalence_high_low, High_low_fitted_microscopy, col = "#ECB176FF", lwd = 3)
lines(PCR_prevalence_low_low, Low_low_fitted_microscopy, col = "darkgrey", lwd = 3)
polygon(x = c(PCR_prevalence_high_high, rev(PCR_prevalence_high_high)), y = c(High_high_credible_upper, rev(High_high_credible_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_high_low, rev(PCR_prevalence_high_low)), y = c(High_low_credible_upper, rev(High_low_credible_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA) 
polygon(x = c(PCR_prevalence_low_low, rev(PCR_prevalence_low_low)), y = c(Low_low_credible_upper, rev(Low_low_credible_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)
legend("topleft", legend = c("High High", "High Low", "Low low"), col = c("#00A600FF", "#ECB176FF", "darkgrey"),
       lty = 1, lwd = 2, pt.cex = 2, cex = 0.8)

# Figure 5A Plotting - Contribution to Transmission in High High Transmission Settings
# Range of High High Survey PCR Prevalences
HH_PCR_Prevalence <- seq(0.03,0.95,0.001) 
# Proportion of individuals with patent infections (equivalent to the sensitivity)
HH_Patent_Percentage <- (High_high_fitted_microscopy/PCR_prevalence_high_high) * 100 
# Proportion of individuals with subpatent infections
HH_Subpatent_Percentage <- (100 - HH_Patent_Percentage) 
# Calculating subpatent contribution to transmission
HH_Subpatent_Contribution20 <- (HH_Subpatent_Percentage) / ((20 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH_Subpatent_Contribution5 <- (HH_Subpatent_Percentage) / ((5 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH_Subpatent_Contribution2 <- (HH_Subpatent_Percentage) / ((2 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), cex = 0)
lines(HH_PCR_Prevalence, HH_Subpatent_Contribution20, xlim = c(0, 1), ylim = c(0, 1))
lines(HH_PCR_Prevalence, HH_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(HH_PCR_Prevalence, HH_Subpatent_Contribution2, xlim = c(0, 1), ylim = c(0, 1))

# Figure 5B Plotting - Contribution to Transmission in High Low Transmission Settings
# Range of High Low Survey PCR Prevalences
HL_PCR_Prevalence <- seq(0.004,0.9,0.001)
# Proportion of individuals with patent infections (equivalent to the sensitivity)
HL_Patent_Percentage <- (High_low_fitted_microscopy/PCR_prevalence_high_low) * 100 
# Proportion of individuals with subpatent infections
HL_Subpatent_Percentage <- (100 - HL_Patent_Percentage) 
# Calculating subpatent contribution to transmission
HL_Subpatent_Contribution20 <- (HL_Subpatent_Percentage) / ((20 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL_Subpatent_Contribution5 <- (HL_Subpatent_Percentage) / ((5 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL_Subpatent_Contribution2 <- (HL_Subpatent_Percentage) / ((2 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), cex = 0)
lines(HL_PCR_Prevalence, HL_Subpatent_Contribution20, xlim = c(0, 1), ylim = c(0, 1))
lines(HL_PCR_Prevalence, HL_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(HL_PCR_Prevalence, HL_Subpatent_Contribution2, xlim = c(0, 1), ylim = c(0, 1))

# Figure 5C Plotting - Contribution to Transmission in Low Low Transmission Settings
# Range of Low Low Survey PCR Prevalences
LL_PCR_Prevalence <- seq(0.01,0.55,0.001)
# Proportion of individuals with patent infections (equivalent to the sensitivity)
LL_Patent_Percentage <- (Low_low_fitted_microscopy/PCR_prevalence_low_low) * 100 
# Proportion of individuals with subpatent infections
LL_Subpatent_Percentage <- (100 - LL_Patent_Percentage) 
# Calculating subpatent contribution to transmission
LL_Subpatent_Contribution20 <- (LL_Subpatent_Percentage) / ((20 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL_Subpatent_Contribution5 <- (LL_Subpatent_Percentage) / ((5 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL_Subpatent_Contribution2 <- (LL_Subpatent_Percentage) / ((2 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), cex = 0)
lines(LL_PCR_Prevalence, LL_Subpatent_Contribution20, xlim = c(0, 1), ylim = c(0, 1))
lines(LL_PCR_Prevalence, LL_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(LL_PCR_Prevalence, LL_Subpatent_Contribution2, xlim = c(0, 1), ylim = c(0, 1))

