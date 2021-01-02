###################################################################################################
##                                                                                               ##
##  Whittaker et al., 2019 : Understanding the Drivers of Submicroscopic Malaria Infection:      ##
##                           Updated Insights from A Systematic Review of Population Surveys     ##
##                                                                                               ##
##    This paper represents an update on the systematic reviews of submicroscopic malaria        ##
##    infection published by:                                                                    ##
##        Okell et al in 2009 (https://academic.oup.com/jid/article/200/10/1509/879741)          ##
##        Okell et al in 2012 (https://www.nature.com/articles/ncomms2241).                      ##
##                                                                                               ##
##    This updated review involves analysis of both data from these previous reviews and new     ##
##    data collected during the updating process.                                                ##
##                                                                                               ##
##    The code below is responsible for the analyses and plotting that produced the following    ##
##    paper figures:                                                                             ##
##        Supplementary Figure 7: Transmission Archetype Stratification Sensitivity Analysis     ##
##                                                                                               ##
##    Any questions, queries, comments, or mistakes, please feel free to get in touch at:        ##
##        charles.whittaker16@imperial.ac.uk :)                                                  ##
##                                                                                               ##
###################################################################################################
library(rjags); library(ssa); library(binom); library(MALDIquant); library(formattable); 
library(tictoc); library(BayesianTools); library(R2jags); library(bayesmix)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/Sub-Patent Malarial Infections/Sub_Patent_Malaria_Analysis/")
source("Functions/Submicroscopic_Analysis_Functions.R")
seed <- 193

# Load in the dataset and subset the data by the survey region's transmission history (African surveys only):
data_frame <- read.csv("Data/SI_Systematic_Review_Results_R_Import.csv")
data_frame <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ] #Remove age disaggregated data
data_frame <- data_frame[order(data_frame$Transmission_Setting_15), ] # orders by transmission history status (puts NAs at end)
data_frame <- data_frame[1:229, ] # removes NAs which are surveys not conducted in Africa for which Trans_Hist data was not available 

##### 20% Threshold ####

# Subsetting the Data by Transmission Archetype
high_high_subset <- data_frame[data_frame$Transmission_Setting_20 == "High_High", ]
high_low_subset <- data_frame[data_frame$Transmission_Setting_20 == "High_Low", ]
low_low_subset <- data_frame[data_frame$Transmission_Setting_20 == "Low_Low", ]

# Specifying the Data in the Format Required for Input into RJAGS
# High High Data
high_high_microscopy_PCR_comparison <- list(prev_pcr = high_high_subset$PCR_N_Positive, ## number positive by PCR,
                                            prev_microscopy = high_high_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                            total_pcr = high_high_subset$PCR_N_Tested, ## number tested by PCR,
                                            total_microscopy = high_high_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                            N = length(high_high_subset$Name)) 

# High Low Data
high_low_microscopy_PCR_comparison <- list(prev_pcr = high_low_subset$PCR_N_Positive, ## number positive by PCR,
                                           prev_microscopy = high_low_subset$Microscopy_N_Positive,## number positive by microscopy,
                                           total_pcr = high_low_subset$PCR_N_Tested, ## number tested by PCR,
                                           total_microscopy = high_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                           N = length(high_low_subset$Name)) 
# Low Low Data
low_low_microscopy_PCR_comparison <- list(prev_pcr = low_low_subset$PCR_N_Positive, ## number positive by PCR,
                                          prev_microscopy = low_low_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                          total_pcr = low_low_subset$PCR_N_Tested, ## number tested by PCR,
                                          total_microscopy = low_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                          N = length(low_low_subset$Name))

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud")
standard_jags_inits <- function(){
  list("beta" = 0.0001, "delt" = 0.001, "taud" = 0.001)
}

# Specifying, updating and iterating the RJAGS model

# High High DATA
High_high_micr_PCR_comp_model <- jags.parallel(data = high_high_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
High_high_micr_PCR_comp_model <- as.mcmc(High_high_micr_PCR_comp_model)
High_high_overall_chain <- rbind(High_high_micr_PCR_comp_model[[1]], High_high_micr_PCR_comp_model[[2]], High_high_micr_PCR_comp_model[[3]], High_high_micr_PCR_comp_model[[4]])
High_high_beta_mean <- mean(as.array(High_high_overall_chain[, 1]))
High_high_delt_mean <- mean(as.array(High_high_overall_chain[, 2]))

# High low DATA
High_low_micr_PCR_comp_model <- jags.parallel(data = high_low_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
High_low_micr_PCR_comp_model <- as.mcmc(High_low_micr_PCR_comp_model)
High_low_overall_chain <- rbind(High_low_micr_PCR_comp_model[[1]], High_low_micr_PCR_comp_model[[2]], High_low_micr_PCR_comp_model[[3]], High_low_micr_PCR_comp_model[[4]])
High_low_beta_mean <- mean(as.array(High_low_overall_chain[, 1]))
High_low_delt_mean <- mean(as.array(High_low_overall_chain[, 2]))

# Low low DATA
Low_low_micr_PCR_comp_model <- jags.parallel(data = low_low_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
Low_low_micr_PCR_comp_model <- as.mcmc(Low_low_micr_PCR_comp_model)
Low_low_overall_chain <- rbind(Low_low_micr_PCR_comp_model[[1]], Low_low_micr_PCR_comp_model[[2]], Low_low_micr_PCR_comp_model[[3]], Low_low_micr_PCR_comp_model[[4]])
Low_low_beta_mean <- mean(as.array(Low_low_overall_chain[, 1]))
Low_low_delt_mean <- mean(as.array(Low_low_overall_chain[, 2]))

# Calculating the Predicted Values of Microscopy Prevalence Given a PCR Prevalence
# Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_high_high <- seq(0.03,0.95,0.001)
logit_PCR_prevalence_high_high <- logit(PCR_prevalence_high_high)
High_high_fitted_logit_microscopy <- High_high_delt_mean - (High_high_beta_mean * mean(logit_PCR_prevalence_high_high)) + ((1 + High_high_beta_mean)  * logit_PCR_prevalence_high_high)
High_high_fitted_microscopy <- expit(High_high_fitted_logit_microscopy)

PCR_prevalence_high_low <- seq(0.004,0.85,0.001)
logit_PCR_prevalence_high_low <- logit(PCR_prevalence_high_low)
High_low_fitted_logit_microscopy <- High_low_delt_mean - (High_low_beta_mean * mean(logit_PCR_prevalence_high_low)) + ((1 + High_low_beta_mean)  * logit_PCR_prevalence_high_low)
High_low_fitted_microscopy <- expit(High_low_fitted_logit_microscopy)

PCR_prevalence_low_low <- seq(0.01,0.55,0.001)
logit_PCR_prevalence_low_low <- logit(PCR_prevalence_low_low)
Low_low_fitted_logit_microscopy <- Low_low_delt_mean - (Low_low_beta_mean * mean(logit_PCR_prevalence_low_low)) + ((1 + Low_low_beta_mean)  * logit_PCR_prevalence_low_low)
Low_low_fitted_microscopy <- expit(Low_low_fitted_logit_microscopy)

# 95% Credible Interval Calculation for Prevalence and Sensitivity
# High_High
High_high_chains <- High_high_micr_PCR_comp_model[[1]]
High_high_pred_mean_dist <- matrix(NA, nrow = nrow(High_high_chains), ncol = length(PCR_prevalence_high_high))
for (i in 1:nrow(High_high_pred_mean_dist)){
  High_high_value_logit_scale <- High_high_chains[i, "delt"] - (High_high_chains[i, "beta"] * mean(logit_PCR_prevalence_high_high)) + ((1 + High_high_chains[i, "beta"])  * logit_PCR_prevalence_high_high)
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
  High_low_value_logit_scale <- High_low_chains[i, "delt"] - (High_low_chains[i, "beta"] * mean(logit_PCR_prevalence_high_low)) + ((1 + High_low_chains[i, "beta"])  * logit_PCR_prevalence_high_low)
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
  Low_low_value_logit_scale <- Low_low_chains[i, "delt"] - (Low_low_chains[i, "beta"] * mean(logit_PCR_prevalence_low_low)) + ((1 + Low_low_chains[i, "beta"])  * logit_PCR_prevalence_low_low)
  Low_low_pred_mean_dist[i, ] <- expit(Low_low_value_logit_scale)
}
Low_low_credible_lower <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
Low_low_credible_upper <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
low_low_sensitivity_upper <- Low_low_credible_upper / PCR_prevalence_low_low
low_low_sensitivity_lower <- Low_low_credible_lower / PCR_prevalence_low_low

# Panel 1 - Historically and Currently High Transmission Settings
par(mfrow = c(4, 3))
plot(high_high_subset$PCR_Prev * 100, high_high_subset$Micro_Prev/high_high_subset$PCR_Prev,  las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#00A600FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_high * 100, High_high_fitted_microscopy/PCR_prevalence_high_high, col = "#00A600FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_high * 100, rev(PCR_prevalence_high_high * 100)), 
        y = c(high_high_sensitivity_upper, rev(high_high_sensitivity_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)

# Panel 2 - Historically High But Currently Low Transmission Settings
plot(high_low_subset$PCR_Prev * 100, high_low_subset$Micro_Prev/high_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#ECB176FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_low * 100, High_low_fitted_microscopy/PCR_prevalence_high_low, col = "#ECB176FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_low * 100, rev(PCR_prevalence_high_low * 100)), 
        y = c(high_low_sensitivity_upper, rev(high_low_sensitivity_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)

# Panel 3 - Historically and Currently Low Transmission Settings
plot(low_low_subset$PCR_Prev * 100, low_low_subset$Micro_Prev/low_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "darkgrey", xlab = "PCR Prevalence", ylab = "Prevalence Ratio")
lines(PCR_prevalence_low_low * 100, Low_low_fitted_microscopy/PCR_prevalence_low_low, col = "darkgrey", lwd = 3)
polygon(x = c(PCR_prevalence_low_low * 100, rev(PCR_prevalence_low_low * 100)), 
        y = c(low_low_sensitivity_upper, rev(low_low_sensitivity_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)

##### 15% Threshold ####

# Subsetting the Data by Transmission Archetype
high_high_subset <- data_frame[data_frame$Transmission_Setting_15 == "High_High", ]
high_low_subset <- data_frame[data_frame$Transmission_Setting_15 == "High_Low", ]
low_low_subset <- data_frame[data_frame$Transmission_Setting_15 == "Low_Low", ]

# Specifying the Data in the Format Required for Input into RJAGS
# High High Data
high_high_microscopy_PCR_comparison <- list(prev_pcr = high_high_subset$PCR_N_Positive, ## number positive by PCR,
                                            prev_microscopy = high_high_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                            total_pcr = high_high_subset$PCR_N_Tested, ## number tested by PCR,
                                            total_microscopy = high_high_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                            N = length(high_high_subset$Name)) 

# High Low Data
high_low_microscopy_PCR_comparison <- list(prev_pcr = high_low_subset$PCR_N_Positive, ## number positive by PCR,
                                           prev_microscopy = high_low_subset$Microscopy_N_Positive,## number positive by microscopy,
                                           total_pcr = high_low_subset$PCR_N_Tested, ## number tested by PCR,
                                           total_microscopy = high_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                           N = length(high_low_subset$Name)) 
# Low Low Data
low_low_microscopy_PCR_comparison <- list(prev_pcr = low_low_subset$PCR_N_Positive, ## number positive by PCR,
                                          prev_microscopy = low_low_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                          total_pcr = low_low_subset$PCR_N_Tested, ## number tested by PCR,
                                          total_microscopy = low_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                          N = length(low_low_subset$Name))

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud")
standard_jags_inits <- function(){
  list("beta" = 0.0001, "delt" = 0.001, "taud" = 0.001)
}

# Specifying, updating and iterating the RJAGS model

# High High DATA
High_high_micr_PCR_comp_model <- jags.parallel(data = high_high_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
High_high_micr_PCR_comp_model <- as.mcmc(High_high_micr_PCR_comp_model)
High_high_overall_chain <- rbind(High_high_micr_PCR_comp_model[[1]], High_high_micr_PCR_comp_model[[2]], High_high_micr_PCR_comp_model[[3]], High_high_micr_PCR_comp_model[[4]])
High_high_beta_mean <- mean(as.array(High_high_overall_chain[, 1]))
High_high_delt_mean <- mean(as.array(High_high_overall_chain[, 2]))

# High low DATA
High_low_micr_PCR_comp_model <- jags.parallel(data = high_low_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
High_low_micr_PCR_comp_model <- as.mcmc(High_low_micr_PCR_comp_model)
High_low_overall_chain <- rbind(High_low_micr_PCR_comp_model[[1]], High_low_micr_PCR_comp_model[[2]], High_low_micr_PCR_comp_model[[3]], High_low_micr_PCR_comp_model[[4]])
High_low_beta_mean <- mean(as.array(High_low_overall_chain[, 1]))
High_low_delt_mean <- mean(as.array(High_low_overall_chain[, 2]))

# Low low DATA
Low_low_micr_PCR_comp_model <- jags.parallel(data = low_low_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
Low_low_micr_PCR_comp_model <- as.mcmc(Low_low_micr_PCR_comp_model)
Low_low_overall_chain <- rbind(Low_low_micr_PCR_comp_model[[1]], Low_low_micr_PCR_comp_model[[2]], Low_low_micr_PCR_comp_model[[3]], Low_low_micr_PCR_comp_model[[4]])
Low_low_beta_mean <- mean(as.array(Low_low_overall_chain[, 1]))
Low_low_delt_mean <- mean(as.array(Low_low_overall_chain[, 2]))

# Calculating the Predicted Values of Microscopy Prevalence Given a PCR Prevalence
# Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_high_high <- seq(0.03,0.95,0.001)
logit_PCR_prevalence_high_high <- logit(PCR_prevalence_high_high)
High_high_fitted_logit_microscopy <- High_high_delt_mean - (High_high_beta_mean * mean(logit_PCR_prevalence_high_high)) + ((1 + High_high_beta_mean)  * logit_PCR_prevalence_high_high)
High_high_fitted_microscopy <- expit(High_high_fitted_logit_microscopy)

PCR_prevalence_high_low <- seq(0.004,0.85,0.001)
logit_PCR_prevalence_high_low <- logit(PCR_prevalence_high_low)
High_low_fitted_logit_microscopy <- High_low_delt_mean - (High_low_beta_mean * mean(logit_PCR_prevalence_high_low)) + ((1 + High_low_beta_mean)  * logit_PCR_prevalence_high_low)
High_low_fitted_microscopy <- expit(High_low_fitted_logit_microscopy)

PCR_prevalence_low_low <- seq(0.01,0.55,0.001)
logit_PCR_prevalence_low_low <- logit(PCR_prevalence_low_low)
Low_low_fitted_logit_microscopy <- Low_low_delt_mean - (Low_low_beta_mean * mean(logit_PCR_prevalence_low_low)) + ((1 + Low_low_beta_mean)  * logit_PCR_prevalence_low_low)
Low_low_fitted_microscopy <- expit(Low_low_fitted_logit_microscopy)

# 95% Credible Interval Calculation for Prevalence and Sensitivity
# High_High
High_high_chains <- High_high_micr_PCR_comp_model[[1]]
High_high_pred_mean_dist <- matrix(NA, nrow = nrow(High_high_chains), ncol = length(PCR_prevalence_high_high))
for (i in 1:nrow(High_high_pred_mean_dist)){
  High_high_value_logit_scale <- High_high_chains[i, "delt"] - (High_high_chains[i, "beta"] * mean(logit_PCR_prevalence_high_high)) + ((1 + High_high_chains[i, "beta"])  * logit_PCR_prevalence_high_high)
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
  High_low_value_logit_scale <- High_low_chains[i, "delt"] - (High_low_chains[i, "beta"] * mean(logit_PCR_prevalence_high_low)) + ((1 + High_low_chains[i, "beta"])  * logit_PCR_prevalence_high_low)
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
  Low_low_value_logit_scale <- Low_low_chains[i, "delt"] - (Low_low_chains[i, "beta"] * mean(logit_PCR_prevalence_low_low)) + ((1 + Low_low_chains[i, "beta"])  * logit_PCR_prevalence_low_low)
  Low_low_pred_mean_dist[i, ] <- expit(Low_low_value_logit_scale)
}
Low_low_credible_lower <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
Low_low_credible_upper <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
low_low_sensitivity_upper <- Low_low_credible_upper / PCR_prevalence_low_low
low_low_sensitivity_lower <- Low_low_credible_lower / PCR_prevalence_low_low

# Panel 1 - Historically and Currently High Transmission Settings
plot(high_high_subset$PCR_Prev * 100, high_high_subset$Micro_Prev/high_high_subset$PCR_Prev,  las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#00A600FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_high * 100, High_high_fitted_microscopy/PCR_prevalence_high_high, col = "#00A600FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_high * 100, rev(PCR_prevalence_high_high * 100)), 
        y = c(high_high_sensitivity_upper, rev(high_high_sensitivity_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)

# Panel 2 - Historically High But Currently Low Transmission Settings
plot(high_low_subset$PCR_Prev * 100, high_low_subset$Micro_Prev/high_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#ECB176FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_low * 100, High_low_fitted_microscopy/PCR_prevalence_high_low, col = "#ECB176FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_low * 100, rev(PCR_prevalence_high_low * 100)), 
        y = c(high_low_sensitivity_upper, rev(high_low_sensitivity_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)

# Panel 3 - Historically and Currently Low Transmission Settings
plot(low_low_subset$PCR_Prev * 100, low_low_subset$Micro_Prev/low_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "darkgrey", xlab = "PCR Prevalence", ylab = "Prevalence Ratio")
lines(PCR_prevalence_low_low * 100, Low_low_fitted_microscopy/PCR_prevalence_low_low, col = "darkgrey", lwd = 3)
polygon(x = c(PCR_prevalence_low_low * 100, rev(PCR_prevalence_low_low * 100)), 
        y = c(low_low_sensitivity_upper, rev(low_low_sensitivity_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)

##### 10% Threshold ####

# Subsetting the Data by Transmission Archetype
high_high_subset <- data_frame[data_frame$Transmission_Setting_10 == "High_High", ]
high_low_subset <- data_frame[data_frame$Transmission_Setting_10 == "High_Low", ]
low_low_subset <- data_frame[data_frame$Transmission_Setting_10 == "Low_Low", ]

# Specifying the Data in the Format Required for Input into RJAGS
# High High Data
high_high_microscopy_PCR_comparison <- list(prev_pcr = high_high_subset$PCR_N_Positive, ## number positive by PCR,
                                            prev_microscopy = high_high_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                            total_pcr = high_high_subset$PCR_N_Tested, ## number tested by PCR,
                                            total_microscopy = high_high_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                            N = length(high_high_subset$Name)) 

# High Low Data
high_low_microscopy_PCR_comparison <- list(prev_pcr = high_low_subset$PCR_N_Positive, ## number positive by PCR,
                                           prev_microscopy = high_low_subset$Microscopy_N_Positive,## number positive by microscopy,
                                           total_pcr = high_low_subset$PCR_N_Tested, ## number tested by PCR,
                                           total_microscopy = high_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                           N = length(high_low_subset$Name)) 
# Low Low Data
low_low_microscopy_PCR_comparison <- list(prev_pcr = low_low_subset$PCR_N_Positive, ## number positive by PCR,
                                          prev_microscopy = low_low_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                          total_pcr = low_low_subset$PCR_N_Tested, ## number tested by PCR,
                                          total_microscopy = low_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                          N = length(low_low_subset$Name))

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud")
standard_jags_inits <- function(){
  list("beta" = 0.0001, "delt" = 0.001, "taud" = 0.001)
}

# Specifying, updating and iterating the RJAGS model

# High High DATA
High_high_micr_PCR_comp_model <- jags.parallel(data = high_high_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
High_high_micr_PCR_comp_model <- as.mcmc(High_high_micr_PCR_comp_model)
High_high_overall_chain <- rbind(High_high_micr_PCR_comp_model[[1]], High_high_micr_PCR_comp_model[[2]], High_high_micr_PCR_comp_model[[3]], High_high_micr_PCR_comp_model[[4]])
High_high_beta_mean <- mean(as.array(High_high_overall_chain[, 1]))
High_high_delt_mean <- mean(as.array(High_high_overall_chain[, 2]))

# High low DATA
High_low_micr_PCR_comp_model <- jags.parallel(data = high_low_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4,n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
High_low_micr_PCR_comp_model <- as.mcmc(High_low_micr_PCR_comp_model)
High_low_overall_chain <- rbind(High_low_micr_PCR_comp_model[[1]], High_low_micr_PCR_comp_model[[2]], High_low_micr_PCR_comp_model[[3]], High_low_micr_PCR_comp_model[[4]])
High_low_beta_mean <- mean(as.array(High_low_overall_chain[, 1]))
High_low_delt_mean <- mean(as.array(High_low_overall_chain[, 2]))

# Low low DATA
Low_low_micr_PCR_comp_model <- jags.parallel(data = low_low_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
Low_low_micr_PCR_comp_model <- as.mcmc(Low_low_micr_PCR_comp_model)
Low_low_overall_chain <- rbind(Low_low_micr_PCR_comp_model[[1]], Low_low_micr_PCR_comp_model[[2]], Low_low_micr_PCR_comp_model[[3]], Low_low_micr_PCR_comp_model[[4]])
Low_low_beta_mean <- mean(as.array(Low_low_overall_chain[, 1]))
Low_low_delt_mean <- mean(as.array(Low_low_overall_chain[, 2]))

# Calculating the Predicted Values of Microscopy Prevalence Given a PCR Prevalence
# Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_high_high <- seq(0.03,0.95,0.001)
logit_PCR_prevalence_high_high <- logit(PCR_prevalence_high_high)
High_high_fitted_logit_microscopy <- High_high_delt_mean - (High_high_beta_mean * mean(logit_PCR_prevalence_high_high)) + ((1 + High_high_beta_mean)  * logit_PCR_prevalence_high_high)
High_high_fitted_microscopy <- expit(High_high_fitted_logit_microscopy)

PCR_prevalence_high_low <- seq(0.004,0.85,0.001)
logit_PCR_prevalence_high_low <- logit(PCR_prevalence_high_low)
High_low_fitted_logit_microscopy <- High_low_delt_mean - (High_low_beta_mean * mean(logit_PCR_prevalence_high_low)) + ((1 + High_low_beta_mean)  * logit_PCR_prevalence_high_low)
High_low_fitted_microscopy <- expit(High_low_fitted_logit_microscopy)

PCR_prevalence_low_low <- seq(0.01,0.55,0.001)
logit_PCR_prevalence_low_low <- logit(PCR_prevalence_low_low)
Low_low_fitted_logit_microscopy <- Low_low_delt_mean - (Low_low_beta_mean * mean(logit_PCR_prevalence_low_low)) + ((1 + Low_low_beta_mean)  * logit_PCR_prevalence_low_low)
Low_low_fitted_microscopy <- expit(Low_low_fitted_logit_microscopy)

# 95% Credible Interval Calculation for Prevalence and Sensitivity
# High_High
High_high_chains <- High_high_micr_PCR_comp_model[[1]]
High_high_pred_mean_dist <- matrix(NA, nrow = nrow(High_high_chains), ncol = length(PCR_prevalence_high_high))
for (i in 1:nrow(High_high_pred_mean_dist)){
  High_high_value_logit_scale <- High_high_chains[i, "delt"] - (High_high_chains[i, "beta"] * mean(logit_PCR_prevalence_high_high)) + ((1 + High_high_chains[i, "beta"])  * logit_PCR_prevalence_high_high)
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
  High_low_value_logit_scale <- High_low_chains[i, "delt"] - (High_low_chains[i, "beta"] * mean(logit_PCR_prevalence_high_low)) + ((1 + High_low_chains[i, "beta"])  * logit_PCR_prevalence_high_low)
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
  Low_low_value_logit_scale <- Low_low_chains[i, "delt"] - (Low_low_chains[i, "beta"] * mean(logit_PCR_prevalence_low_low)) + ((1 + Low_low_chains[i, "beta"])  * logit_PCR_prevalence_low_low)
  Low_low_pred_mean_dist[i, ] <- expit(Low_low_value_logit_scale)
}
Low_low_credible_lower <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
Low_low_credible_upper <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
low_low_sensitivity_upper <- Low_low_credible_upper / PCR_prevalence_low_low
low_low_sensitivity_lower <- Low_low_credible_lower / PCR_prevalence_low_low

# Panel 1 - Historically and Currently High Transmission Settings
plot(high_high_subset$PCR_Prev * 100, high_high_subset$Micro_Prev/high_high_subset$PCR_Prev,  las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#00A600FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_high * 100, High_high_fitted_microscopy/PCR_prevalence_high_high, col = "#00A600FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_high * 100, rev(PCR_prevalence_high_high * 100)), 
        y = c(high_high_sensitivity_upper, rev(high_high_sensitivity_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)

# Panel 2 - Historically High But Currently Low Transmission Settings
plot(high_low_subset$PCR_Prev * 100, high_low_subset$Micro_Prev/high_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#ECB176FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_low * 100, High_low_fitted_microscopy/PCR_prevalence_high_low, col = "#ECB176FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_low * 100, rev(PCR_prevalence_high_low * 100)), 
        y = c(high_low_sensitivity_upper, rev(high_low_sensitivity_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)

# Panel 3 - Historically and Currently Low Transmission Settings
plot(low_low_subset$PCR_Prev * 100, low_low_subset$Micro_Prev/low_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "darkgrey", xlab = "PCR Prevalence", ylab = "Prevalence Ratio")
lines(PCR_prevalence_low_low * 100, Low_low_fitted_microscopy/PCR_prevalence_low_low, col = "darkgrey", lwd = 3)
polygon(x = c(PCR_prevalence_low_low * 100, rev(PCR_prevalence_low_low * 100)), 
        y = c(low_low_sensitivity_upper, rev(low_low_sensitivity_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)

##### 5% Threshold ####

# Subsetting the Data by Transmission Archetype
high_high_subset <- data_frame[data_frame$Transmission_Setting_5 == "High_High", ]
high_low_subset <- data_frame[data_frame$Transmission_Setting_5 == "High_Low", ]
low_low_subset <- data_frame[data_frame$Transmission_Setting_5 == "Low_Low", ]

# Specifying the Data in the Format Required for Input into RJAGS
# High High Data
high_high_microscopy_PCR_comparison <- list(prev_pcr = high_high_subset$PCR_N_Positive, ## number positive by PCR,
                                            prev_microscopy = high_high_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                            total_pcr = high_high_subset$PCR_N_Tested, ## number tested by PCR,
                                            total_microscopy = high_high_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                            N = length(high_high_subset$Name)) 

# High Low Data
high_low_microscopy_PCR_comparison <- list(prev_pcr = high_low_subset$PCR_N_Positive, ## number positive by PCR,
                                           prev_microscopy = high_low_subset$Microscopy_N_Positive,## number positive by microscopy,
                                           total_pcr = high_low_subset$PCR_N_Tested, ## number tested by PCR,
                                           total_microscopy = high_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                           N = length(high_low_subset$Name)) 
# Low Low Data
low_low_microscopy_PCR_comparison <- list(prev_pcr = low_low_subset$PCR_N_Positive, ## number positive by PCR,
                                          prev_microscopy = low_low_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                          total_pcr = low_low_subset$PCR_N_Tested, ## number tested by PCR,
                                          total_microscopy = low_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                          N = length(low_low_subset$Name))

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud")
standard_jags_inits <- function(){
  list("beta" = 0.0001, "delt" = 0.001, "taud" = 0.001)
}

# Specifying, updating and iterating the RJAGS model

# High High DATA
High_high_micr_PCR_comp_model <- jags.parallel(data = high_high_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
High_high_micr_PCR_comp_model <- as.mcmc(High_high_micr_PCR_comp_model)
High_high_overall_chain <- rbind(High_high_micr_PCR_comp_model[[1]], High_high_micr_PCR_comp_model[[2]], High_high_micr_PCR_comp_model[[3]], High_high_micr_PCR_comp_model[[4]])
High_high_beta_mean <- mean(as.array(High_high_overall_chain[, 1]))
High_high_delt_mean <- mean(as.array(High_high_overall_chain[, 2]))

# High low DATA
High_low_micr_PCR_comp_model <- jags.parallel(data = high_low_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
High_low_micr_PCR_comp_model <- as.mcmc(High_low_micr_PCR_comp_model)
High_low_overall_chain <- rbind(High_low_micr_PCR_comp_model[[1]], High_low_micr_PCR_comp_model[[2]], High_low_micr_PCR_comp_model[[3]], High_low_micr_PCR_comp_model[[4]])
High_low_beta_mean <- mean(as.array(High_low_overall_chain[, 1]))
High_low_delt_mean <- mean(as.array(High_low_overall_chain[, 2]))

# Low low DATA
Low_low_micr_PCR_comp_model <- jags.parallel(data = low_low_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = params, jags.seed = seed, model.file = "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 10000, n.burnin = 5000, n.thin = 25, DIC = TRUE)
Low_low_micr_PCR_comp_model <- as.mcmc(Low_low_micr_PCR_comp_model)
Low_low_overall_chain <- rbind(Low_low_micr_PCR_comp_model[[1]], Low_low_micr_PCR_comp_model[[2]], Low_low_micr_PCR_comp_model[[3]], Low_low_micr_PCR_comp_model[[4]])
Low_low_beta_mean <- mean(as.array(Low_low_overall_chain[, 1]))
Low_low_delt_mean <- mean(as.array(Low_low_overall_chain[, 2]))

# Calculating the Predicted Values of Microscopy Prevalence Given a PCR Prevalence
# Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_high_high <- seq(0.03,0.95,0.001)
logit_PCR_prevalence_high_high <- logit(PCR_prevalence_high_high)
High_high_fitted_logit_microscopy <- High_high_delt_mean - (High_high_beta_mean * mean(logit_PCR_prevalence_high_high)) + ((1 + High_high_beta_mean)  * logit_PCR_prevalence_high_high)
High_high_fitted_microscopy <- expit(High_high_fitted_logit_microscopy)

PCR_prevalence_high_low <- seq(0.004,0.85,0.001)
logit_PCR_prevalence_high_low <- logit(PCR_prevalence_high_low)
High_low_fitted_logit_microscopy <- High_low_delt_mean - (High_low_beta_mean * mean(logit_PCR_prevalence_high_low)) + ((1 + High_low_beta_mean)  * logit_PCR_prevalence_high_low)
High_low_fitted_microscopy <- expit(High_low_fitted_logit_microscopy)

PCR_prevalence_low_low <- seq(0.01,0.55,0.001)
logit_PCR_prevalence_low_low <- logit(PCR_prevalence_low_low)
Low_low_fitted_logit_microscopy <- Low_low_delt_mean - (Low_low_beta_mean * mean(logit_PCR_prevalence_low_low)) + ((1 + Low_low_beta_mean)  * logit_PCR_prevalence_low_low)
Low_low_fitted_microscopy <- expit(Low_low_fitted_logit_microscopy)

# 95% Credible Interval Calculation for Prevalence and Sensitivity
# High_High
High_high_chains <- High_high_micr_PCR_comp_model[[1]]
High_high_pred_mean_dist <- matrix(NA, nrow = nrow(High_high_chains), ncol = length(PCR_prevalence_high_high))
for (i in 1:nrow(High_high_pred_mean_dist)){
  High_high_value_logit_scale <- High_high_chains[i, "delt"] - (High_high_chains[i, "beta"] * mean(logit_PCR_prevalence_high_high)) + ((1 + High_high_chains[i, "beta"])  * logit_PCR_prevalence_high_high)
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
  High_low_value_logit_scale <- High_low_chains[i, "delt"] - (High_low_chains[i, "beta"] * mean(logit_PCR_prevalence_high_low)) + ((1 + High_low_chains[i, "beta"])  * logit_PCR_prevalence_high_low)
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
  Low_low_value_logit_scale <- Low_low_chains[i, "delt"] - (Low_low_chains[i, "beta"] * mean(logit_PCR_prevalence_low_low)) + ((1 + Low_low_chains[i, "beta"])  * logit_PCR_prevalence_low_low)
  Low_low_pred_mean_dist[i, ] <- expit(Low_low_value_logit_scale)
}
Low_low_credible_lower <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
Low_low_credible_upper <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
low_low_sensitivity_upper <- Low_low_credible_upper / PCR_prevalence_low_low
low_low_sensitivity_lower <- Low_low_credible_lower / PCR_prevalence_low_low

# Panel 1 - Historically and Currently High Transmission Settings
plot(high_high_subset$PCR_Prev * 100, high_high_subset$Micro_Prev/high_high_subset$PCR_Prev,  las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#00A600FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_high * 100, High_high_fitted_microscopy/PCR_prevalence_high_high, col = "#00A600FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_high * 100, rev(PCR_prevalence_high_high * 100)), 
        y = c(high_high_sensitivity_upper, rev(high_high_sensitivity_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)

# Panel 2 - Historically High But Currently Low Transmission Settings
plot(high_low_subset$PCR_Prev * 100, high_low_subset$Micro_Prev/high_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#ECB176FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_low * 100, High_low_fitted_microscopy/PCR_prevalence_high_low, col = "#ECB176FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_low * 100, rev(PCR_prevalence_high_low * 100)), 
        y = c(high_low_sensitivity_upper, rev(high_low_sensitivity_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)

# Panel 3 - Historically and Currently Low Transmission Settings
plot(low_low_subset$PCR_Prev * 100, low_low_subset$Micro_Prev/low_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "darkgrey", xlab = "PCR Prevalence", ylab = "Prevalence Ratio")
lines(PCR_prevalence_low_low * 100, Low_low_fitted_microscopy/PCR_prevalence_low_low, col = "darkgrey", lwd = 3)
polygon(x = c(PCR_prevalence_low_low * 100, rev(PCR_prevalence_low_low * 100)), 
        y = c(low_low_sensitivity_upper, rev(low_low_sensitivity_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)


# ANOVA and TukeyHSD
# Load in the dataset and subset the data by the survey region's transmission history (African surveys only):
data_frame <- read.csv("Data/SI_Systematic_Review_Results_R_Import.csv")
data_frame <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ] #Remove age disaggregated data
data_frame <- data_frame[order(data_frame$Transmission_Setting_15), ] # orders by transmission history status (puts NAs at end)
data_frame <- data_frame[1:229, ]
data_frame$Sensitivity <- data_frame$Micro_Prev/data_frame$PCR_Prev
variance <- data_frame$PCR_N_Tested * data_frame$PCR_Prev * (1 - data_frame$PCR_Prev)

weighted_20_model <- lm(Sensitivity ~ Transmission_Setting_20, data = data_frame, na.action = na.omit, weights = 1/variance)
summary(weighted_20_model)
summary(aov(weighted_20_model))
TukeyHSD(aov(weighted_20_model))

weighted_15_model <- lm(Sensitivity ~ Transmission_Setting_15, data = data_frame, na.action = na.omit, weights = 1/variance)
summary(weighted_15_model)
summary(aov(weighted_15_model))
TukeyHSD(aov(weighted_15_model))

weighted_10_model <- lm(Sensitivity ~ Transmission_Setting_10, data = data_frame, na.action = na.omit, weights = 1/variance)
summary(weighted_10_model)
summary(aov(weighted_10_model))
TukeyHSD(aov(weighted_10_model))

weighted_5_model <- lm(Sensitivity ~ Transmission_Setting_5, data = data_frame, na.action = na.omit, weights = 1/variance)
summary(weighted_5_model)
summary(aov(weighted_5_model))
TukeyHSD(aov(weighted_5_model))

# Including PCR Prevalence 
weighted_20_model <- lm(Sensitivity ~ PCR_Prev + Transmission_Setting_20, data = data_frame, na.action = na.omit, weights = 1/variance)
summary(weighted_20_model)
summary(aov(weighted_20_model))

weighted_15_model <- lm(Sensitivity ~ PCR_Prev + Transmission_Setting_15, data = data_frame, na.action = na.omit, weights = 1/variance)
summary(weighted_15_model)
summary(aov(weighted_15_model))

weighted_10_model <- lm(Sensitivity ~ PCR_Prev + Transmission_Setting_10, data = data_frame, na.action = na.omit, weights = 1/variance)
summary(weighted_10_model)
summary(aov(weighted_10_model))

weighted_5_model <- lm(Sensitivity ~ PCR_Prev + Transmission_Setting_5, data = data_frame, na.action = na.omit, weights = 1/variance)
summary(weighted_5_model)
summary(aov(weighted_10_model))

