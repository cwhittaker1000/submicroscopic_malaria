# Whittaker et al., 2019 : Variation in the Prevalence of Submicroscopic Malaria Infections: Historical Transmission 
#                          Intensity and Age as Key Determinants

# This paper represents an update on the systematic reviews of submicroscopic malaria infections  published
# by Okell et al in 2009 (https://academic.oup.com/jid/article/200/10/1509/879741) and by Okell et al in 2012
# (https://www.nature.com/articles/ncomms2241). In it, statistical analyses of both data from these previous 
# reviews and new data collected during the updating process are undertaken and results displayed graphically.

# The code below is responsible for the analyses and plotting that produced Figure 1 of the paper. Any questions, 
# queries, comments, or mistakes, please feel free to get in touch at charles.whittaker16@imperial.ac.uk :) 

# Analyses Responsible for Producing Figure 3 In the Paper

# Access and load various required packages for the analyses
library(rjags); library(ssa); library(binom);

# Load in the dataset and subset the data by global region the survey was carried out in:
data_frame <- Submicroscopic_Review_Data_R_Import
Asia_Oceania <- data_frame[data_frame$Global_Region == "Asia&Oceania" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
South_America <- data_frame[data_frame$Global_Region == "South America" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa" & data_frame$Full_Or_Age_Disagg_Data == 2, ]

# Specifying the Data in the Format Required for Input into RJAGS
# Asian data
Asia_microscopy_PCR_comparison <- list(prev_pcr = Asia_Oceania$PCR_N_Positive, ## number positive by PCR,
                                       prev_microscopy = Asia_Oceania$Microscopy_N_Positive, ## number positive by microscopy,
                                       total_pcr = Asia_Oceania$PCR_N_Tested, ## number tested by PCR,
                                       total_microscopy = Asia_Oceania$Microscopy_N_Tested, ## number tested by microscopy,
                                       N = 71)               # Total sample data overall

# East African data
East_Africa_microscopy_PCR_comparison <- list(prev_pcr = East_Africa$PCR_N_Positive, ## number positive by PCR,
                                              prev_microscopy = East_Africa$Microscopy_N_Positive,## number positive by microscopy,
                                              total_pcr = East_Africa$PCR_N_Tested, ## number tested by PCR,
                                              total_microscopy = East_Africa$Microscopy_N_Tested, ## number tested by microscopy,
                                              N = 121)            # Total sample data overall

# South American data
South_America_microscopy_PCR_comparison <- list(prev_pcr = South_America$PCR_N_Positive, ## number positive by PCR,
                                                prev_microscopy = South_America$Microscopy_N_Positive, ## number positive by microscopy,
                                                total_pcr = South_America$PCR_N_Tested, ## number tested by PCR,
                                                total_microscopy = South_America$Microscopy_N_Tested, ## number tested by microscopy,
                                                N = 27) 

# West African data
West_Africa_microscopy_PCR_comparison <- list(prev_pcr = West_Africa$PCR_N_Positive, ## number positive by PCR,
                                              prev_microscopy = West_Africa$Microscopy_N_Positive, ## number positive by microscopy,
                                              total_pcr = West_Africa$PCR_N_Tested, ## number tested by PCR,
                                              total_microscopy = West_Africa$Microscopy_N_Tested, ## number tested by microscopy,
                                              N = 63) 

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud")
inits <- list(beta = 0.0001, delt = 0.0001, taud = 0.0001)

# Specifying and initialising the RJAGS model- creates a JAGS model object
Asia_micr_PCR_comp_model <- jags.model('Bayesian_Logit_Linear_Model.txt',   # OLD DATA
                                       data = Asia_microscopy_PCR_comparison, n.chains = 1, inits = inits)

East_Africa_micr_PCR_comp_model <- jags.model('Bayesian_Logit_Linear_Model.txt',   # NEW DATA
                                              data = East_Africa_microscopy_PCR_comparison, n.chains = 1, inits = inits)

South_America_micr_PCR_comp_model <- jags.model('Bayesian_Logit_Linear_Model.txt',   # ALL DATA
                                                data = South_America_microscopy_PCR_comparison, n.chains = 1, inits = inits)

West_Africa_micr_PCR_comp_model <- jags.model('Bayesian_Logit_Linear_Model.txt',   # ALL DATA
                                              data = West_Africa_microscopy_PCR_comparison, n.chains = 1, inits = inits)

# Running, updating and iterating the RJAGS model

# ASIA DATA
update(Asia_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
Asia_micr_PCR_comp_model <- coda.samples(Asia_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating

# Exploring MCMC output
summary(Asia_micr_PCR_comp_model)
plot(Asia_micr_PCR_comp_model)
autocorr.plot(Asia_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
Asia_beta_mean <- mean(as.array(Asia_micr_PCR_comp_model[, 1]))
Asia_delt_mean <- mean(as.array(Asia_micr_PCR_comp_model[, 2]))

# East Africa DATA
update(East_Africa_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
East_Africa_micr_PCR_comp_model <- coda.samples(East_Africa_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating

# Exploring MCMC output
summary(East_Africa_micr_PCR_comp_model)
plot(East_Africa_micr_PCR_comp_model)
autocorr.plot(East_Africa_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
East_Africa_beta_mean <- mean(as.array(East_Africa_micr_PCR_comp_model[, 1]))
East_Africa_delt_mean <- mean(as.array(East_Africa_micr_PCR_comp_model[, 2]))

# South America DATA
update(South_America_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
South_America_micr_PCR_comp_model <- coda.samples(South_America_micr_PCR_comp_model, params, n.iter = 200000, thin = 70) # Model updating

# Exploring MCMC output
summary(South_America_micr_PCR_comp_model)
plot(South_America_micr_PCR_comp_model)
autocorr.plot(South_America_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
South_America_beta_mean <- mean(as.array(South_America_micr_PCR_comp_model[, 1]))
South_America_delt_mean <- mean(as.array(South_America_micr_PCR_comp_model[, 2]))

# WEST AFRICA DATA
update(West_Africa_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
West_Africa_micr_PCR_comp_model <- coda.samples(West_Africa_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating

# Exploring MCMC output
summary(West_Africa_micr_PCR_comp_model)
plot(West_Africa_micr_PCR_comp_model)
autocorr.plot(West_Africa_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
West_Africa_beta_mean <- mean(as.array(West_Africa_micr_PCR_comp_model[, 1]))
West_Africa_delt_mean <- mean(as.array(West_Africa_micr_PCR_comp_model[, 2]))

# Calculating the Predicted Values of Microscopy Prevalence Given a PCR Prevalence
  # Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_Asia <- seq(0,0.8,0.001)
Asia_logit_PCR_prevalence <- logit(PCR_prevalence_Asia)
Asia_fitted_logit_microscopy <- Asia_delt_mean + (1 + Asia_beta_mean)  * Asia_logit_PCR_prevalence
Asia_fitted_microscopy <- expit(Asia_fitted_logit_microscopy)

PCR_prevalence_East_Africa <- seq(0,0.92,0.001)
East_Africa_logit_PCR_prevalence <- logit(PCR_prevalence_East_Africa)
East_Africa_fitted_logit_microscopy <- East_Africa_delt_mean + (1 + East_Africa_beta_mean) * East_Africa_logit_PCR_prevalence
East_Africa_fitted_microscopy <- expit(East_Africa_fitted_logit_microscopy)

PCR_prevalence_West_Africa <- seq(0.008,0.97,0.001)
West_Africa_logit_PCR_prevalence <- logit(PCR_prevalence_West_Africa)
West_Africa_fitted_logit_microscopy <- West_Africa_delt_mean + (1 + West_Africa_beta_mean) * West_Africa_logit_PCR_prevalence
West_Africa_fitted_microscopy <- expit(West_Africa_fitted_logit_microscopy)

PCR_prevalence_South_America <- seq(0,0.47,0.001)
South_America_logit_PCR_prevalence <- logit(PCR_prevalence_South_America)
South_America_fitted_logit_microscopy <- South_America_delt_mean + (1 + South_America_beta_mean)  * South_America_logit_PCR_prevalence
South_America_fitted_microscopy <- expit(South_America_fitted_logit_microscopy)

# 95% Credible Interval Calculation for Prevalence and Sensitivity
  # Asia
Asia_values <- seq(0, 0.8, 0.001)
Asia_chains <- Asia_micr_PCR_comp_model[[1]]
Asia_pred_mean_dist <- matrix(NA, nrow = nrow(Asia_chains), ncol = length(Asia_values))
for (i in 1:nrow(Asia_pred_mean_dist)){
  Asia_value_logit_scale <- Asia_chains[i, "delt"] + 
    (1 + Asia_chains[i, "beta"]) * logit(Asia_values)
  Asia_pred_mean_dist[i, ] <- expit(Asia_value_logit_scale)
}
Asia_credible_lower <- apply(Asia_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
Asia_credible_upper <- apply(Asia_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
Asia_sensitivity_upper <- Asia_credible_upper / Asia_values
Asia_sensitivity_lower <- Asia_credible_lower / Asia_values

  # East Africa
East_Africa_values <- seq(0, 0.92, 0.001)
East_Africa_chains <- East_Africa_micr_PCR_comp_model[[1]]
East_Africa_pred_mean_dist <- matrix(NA, nrow = nrow(East_Africa_chains), ncol = length(East_Africa_values))
for (i in 1:nrow(East_Africa_pred_mean_dist)){
  East_Africa_value_logit_scale <- East_Africa_chains[i, "delt"] + 
    (1 + East_Africa_chains[i, "beta"]) * logit(East_Africa_values)
  East_Africa_pred_mean_dist[i, ] <- expit(East_Africa_value_logit_scale)
}
East_Africa_credible_lower <- apply(East_Africa_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
East_Africa_credible_upper <- apply(East_Africa_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
East_Africa_sensitivity_upper <- East_Africa_credible_upper / East_Africa_values
East_Africa_sensitivity_lower <- East_Africa_credible_lower / East_Africa_values

  # South America
South_America_values <- seq(0, 0.47, 0.001)
South_America_chains <- South_America_micr_PCR_comp_model[[1]]
South_America_pred_mean_dist <- matrix(NA, nrow = nrow(South_America_chains), ncol = length(South_America_values))
for (i in 1:nrow(South_America_pred_mean_dist)){
  South_America_value_logit_scale <- South_America_chains[i, "delt"] + 
    (1 + South_America_chains[i, "beta"]) * logit(South_America_values)
  South_America_pred_mean_dist[i, ] <- expit(South_America_value_logit_scale)
}
South_America_credible_lower <- apply(South_America_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
South_America_credible_upper <- apply(South_America_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
South_America_sensitivity_upper <- South_America_credible_upper / South_America_values
South_America_sensitivity_lower <- South_America_credible_lower / South_America_values

  # West Africa
West_Africa_values <- seq(0.008, 0.97, 0.001)
West_Africa_chains <- West_Africa_micr_PCR_comp_model[[1]]
West_Africa_pred_mean_dist <- matrix(NA, nrow = nrow(West_Africa_chains), ncol = length(West_Africa_values))
for (i in 1:nrow(West_Africa_pred_mean_dist)){
  West_Africa_value_logit_scale <- West_Africa_chains[i, "delt"] + 
    (1 + West_Africa_chains[i, "beta"]) * logit(West_Africa_values)
  West_Africa_pred_mean_dist[i, ] <- expit(West_Africa_value_logit_scale)
}
West_Africa_credible_lower <- apply(West_Africa_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
West_Africa_credible_upper <- apply(West_Africa_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
West_Africa_sensitivity_upper <- West_Africa_credible_upper / West_Africa_values
West_Africa_sensitivity_lower <- West_Africa_credible_lower / West_Africa_values

# Figure 3A Plotting - Microscopy Prevalence Against PCR Prevalence for West Africa - Data & Modelled Relationship
plot(West_Africa$PCR_Prev, West_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "red",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence_West_Africa, West_Africa_fitted_microscopy, col = "red", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 2, lty = 2)
polygon(x = c(West_Africa_values, rev(West_Africa_values)), 
        y = c(West_Africa_credible_upper, rev(West_Africa_credible_lower)), 
        col = adjustcolor("red", alpha.f = 0.5), border = NA)

# Figure 3B Plotting - Microscopy Prevalence Against PCR Prevalence for East Africa - Data & Modelled Relationship
plot(East_Africa$PCR_Prev, East_Africa$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "green3",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence_East_Africa, East_Africa_fitted_microscopy, col = "green3", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 2, lty = 2)
polygon(x = c(East_Africa_values, rev(East_Africa_values)), 
        y = c(East_Africa_credible_upper, rev(East_Africa_credible_lower)), 
        col = adjustcolor("green3", alpha.f = 0.5), border = NA)

# Figure 3D Plotting - Microscopy Prevalence Against PCR Prevalence for South America - Data & Modelled Relationship
plot(South_America$PCR_Prev, South_America$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "blue",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence_South_America, South_America_fitted_microscopy, col = "blue", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 2, lty = 2)
polygon(x = c(South_America_values, rev(South_America_values)), 
        y = c(South_America_credible_upper, rev(South_America_credible_lower)), 
        col = adjustcolor("blue", alpha.f = 0.5), border = NA)

# Figure 3D Plotting - Microscopy Prevalence Against PCR Prevalence for Asia & Oceania - Data & Modelled Relationship
plot(Asia_Oceania$PCR_Prev, Asia_Oceania$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "black",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence_Asia, Asia_fitted_microscopy, col = "black", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 2, lty = 2)
polygon(x = c(Asia_values, rev(Asia_values)), 
        y = c(Asia_credible_upper, rev(Asia_credible_lower)), 
        col = adjustcolor("black", alpha.f = 0.5), border = NA)

# Figure 3E Plotting - Microscopy Sensitivity Against PCR Prevalence for Asia & Oceania - Modelled Relationship
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), pch = 1, xlab = "PCR Prevalence", ylab = "Sensitivity", cex = 0)
lines(PCR_prevalence_East_Africa, East_Africa_fitted_microscopy/PCR_prevalence_East_Africa, col = adjustcolor("green3"), lwd = 3)
lines(PCR_prevalence_West_Africa, West_Africa_fitted_microscopy/PCR_prevalence_West_Africa, col = adjustcolor("red"), lwd = 3)
lines(PCR_prevalence_Asia, Asia_fitted_microscopy/PCR_prevalence_Asia, col = adjustcolor("black"), lwd = 3)
lines(PCR_prevalence_South_America, South_America_fitted_microscopy/PCR_prevalence_South_America, col = adjustcolor("blue"), lwd = 3)
polygon(x = c(West_Africa_values, rev(West_Africa_values)), 
        y = c(West_Africa_sensitivity_upper, rev(West_Africa_sensitivity_lower)), 
        col = adjustcolor("red", alpha.f = 0.5), border = NA)
polygon(x = c(South_America_values, rev(South_America_values)), 
        y = c(South_America_sensitivity_upper, rev(South_America_sensitivity_lower)), 
        col = adjustcolor("blue", alpha.f = 0.5), border = NA)
polygon(x = c(East_Africa_values, rev(East_Africa_values)), 
        y = c(East_Africa_sensitivity_upper, rev(East_Africa_sensitivity_lower)), 
        col = adjustcolor("green3", alpha.f = 0.5), border = NA)
polygon(x = c(Asia_values, rev(Asia_values)), 
        y = c(Asia_sensitivity_upper, rev(Asia_sensitivity_lower)), 
        col = adjustcolor("black", alpha.f = 0.5), border = NA)

# Statistical Tests Carried Out On The Data
# ANOVA - Testing for Differences in Means
data_frame_ANOVA <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ] # all non age-disaggregated data
data_frame_ANOVA$Sensitivity <- (data_frame_ANOVA$Microscopy_N_Positive/data_frame_ANOVA$Microscopy_N_Tested)/
                                (data_frame_ANOVA$PCR_N_Positive/data_frame_ANOVA$PCR_N_Tested)
ANOVA_object <- aov(Sensitivity ~ Global_Region, data = data_frame_ANOVA)
summary(ANOVA_object)
TukeyHSD(ANOVA_object)


