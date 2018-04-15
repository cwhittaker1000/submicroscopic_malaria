# Access relevant packages
library(rjags)
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

# Load in the dataset
# data_frame <- Data.4th.February
# data_frame <- Data.4th.February.Dec.Rounded
data_frame <- Zero.Tester.Data.4th.February.Dec.Rounded

# Manipulate data_frame in prep for analysis
data_frame$Old_or_New <- as.factor(data_frame$Old_or_New)
data_frame$Trans_Hist <- data_frame$Transmission_Setting_History
data_frame$Trans_Hist <- as.factor(data_frame$Trans_Hist)
data_frame$Trans_Now <- data_frame$Transmission_Setting_Current
data_frame$Trans_Now <- as.factor(data_frame$Trans_Now)

#Remove old age disaggregated data
data_frame <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 1 | data_frame$Full_Or_Age_Disagg_Data == 2, ]
# Subsetting data_frame, do if 0s NOT removed
data_frame <- data_frame[order(data_frame$Trans_Hist), ]
data_frame <- data_frame[1:213, ]

# OPTIONAL- REMOVE THOSE WHICH WERE INITIALLY ZERO NOT NECESSARY ANYMORE
data_frame <- data_frame[data_frame$Was_Initially_Zero. == "N", ]

#For when subsetting data with < 40% PCR Prevalence
data_frame <- data_frame[(data_frame$Full_Or_Age_Disagg_Data == 1 | data_frame$Full_Or_Age_Disagg_Data == 2)
                         & data_frame$PCR_Prev < 0.4, ]
data_frame <- data_frame[order(data_frame$Trans_Hist), ]
data_frame <- data_frame[1:139, ]

# Add columns for composite components
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

# Create individual datasets for each composite component
high_high_subset <- data_frame[data_frame$high_high == 1,]
high_low_subset <- data_frame[data_frame$high_low == 1,]
low_low_subset <- data_frame[data_frame$low_low == 1,]

## Specifying the data in the format required for input into RJAGS

#### NOTE THAT THERE ARE CURRENTLY DECIMALS IN THE DATA WHERE I USED LUCY'S TRICK- UNSURE HOW TO INCORPORATE
#### THESE CURRENTLY, MUST CHECK WITH HER. 

# High high data
high_high_microscopy_PCR_comparison <- list(prev_pcr = high_high_subset$PCR_N_Positive, ## number positive by PCR,
                                       prev_microscopy = high_high_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                       total_pcr = high_high_subset$PCR_N_Tested, ## number tested by PCR,
                                       total_microscopy = high_high_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                       N = 95) 
                                        #44 if not including > 0.4) 

# High low data
high_low_microscopy_PCR_comparison <- list(prev_pcr = high_low_subset$PCR_N_Positive, ## number positive by PCR,
                                              prev_microscopy = high_low_subset$Microscopy_N_Positive,## number positive by microscopy,
                                              total_pcr = high_low_subset$PCR_N_Tested, ## number tested by PCR,
                                              total_microscopy = high_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                              N = 88) 
                                                #69 if not including > 0.4) 

# Low low data
low_low_microscopy_PCR_comparison <- list(prev_pcr = low_low_subset$PCR_N_Positive, ## number positive by PCR,
                                                prev_microscopy = low_low_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                                total_pcr = low_low_subset$PCR_N_Tested, ## number tested by PCR,
                                                total_microscopy = low_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                                N = 24) 
                                                    #28 if not including > 0.4 

# Specifying the parameters of interest that RJAGS will track and output, and the 
# initial values to start the chain with
params <- c("delt", "beta", "taud")
inits <- list(beta=0.0001,delt=0.0001,taud=0.0001)

# Specifying and initialising the RJAGS model- creates a JAGS model object
High_high_micr_PCR_comp_model <- jags.model('1_Reversed_Basic_Model_Old_And_New.txt',   # OLD DATA
                                       data = high_high_microscopy_PCR_comparison, n.chains = 1, inits = inits)

High_low_micr_PCR_comp_model <- jags.model('1_Reversed_Basic_Model_Old_And_New.txt',   # NEW DATA
                                              data = high_low_microscopy_PCR_comparison, n.chains = 1, inits = inits)

Low_low_micr_PCR_comp_model <- jags.model('1_Reversed_Basic_Model_Old_And_New.txt',   # ALL DATA
                                                data = low_low_microscopy_PCR_comparison, n.chains = 1, inits = inits)

# Running, updating and iterating the RJAGS model

# High high DATA
update(High_high_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
High_high_micr_PCR_comp_model <- coda.samples(High_high_micr_PCR_comp_model, params, 
                                         n.iter = 35000, thin = 70) # Model updating
summary(High_high_micr_PCR_comp_model)
plot(High_high_micr_PCR_comp_model)
autocorr.plot(High_high_micr_PCR_comp_model)

High_high_beta_mean <- mean(as.array(High_high_micr_PCR_comp_model[, 1]))
High_high_delt_mean <- mean(as.array(High_high_micr_PCR_comp_model[, 2]))

# High low DATA
update(High_low_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
High_low_micr_PCR_comp_model <- coda.samples(High_low_micr_PCR_comp_model, params, 
                                                n.iter = 35000, thin = 70) # Model updating
summary(High_low_micr_PCR_comp_model)
plot(High_low_micr_PCR_comp_model)
autocorr.plot(High_low_micr_PCR_comp_model)

High_low_beta_mean <- mean(as.array(High_low_micr_PCR_comp_model[, 1]))
High_low_delt_mean <- mean(as.array(High_low_micr_PCR_comp_model[, 2]))

# Low low DATA
update(Low_low_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
Low_low_micr_PCR_comp_model <- coda.samples(Low_low_micr_PCR_comp_model, params, 
                                                  n.iter = 35000, thin = 70) # Model updating
summary(Low_low_micr_PCR_comp_model)
plot(Low_low_micr_PCR_comp_model)
autocorr.plot(Low_low_micr_PCR_comp_model)

Low_low_beta_mean <- mean(as.array(Low_low_micr_PCR_comp_model[, 1]))
Low_low_delt_mean <- mean(as.array(Low_low_micr_PCR_comp_model[, 2]))

# Specifying PCR prevalence to plot against
LM_prevalence_high_high <- seq(0,0.95,0.001)
logit_LM_prevalence_high_high <- logit(LM_prevalence_high_high)

LM_prevalence_high_low <- seq(0,0.9,0.001)
logit_LM_prevalence_high_low <- logit(LM_prevalence_high_low)

LM_prevalence_low_low <- seq(0,0.55,0.001)
logit_LM_prevalence_low_low <- logit(LM_prevalence_low_low)

# Calculating the fitted values specifying the best fit line on the logit scale
# DON'T FORGET IT'S delta' + (1 + BETA) * logit_PCR_prevalence
High_high_fitted_logit_PCR <- High_high_delt_mean + (1 + High_high_beta_mean)  * logit_LM_prevalence_high_high
High_low_fitted_logit_PCR <- High_low_delt_mean + (1 + High_low_beta_mean) * logit_LM_prevalence_high_low
Low_low_fitted_logit_PCR <- Low_low_delt_mean + (1 + Low_low_beta_mean)  * logit_LM_prevalence_low_low

# Converting them to the natural scale
High_high_fitted_PCR <- expit(High_high_fitted_logit_PCR)
High_low_fitted_PCR  <- expit(High_low_fitted_logit_PCR)
Low_low_fitted_PCR  <- expit(Low_low_fitted_logit_PCR)

# 95% Credible Intervals

# High_high
High_high_values <- seq(0, 0.95, 0.001)
High_high_chains <- High_high_micr_PCR_comp_model[[1]]
High_high_pred_mean_dist <- matrix(NA, nrow = nrow(High_high_chains), ncol = length(High_high_values))
for (i in 1:nrow(High_high_pred_mean_dist)){
  High_high_value_logit_scale <- High_high_chains[i, "delt"] + 
    (1 + High_high_chains[i, "beta"]) * logit(High_high_values)
  High_high_pred_mean_dist[i, ] <- expit(High_high_value_logit_scale)
}
High_high_credible_lower <- apply(High_high_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
High_high_credible_upper <- apply(High_high_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

# High_low
High_low_values <- seq(0, 0.9, 0.001)
High_low_chains <- High_low_micr_PCR_comp_model[[1]]
High_low_pred_mean_dist <- matrix(NA, nrow = nrow(High_low_chains), ncol = length(High_low_values))
for (i in 1:nrow(High_low_pred_mean_dist)){
  High_low_value_logit_scale <- High_low_chains[i, "delt"] + 
    (1 + High_low_chains[i, "beta"]) * logit(High_low_values)
  High_low_pred_mean_dist[i, ] <- expit(High_low_value_logit_scale)
}
High_low_credible_lower <- apply(High_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
High_low_credible_upper <- apply(High_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

# Low_low
Low_low_values <- seq(0, 0.55, 0.001)
Low_low_chains <- Low_low_micr_PCR_comp_model[[1]]
Low_low_pred_mean_dist <- matrix(NA, nrow = nrow(Low_low_chains), ncol = length(Low_low_values))
for (i in 1:nrow(Low_low_pred_mean_dist)){
  Low_low_value_logit_scale <- Low_low_chains[i, "delt"] + 
    (1 + Low_low_chains[i, "beta"]) * logit(Low_low_values)
  Low_low_pred_mean_dist[i, ] <- expit(Low_low_value_logit_scale)
}
Low_low_credible_lower <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
Low_low_credible_upper <- apply(Low_low_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

# Plotting the results
plot(high_high_subset$PCR_Prev, high_high_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#00A600FF",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(high_low_subset$PCR_Prev, high_low_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#ECB176FF")
points(low_low_subset$PCR_Prev, low_low_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "darkgrey")

lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)
lines(High_high_fitted_PCR, LM_prevalence_high_high,  col = "#00A600FF", lwd = 3)
lines(High_low_fitted_PCR, LM_prevalence_high_low,  col = "#ECB176FF", lwd = 3)
lines(Low_low_fitted_PCR, LM_prevalence_low_low,  col = "darkgrey", lwd = 3)

polygon(y = c(High_high_values, rev(High_high_values)), 
        x = c(High_high_credible_upper, rev(High_high_credible_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)

polygon(y = c(High_low_values, rev(High_low_values)), 
        x = c(High_low_credible_upper, rev(High_low_credible_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)

polygon(y = c(Low_low_values, rev(Low_low_values)), 
        x = c(Low_low_credible_upper, rev(Low_low_credible_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)

legend("topleft", 
       legend = c("High High", "High Low", "Low low"), 
       col = c("#00A600FF", "#ECB176FF", "darkgrey"),
       lty = 1,
       lwd = 2,
       pt.cex = 2, 
       cex = 0.8)

# Plotting the results separately

# High High Transmission Settings
plot(high_high_subset$PCR_Prev, high_high_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#00A600FF",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")

lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)
lines(PCR_prevalence_high_high, High_high_fitted_microscopy, col = "#00A600FF", lwd = 3)

polygon(x = c(PCR_prevalence_high_high, rev(PCR_prevalence_high_high)), 
        y = c(High_high_credible_upper, rev(High_high_credible_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)

# High Low Transmission Settings
plot(high_low_subset$PCR_Prev, high_low_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#ECB176FF",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")

lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)
lines(PCR_prevalence_high_low, High_low_fitted_microscopy, col = "#ECB176FF", lwd = 3)

polygon(x = c(PCR_prevalence_high_low, rev(PCR_prevalence_high_low)), 
        y = c(High_low_credible_upper, rev(High_low_credible_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)

# Low Low Transmission Settings
plot(low_low_subset$PCR_Prev, low_low_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "darkgrey",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")

lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)
lines(PCR_prevalence_low_low, Low_low_fitted_microscopy, col = "darkgrey", lwd = 3)

polygon(x = c(PCR_prevalence_low_low, rev(PCR_prevalence_low_low)), 
        y = c(Low_low_credible_upper, rev(Low_low_credible_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)


## Plotting the sensitivities

#High High Transmission Settings 

plot(high_high_subset$PCR_Prev, high_high_subset$Micro_Prev/high_high_subset$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#00A600FF",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(High_high_fitted_PCR, High_high_values/High_high_fitted_PCR, col = "#00A600FF", lwd = 3)

high_high_data_values <- seq(0, 0.95, 0.001)
high_high_data_chains_sens <- High_high_micr_PCR_comp_model[[1]]
high_high_data_pred_mean_dist_sens <- matrix(NA, nrow = nrow(high_high_data_chains_sens), ncol = length(high_high_data_values))
for (i in 1:nrow(high_high_data_pred_mean_dist_sens)){
  high_high_data_value_logit_scale <- high_high_data_chains_sens[i, "delt"] + 
    (1 + high_high_data_chains_sens[i, "beta"]) * logit(high_high_data_values)
  high_high_data_pred_mean_dist_sens[i, ] <- expit(high_high_data_value_logit_scale)
}
high_high_data_credible_lower_sens <- apply(high_high_data_pred_mean_dist_sens, MARGIN = 2, quantile, prob = 0.025)
high_high_data_credible_upper_sens <- apply(high_high_data_pred_mean_dist_sens, MARGIN = 2, quantile, prob = 0.975)

# Actually = lower because doing the 1 -
high_high_data_sensitivity_upper <- high_high_data_credible_upper_sens / high_high_data_values

# Actually = upper because doing the 1 -
high_high_data_sensitivity_lower <- high_high_data_credible_lower_sens / high_high_data_values

polygon(x = c(high_high_data_values, rev(high_high_data_values)), 
        y = c(high_high_data_sensitivity_upper, rev(high_high_data_sensitivity_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)

# High Low Transmission Settings

plot(high_low_subset$PCR_Prev, high_low_subset$Micro_Prev/high_low_subset$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#ECB176FF",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(High_low_fitted_PCR, High_low_values/High_low_fitted_PCR, col = "#ECB176FF", lwd = 3)

high_low_data_values <- seq(0, 0.9, 0.001)
high_low_data_chains_sens <- High_low_micr_PCR_comp_model[[1]]
high_low_data_pred_mean_dist_sens <- matrix(NA, nrow = nrow(high_low_data_chains_sens), ncol = length(high_low_data_values))
for (i in 1:nrow(high_low_data_pred_mean_dist_sens)){
  high_low_data_value_logit_scale <- high_low_data_chains_sens[i, "delt"] + 
    (1 + high_low_data_chains_sens[i, "beta"]) * logit(high_low_data_values)
  high_low_data_pred_mean_dist_sens[i, ] <- expit(high_low_data_value_logit_scale)
}
high_low_data_credible_lower_sens <- apply(high_low_data_pred_mean_dist_sens, MARGIN = 2, quantile, prob = 0.025)
high_low_data_credible_upper_sens <- apply(high_low_data_pred_mean_dist_sens, MARGIN = 2, quantile, prob = 0.975)

# Actually = lower because doing the 1 -
high_low_data_sensitivity_upper <- high_low_data_credible_upper_sens / high_low_data_values

# Actually = upper because doing the 1 -
high_low_data_sensitivity_lower <- high_low_data_credible_lower_sens / high_low_data_values

polygon(x = c(high_low_data_values, rev(high_low_data_values)), 
        y = c(high_low_data_sensitivity_upper, rev(high_low_data_sensitivity_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)

# Low Low Transmission Settings

plot(low_low_subset$PCR_Prev, low_low_subset$Micro_Prev/low_low_subset$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "darkgrey",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(Low_low_fitted_PCR, Low_low_values/Low_low_fitted_PCR, col = "darkgrey", lwd = 3)

low_low_data_values <- seq(0, 0.55, 0.001)
low_low_data_chains_sens <- Low_low_micr_PCR_comp_model[[1]]
low_low_data_pred_mean_dist_sens <- matrix(NA, nrow = nrow(low_low_data_chains_sens), ncol = length(low_low_data_values))
for (i in 1:nrow(low_low_data_pred_mean_dist_sens)){
  low_low_data_value_logit_scale <- low_low_data_chains_sens[i, "delt"] + 
    (1 + low_low_data_chains_sens[i, "beta"]) * logit(low_low_data_values)
  low_low_data_pred_mean_dist_sens[i, ] <- expit(low_low_data_value_logit_scale)
}
low_low_data_credible_lower_sens <- apply(low_low_data_pred_mean_dist_sens, MARGIN = 2, quantile, prob = 0.025)
low_low_data_credible_upper_sens <- apply(low_low_data_pred_mean_dist_sens, MARGIN = 2, quantile, prob = 0.975)

# Actually = lower because doing the 1 -
low_low_data_sensitivity_upper <- low_low_data_values/low_low_data_credible_upper_sens

# Actually = upper because doing the 1 -
low_low_data_sensitivity_lower <- low_low_data_values/low_low_data_credible_lower_sens 

polygon(x = c(low_low_data_values, rev(low_low_data_values)), 
        y = c(low_low_data_sensitivity_upper, rev(low_low_data_sensitivity_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)

plot(high_high_subset$PCR_Prev, high_high_subset$Micro_Prev/high_high_subset$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#00A600FF",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(High_high_fitted_PCR, High_high_values/High_high_fitted_PCR, col = "#00A600FF", lwd = 3)
points(high_low_subset$PCR_Prev, high_low_subset$Micro_Prev/high_low_subset$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#ECB176FF",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(High_low_fitted_PCR, High_low_values/High_low_fitted_PCR, col = "#ECB176FF", lwd = 3)
points(low_low_subset$PCR_Prev, low_low_subset$Micro_Prev/low_low_subset$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "darkgrey",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(Low_low_fitted_PCR, Low_low_values/Low_low_fitted_PCR, col = "darkgrey", lwd = 3)

low_low_sensitivity <- low_low_subset$Micro_Prev/low_low_subset$PCR_Prev
high_low_sensitivity <- high_low_subset$Micro_Prev/high_low_subset$PCR_Prev
high_high_sensitivity <- high_high_subset$Micro_Prev/high_high_subset$PCR_Prev

boxplot(low_low_sensitivity, high_low_sensitivity, high_high_sensitivity)
