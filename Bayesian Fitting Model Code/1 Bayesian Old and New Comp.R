# Whittaker et al., 2019 Submicroscopic Malaria Systematic Review
# Analyses Responsible for Producing Figure 1 In the Paper

# Access and load various required packages for the analyses
library(rjags)
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)
library(binom)

# Load in the dataset and subset the data into:
#        1) Data from the previous review, Okell et al., 2012 --> old_data
#        2) Data from this review, Whittaker et al., 2019 --> new_data
#        3) A combination of the data from both reviews --> full_data
data_frame <- Whittaker.et.al.R.Import
old_data <- data_frame[data_frame$Old_or_New == "Old" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
new_data <- data_frame[data_frame$Old_or_New == "New" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]

# Specifying the Data in the Format Required for Input into RJAGS

# Data from the old systematic review
old_microscopy_PCR_comparison <- list(prev_pcr = old_data$PCR_N_Positive, ## number positive by PCR,
                                  prev_microscopy = old_data$Microscopy_N_Positive, ## number positive by microscopy,
                                  total_pcr = old_data$PCR_N_Tested, ## number tested by PCR,
                                  total_microscopy = old_data$Microscopy_N_Tested, ## number tested by microscopy,
                                  N = 105) # Total sample data overall
  
# Data from the new systematic review
new_microscopy_PCR_comparison <- list(prev_pcr = new_data$PCR_N_Positive, ## number positive by PCR,
                                        prev_microscopy = new_data$Microscopy_N_Positive,## number positive by microscopy,
                                        total_pcr = new_data$PCR_N_Tested, ## number tested by PCR,
                                        total_microscopy = new_data$Microscopy_N_Tested, ## number tested by microscopy,
                                      N = 177) # Total sample data overall

# Combining old and new data together
full_microscopy_PCR_comparison <- list(prev_pcr = full_data$PCR_N_Positive, ## number positive by PCR,
                                        prev_microscopy = full_data$Microscopy_N_Positive, ## number positive by microscopy,
                                        total_pcr = full_data$PCR_N_Tested, ## number tested by PCR,
                                        total_microscopy = full_data$Microscopy_N_Tested, ## number tested by microscopy
                                        N = 282) # Total sample data overall

# Specifying the parameters of interest that RJAGS will track and output, and the 
# initial values to start the chain with
params <- c("delt", "beta", "taud", "gamma")
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
old_micr_PCR_comp_model <- coda.samples(old_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating

summary(old_micr_PCR_comp_model)
plot(old_micr_PCR_comp_model)
autocorr.plot(old_micr_PCR_comp_model)

old_beta_mean <- mean(as.array(old_micr_PCR_comp_model[, 1]))
old_delt_mean <- mean(as.array(old_micr_PCR_comp_model[, 2]))

# NEW DATA
update(new_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
new_micr_PCR_comp_model <- coda.samples(new_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating
summary(new_micr_PCR_comp_model)
plot(new_micr_PCR_comp_model)
autocorr.plot(new_micr_PCR_comp_model)

new_beta_mean <- mean(as.array(new_micr_PCR_comp_model[, 1]))
new_delt_mean <- mean(as.array(new_micr_PCR_comp_model[, 2]))

# FULL DATA
update(full_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
full_micr_PCR_comp_model <- coda.samples(full_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating
summary(full_micr_PCR_comp_model)
plot(full_micr_PCR_comp_model)
autocorr.plot(full_micr_PCR_comp_model)

full_beta_mean <- mean(as.array(full_micr_PCR_comp_model[, 1]))
full_delt_mean <- mean(as.array(full_micr_PCR_comp_model[, 2]))


# Specifying PCR prevalence to plot against
PCR_prevalence_old <- seq(0,1,0.001)
old_logit_PCR_prevalence <- logit(PCR_prevalence_old)

PCR_prevalence_new <- seq(0,1,0.001)
new_logit_PCR_prevalence <- logit(PCR_prevalence_new)

PCR_prevalence_full <- seq(0,0.96,0.001)
full_logit_PCR_prevalence <- logit(PCR_prevalence_full)

# Calculating the fitted values specifying the best fit line on the logit scale
# DON'T FORGET IT'S delta' + (1 + BETA) * logit_PCR_prevalence
old_fitted_logit_microscopy <- old_delt_mean + (1 + old_beta_mean)  * old_logit_PCR_prevalence
new_fitted_logit_microscopy <- new_delt_mean + (1 + new_beta_mean) * new_logit_PCR_prevalence
full_fitted_logit_microscopy <- full_delt_mean + (1 + full_beta_mean)  * full_logit_PCR_prevalence

plot(full_logit_PCR_prevalence, (1/(1 + full_beta_mean)) * full_fitted_logit_microscopy - (1/(1 + full_beta_mean)) * full_delt_mean)


# Converting them to the natural scale
old_fitted_microscopy <- expit(old_fitted_logit_microscopy)
new_fitted_microscopy <- expit(new_fitted_logit_microscopy)
full_fitted_microscopy <- expit(full_fitted_logit_microscopy)

# 95% credible interval plotting

# Old Data
old_data_values <- seq(0, 1, 0.001)
old_data_chains <- old_micr_PCR_comp_model[[1]]
old_data_pred_mean_dist <- matrix(NA, nrow = nrow(old_data_chains), ncol = length(old_data_values))
for (i in 1:nrow(old_data_pred_mean_dist)){
  old_data_value_logit_scale <- old_data_chains[i, "delt"] + 
    (1 + old_data_chains[i, "beta"]) * logit(old_data_values)
  old_data_pred_mean_dist[i, ] <- expit(old_data_value_logit_scale)
}
old_data_credible_lower <- apply(old_data_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
old_data_credible_upper <- apply(old_data_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)


# New Data
new_data_values <- seq(0, 1, 0.001)
new_data_chains <- new_micr_PCR_comp_model[[1]]
new_data_pred_mean_dist <- matrix(NA, nrow = nrow(new_data_chains), ncol = length(new_data_values))
for (i in 1:nrow(new_data_pred_mean_dist)){
  new_data_value_logit_scale <- new_data_chains[i, "delt"] + 
    (1 + new_data_chains[i, "beta"]) * logit(new_data_values)
  new_data_pred_mean_dist[i, ] <- expit(new_data_value_logit_scale)
}
new_data_credible_lower <- apply(new_data_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
new_data_credible_upper <- apply(new_data_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

# Full Data
full_data_values <- seq(0, 0.96, 0.001)
full_data_chains <- full_micr_PCR_comp_model[[1]]
full_data_pred_mean_dist <- matrix(NA, nrow = nrow(full_data_chains), ncol = length(full_data_values))
for (i in 1:nrow(full_data_pred_mean_dist)){
  full_data_value_logit_scale <- full_data_chains[i, "delt"] + 
    (1 + full_data_chains[i, "beta"]) * logit(full_data_values)
  full_data_pred_mean_dist[i, ] <- expit(full_data_value_logit_scale)
}
full_data_credible_lower <- apply(full_data_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
full_data_credible_upper <- apply(full_data_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)

# Plotting the results
plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#FF5900",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#019CD9")
lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)
lines(PCR_prevalence_new, new_fitted_microscopy, col = "#019CD9", lwd = 3)
lines(PCR_prevalence_old, old_fitted_microscopy, col = "#FF5900", lwd = 3)
lines(PCR_prevalence_full, full_fitted_microscopy, col = "#FEC93C", lwd = 3)

# Testing that backcalculation of what the PCR predicting equation would look like
LM_prevalence_tester <- seq(0,0.96,0.001)
PCR_Predicted_logit <- (1.075 + 0.973 * logit(LM_prevalence_tester))
PCR_Predicted_normal <- expit(PCR_Predicted_logit)

lines(PCR_Predicted_normal, LM_prevalence_tester, col = "black", lwd = 3)

polygon(x = c(old_data_values, rev(old_data_values)), 
        y = c(old_data_credible_upper, rev(old_data_credible_lower)), 
        col = adjustcolor("#04E824", alpha.f = 0.5), border = NA)
# optional lines: lines(old_data_values, old_data_credible_lower, col = "blue", lwd = 1)
# optional lines: lines(old_data_values, old_data_credible_upper, col = "blue", lwd = 1)

polygon(x = c(new_data_values, rev(new_data_values)), 
        y = c(new_data_credible_upper, rev(new_data_credible_lower)), 
        col = adjustcolor("#FF5900", alpha.f = 0.5), border = NA)
# optional lines: lines(new_data_values, new_data_credible_lower, col = "orange", lwd = 1)
# optional lines: lines(new_data_values, new_data_credible_upper, col = "orange", lwd = 1)

polygon(x = c(full_data_values, rev(full_data_values)), 
        y = c(full_data_credible_upper, rev(full_data_credible_lower)), 
        col = adjustcolor("#FEC93C", alpha.f = 0.5), border = NA)
# optional lines: lines(full_data_values, full_data_credible_lower, col = "black", lwd = 1)
# optional lines: lines(full_data_values, full_data_credible_upper, col = "black", lwd = 1)

legend("topleft", 
       legend = c("Old Data", "New Data",  "Old and New Data"), 
       col = c("#019CD9", "#FF5900", "#FEC93C"),
       pch = 20,
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
plot(full_data$PCR_Prev, full_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "grey",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence_full, full_fitted_microscopy, col = "black", lwd = 3)
lines(seq(0,1,0.01), seq(0,1,0.01), col = "black", lwd = 3)

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
       full_data$Micro_Prev, length=0, angle=90, code=1, col = "grey")
arrows(full_data$PCR_Prev, full_data_Lower_Micro_Confint, full_data$PCR_Prev, 
       full_data_Upper_Micro_Confint, length=0, angle=90, code=1, col = "grey")

# Plotting the sensitivities 

# Full data
plot(full_data$PCR_Prev, full_data$Micro_Prev/full_data$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "grey",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence")
lines(PCR_prevalence_full, full_fitted_microscopy/PCR_prevalence_full, col = "#383838", lwd = 3)

full_data_values <- seq(0, 0.96, 0.001)
full_data_chains_sens <- full_micr_PCR_comp_model[[1]]
full_data_pred_mean_dist_sens <- matrix(NA, nrow = nrow(full_data_chains_sens), ncol = length(full_data_values))
for (i in 1:nrow(full_data_pred_mean_dist_sens)){
  full_data_value_logit_scale <- full_data_chains_sens[i, "delt"] + 
    (1 + full_data_chains_sens[i, "beta"]) * logit(full_data_values)
  full_data_pred_mean_dist_sens[i, ] <- expit(full_data_value_logit_scale)
}
full_data_credible_lower_sens <- apply(full_data_pred_mean_dist_sens, MARGIN = 2, quantile, prob = 0.025)
full_data_credible_upper_sens <- apply(full_data_pred_mean_dist_sens, MARGIN = 2, quantile, prob = 0.975)

# Actually = lower because doing the 1 -
full_data_sensitivity_upper <- full_data_credible_upper_sens / full_data_values

# Actually = upper because doing the 1 -
full_data_sensitivity_lower <- full_data_credible_lower_sens / full_data_values

polygon(x = c(full_data_values, rev(full_data_values)), 
        y = c(full_data_sensitivity_upper, rev(full_data_sensitivity_lower)), 
        col = adjustcolor("#383838", alpha.f = 0.5), border = NA)


### PLOTTING THE PROPORTION SUBPATENT FOR EVERY 5% PREVALENCE

PCR_prevalence_full_barplot <- seq(0.1, 1, 0.2)
full_logit_PCR_prevalence_barplot <- logit(PCR_prevalence_full_barplot)
full_fitted_logit_microscopy_barplot <- full_delt_mean + (1 + full_beta_mean)  * full_logit_PCR_prevalence_barplot
full_fitted_microscopy_prevalence <- expit(full_fitted_logit_microscopy_barplot)

full_data_sensitivity <- full_fitted_microscopy_prevalence / PCR_prevalence_full_barplot
full_data_proportion_subpatent <- 1 - full_data_sensitivity

# for use with 10 bars
bloop <- barplot(full_data_proportion_subpatent, ylim = c(0,1), names.arg = c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60",
                                                           "60-70", "70-80", "80-90", "90-100"), 
        xlab = "PCR Prevalence", ylab = "Proportion of Population \n Submicroscopic", 
        col = adjustcolor("#383838", alpha.f = 0.9))

# for use with 5 bars
bloop <- barplot(full_data_proportion_subpatent, ylim = c(0,1), names.arg = c("0-20", "20-40", "40-60",
                                                                              "60-80", "80-100"),  
                 xlab = "PCR Prevalence", ylab = "Proportion of Population \n Submicroscopic", 
                 col = adjustcolor("#383838", alpha.f = 0.9))


# Upper and Lower Confidence Intervals

full_data_chains_prop <- full_micr_PCR_comp_model[[1]]
full_data_pred_mean_dist_prop <- matrix(NA, nrow = nrow(full_data_chains_prop), ncol = length(PCR_prevalence_full_barplot))
for (i in 1:nrow(full_data_pred_mean_dist_prop)){
  full_data_value_logit_scale <- full_data_chains_prop[i, "delt"] + 
    (1 + full_data_chains_prop[i, "beta"]) * logit(PCR_prevalence_full_barplot)
  full_data_pred_mean_dist_prop[i, ] <- expit(full_data_value_logit_scale)
}
full_data_credible_lower_prop <- apply(full_data_pred_mean_dist_prop, MARGIN = 2, quantile, prob = 0.025)
full_data_credible_upper_prop <- apply(full_data_pred_mean_dist_prop, MARGIN = 2, quantile, prob = 0.975)

# Actually = lower because doing the 1 -
full_data_sensitivity_upper <- full_data_credible_upper_prop / PCR_prevalence_full_barplot
full_data_proportion_subpatent_upper <- 1 - full_data_sensitivity_upper

# Actually = upper because doing the 1 -
full_data_sensitivity_lower <- full_data_credible_lower_prop / PCR_prevalence_full_barplot
full_data_proportion_subpatent_lower <- 1 - full_data_sensitivity_lower

# Plotting the confidence intervals
arrows(bloop, full_data_proportion_subpatent_upper, bloop, 
       full_data_proportion_subpatent_lower, length=0.1, angle=90, code=3, col = "black")



sens <- full_data$Micro_Prev/full_data$PCR_Prev
mean(sens)

sd <- sd(sens)
stderr <- sd/sqrt(length(sens))

mean(sens) - 1.96 * stderr
mean(sens) + 1.96 * stderr


hist(full_data$PCR_Prev, breaks = 20)

zero_five <- full_data[full_data$PCR_Prev < 0.05, ]
five_ten <- full_data[full_data$PCR_Prev >= 0.05 & full_data$PCR_Prev < 0.1, ]
ten_fifteen <- full_data[full_data$PCR_Prev >= 0.1 & full_data$PCR_Prev < 0.15, ]
fifteen_twenty <- full_data[full_data$PCR_Prev >= 0.15 & full_data$PCR_Prev < 0.2, ]
twnty_thirty <- full_data[full_data$PCR_Prev >= 0.2 & full_data$PCR_Prev < 0.30, ]
thirty_fourty <- full_data[full_data$PCR_Prev >= 0.3 & full_data$PCR_Prev < 0.40, ] 
fourty_sixty <- full_data[full_data$PCR_Prev >= 0.4 & full_data$PCR_Prev < 0.6, ]
sixty_hundred <- full_data[full_data$PCR_Prev >= 0.6 & full_data$PCR_Prev < 1, ]

eins <- mean(zero_five$Micro_Prev/zero_five$PCR_Prev)
zwei <- mean(five_ten$Micro_Prev/five_ten$PCR_Prev)
drei <- mean(ten_fifteen$Micro_Prev/ten_fifteen$PCR_Prev)
vier <- mean(fifteen_twenty$Micro_Prev/fifteen_twenty$PCR_Prev)
funf <- mean(twnty_thirty$Micro_Prev/twnty_thirty$PCR_Prev)
sechs <- mean(thirty_fourty$Micro_Prev/thirty_fourty$PCR_Prev)
sieben <- mean(fourty_sixty$Micro_Prev/fourty_sixty$PCR_Prev)
acht <- mean(sixty_hundred$Micro_Prev/sixty_hundred$PCR_Prev)

barplot(c(eins, zwei, drei, vier, funf, sechs, sieben, acht))


sens <- (full_data$Microscopy_N_Positive/full_data$Microscopy_N_Tested)/(full_data$PCR_N_Positive/full_data$PCR_N_Tested)

mean(sens) - 1.96 * sd(sens)/sqrt(length(sens))
mean(sens)
mean(sens) + 1.96 * sd(sens)/sqrt(length(sens))


# Calculating posterior densities for Old and New data parameter estimates and comparing them graphically
old_beta <- density(as.array(old_micr_PCR_comp_model[, 1]))
plot(old_beta, col = "blue", ylim = c(0, 13), xlim = c(-0.13, 0.24), main = "Beta")
polygon(old_beta, col = adjustcolor("blue", 0.4), border = "blue", lwd = 2)
new_beta <- density(as.array(new_micr_PCR_comp_model[, 1]))
polygon(new_beta, col = adjustcolor("orange", 0.4), border = "orange", lwd = 2)
full_beta <- density(as.array(full_micr_PCR_comp_model[, 1]))
polygon(full_beta, col =  adjustcolor("green", 0.4), border = "green", lwd = 2)
legend(-0.10, 10, legend=c("Old", "New", "Full"), col=c("blue", "orange", "green"), lty=1, cex=1, lwd = 2)

old_delt <- density(as.array(old_micr_PCR_comp_model[, 2]))
plot(old_delt, col = "blue", ylim = c(0, 8), xlim = c(-1.6, -0.8), main = "Delta")
polygon(old_delt, col = adjustcolor("blue", 0.4), border = "blue", lwd = 2)
new_delt <- density(as.array(new_micr_PCR_comp_model[, 2]))
polygon(new_delt, col = adjustcolor("orange", 0.4), border = "orange", lwd = 2)
full_delt <- density(as.array(full_micr_PCR_comp_model[, 2]))
polygon(full_delt, col =  adjustcolor("green", 0.4), border = "green", lwd = 2)
legend(-1.3, 10, legend=c("Old", "New", "Full"), col=c("blue", "orange", "green"), lty=1, cex=1, lwd = 2)

old_taud <- density(as.array(old_micr_PCR_comp_model[, 3]))
plot(old_taud, col = "blue", ylim = c(0, 2.5), xlim = c(0, 10), main = "tauda")
polygon(old_taud, col = adjustcolor("blue", 0.4), border = "blue", lwd = 2)
new_taud <- density(as.array(new_micr_PCR_comp_model[, 3]))
polygon(new_taud, col = adjustcolor("orange", 0.4), border = "orange", lwd = 2)
full_taud <- density(as.array(full_micr_PCR_comp_model[, 3]))
polygon(full_taud, col =  adjustcolor("green", 0.4), border = "green", lwd = 2)
legend(-1.3, 10, legend=c("Old", "New", "Full"), col=c("blue", "orange", "green"), lty=1, cex=1, lwd = 2)


