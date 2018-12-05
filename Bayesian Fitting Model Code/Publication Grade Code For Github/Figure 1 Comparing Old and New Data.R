# Whittaker et al., 2019 : Variation in the Prevalence of Submicroscopic Malaria Infections: Historical Transmission 
#                          Intensity and Age as Key Determinants

# This paper represents an update on the systematic reviews of submicroscopic malaria infections  published
# by Okell et al in 2009 (https://academic.oup.com/jid/article/200/10/1509/879741) and by Okell et al in 2012
# (https://www.nature.com/articles/ncomms2241). In it, statistical analyses of both data from these previous 
# reviews and new data collected during the updating process are undertaken and results displayed graphically.

# The code below is responsible for the analyses and plotting that produced Figure 1 of the paper. Any questions, 
# queries, comments, or mistakes, please feel free to get in touch at charles.whittaker16@imperial.ac.uk :) 

# Analyses Responsible for Producing Figure 1 In the Paper

# Access and load various required packages for the analyses
library(rjags); library(moments); library(gtools); library(nortest); library(ggplot2); library(LMest); library(ssa); library(binom)

# Load in the dataset and subset the data into:
#        1) Data from the previous review, Okell et al., 2012 --> old_data
#        2) Data from this review, Whittaker et al., 2019 --> new_data
#        3) A combination of the data from both reviews --> full_data
data_frame <- Whittaker.et.al.R.Import
old_data <- data_frame[data_frame$Old_or_New == "Old" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
new_data <- data_frame[data_frame$Old_or_New == "New" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]

# Specifying the Data in the Format Required for Input into RJAGS
full_microscopy_PCR_comparison <- list(prev_pcr = full_data$PCR_N_Positive, ## number positive by PCR,
                                       prev_microscopy = full_data$Microscopy_N_Positive, ## number positive by microscopy,
                                       total_pcr = full_data$PCR_N_Tested, ## number tested by PCR,
                                       total_microscopy = full_data$Microscopy_N_Tested, ## number tested by microscopy
                                       N = 282) # Total sample data overall

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud", "gamma")
inits <- list(beta=0.0001,delt=0.0001,taud=0.0001)

# Specifying and initialising the RJAGS model- creates a JAGS model object
full_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # ALL DATA
                                       data = full_microscopy_PCR_comparison, n.chains = 1, inits = inits)

# Running, updating and iterating the RJAGS model
update(full_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, the burn in
full_micr_PCR_comp_model <- coda.samples(full_micr_PCR_comp_model, params, n.iter = 100000, thin = 70) # Model updating

# Exploring MCMC output
summary(full_micr_PCR_comp_model)
plot(full_micr_PCR_comp_model)
autocorr.plot(full_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
full_beta_mean <- mean(as.array(full_micr_PCR_comp_model[, 1])) 
full_delt_mean <- mean(as.array(full_micr_PCR_comp_model[, 2])) 

# Calculating the Predicted Values of Microscopy Prevalence Given a PCR Prevalence
  # Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_full <- seq(0, 0.96, 0.001)
full_logit_PCR_prevalence <- logit(PCR_prevalence_full)
full_fitted_logit_microscopy <- full_delt_mean + (1 + full_beta_mean)  * full_logit_PCR_prevalence
full_fitted_microscopy <- expit(full_fitted_logit_microscopy) # Converting logit to natural scale

# 95% Credible Interval Calculation for Prevalence and Sensitivity
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
full_data_sensitivity_credible_upper <- full_data_credible_upper / full_data_values
full_data_sensitivity_credible_lower <- full_data_credible_lower / full_data_values

# Figure 1A Plotting - Microscopy Prevalence Against PCR Prevalence - Data & Modelled Relationship
plot(new_data$PCR_Prev, new_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#C6D3F3", xlab = "PCR Prevalence", 
     ylab = "LM Prevalence", cex = 1.5)
points(old_data$PCR_Prev, old_data$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#6BC24E", cex = 1.5)
lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)
polygon(x = c(full_data_values, rev(full_data_values)), y = c(full_data_credible_upper, rev(full_data_credible_lower)), 
        col = adjustcolor("#F44A96", alpha.f = 0.5), border = NA)
lines(PCR_prevalence_full, full_fitted_microscopy, col = "#F44A96", lwd = 3)

legend("topleft", legend = c("Old Data", "New Data",  "Old and New Data"), col = c("#6BC24E", "#C6D3F3", "#F44A96"),
       pch = 20, pt.cex = 2, cex = 0.8)

# Figure 1B Plotting - Microscopy Sensitivity Against PCR Prevalence - Raw Data & Modelled Relationship
plot(full_data$PCR_Prev, full_data$Micro_Prev/full_data$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#BDBDBD",
     xlab = "PCR Prevalence", ylab = "Slide Prevalence", cex = 1.5)
lines(PCR_prevalence_full, full_fitted_microscopy/PCR_prevalence_full, col = "#F44A96", lwd = 3)
polygon(x = c(full_data_values, rev(full_data_values)), y = c(full_data_sensitivity_credible_upper, rev(full_data_sensitivity_credible_lower)), 
        col = adjustcolor("#F44A96", alpha.f = 0.5), border = NA)

# Figure 1C Plotting - Microscopy Sensitivity Against PCR Prevalence - Data Binned Into 9 Categories With Identical # of Studies
ordered_PCR_prevalences <- order(full_data$PCR_Prev)
full_data_bar_plotting <- full_data[ordered_PCR_prevalences, ]
size_of_group <- round(length(full_data_bar_plotting$Name)/9)
breaks <- full_data_bar_plotting$PCR_Prev[seq(1, length(full_data_bar_plotting$Name), size_of_group)]
breaks[1] <- -0.00001
breaks[10] <- 0.98
full_data_bar_plotting$Groups <- cut(x = full_data_bar_plotting$PCR_Prev, breaks = breaks)
full_data_bar_plotting$Sensitivity <- full_data_bar_plotting$Micro_Prev/full_data_bar_plotting$PCR_Prev
Sens_by_group <- tapply(full_data_bar_plotting$Sensitivity, full_data_bar_plotting$Groups, mean)
names_for_subsetting <- names(Sens_by_group)

lower_ci <- c()
upper_ci <- c()
for (i in 1:length(names_for_subsetting)) {
  sd <- sd(full_data_bar_plotting$Sensitivity[full_data_bar_plotting$Groups == names_for_subsetting[i]])
  stderr <- sd/sqrt(length(full_data_bar_plotting$Sensitivity[full_data_bar_plotting$Groups == names_for_subsetting[i]]))
  low_ci <- Sens_by_group[i] - 1.96 * stderr
  high_ci <- Sens_by_group[i] + 1.96 * stderr 
  lower_ci <- c(lower_ci, low_ci)
  upper_ci <- c(upper_ci, high_ci)
}

plot(0, 0, ylim = c(0, 1), xlim = c(0, 1), xlab = "PCR Prevalence (%)", ylab = "Microscopy Sensitivity", cex = 0)
for (i in 1:length(Sens_by_group)){
  rect(xleft = breaks[i], ybottom = 0, xright = breaks[i + 1], ytop = Sens_by_group[i], col = "#D3ECF3")
  arrows(x0 = (breaks[i] + breaks[i + 1])/2, y0 = lower_ci[i], 
         x1 = (breaks[i] + breaks[i + 1])/2, y1 = upper_ci[i], 
         col=1, angle=90, code=3, length = 0.05)
}
