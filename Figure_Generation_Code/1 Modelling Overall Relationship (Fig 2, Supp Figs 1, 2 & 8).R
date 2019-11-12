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
##        Figure 2: Modelling the Overall Relationship Between LM & PCR Prevalence               ##
##        Supplementary Figure 1: Modelling the Overall Relationship for Old & New Data          ##
##                                Separately                                                     ##
##        Suppplementary Figure 2: Comparing Observed and Model Predictions of LM Prevalence     ##
##        Supplementary Figure 9: MCMC Output from JAGS Model Fitting to Overall Data            ##
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

# Load in the dataset and subset the data into data from the previous review (old_data), data from 
# this review (new_data) all data together (full_data). Note: non-age disaggregated data is assigned coding 2)
data_frame <- read.csv("Data/Submicroscopic_Review_Data_R_Import.csv")
old_data <- data_frame[data_frame$Old_or_New == "Old" & data_frame$Full_Or_Age_Disagg_Data == 2, ] 
new_data <- data_frame[data_frame$Old_or_New == "New" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]


###################################################################################################
##                                                                                               ##
##         Fig 2 & Supp Fig 8: Full Data - Running the Bayesian Log-Linear Regression            ##
##                                                                                               ##
##    This section of the code runs the Bayesian Log-Linear Regression, with the model           ##
##    implemented in the statistical programming software JAGS, more info available here:        ##
##    (http://mcmc-jags.sourceforge.net/). The output of fitting the model to the collated       ##
##    malaria survey data are then processed and used to generate model predictions.             ##
##                                                                                               ##
###################################################################################################
params <- c("delt", "beta", "taud")
initial_values_function <- function(){
  list("beta" = 0.0001, "delt" = 0.001, "taud" = 0.001)
}
model_file <- "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt"
full_data_output <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, full_data)
  
# Supplementary Figure 8 Materials
full_data_param_table <- param_table(full_data_output, params)
plot(full_data_output)

# Saving the MCMC estimates of parameter mean values
full_overall_chain <- rbind(full_data_output[[1]], full_data_output[[2]], full_data_output[[3]], full_data_output[[4]])
full_beta_mean <- mean(as.array(full_overall_chain[, "beta"])) 
full_delt_mean <- mean(as.array(full_overall_chain[, "delt"])) 

# Calculating the Predicted Values and 95% CI of Microscopy Prevalence Given a PCR Prevalence
  # Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_full <- seq(0.001, 0.96, 0.001)
full_fitted_microscopy <- mean_output(PCR_prevalence_full, full_beta_mean, full_delt_mean)
credibles <- credible_intervals(PCR_prevalence_full, full_overall_chain)
full_data_credible_lower <- credibles$credible_lower
full_data_credible_upper <- credibles$credible_upper
full_data_sensitivity_credible_lower <- credibles$sensitivity_lower
full_data_sensitivity_credible_upper <- credibles$sensitivity_upper

# Figure 1A Plotting - Microscopy Prevalence Against PCR Prevalence - Data & Modelled Relationship
par(mfrow = c(1, 1)); xlim <- c(0, 100); ylim <- c(0, 100); 
plot(new_data$PCR_Prev * 100, new_data$Micro_Prev * 100, xlim = xlim, ylim = ylim, pch = 20, col = "#C6D3F3", 
     xlab = "PCR Prevalence (%)", ylab = "Microscopy Prevalence (%)", cex = 2, las = 1)
points(old_data$PCR_Prev * 100, old_data$Micro_Prev * 100, pch = 20, col = "#6BC24E", cex = 2)
lines(seq(0, 100, 0.01), seq(0, 100, 0.01), lwd = 2, lty = 2)
polygon(x = c(PCR_prevalence_full * 100, rev(PCR_prevalence_full * 100)), 
        y = c(full_data_credible_upper * 100, rev(full_data_credible_lower * 100)), 
        col = adjustcolor("#F44A96", alpha.f = 0.5), border = NA)
lines(PCR_prevalence_full * 100, full_fitted_microscopy * 100, col = "#F44A96", lwd = 3)
legend("topleft", legend = c("Old Data", "New Data",  "Old and New Data"), 
       col = c("#6BC24E", "#C6D3F3", "#F44A96"), pch = 20, pt.cex = 2, cex = 1.1)

# Figure 1B Plotting - Microscopy Sensitivity Against PCR Prevalence - Raw Data & Modelled Relationship
plot(full_data$PCR_Prev * 100, full_data$Micro_Prev/full_data$PCR_Prev, xlim = xlim, ylim = c(0, 1), pch = 20, col = "#BDBDBD",
     xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio", cex = 1.5, las = 1)
lines(PCR_prevalence_full * 100, (full_fitted_microscopy * 100)/(PCR_prevalence_full * 100), col = "#F44A96", lwd = 3)
polygon(x = c(PCR_prevalence_full * 100, rev(PCR_prevalence_full * 100)), 
        y = c(full_data_sensitivity_credible_upper, rev(full_data_sensitivity_credible_lower)), 
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

plot(0, 0, ylim = c(0, 1), xlim = c(0, 100), xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio", cex = 0, las = 1)
for (i in 1:length(Sens_by_group)){
  rect(xleft = breaks[i] * 100, ybottom = 0, xright = breaks[i + 1] * 100, ytop = Sens_by_group[i], col = "#D3ECF3")
  arrows(x0 = 100 * (breaks[i] + breaks[i + 1])/2, y0 = lower_ci[i],
         x1 = 100 * (breaks[i] + breaks[i + 1])/2, y1 = upper_ci[i],
         col=1, angle=90, code=3, length = 0.05)
}



###################################################################################################
##                                                                                               ##
##           Supp Fig 1: Old & New Data - Running the Bayesian Log-Linear Regression             ##
##                                                                                               ##
##    This section of the code runs the Bayesian Log-Linear Regression, with the model           ##
##    implemented in the statistical programming software JAGS, more info available here:        ##
##    (http://mcmc-jags.sourceforge.net/). The output of fitting the model to the collated       ##
##    malaria survey data are then processed and used to generate model predictions.             ##
##                                                                                               ##
###################################################################################################
old_data_output <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, old_data)
old_chain_overall <- rbind(old_data_output[[1]], old_data_output[[2]], old_data_output[[3]], old_data_output[[4]])
old_beta_mean <- mean(as.array(old_chain_overall[, "beta"])) 
old_delt_mean <- mean(as.array(old_chain_overall[, "delt"])) 
old_fitted_microscopy <- mean_output(PCR_prevalence_full, old_beta_mean, old_delt_mean)
old_credibles <- credible_intervals(PCR_prevalence_full, old_chain_overall)
old_data_credible_lower <- old_credibles$credible_lower
old_data_credible_upper <- old_credibles$credible_upper
old_data_sensitivity_credible_lower <- old_credibles$sensitivity_lower
old_data_sensitivity_credible_upper <- old_credibles$sensitivity_upper

new_data_output <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, new_data)
new_chain_overall <- rbind(new_data_output[[1]], new_data_output[[2]], new_data_output[[3]], new_data_output[[4]])
new_beta_mean <- mean(as.array(new_chain_overall[, "beta"])) 
new_delt_mean <- mean(as.array(new_chain_overall[, "delt"])) 
new_fitted_microscopy <- mean_output(PCR_prevalence_full, new_beta_mean, new_delt_mean)
new_credibles <- credible_intervals(PCR_prevalence_full, new_chain_overall)
new_data_credible_lower <- new_credibles$credible_lower
new_data_credible_upper <- new_credibles$credible_upper
new_data_sensitivity_credible_lower <- new_credibles$sensitivity_lower
new_data_sensitivity_credible_upper <- new_credibles$sensitivity_upper

# Supplementary Figure 1A Plotting - Microscopy Prevalence Against PCR Prevalence for Old and New Data Separately - Data & Modelled Relationship
par(mfrow = c(1, 2))
plot(new_data$PCR_Prev * 100, new_data$Micro_Prev * 100, xlim = xlim, ylim = ylim, pch = 20, col = "#C6D3F3", 
     xlab = "PCR Prevalence (%)", ylab = "Microscopy Prevalence (%)", cex = 1.5, las = 1, cex.lab = 1.1, cex.axis = 1.1)
points(old_data$PCR_Prev * 100, old_data$Micro_Prev * 100, xlim = xlim, ylim = ylim, pch = 20, col = "#6BC24E", cex = 1.5)
lines(seq(0, 100, 0.01), seq(0, 100, 0.01), lwd = 2, lty = 2)
polygon(x = c(PCR_prevalence_full * 100, rev(PCR_prevalence_full * 100)), 
        y = c(new_data_credible_upper * 100, rev(new_data_credible_lower * 100)), 
        col = adjustcolor("#C6D3F3", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_full * 100, rev(PCR_prevalence_full*100)), 
        y = c(old_data_credible_upper * 100, rev(old_data_credible_lower * 100)), 
        col = adjustcolor("#6BC24E", alpha.f = 0.5), border = NA)
lines(PCR_prevalence_full * 100, old_fitted_microscopy * 100, col = "#6BC24E", lwd = 3)
lines(PCR_prevalence_full * 100, new_fitted_microscopy * 100, col = "#C6D3F3", lwd = 3)
legend("topleft", legend = c("Old Data", "New Data"), col = c("#6BC24E", "#C6D3F3"), pch = 20, pt.cex = 2, cex = 0.8)

# Supplementary Figure 1B Plotting - Microscopy Sensitivity Against PCR Prevalence for Old and New Data Separately - Data & Modelled Relationship
plot(new_data$PCR_Prev * 100, new_data$Micro_Prev/new_data$PCR_Prev, xlim = xlim, ylim = c(0, 1), pch = 20, col = "#C6D3F3", 
     xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio", cex = 1.5, las = 1, las = 1, cex.lab = 1.1, cex.axis = 1.1)
points(old_data$PCR_Prev*100, old_data$Micro_Prev/old_data$PCR_Prev, xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#6BC24E", cex = 1.5)
polygon(x = c(PCR_prevalence_full * 100, rev(PCR_prevalence_full * 100)), 
        y = c(new_data_sensitivity_credible_upper, rev(new_data_sensitivity_credible_lower)), 
        col = adjustcolor("#C6D3F3", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_full*100, rev(PCR_prevalence_full*100)), 
        y = c(old_data_sensitivity_credible_upper, rev(old_data_sensitivity_credible_lower)), 
        col = adjustcolor("#6BC24E", alpha.f = 0.5), border = NA)
lines(PCR_prevalence_full * 100, old_fitted_microscopy/PCR_prevalence_full, col = "#6BC24E", lwd = 3)
lines(PCR_prevalence_full * 100, new_fitted_microscopy/PCR_prevalence_full, col = "#C6D3F3", lwd = 3)
legend("topleft", legend = c("Old Data", "New Data"), col = c("#6BC24E", "#C6D3F3"),
       pch = 20, pt.cex = 2, cex = 0.8)

# Supplementary Figure 2 Plotting - Predicted vs Observed Microscopy Prevalence
observed_PCR_prevalences_full <- full_data$PCR_Prev
matched_PCR_indices_full <-  match.closest(observed_PCR_prevalences_full, PCR_prevalence_full)
logit_matched_PCR <- logit(PCR_prevalence_full[matched_PCR_indices_full])
predicted_LM_logit <- (full_delt_mean - (full_beta_mean * mean(logit_matched_PCR))) + ((1 + full_beta_mean)  * logit_matched_PCR)
predicted_LM <- expit(predicted_LM_logit) 
plot(full_data$Micro_Prev, predicted_LM, xlim = c(0, 1), ylim = c(0, 1))
cor(full_data$Micro_Prev, predicted_LM)

observed_PCR_prevalences_old <- old_data$PCR_Prev
matched_PCR_indices_old <-  match.closest(observed_PCR_prevalences_old, PCR_prevalence_full)
logit_matched_PCR_old <- logit(PCR_prevalence_full[matched_PCR_indices_old])
predicted_LM_logit_old <- (full_delt_mean - (full_beta_mean * mean(logit_matched_PCR_old))) + ((1 + full_beta_mean)  * logit_matched_PCR_old)
predicted_LM_old <- expit(predicted_LM_logit_old) 
plot(old_data$Micro_Prev, predicted_LM_old, xlim = c(0, 1), ylim = c(0, 1), col = "#C6D3F3", pch = 20)
cor(old_data$Micro_Prev, predicted_LM_old)

observed_PCR_prevalences_new <- new_data$PCR_Prev
matched_PCR_indices_new <-  match.closest(observed_PCR_prevalences_new, PCR_prevalence_full)
logit_matched_PCR_new <- logit(PCR_prevalence_full[matched_PCR_indices_new])
predicted_LM_logit_new <- (full_delt_mean - (full_beta_mean * mean(logit_matched_PCR_new))) + ((1 + full_beta_mean)  * logit_matched_PCR_new)
predicted_LM_new <- expit(predicted_LM_logit_new) 
plot(new_data$Micro_Prev, predicted_LM_new, xlim = c(0, 1), ylim = c(0, 1), col = "#6BC24E", pch = 20)
cor(new_data$Micro_Prev, predicted_LM_new)

par(mfrow = c(1, 1))
plot(old_data$Micro_Prev*100, predicted_LM_old*100, xlim = c(0, 100), ylim = c(0, 100), col = "#C6D3F3", pch = 1, lwd = 2, ylab = "Predicted PCR Prevalence", xlab = "Observed PCR Prevalence", las = 1)
points(new_data$Micro_Prev*100, predicted_LM_new*100, xlim = c(0, 100), ylim = c(0, 100), pch = 1, lwd = 2, col = "#6BC24E")
lines(seq(0,100,0.01), seq(0,100,0.01), lwd = 2, lty = 2)
cor(full_data$Micro_Prev, predicted_LM)


###################################################################################################
##                                                                                               ##
##         Calculating average prevalence ratio and associated uncertainty overall               ##
##                                                                                               ##
###################################################################################################
LM_prev <- full_data$Micro_Prev
PCR_prev <- full_data$PCR_Prev
prev_ratio <- 100 * (LM_prev/PCR_prev) 
mean_prev_ratio <- mean(prev_ratio)
std_err <- sd(prev_ratio)/sqrt(length(prev_ratio))
mean_prev_ratio - 1.96 * std_err
mean_prev_ratio + 1.96 * std_err

