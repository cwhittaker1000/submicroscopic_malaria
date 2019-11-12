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
##        Figure 5: Modelling the Relationship Between LM & PCR Prevalence for Different         ## 
##                  Global Retions                                                               ##
##        Supplementary Figure 10: MCMC Output from JAGS Model Fitting to Global Region Data     ##
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

# Load in the dataset and subset the data by global region the survey was carried out in:
data_frame <- read.csv("Data/Submicroscopic_Review_Data_R_Import.csv")
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data$Sensitivity <- full_data$Micro_Prev/full_data$PCR_Prev
Asia_Oceania <- data_frame[data_frame$Global_Region == "Asia&Oceania" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
South_America <- data_frame[data_frame$Global_Region == "South America" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa" & data_frame$Full_Or_Age_Disagg_Data == 2, ]

###################################################################################################
##                                                                                               ##
##  Fig 5 & Supp Fig 10: Global Region Disaggregated Data - Running the Log-Linear Regression    ##
##                                                                                               ##
##    This section of the code runs the Bayesian Log-Linear Regression, with the model           ##
##    implemented in the statistical programming software JAGS, more info available here:        ##
##    (http://mcmc-jags.sourceforge.net/). The output of fitting the model to the collated       ##
##    malaria survey data are then processed and used to generate model predictions.             ##
##                                                                                               ##
###################################################################################################
# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud")
initial_values_function <- function(){
  list("beta" = 0.0001, "delt" = 0.001, "taud" = 0.001)
}
model_file <- "JAGS_Model/LM_Standard_Bayesian_Logit_Linear_Model.txt"

# Running the JAGS models for each dataset
Asia_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, Asia_Oceania)
East_Africa_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, East_Africa)
South_America_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, South_America)
West_Africa_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, West_Africa)

# Supplementary Figure 10 - MCMC Output and Parameter Tables
Asia_param_table <- param_table(Asia_model, params)
plot(Asia_model, col = c("black"), las = 1)
East_africa_param_table <- param_table(East_Africa_model, params)
plot(East_Africa_model, col = c("green"), las = 1)
South_America_param_table <- param_table(South_America_model, params)
plot(South_America_model, col = c("blue"), las = 1)
West_Africa_param_table <- param_table(West_Africa_model, params)
plot(West_Africa_model, col = c("red"), las = 1)

# Processing the Output from the JAGS Models
    # Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
Asia_overall_chain <- rbind(Asia_model[[1]], Asia_model[[2]], Asia_model[[3]], Asia_model[[4]])
Asia_beta_mean <- mean(as.array(Asia_overall_chain[, "beta"]))
Asia_delt_mean <- mean(as.array(Asia_overall_chain[, "delt"]))
PCR_prevalence_Asia <- seq(0.001,0.8,0.001)
Asia_fitted_microscopy <- mean_output(PCR_prevalence_Asia, Asia_beta_mean, Asia_delt_mean)
Asia_credibles <- credible_intervals(PCR_prevalence_Asia, Asia_overall_chain)
Asia_credible_lower <- Asia_credibles$credible_lower
Asia_credible_upper <- Asia_credibles$credible_upper
Asia_prev_ratio_lower <- Asia_credibles$sensitivity_lower
Asia_prev_ratio_upper <- Asia_credibles$sensitivity_upper

East_Africa_overall_chain <- rbind(East_Africa_model[[1]], East_Africa_model[[2]], East_Africa_model[[3]], East_Africa_model[[4]])
East_Africa_beta_mean <- mean(as.array(East_Africa_overall_chain[, "beta"]))
East_Africa_delt_mean <- mean(as.array(East_Africa_overall_chain[, "delt"]))
PCR_prevalence_East_Africa <- seq(0.001,0.92,0.001)
East_Africa_fitted_microscopy <- mean_output(PCR_prevalence_East_Africa, East_Africa_beta_mean, East_Africa_delt_mean)
East_Africa_credibles <- credible_intervals(PCR_prevalence_East_Africa, East_Africa_overall_chain)
East_Africa_credible_lower <- East_Africa_credibles$credible_lower
East_Africa_credible_upper <- East_Africa_credibles$credible_upper
East_Africa_prev_ratio_lower <- East_Africa_credibles$sensitivity_lower
East_Africa_prev_ratio_upper <- East_Africa_credibles$sensitivity_upper

South_America_overall_chain <- rbind(South_America_model[[1]], South_America_model[[2]], South_America_model[[3]], South_America_model[[4]])
South_America_beta_mean <- mean(as.array(South_America_overall_chain[, "beta"]))
South_America_delt_mean <- mean(as.array(South_America_overall_chain[, "delt"]))
PCR_prevalence_South_America <- seq(0.001,0.47,0.001)
South_America_fitted_microscopy <- mean_output(PCR_prevalence_South_America, South_America_beta_mean, South_America_delt_mean)
South_America_credibles <- credible_intervals(PCR_prevalence_South_America, South_America_overall_chain)
South_America_credible_lower <- South_America_credibles$credible_lower
South_America_credible_upper <- South_America_credibles$credible_upper
South_America_prev_ratio_lower <- South_America_credibles$sensitivity_lower
South_America_prev_ratio_upper <- South_America_credibles$sensitivity_upper

West_Africa_overall_chain <- rbind(West_Africa_model[[1]], West_Africa_model[[2]], West_Africa_model[[3]], West_Africa_model[[4]])
West_Africa_beta_mean <- mean(as.array(West_Africa_model[, "beta"]))
West_Africa_delt_mean <- mean(as.array(West_Africa_model[, "delt"]))
PCR_prevalence_West_Africa <- seq(0.008,0.97,0.001)
West_Africa_fitted_microscopy <- mean_output(PCR_prevalence_West_Africa, West_Africa_beta_mean, West_Africa_delt_mean)
West_Africa_credibles <- credible_intervals(PCR_prevalence_West_Africa, West_Africa_overall_chain)
West_Africa_credible_lower <- West_Africa_credibles$credible_lower
West_Africa_credible_upper <- West_Africa_credibles$credible_upper
West_Africa_prev_ratio_lower <- West_Africa_credibles$sensitivity_lower
West_Africa_prev_ratio_upper <- West_Africa_credibles$sensitivity_upper


# Figure 3A Plotting - Microscopy Prevalence Against PCR Prevalence for West Africa - Data & Modelled Relationship
par(mfrow = c(1, 1))
plot(West_Africa$PCR_Prev * 100, West_Africa$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "red",
     xlab = "PCR Prevalence (%)", ylab = "Slide Prevalence (%)", las = 1, cex = 1.5)
lines(PCR_prevalence_West_Africa * 100, West_Africa_fitted_microscopy * 100, col = "red", lwd = 3)
lines(seq(0,100,0.01), seq(0,100,0.01), col = "black", lwd = 2, lty = 2)
polygon(x = c(PCR_prevalence_West_Africa * 100, rev(PCR_prevalence_West_Africa * 100)), 
        y = c(West_Africa_credible_upper * 100, rev(West_Africa_credible_lower * 100)), 
        col = adjustcolor("red", alpha.f = 0.5), border = NA)

# Figure 3B Plotting - Microscopy Prevalence Against PCR Prevalence for East Africa - Data & Modelled Relationship
plot(East_Africa$PCR_Prev * 100, East_Africa$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "green3",
     xlab = "PCR Prevalence (%)", ylab = "Slide Prevalence (%)", las = 1, cex = 1.5)
lines(PCR_prevalence_East_Africa * 100, East_Africa_fitted_microscopy * 100, col = "green3", lwd = 3)
lines(seq(0,100,0.01), seq(0,100,0.01), col = "black", lwd = 2, lty = 2)
polygon(x = c(PCR_prevalence_East_Africa * 100, rev(PCR_prevalence_East_Africa * 100)), 
        y = c(East_Africa_credible_upper * 100, rev(East_Africa_credible_lower * 100)), 
        col = adjustcolor("green3", alpha.f = 0.5), border = NA)

# Figure 3D Plotting - Microscopy Prevalence Against PCR Prevalence for South America - Data & Modelled Relationship
plot(South_America$PCR_Prev * 100, South_America$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "blue",
     xlab = "PCR Prevalence (%)", ylab = "Slide Prevalence (%)", las = 1)
lines(PCR_prevalence_South_America * 100, South_America_fitted_microscopy * 100, col = "blue", lwd = 3)
lines(seq(0,100,0.01), seq(0,100,0.01), col = "black", lwd = 2, lty = 2)
polygon(x = c(PCR_prevalence_South_America * 100, rev(PCR_prevalence_South_America * 100)), 
        y = c(South_America_credible_upper * 100, rev(South_America_credible_lower * 100)), 
        col = adjustcolor("blue", alpha.f = 0.5), border = NA)

# Figure 3D Plotting - Microscopy Prevalence Against PCR Prevalence for Asia & Oceania - Data & Modelled Relationship
plot(Asia_Oceania$PCR_Prev * 100, Asia_Oceania$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "black",
     xlab = "PCR Prevalence (%)", ylab = "Slide Prevalence (%)", las = 1)
lines(PCR_prevalence_Asia * 100, Asia_fitted_microscopy * 100, col = "black", lwd = 3)
lines(seq(0,100,0.01), seq(0,100,0.01), col = "black", lwd = 2, lty = 2)
polygon(x = c(PCR_prevalence_Asia * 100, rev(PCR_prevalence_Asia * 100)), 
        y = c(Asia_credible_upper * 100, rev(Asia_credible_lower * 100)), 
        col = adjustcolor("black", alpha.f = 0.5), border = NA)

# Figure 3E Plotting - Microscopy Sensitivity Against PCR Prevalence for Asia & Oceania - Modelled Relationship
par(mfrow = c(1, 1))
plot(0, 0, xlim = c(0, 100), ylim = c(0, 1), pch = 1, xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio", cex = 0, las = 1)
lines(PCR_prevalence_East_Africa * 100, East_Africa_fitted_microscopy/PCR_prevalence_East_Africa, col = adjustcolor("green3"), lwd = 3)
lines(PCR_prevalence_West_Africa * 100, West_Africa_fitted_microscopy/PCR_prevalence_West_Africa, col = adjustcolor("red"), lwd = 3)
lines(PCR_prevalence_Asia * 100, Asia_fitted_microscopy/PCR_prevalence_Asia, col = adjustcolor("black"), lwd = 3)
lines(PCR_prevalence_South_America * 100, South_America_fitted_microscopy/PCR_prevalence_South_America, col = adjustcolor("blue"), lwd = 3)
polygon(x = c(PCR_prevalence_West_Africa * 100, rev(PCR_prevalence_West_Africa * 100)), 
        y = c(West_Africa_prev_ratio_upper, rev(West_Africa_prev_ratio_lower)), 
        col = adjustcolor("red", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_South_America * 100, rev(PCR_prevalence_South_America * 100)), 
        y = c(South_America_prev_ratio_upper, rev(South_America_prev_ratio_lower)), 
        col = adjustcolor("blue", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_East_Africa * 100, rev(PCR_prevalence_East_Africa * 100)), 
        y = c(East_Africa_prev_ratio_upper, rev(East_Africa_prev_ratio_lower)), 
        col = adjustcolor("green3", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_Asia * 100, rev(PCR_prevalence_Asia * 100)), 
        y = c(Asia_prev_ratio_upper, rev(Asia_prev_ratio_lower)), 
        col = adjustcolor("black", alpha.f = 0.5), border = NA)

# Statistical Tests Carried Out On The Data
# ANOVA - Testing for Differences in Means
variance <- full_data$PCR_N_Tested * full_data$PCR_Prev * (1 - full_data$PCR_Prev)
stdev <- sqrt(variance)
weighted_global_region_model <- lm(Sensitivity ~ Global_Region, data = full_data, na.action = na.omit, weights = 1/stdev) # similar results with 1/variance
summary(weighted_global_region_model)
ANOVA_object <- aov(weighted_global_region_model)
summary(ANOVA_object)
TukeyHSD(ANOVA_object)


