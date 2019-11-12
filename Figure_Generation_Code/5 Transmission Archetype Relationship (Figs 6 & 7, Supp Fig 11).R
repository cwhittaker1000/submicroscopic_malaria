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
##        Figure 5A: Africa Map With Surveyed Admin Units Coloured                               ##                    
##        Figure 5B: Modelling the Influence of Transmission Archetype on the Prevalence Ratio   ##
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
data_frame <- read.csv("Data/Submicroscopic_Review_Data_R_Import.csv")
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data$Sensitivity <- full_data$Micro_Prev/full_data$PCR_Prev
data_frame <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ] #Remove age disaggregated data
data_frame <- data_frame[order(data_frame$Transmission_Setting_15), ] # orders by transmission history status (puts NAs at end)
data_frame <- data_frame[1:167, ] # removes NAs which are surveys not conducted in Africa for which Trans_Hist data was not available 

# Subsetting the Data by Transmission Archetype
high_high_subset <- data_frame[data_frame$Transmission_Setting_15 == "High_High", ]
high_low_subset <- data_frame[data_frame$Transmission_Setting_15 == "High_Low", ]
low_low_subset <- data_frame[data_frame$Transmission_Setting_15 == "Low_Low", ]

###################################################################################################
##                                                                                               ##
##  Fig 6 & Supp Fig 11: Trans. Archetype Disaggregated Data - Running the Log-Linear Regression ##
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
High_high_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, high_high_subset)
High_low_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, high_low_subset)
Low_low_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, low_low_subset)

# Supplementary Figure 11 - MCMC Output and Parameter Tables
high_high_param_table <- param_table(High_high_model, params)
plot(High_high_model, col = c("#00A600FF"), las = 1)
high_low_param_table <- param_table(High_low_model, params)
plot(High_low_model, col = c("#ECB176FF"), las = 1)
low_low_param_table <- param_table(Low_low_model, params)
plot(Low_low_model, col = c("#ECB176FF"), las = 1)

# Processing the Output from the JAGS Models
    # Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
High_high_overall_chain <- rbind(High_high_model[[1]], High_high_model[[2]], High_high_model[[3]], High_high_model[[4]])
High_high_beta_mean <- mean(as.array(High_high_overall_chain[, "beta"]))
High_high_delt_mean <- mean(as.array(High_high_overall_chain[, "delt"]))
PCR_prevalence_high_high <- seq(0.03,0.95,0.001)
High_high_fitted_microscopy <- mean_output(PCR_prevalence_high_high, High_high_beta_mean, High_high_delt_mean)
High_high_credibles <- credible_intervals(PCR_prevalence_high_high, High_high_overall_chain)
High_high_credible_lower <- High_high_credibles$credible_lower
High_high_credible_upper <- High_high_credibles$credible_upper
High_high_prev_ratio_lower <- High_high_credibles$sensitivity_lower
High_high_prev_ratio_upper <- High_high_credibles$sensitivity_upper
  

High_low_overall_chain <- rbind(High_low_model[[1]], High_low_model[[2]], High_low_model[[3]], High_low_model[[4]])
High_low_beta_mean <- mean(as.array(High_low_overall_chain[, "beta"])) 
High_low_delt_mean <- mean(as.array(High_low_overall_chain[, "delt"]))
PCR_prevalence_high_low <- seq(0.004,0.85,0.001)
High_low_fitted_microscopy <- mean_output(PCR_prevalence_high_low, High_low_beta_mean, High_low_delt_mean)
High_low_credibles <- credible_intervals(PCR_prevalence_high_low, High_low_overall_chain)
High_low_credible_lower <- High_low_credibles$credible_lower
High_low_credible_upper <- High_low_credibles$credible_upper
High_low_prev_ratio_lower <- High_low_credibles$sensitivity_lower
High_low_prev_ratio_upper <- High_low_credibles$sensitivity_upper
  
Low_low_overall_chain <- rbind(Low_low_model[[1]], Low_low_model[[2]], Low_low_model[[3]], Low_low_model[[4]])
Low_low_beta_mean <- mean(as.array(Low_low_overall_chain[, "beta"]))
Low_low_delt_mean <- mean(as.array(Low_low_overall_chain[, "delt"]))
PCR_prevalence_low_low <- seq(0.01,0.55,0.001)
Low_low_fitted_microscopy <- mean_output(PCR_prevalence_low_low, Low_low_beta_mean, Low_low_delt_mean)
Low_low_credibles <- credible_intervals(PCR_prevalence_low_low, Low_low_overall_chain)
Low_low_credible_lower <- Low_low_credibles$credible_lower
Low_low_credible_upper <- Low_low_credibles$credible_upper
Low_low_prev_ratio_lower <- Low_low_credibles$sensitivity_lower
Low_low_prev_ratio_upper <- Low_low_credibles$sensitivity_upper


# Figure 6B Plotting - Sensitivity Against PCR Prevalence for All 3 Transmission Settings - Data & Modelled Relationship
# Panel 1 - Historically and Currently High Transmission Settings
plot(high_high_subset$PCR_Prev * 100, high_high_subset$Micro_Prev/high_high_subset$PCR_Prev,  las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#00A600FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_high * 100, High_high_fitted_microscopy/PCR_prevalence_high_high, col = "#00A600FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_high * 100, rev(PCR_prevalence_high_high * 100)), 
        y = c(High_high_prev_ratio_upper, rev(High_high_prev_ratio_lower)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)
# Panel 2 - Historically High But Currently Low Transmission Settings
plot(high_low_subset$PCR_Prev * 100, high_low_subset$Micro_Prev/high_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#ECB176FF", xlab = "PCR Prevalence (%)", ylab = "Prevalence Ratio")
lines(PCR_prevalence_high_low * 100, High_low_fitted_microscopy/PCR_prevalence_high_low, col = "#ECB176FF", lwd = 3)
polygon(x = c(PCR_prevalence_high_low * 100, rev(PCR_prevalence_high_low * 100)), 
        y = c(High_low_prev_ratio_upper, rev(High_low_prev_ratio_lower)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA)
# Panel 3 - Historically and Currently Low Transmission Settings
plot(low_low_subset$PCR_Prev * 100, low_low_subset$Micro_Prev/low_low_subset$PCR_Prev, las = 1, cex = 2,
     xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "darkgrey", xlab = "PCR Prevalence", ylab = "Prevalence Ratio")
lines(PCR_prevalence_low_low * 100, Low_low_fitted_microscopy/PCR_prevalence_low_low, col = "darkgrey", lwd = 3)
polygon(x = c(PCR_prevalence_low_low * 100, rev(PCR_prevalence_low_low * 100)), 
        y = c(Low_low_prev_ratio_upper, rev(Low_low_prev_ratio_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)

# Figure 6C Plotting - Microscopy Prevalence Against PCR Prevalence for All 3 Transmission Settings - Data & Modelled Relationship
par(mfrow = c(1, 1))
plot(high_high_subset$PCR_Prev * 100, high_high_subset$Micro_Prev * 100, xlim = c(0, 20), ylim = c(0, 20), pch = 20, col = "#00A600FF", las = 1, cex = 1.5,
     xlab = "PCR Prevalence (%)", ylab = "LM Prevalence (%)") # change xlim and ylim to capture the main figure panel xlim = c(0, 1), ylim = c(0, 1) and the inset (0, 0.2)
points(high_low_subset$PCR_Prev * 100, high_low_subset$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "#ECB176FF", cex = 1.5)
points(low_low_subset$PCR_Prev * 100, low_low_subset$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "darkgrey", cex = 1.5)
lines(seq(0,100,0.01), seq(0,100,0.01), lwd = 2, lty = 2)
lines(PCR_prevalence_high_high * 100, High_high_fitted_microscopy * 100, col = "#00A600FF", lwd = 3)
lines(PCR_prevalence_high_low * 100, High_low_fitted_microscopy * 100, col = "#ECB176FF", lwd = 3)
lines(PCR_prevalence_low_low * 100, Low_low_fitted_microscopy * 100, col = "darkgrey", lwd = 3)
polygon(x = c(PCR_prevalence_high_high * 100, rev(PCR_prevalence_high_high * 100)), y = c(High_high_credible_upper * 100, rev(High_high_credible_lower * 100)), 
        col = adjustcolor("#00A600FF", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_high_low * 100, rev(PCR_prevalence_high_low * 100)), y = c(High_low_credible_upper * 100, rev(High_low_credible_lower * 100)), 
        col = adjustcolor("#ECB176FF", alpha.f = 0.5), border = NA) 
polygon(x = c(PCR_prevalence_low_low * 100, rev(PCR_prevalence_low_low * 100)), y = c(Low_low_credible_upper * 100, rev(Low_low_credible_lower * 100)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)
legend("topleft", legend = c("High High", "High Low", "Low low"), col = c("#00A600FF", "#ECB176FF", "darkgrey"),
       lty = 1, lwd = 2, pt.cex = 2, cex = 0.8)


# Figure 7A Plotting - Contribution to Transmission in High High Transmission Settings
HH_PCR_Prevalence <- seq(0.03,0.95,0.001) 
HH_Patent_Percentage <- (High_high_fitted_microscopy/PCR_prevalence_high_high) * 100 
HH_Subpatent_Percentage <- (100 - HH_Patent_Percentage) 
HH_Subpatent_Contribution20 <- (HH_Subpatent_Percentage) / ((20 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH_Subpatent_Contribution5 <- (HH_Subpatent_Percentage) / ((5 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH_Subpatent_Contribution2 <- (HH_Subpatent_Percentage) / ((2 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), cex = 0)
lines(HH_PCR_Prevalence, HH_Subpatent_Contribution20, xlim = c(0, 1), ylim = c(0, 1))
lines(HH_PCR_Prevalence, HH_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(HH_PCR_Prevalence, HH_Subpatent_Contribution2, xlim = c(0, 1), ylim = c(0, 1))

# Figure 7B Plotting - Contribution to Transmission in High Low Transmission Settings
HL_PCR_Prevalence <- seq(0.004,0.85,0.001)
HL_Patent_Percentage <- (High_low_fitted_microscopy/PCR_prevalence_high_low) * 100 
HL_Subpatent_Percentage <- (100 - HL_Patent_Percentage) 
HL_Subpatent_Contribution20 <- (HL_Subpatent_Percentage) / ((20 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL_Subpatent_Contribution5 <- (HL_Subpatent_Percentage) / ((5 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL_Subpatent_Contribution2 <- (HL_Subpatent_Percentage) / ((2 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), cex = 0)
lines(HL_PCR_Prevalence, HL_Subpatent_Contribution20, xlim = c(0, 1), ylim = c(0, 1))
lines(HL_PCR_Prevalence, HL_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(HL_PCR_Prevalence, HL_Subpatent_Contribution2, xlim = c(0, 1), ylim = c(0, 1))

# Figure 7C Plotting - Contribution to Transmission in Low Low Transmission Settings
LL_PCR_Prevalence <- seq(0.01,0.55,0.001)
LL_Patent_Percentage <- (Low_low_fitted_microscopy/PCR_prevalence_low_low) * 100 
LL_Subpatent_Percentage <- (100 - LL_Patent_Percentage) 
LL_Subpatent_Contribution20 <- (LL_Subpatent_Percentage) / ((20 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL_Subpatent_Contribution5 <- (LL_Subpatent_Percentage) / ((5 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL_Subpatent_Contribution2 <- (LL_Subpatent_Percentage) / ((2 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), cex = 0)
lines(LL_PCR_Prevalence, LL_Subpatent_Contribution20, xlim = c(0, 1), ylim = c(0, 1))
lines(LL_PCR_Prevalence, LL_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(LL_PCR_Prevalence, LL_Subpatent_Contribution2, xlim = c(0, 1), ylim = c(0, 1))


