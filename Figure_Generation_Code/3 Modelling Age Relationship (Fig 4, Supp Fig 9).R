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
##        Figure 4: Modelling the Relationship Between LM & PCR Prevalence for Different         ## 
##                  Age Groups                                                                   ##
##        Supplementary Figure 9: MCMC Output from JAGS Model Fitting to Age Group Data          ##
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

# Load in the dataset and subset the data into children aged 0-5 years (young_children), old children, 
# aged 5-15 years (old_children) and data from adults, aged 15+ years (adults)
data_frame <- read.csv("Data/Submicroscopic_Review_Data_R_Import.csv")
young_children <- data_frame[data_frame$Age_Group == "0-5" , ]
old_children <- data_frame[data_frame$Age_Group ==  "5-15years", ]
adults <- data_frame[data_frame$Age_Group == "15+", ]

###################################################################################################
##                                                                                               ##
##    Fig 3 & Supp Fig 9: Age Disaggregated Data - Running the Bayesian Log-Linear Regression    ##
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
young_child_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, young_children)
old_child_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, old_children)
adult_model <- run_rJAGS_model(45000, 4, model_file, params, initial_values_function, adults)

# Supplementary Figure 9 - MCMC Output and Parameter Tables
young_children_param_table <- param_table(young_child_model, params)
plot(young_child_model, col = c("yellow"), las = 1)
old_children_param_table <- param_table(old_child_model, params)
plot(old_child_model, col = c("orange"), las = 1)
old_children_param_table <- param_table(adult_model, params)
plot(old_child_model, col = c("red"), las = 1)

# Processing the Output from the JAGS Models
    # Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
young_child <- rbind(young_child_model[[1]], young_child_model[[2]], young_child_model[[3]], young_child_model[[4]])
young_child_beta_mean <- mean(as.array(young_child[, "beta"]))
young_child_delt_mean <- mean(as.array(young_child[, "delt"]))
PCR_prevalence_young <- seq(0.005, 0.9, 0.001)
young_child_fitted_microscopy <- mean_output(PCR_prevalence_young, young_child_beta_mean, young_child_delt_mean)
young_credibles <- credible_intervals(PCR_prevalence_young, young_child)
young_children_credible_lower <- young_credibles$credible_lower
young_children_credible_upper <- young_credibles$credible_upper
young_prev_ratio_lower <- young_credibles$sensitivity_lower
young_prev_ratio_upper <- young_credibles$sensitivity_upper

old_child <- rbind(old_child_model[[1]], old_child_model[[2]], old_child_model[[3]], old_child_model[[4]])
old_child_beta_mean <- mean(as.array(old_child[, "beta"]))
old_child_delt_mean <- mean(as.array(old_child[, "delt"]))
PCR_prevalence_old <- seq(0.04, 0.97, 0.001)
old_child_fitted_microscopy <- mean_output(PCR_prevalence_old, old_child_beta_mean, old_child_delt_mean)
old_credibles <- credible_intervals(PCR_prevalence_old, old_child)
old_children_credible_lower <- old_credibles$credible_lower
old_children_credible_upper <- old_credibles$credible_upper
old_prev_ratio_lower <- old_credibles$sensitivity_lower
old_prev_ratio_upper <- old_credibles$sensitivity_upper

adult <- rbind(adult_model[[1]], adult_model[[2]], adult_model[[3]], adult_model[[4]])
adult_beta_mean <- mean(as.array(adult[, "beta"]))
adult_delt_mean <- mean(as.array(adult[, "delt"]))
PCR_prevalence_adults <- seq(0.02, 0.8, 0.001)
adult_fitted_microscopy <- mean_output(PCR_prevalence_adults, adult_beta_mean, adult_delt_mean)
adult_credibles <- credible_intervals(PCR_prevalence_adults, adult)
adult_credible_lower <- adult_credibles$credible_lower
adult_credible_upper <- adult_credibles$credible_upper
adult_prev_ratio_lower <- adult_credibles$sensitivity_lower
adult_prev_ratio_upper <- adult_credibles$sensitivity_upper

# Figure 2A Plotting - Microscopy Prevalence Against PCR Prevalence for All 3 Age Groups- Data & Modelled Relationship
mat <- matrix(c(1, 1, 1,
                1, 1, 1,
                1, 1, 1,
                1, 1, 1,
                2, 3, 4,
                2, 3, 4), nrow = 6, byrow = TRUE)
par(mar = c(6, 5, 4, 1))
layout(mat)
plot(young_children$PCR_Prev * 100, young_children$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "#EAE600",
     xlab = "PCR Prevalence (%)", ylab = "LM Prevalence (%)", las = 1, cex = 2, cex.lab = 1.75, cex.axis = 1.75)
points(old_children$PCR_Prev * 100, old_children$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "orange", cex = 1.5)
points(adults$PCR_Prev * 100, adults$Micro_Prev * 100, xlim = c(0, 100), ylim = c(0, 100), pch = 20, col = "red", cex = 1.5)
lines(seq(0,100,0.01), seq(0,100,0.01), lwd = 2, lty = 2)
lines(PCR_prevalence_young * 100, young_child_fitted_microscopy * 100, col = "#EAE600", lwd = 4)
lines(PCR_prevalence_old * 100, old_child_fitted_microscopy * 100, col = "orange", lwd = 4)
lines(PCR_prevalence_adults * 100, adult_fitted_microscopy * 100, col = "red", lwd = 4)
polygon(x = c(PCR_prevalence_young * 100, rev(PCR_prevalence_young * 100)), 
        y = c(young_children_credible_upper * 100, rev(young_children_credible_lower * 100)), 
        col = adjustcolor("#EAE600", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_old * 100, rev(PCR_prevalence_old * 100)), 
        y = c(old_children_credible_upper * 100, rev(old_children_credible_lower * 100)), 
        col = adjustcolor("orange", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_adults * 100, rev(PCR_prevalence_adults * 100)), 
        y = c(adult_credible_upper * 100, rev(adult_credible_lower * 100)), 
        col = adjustcolor("red", alpha.f = 0.5), border = NA)
legend("topleft", legend = c("Young Children", "Older Children", "Adults"), col = c("#EAE600", "orange", "red"),
       lty = 1, lwd = 2, pt.cex = 2, cex = 1.5)

# Figure 2B Panel 1 Plotting - Microscopy prev_ratio Against PCR Prevalence for YOUNG CHILDREN - Raw Data & Modelled 
par(mar = c(4, 5, 1, 0))
plot(young_children$PCR_Prev * 100, young_children$Micro_Prev/young_children$PCR_Prev, xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "#EAE600",
     xlab = "PCR Prevalence (%)", las = 1, cex = 2, ylab = "Prevalence Ratio",  cex.lab = 1.5, cex.axis = 1.5)
lines(PCR_prevalence_young * 100, young_child_fitted_microscopy/PCR_prevalence_young, col = "#EAE600", lwd = 3)
polygon(x = c(PCR_prevalence_young * 100, rev(PCR_prevalence_young * 100)), 
        y = c(young_prev_ratio_upper, rev(young_prev_ratio_lower)), 
        col = adjustcolor("#EAE600", alpha.f = 0.5), border = NA)

# Figure 2B Panel 2 Plotting - Microscopy prev_ratio Against PCR Prevalence for OLD CHILDREN - Raw Data & Modelled Relationship
par(mar = c(4, 4, 1, 0.5))
plot(old_children$PCR_Prev * 100, old_children$Micro_Prev/old_children$PCR_Prev, xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "orange",
     xlab = "PCR Prevalence (%)", las = 1, cex = 2, yaxt = "n", ylab = "",  cex.lab = 1.5, cex.axis = 1.5)
lines(PCR_prevalence_old * 100, old_child_fitted_microscopy/PCR_prevalence_old, col = "orange", lwd = 3)
polygon(x = c(PCR_prevalence_old * 100, rev(PCR_prevalence_old * 100)), 
        y = c(old_prev_ratio_upper, rev(old_prev_ratio_lower)), 
        col = adjustcolor("orange", alpha.f = 0.5), border = NA)

# Figure 2B Panel 3 Plotting - Microscopy prev_ratio Against PCR Prevalence for ADULTS - Raw Data & Modelled Relationship
par(mar = c(4, 4, 1, 1))
plot(adults$PCR_Prev * 100, adults$Micro_Prev/adults$PCR_Prev, xlim = c(0, 100), ylim = c(0, 1), pch = 20, col = "red",
     xlab = "PCR Prevalence (%)", las = 1, cex = 2, yaxt = "n", ylab = "",  cex.lab = 1.5, cex.axis = 1.5)
lines(PCR_prevalence_adults * 100, adult_fitted_microscopy/PCR_prevalence_adults, col = "red", lwd = 3)
polygon(x = c(PCR_prevalence_adults * 100, rev(PCR_prevalence_adults * 100)), 
        y = c(adult_prev_ratio_upper, rev(adult_prev_ratio_lower)), 
        col = adjustcolor("red", alpha.f = 0.5), border = NA)

# Statistical Tests Carried Out On The Data
# ANOVA - Testing for Differences in Means
data_frame_ANOVA <- data_frame[(data_frame$Age_Group == "0-5") | (data_frame$Age_Group ==  "5-15years") | (data_frame$Age_Group ==  "15+"), ] # all age-disaggregated data
data_frame_ANOVA$prev_ratio <- (data_frame_ANOVA$Microscopy_N_Positive/data_frame_ANOVA$Microscopy_N_Tested)/
                                (data_frame_ANOVA$PCR_N_Positive/data_frame_ANOVA$PCR_N_Tested)
ANOVA_object <- aov(prev_ratio ~ Age_Group, data = data_frame_ANOVA)
summary(ANOVA_object)
TukeyHSD(ANOVA_object)

