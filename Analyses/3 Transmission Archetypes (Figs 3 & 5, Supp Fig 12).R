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
library(tictoc); library(BayesianTools); library(R2jags); library(bayesmix); library(dplyr);
library(ggplot2); library(cowplot)

setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/Sub-Patent Malarial Infections/Sub_Patent_Malaria_Analysis/")
source("Functions/Submicroscopic_Analysis_Functions.R")
seed <- 193
fresh_run <- FALSE

# Load in the dataset and subset the data by the survey region's transmission history (African surveys only):
data_frame <- read.csv("Data/SI_Systematic_Review_Results_R_Import.csv")
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data <- full_data[order(full_data$Transmission_Setting_15), ]
full_data <- full_data[1:229, ] # removes NAs which are surveys not conducted in Africa for which Trans_Hist data was not available 
full_data <- full_data %>%
  mutate(PCR_Prev = 100 * (PCR_N_Positive/PCR_N_Tested)) %>%
  mutate(LM_Prev = 100 * (Microscopy_N_Positive/Microscopy_N_Tested)) %>%
  mutate(Prev_Ratio = LM_Prev/PCR_Prev)

# Subsetting the Data by Transmission Archetype
high_high_subset <- full_data[full_data$Transmission_Setting_15 == "High_High", ]
high_low_subset <- full_data[full_data$Transmission_Setting_15 == "High_Low", ]
low_low_subset <- full_data[full_data$Transmission_Setting_15 == "Low_Low", ]

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
if (fresh_run == TRUE) {
  High_high_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, high_high_subset)
  High_low_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, high_low_subset)
  Low_low_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, low_low_subset)
  saveRDS(High_high_model, "Outputs/High_High_MCMC_Output.rds")
  saveRDS(High_low_model, "Outputs/High_Low_MCMC_Output.rds")
  saveRDS(Low_low_model, "Outputs/Low_Low_MCMC_Output.rds")
} else {
  High_high_model <- readRDS("Outputs/High_High_MCMC_Output.rds")
  High_low_model <- readRDS("Outputs/High_Low_MCMC_Output.rds")
  Low_low_model <- readRDS("Outputs/Low_Low_MCMC_Output.rds")
}


# Supplementary Figure 11 - MCMC Output and Parameter Tables
pdf("Figures/Supplementary/Supp Figure 12 - Transmission Archetype MCMC Output/Supp Figure 12 - High High MCMC Output.pdf", width = 7.33, height = 7.51)
high_high_param_table <- param_table(High_high_model, params)
plot(High_high_model, col = c("#00A600FF"), las = 1)
dev.off()

pdf("Figures/Supplementary/Supp Figure 12 - Transmission Archetype MCMC Output/Supp Figure 12 - High Low MCMC Output.pdf", width = 7.33, height = 7.51)
high_low_param_table <- param_table(High_low_model, params)
plot(High_low_model, col = c("#ECB176FF"), las = 1)
dev.off()

pdf("Figures/Supplementary/Supp Figure 12 - Transmission Archetype MCMC Output/Supp Figure 12 - Low Low MCMC Output.pdf", width = 7.33, height = 7.51)
low_low_param_table <- param_table(Low_low_model, params)
plot(Low_low_model, col = c("dark grey"), las = 1)
dev.off()

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
High_High_plotting <- data.frame(LM_Prev = 100 * High_high_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_high_high, 
                            lower = 100 * High_high_credible_lower, upper = 100 * High_high_credible_upper,
                            Prev_Ratio = (High_high_fitted_microscopy/PCR_prevalence_high_high),
                            lower_pr = High_high_prev_ratio_lower, upper_pr = High_high_prev_ratio_upper)
  

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
High_Low_plotting <- data.frame(LM_Prev = 100 * High_low_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_high_low, 
                                lower = 100 * High_low_credible_lower, upper = 100 * High_low_credible_upper,
                                Prev_Ratio = (High_low_fitted_microscopy/PCR_prevalence_high_low),
                                lower_pr = High_low_prev_ratio_lower, upper_pr = High_low_prev_ratio_upper)
  
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
Low_Low_plotting <- data.frame(LM_Prev = 100 * Low_low_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_low_low, 
                                lower = 100 * Low_low_credible_lower, upper = 100 * Low_low_credible_upper,
                                Prev_Ratio = (Low_low_fitted_microscopy/PCR_prevalence_low_low),
                                lower_pr = Low_low_prev_ratio_lower, upper_pr = Low_low_prev_ratio_upper)


a <- ggplot(high_high_subset, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "#00A600FF") +
  geom_ribbon(data = High_High_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#00A600FF", linetype = 0) + 
  geom_line(data = High_High_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "#00A600FF", size = 2) +
  labs(y = "Prevalence Ratio", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.15, 1, 0.7),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

b <- ggplot(high_low_subset, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "#ECB176FF") +
  geom_ribbon(data = High_Low_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#ECB176FF", linetype = 0) + 
  geom_line(data = High_Low_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "#ECB176FF", size = 2) +
  labs(y = "", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.2, 1, 0),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

c <- ggplot(low_low_subset, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "darkgrey") +
  geom_ribbon(data = Low_Low_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "darkgrey", linetype = 0) + 
  geom_line(data = Low_Low_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "darkgrey", size = 2) +
  labs(y = "", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.15, 1, 0),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

d <- ggplot(low_low_subset, aes(x = PCR_Prev, y = LM_Prev)) +
  geom_point(size = 2, col = "darkgrey") +
  geom_point(data = high_low_subset, aes(x = PCR_Prev, y = LM_Prev), col = "#ECB176FF", size = 2) +
  geom_point(data = high_high_subset, aes(x = PCR_Prev, y = LM_Prev), col = "#00A600FF", size = 2) +
  geom_ribbon(data = Low_Low_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "darkgrey", linetype = 0) +
  geom_ribbon(data = High_Low_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "#ECB176FF", linetype = 0) + 
  geom_ribbon(data = High_High_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "#00A600FF", linetype = 0) +
  geom_line(data = Low_Low_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "darkgrey", size = 2) +
  geom_line(data = High_Low_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#ECB176FF", size = 2) +
  geom_line(data = High_High_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#00A600FF", size = 2) +
  geom_segment(x = 0, y = 0, xend = 100, yend = 100, colour = "black", size = 0.5, linetype = 2) +
  labs(y = "Microscopy\nPrevalence (%)", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.2, 2, 1, 2),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_rect(xmin = 0, xmax = 20, ymin = 0, ymax = 20, color = "black", fill = NA, size = 0.25) +
  geom_segment(x = 0, y = 20, xend = -3.3, yend = 43, colour = "black", size = 0.25, linetype = 1) +
  geom_segment(x = 20, y = 20, xend = 33.8, yend = 42, colour = "black", size = 0.25, linetype = 1)
  

e <- ggplot(low_low_subset, aes(x = PCR_Prev, y = LM_Prev)) +
  geom_point(size = 2, col = "darkgrey") +
  geom_point(data = high_low_subset, aes(x = PCR_Prev, y = LM_Prev), col = "#ECB176FF", size = 2) +
  geom_point(data = high_high_subset, aes(x = PCR_Prev, y = LM_Prev), col = "#00A600FF", size = 2) +
  geom_ribbon(data = Low_Low_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "darkgrey", linetype = 0) +
  geom_ribbon(data = High_Low_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "#ECB176FF", linetype = 0) + 
  geom_ribbon(data = High_High_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "#00A600FF", linetype = 0) +
  geom_line(data = Low_Low_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "darkgrey", size = 2) +
  geom_line(data = High_Low_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#ECB176FF", size = 2) +
  geom_line(data = High_High_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#00A600FF", size = 2) +
  geom_segment(x = 0, y = 0, xend = 100, yend = 100, colour = "black", size = 0.5, linetype = 2) +
  labs(y = "", x = "") + 
  scale_y_continuous(limits = c(0, 20)) +
  scale_x_continuous(limits = c(0, 20)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_blank(), axis.ticks.length = unit(0, "mm"), axis.title = element_text(size = 12,face="bold"),
        plot.margin = unit(c(-1, -1, -5, -5),"mm"), legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

combined <- plot_grid(a, b, c, label_size = 30, ncol = 3, rel_widths = c(1.25, 1.05, 1.03))
combined_two <- plot_grid(combined, d, nrow = 2, rel_heights = c(0.9, 1.2))
combined_two +
  draw_plot_label(
    c("B", "C"), 
    c(0, 0),
    c(1, 0.55),
    size = 30) +
  draw_plot(e, x = 0.175, y = 0.30, width = 0.25, height = 0.25)
ggsave("Figures/Figure 3 - Transmission Archetype/Figure_3.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 10.21, height = 8.20, units = c("in", "cm", "mm"),
       dpi = 300, useDingbats = FALSE)
  
# Statistical Tests Carried Out On The Data
# ANOVA - Testing for Differences in Means
variance <- full_data$PCR_N_Tested * (full_data$PCR_Prev/100) * (1 - full_data$PCR_Prev/100)
stdev <- sqrt(variance)
weighted_arch_region_model <- lm(Prev_Ratio ~ Transmission_Setting_15 + PCR_Prev, data = full_data, na.action = na.omit, weights = 1/variance) # similar results with 1/variance
summary(weighted_arch_region_model)
ANOVA_object <- aov(weighted_arch_region_model)
summary(ANOVA_object)
TukeyHSD(ANOVA_object, which = "Transmission_Setting_15")

variance <- full_data$PCR_N_Tested * (full_data$PCR_Prev/100) * (1 - full_data$PCR_Prev/100)
stdev <- sqrt(variance)
weighted_arch_region_model <- lm(Prev_Ratio ~ Hist_Trans + Curr_Trans + PCR_Prev, data = full_data, na.action = na.omit, weights = 1/variance) # similar results with 1/variance
summary(weighted_arch_region_model)
ANOVA_object <- aov(weighted_arch_region_model)
summary(ANOVA_object)
TukeyHSD(ANOVA_object, which = "Transmission_Setting_15")

arch_region_model <- lm(Prev_Ratio ~ Transmission_Setting_15 + PCR_Prev, data = full_data, na.action = na.omit) # similar results with 1/variance
summary(arch_region_model)
ANOVA_object <- aov(arch_region_model)
summary(ANOVA_object)
TukeyHSD(ANOVA_object, which = "Transmission_Setting_15")

# Figure 5A Plotting - Contribution to Transmission in High High Transmission Settings
HH_PCR_Prevalence <- seq(0.03,0.95,0.001) 
HH_zeroes <- rep(0, length(HH_PCR_Prevalence))
HH_Patent_Percentage <- (High_high_fitted_microscopy/PCR_prevalence_high_high) * 100 
HH_Subpatent_Percentage <- (100 - HH_Patent_Percentage) 
HH_Subpatent_Contribution20 <- 100 * (HH_Subpatent_Percentage) / ((20 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH_Subpatent_Contribution5 <- 100 * (HH_Subpatent_Percentage) / ((5 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH_Subpatent_Contribution2 <- 100 * (HH_Subpatent_Percentage) / ((2 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH <- data.frame(Zeroes = HH_zeroes, Least = HH_Subpatent_Contribution20, Mid = HH_Subpatent_Contribution5, 
                 Most = HH_Subpatent_Contribution2, PCR_Prev = HH_PCR_Prevalence)

# Figure 5B Plotting - Contribution to Transmission in High Low Transmission Settings
HL_PCR_Prevalence <- seq(0.004,0.85,0.001)
HL_zeroes <- rep(0, length(HL_PCR_Prevalence))
HL_Patent_Percentage <- (High_low_fitted_microscopy/PCR_prevalence_high_low) * 100 
HL_Subpatent_Percentage <- (100 - HL_Patent_Percentage) 
HL_Subpatent_Contribution20 <- 100 * (HL_Subpatent_Percentage) / ((20 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL_Subpatent_Contribution5 <- 100 * (HL_Subpatent_Percentage) / ((5 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL_Subpatent_Contribution2 <- 100 * (HL_Subpatent_Percentage) / ((2 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL <- data.frame(Zeroes = HL_zeroes, Least = HL_Subpatent_Contribution20, Mid = HL_Subpatent_Contribution5, 
                 Most = HL_Subpatent_Contribution2, PCR_Prev = HL_PCR_Prevalence)

# Figure 5C Plotting - Contribution to Transmission in Low Low Transmission Settings
LL_PCR_Prevalence <- seq(0.01,0.55,0.001)
LL_zeroes <- rep(0, length(LL_PCR_Prevalence))
LL_Patent_Percentage <- (Low_low_fitted_microscopy/PCR_prevalence_low_low) * 100 
LL_Subpatent_Percentage <- (100 - LL_Patent_Percentage) 
LL_Subpatent_Contribution20 <- 100 * (LL_Subpatent_Percentage) / ((20 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL_Subpatent_Contribution5 <- 100 * (LL_Subpatent_Percentage) / ((5 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL_Subpatent_Contribution2 <- 100 * (LL_Subpatent_Percentage) / ((2 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL <- data.frame(Zeroes = LL_zeroes, Least = LL_Subpatent_Contribution20, Mid = LL_Subpatent_Contribution5, 
                 Most = LL_Subpatent_Contribution2, PCR_Prev = LL_PCR_Prevalence)

a <- ggplot(HH, aes(x = Zeroes, y = PCR_Prev)) +
  geom_ribbon(data = HH, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Most), fill = "#CCE4CC", linetype = 0) +
  geom_ribbon(data = HH, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Mid), fill = "#96CC9A", linetype = 0) +
  geom_ribbon(data = HH, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Least), fill = "#2DA440", linetype = 0) +
  labs(y = "Proportion Transmission Attributable\n to Submicroscopic Infections", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) + 
  theme_bw() +
  theme(legend.position = c(0.13, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 12, vjust = +6),
        plot.margin = unit(c(1.5, 0.1, 1, 1),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

b <- ggplot(HL, aes(x = Zeroes, y = PCR_Prev)) +
  geom_ribbon(data = HL, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Most), fill = "#FBEBDC", linetype = 0) +
  geom_ribbon(data = HL, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Mid), fill = "#F5D6B7", linetype = 0) +
  geom_ribbon(data = HL, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Least), fill = "#EDB175", linetype = 0) +
  labs(y = "", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 90)) + 
  theme_bw() +
  theme(legend.position = c(0.13, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(1.5, 0.1, 1, 0.1),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


c <- ggplot(LL, aes(x = Zeroes, y = PCR_Prev)) +
  geom_ribbon(data = LL, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Most), fill = "#E8E8EA", linetype = 0) +
  geom_ribbon(data = LL, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Mid), fill = "#D3D1D2", linetype = 0) +
  geom_ribbon(data = LL, aes(x = 100 * PCR_Prev, ymin = Zeroes, ymax = Least), fill = "#A5A5A7", linetype = 0) +
  labs(y = "", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 60)) + 
  theme_bw() +
  theme(legend.position = c(0.13, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(1.5, 1.5, 1, 0.1),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#12.23 and 4.6
combined <- plot_grid(a, b, c, labels = c('', '', ''), label_size = 30, ncol = 3, rel_widths = c(1, 0.8, 0.75)) +
  draw_plot_label(
    c("A", "B", "C"), 
    c(0.083, 0.41, 0.71),
    c(1, 1, 1),
    size = 30) +
  draw_plot_label(
    c("Historically High, Currently High", "Historically High, Currently Low", "Historically Low, Currently Low"), 
    c(0.03, 0.36, 0.66),
    c(0.96, 0.96, 0.96),
    size = 12)
combined
ggsave("Figures/Figure 5 - Contribution To Transmission/Figure_5.pdf", width = 12.23, height = 4.6, plot = last_plot(), device = NULL, path = NULL,
       scale = 1, units = c("in", "cm", "mm"),
       dpi = 300)
