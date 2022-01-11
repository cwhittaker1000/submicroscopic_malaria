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
library(tictoc); library(BayesianTools); library(R2jags); library(bayesmix); library(ggplot2);
library(cowplot); library(dplyr)
source("Functions/Submicroscopic_Analysis_Functions.R")
seed <- 193
fresh_run <- FALSE

# Load in the dataset and subset the data by global region the survey was carried out in:
data_frame <- read.csv("Data/SI_Systematic_Review_Results_R_Import.csv")
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data <- full_data %>%
  mutate(PCR_Prev = 100 * (PCR_N_Positive/PCR_N_Tested)) %>%
  mutate(LM_Prev = 100 * (Microscopy_N_Positive/Microscopy_N_Tested)) %>%
  mutate(Prev_Ratio = LM_Prev/PCR_Prev)
Asia_Oceania <- full_data[full_data$Global_Region == "Asia&Oceania" & full_data$Full_Or_Age_Disagg_Data == 2, ]
East_Africa <- full_data[full_data$Global_Region == "East Africa" & full_data$Full_Or_Age_Disagg_Data == 2, ]
South_America <- full_data[full_data$Global_Region == "South America" & full_data$Full_Or_Age_Disagg_Data == 2, ]
West_Africa <- full_data[full_data$Global_Region == "West Africa" & full_data$Full_Or_Age_Disagg_Data == 2, ]

###################################################################################################
##                                                                                               ##
##  Fig 2 & Supp Fig 11: Global Region Disaggregated Data - Running the Log-Linear Regression    ##
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
  Asia_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, Asia_Oceania)
  East_Africa_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, East_Africa)
  South_America_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, South_America)
  West_Africa_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, West_Africa)
  saveRDS(Asia_model, "Outputs/Asia_MCMC_Output.rds")
  saveRDS(East_Africa_model, "Outputs/East_Africa_MCMC_Output.rds")
  saveRDS(South_America_model, "Outputs/South_America_MCMC_Output.rds")
  saveRDS(West_Africa_model, "Outputs/West_Africa_MCMC_Output.rds")
} else {
  Asia_model <- readRDS("Outputs/Asia_MCMC_Output.rds")
  East_Africa_model <- readRDS("Outputs/East_Africa_MCMC_Output.rds")
  South_America_model <- readRDS("Outputs/South_America_MCMC_Output.rds")
  West_Africa_model <- readRDS("Outputs/West_Africa_MCMC_Output.rds")
}

# Supplementary Figure 11 - MCMC Output and Parameter Tables
pdf("Figures/Supplementary/Supp Figure 11 - Global Region MCMC Output/Supp Figure 9 - Asia MCMC Output.pdf", width = 7.33, height = 7.51)
plot(Asia_model, col = c("black"), las = 1)
dev.off()
Asia_param_table <- param_table(Asia_model, params)

pdf("Figures/Supplementary/Supp Figure 11 - Global Region MCMC Output/Supp Figure 9 - East Africa MCMC Output.pdf", width = 7.33, height = 7.51)
plot(East_Africa_model, col = c("green"), las = 1)
dev.off()
East_africa_param_table <- param_table(East_Africa_model, params)

pdf("Figures/Supplementary/Supp Figure 11 - Global Region MCMC Output/Supp Figure 9 - South America MCMC Output.pdf", width = 7.33, height = 7.51)
plot(South_America_model, col = c("blue"), las = 1)
dev.off()
South_America_param_table <- param_table(South_America_model, params)

pdf("Figures/Supplementary/Supp Figure 11 - Global Region MCMC Output/Supp Figure 9 - West Africa MCMC Output.pdf", width = 7.33, height = 7.51)
plot(West_Africa_model, col = c("red"), las = 1)
dev.off()
West_Africa_param_table <- param_table(West_Africa_model, params)

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
Asia_plotting <- data.frame(LM_Prev = 100 * Asia_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_Asia, 
                            lower = 100 * Asia_credible_lower, upper = 100 * Asia_credible_upper,
                            prev_ratio = (Asia_fitted_microscopy/PCR_prevalence_Asia),
                            lower_pr = Asia_prev_ratio_lower, upper_pr = Asia_prev_ratio_upper)

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
East_Africa_plotting <- data.frame(LM_Prev = 100 * East_Africa_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_East_Africa, 
                                   lower = 100 * East_Africa_credible_lower, upper = 100 * East_Africa_credible_upper,
                                   prev_ratio = (East_Africa_fitted_microscopy/PCR_prevalence_East_Africa),
                                   lower_pr = East_Africa_prev_ratio_lower, upper_pr = East_Africa_prev_ratio_upper)


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
South_America_plotting <- data.frame(LM_Prev = 100 * South_America_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_South_America, 
                                    lower = 100 * South_America_credible_lower, upper = 100 * South_America_credible_upper,
                                    prev_ratio = (South_America_fitted_microscopy/PCR_prevalence_South_America),
                                    lower_pr = South_America_prev_ratio_lower, upper_pr = South_America_prev_ratio_upper)

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
West_Africa_plotting <- data.frame(LM_Prev = 100 * West_Africa_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_West_Africa, 
                                   lower = 100 * West_Africa_credible_lower, upper = 100 * West_Africa_credible_upper,
                                   prev_ratio = (West_Africa_fitted_microscopy/PCR_prevalence_West_Africa),
                                   lower_pr = West_Africa_prev_ratio_lower, upper_pr = West_Africa_prev_ratio_upper)

a <- ggplot(Asia_Oceania, aes(x = PCR_Prev, y = LM_Prev)) +
  geom_point(size = 2, col = "#63AD4A") +
  geom_ribbon(data = Asia_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "#63AD4A", linetype = 0) + 
  geom_line(data = Asia_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#63AD4A", size = 1) +
  labs(y = "Microscopy \nPrevalence (%)", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  geom_segment(x = 0, y = 0, xend = 100, yend = 100, colour = "black", size = 0.5, linetype = 2) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.7, 1, 0.7),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(alpha = FALSE, colour = guide_legend(override.aes = list(size = 6)))

b <- ggplot(West_Africa, aes(x = PCR_Prev, y = LM_Prev)) +
  geom_point(size = 2, col = "#5299D3") +
  geom_ribbon(data = West_Africa_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "#5299D3", linetype = 0) + 
  geom_line(data = West_Africa_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#5299D3", size = 1) +
  labs(y = "Microscopy \nPrevalence (%)", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  geom_segment(x = 0, y = 0, xend = 100, yend = 100, colour = "black", size = 0.5, linetype = 2) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.7, 1, 0.7),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(alpha = FALSE, colour = guide_legend(override.aes = list(size = 6)))

c <- ggplot(East_Africa, aes(x = PCR_Prev, y = LM_Prev)) +
  geom_point(size = 2, col = "#DD954D") +
  geom_ribbon(data = East_Africa_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "#DD954D", linetype = 0) + 
  geom_line(data = East_Africa_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#DD954D", size = 1) +
  labs(y = "Microscopy \nPrevalence (%)", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  geom_segment(x = 0, y = 0, xend = 100, yend = 100, colour = "black", size = 0.5, linetype = 2) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.7, 1, 0.7),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(alpha = FALSE, colour = guide_legend(override.aes = list(size = 6)))


d <- ggplot(South_America, aes(x = PCR_Prev, y = LM_Prev)) +
  geom_point(size = 2, col = "#DBC453") +
  geom_ribbon(data = South_America_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "#DBC453", linetype = 0) + 
  geom_line(data = South_America_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#DBC453", size = 1) +
  labs(y = "Microscopy \nPrevalence (%)", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  geom_segment(x = 0, y = 0, xend = 100, yend = 100, colour = "black", size = 0.5, linetype = 2) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.7, 1, 0.7),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(alpha = FALSE, colour = guide_legend(override.aes = list(size = 6)))

e <- ggplot(data = South_America_plotting, aes(x = PCR_Prev, ymin = lower)) +
  geom_ribbon(data = South_America_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#DBC453", linetype = 0) + 
  geom_line(data = South_America_plotting, aes(x = PCR_Prev, y = prev_ratio), col = "#DBC453", size = 1) +
  geom_ribbon(data = East_Africa_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#DD954D", linetype = 0) + 
  geom_line(data = East_Africa_plotting, aes(x = PCR_Prev, y = prev_ratio), col = "#DD954D", size = 1) +
  geom_ribbon(data = West_Africa_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#5299D3", linetype = 0) + 
  geom_line(data = West_Africa_plotting, aes(x = PCR_Prev, y = prev_ratio), col = "#5299D3", size = 1) +
  geom_ribbon(data = Asia_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#63AD4A", linetype = 0) + 
  geom_line(data = Asia_plotting, aes(x = PCR_Prev, y = prev_ratio), col = "#63AD4A", size = 1) +
  labs(y = "Prevalence Ratio", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() +
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.7, 1, 0.7),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  guides(alpha = FALSE, colour = guide_legend(override.aes = list(size = 6)))

top_row <- plot_grid(a, b, labels = c('', ''), label_size = 30, nrow = 1)
middle_set <- plot_grid(c, d, labels = c('', ''), label_size = 30, ncol = 1)
final <- plot_grid(e, labels = c(''), label_size = 30, ncol = 1)
combined <- plot_grid(middle_set, final, label_size = 30, ncol = 2)
plot_grid(top_row, combined, label_size = 30, ncol = 1, rel_heights = c(1, 2))  +
  draw_plot_label(
    c("A", "B", "C", "D", "E"), 
    c(0, 0.5, 0, 0, 0.5),
    c(1, 1, 0.7, 0.38, 0.7),
    size = 30)
ggsave("Figures/Figure 2 - Global Regions/Figure_2.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 10, height = 9, units = c("in", "cm", "mm"),
       dpi = 300, useDingbats = FALSE)

# Statistical Tests Carried Out On The Data
# ANOVA - Testing for Differences in Means
variance <- full_data$PCR_N_Tested * (full_data$PCR_Prev/100) * (1 - full_data$PCR_Prev/100)
stdev <- sqrt(variance)
weighted_global_region_model <- lm(Prev_Ratio ~ Global_Region + PCR_Prev, data = full_data, na.action = na.omit, weights = 1/variance) # similar results with 1/variance
summary(weighted_global_region_model)
ANOVA_object <- aov(weighted_global_region_model)
summary(ANOVA_object)
TukeyHSD(ANOVA_object, which = "Global_Region")

global_region_model <- lm(Prev_Ratio ~ Global_Region + PCR_Prev, data = full_data, na.action = na.omit) # similar results with 1/variance
summary(global_region_model)
ANOVA_object <- aov(global_region_model)
summary(ANOVA_object)
TukeyHSD(ANOVA_object, which = "Global_Region")



