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
library(rjags); library(ssa); library(binom); library(MALDIquant); library(formattable); library(dplyr)
library(tictoc); library(BayesianTools); library(R2jags); library(bayesmix); library(cowplot); library(extrafont)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/Sub-Patent Malarial Infections/Sub_Patent_Malaria_Analysis/")
source("Functions/Submicroscopic_Analysis_Functions.R")
seed <- 193
theme_set(theme_grey())
fresh_run <- FALSE

# Load in the dataset and subset the data into data from the previous review (old_data), data from 
# this review (new_data) all data together (full_data). Note: non-age disaggregated data is assigned coding 2)
data_frame <- read.csv("Data/SI_Systematic_Review_Results_R_Import.csv")
old_data <- data_frame[data_frame$Old_or_New == "Old" & data_frame$Full_Or_Age_Disagg_Data == 2, ] 
new_data <- data_frame[data_frame$Old_or_New == "New" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data <- full_data %>%
  mutate(PCR_Prev = 100 * (PCR_N_Positive/PCR_N_Tested)) %>%
  mutate(LM_Prev = 100 * (Microscopy_N_Positive/Microscopy_N_Tested)) %>%
  mutate(Prev_Ratio = LM_Prev/PCR_Prev)

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
if (fresh_run) {
  full_data_output <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, full_data)
  saveRDS(full_data_output, "Outputs/Full_Data_MCMC_Output.Rds")
} else {
  full_data_output <- readRDS("Outputs/Full_Data_MCMC_Output.Rds")
}

# Supplementary Figure 8 Materials
pdf("Figures/Supplementary/Supp Figure 8 - Full Model MCMC Output/Supp Figure 8 - Full Model MCMC Output.pdf", width = 7.33, height = 7.51)
full_data_param_table <- param_table(full_data_output, params)
plot(full_data_output)
dev.off()

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
geog <- full_data %>%
  group_by(Global_Region) %>%
  summarise(count = n(), total_sample_size = sum(PCR_N_Tested), mean_sample_size = mean(PCR_N_Tested), mean = mean(Prev_Ratio), sd = sd(Prev_Ratio), se = sd/sqrt(n()))

full_data <- full_data %>% 
  left_join(geog, by = "Global_Region") %>%
  mutate(weight = PCR_N_Tested/total_sample_size) %>%
  mutate(weight_for_plotting = PCR_N_Tested/mean_sample_size) %>%
  mutate(weight_for_plotting = ifelse(weight_for_plotting > 3, 3, weight_for_plotting)) %>%
  mutate(weight_for_plotting = ifelse(weight_for_plotting < 1, 1.5, weight_for_plotting))

summarised <- full_data %>% 
  mutate(weighted_prev_ratio = weight * Prev_Ratio) %>%
  group_by(Global_Region) %>%
  summarise(mean_prev_ratio = sum(weighted_prev_ratio))

full_data <- full_data %>% 
  left_join(summarised, by = "Global_Region")

fitted_model <- data.frame(PCR_Prev = PCR_prevalence_full * 100, LM_Prev = full_fitted_microscopy * 100,
                           Prev_Ratio = full_fitted_microscopy/PCR_prevalence_full,
                           lower = 100 * full_data_credible_lower, upper = 100 * full_data_credible_upper,
                           lower_pr = full_data_sensitivity_credible_lower,
                           upper_pr = full_data_sensitivity_credible_upper)

full_data$Global_Region <- factor(full_data$Global_Region, levels = c("South America", "East Africa", "West Africa", "Asia&Oceania")) 

a <- ggplot(full_data, aes(x = PCR_Prev, y = LM_Prev, col = Global_Region)) +
        geom_point(size = 2) +
        geom_ribbon(data = fitted_model, aes(x = PCR_Prev, ymin = lower, ymax = upper, alpha = 0.01, border = NULL), col = "white", fill = "#F7A82A", linetype = 0) + 
        geom_line(data = fitted_model, aes(x = PCR_Prev, y = LM_Prev), col = "#F7A82A", size = 1) +
        labs(y = "Microscopy\nPrevalence (%)", x = "PCR Prevalence (%)") + 
        scale_y_continuous(limits = c(0, 100)) +
        geom_segment(x = 0, y = 0, xend = 100, yend = 100, colour = "black", size = 0.5, linetype = 2) +
        scale_colour_manual(breaks = c("Asia&Oceania", "West Africa", "East Africa", "South America"),
                            labels = c("Asia & Oceania", "West Africa", "East Africa", "South America"),
                            values = c("#DBC453", "#DD954D", "#5299D3", "#63AD4A")) +
        theme_bw() +
        theme(legend.position = c(0.13, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
              axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
              plot.margin = unit(c(0.5, 0.7, 1, 0.7),"cm"), legend.title = element_blank(),
              legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
              legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        guides(alpha = FALSE, colour = guide_legend(override.aes = list(size = 6)))

b <- ggplot(full_data, aes(x = PCR_Prev, y = Prev_Ratio)) +
        geom_point(size = 2, aes(col = full_data$Global_Region)) +
        geom_line(data = fitted_model, aes(x = PCR_Prev, y = Prev_Ratio), col = "#F7A82A", size = 1) +
        geom_ribbon(data = fitted_model, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, alpha = 0.01, border = NULL), fill = "#F7A82A", linetype = 0) +
        labs(x = "PCR Prevalence (%)", y = "Prevalence Ratio") +
        scale_color_manual(values = c("#DBC453", "#DD954D", "#5299D3", "#63AD4A")) +
        scale_y_continuous(limits = c(0, 1)) +
        theme_bw() +
        theme(legend.position = "none", axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
              axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
              plot.margin = unit(c(0.5,0.32,1.5,0.5), "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

full_data$Global_Region <- factor(full_data$Global_Region, levels = c("South America", "East Africa", "West Africa", "Asia&Oceania")) 

c <- ggplot(full_data, aes(x = Global_Region, y = Prev_Ratio, col = Global_Region)) + 
        geom_boxplot(outlier.alpha = 0) +
        labs(x = "", y = "Prevalence Ratio") + 
        coord_flip() +
        geom_jitter(shape = 16, position = position_jitter(0.2), size = full_data$weight_for_plotting, alpha = 1) + 
        scale_x_discrete(labels = c("Asia&Oceania" = paste0("Asia & Oceania\n n = (", geog$count[1], ")"), 
                                    "West Africa" = paste0("West Africa\n n = (", geog$count[4], ")"),
                                    "East Africa" = paste0("East Africa\n n = (", geog$count[2], ")"),
                                    "South America" = paste0("South America\n n = (", geog$count[3], ")")))  +
        scale_color_manual(values = c("#DBC453", "#DD954D", "#5299D3", "#63AD4A")) +
        theme_bw() + 
        theme(legend.position = "none", axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face = "bold"),
              axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +5),
              plot.margin = unit(c(0.5,1,1.5,0),"cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
        geom_errorbar(aes(ymin = mean_prev_ratio, ymax = mean_prev_ratio, size = 3), linetype = 1, width = 0.75)


bottom_row <- plot_grid(b, c, labels = c('', ''), label_size = 30, rel_widths = c(1, 1.5))
plot_grid(a, bottom_row, labels = c('', ''), label_size = 30, ncol = 1) +
  draw_plot_label(c("A", "B", "C"), 
                  c(0, 0, 0.385),
                  c(1, 0.53, 0.53), 
                  size = 30)
ggsave("Figures/Figure 1 - Overall/Figure_1.pdf", plot = last_plot(), device = NULL, path = NULL,
       scale = 1, width = 9.35, height = 9, units = c("in", "cm", "mm"),
       dpi = 300)

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
if (fresh_run) {
  old_data_output <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, old_data)
  new_data_output <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, new_data)
  saveRDS(old_data_output, "Outputs/Old_Data_MCMC_Output.Rds")
  saveRDS(new_data_output, "Outputs/New_Data_MCMC_Output.Rds")
} else {
  old_data_output <- readRDS("Outputs/Old_Data_MCMC_Output.Rds")
  new_data_output <- readRDS("Outputs/New_Data_MCMC_Output.Rds")
}

old_chain_overall <- rbind(old_data_output[[1]], old_data_output[[2]], old_data_output[[3]], old_data_output[[4]])
old_beta_mean <- mean(as.array(old_chain_overall[, "beta"])) 
old_delt_mean <- mean(as.array(old_chain_overall[, "delt"])) 
old_fitted_microscopy <- mean_output(PCR_prevalence_full, old_beta_mean, old_delt_mean)
old_credibles <- credible_intervals(PCR_prevalence_full, old_chain_overall)
old_data_credible_lower <- old_credibles$credible_lower
old_data_credible_upper <- old_credibles$credible_upper
old_data_sensitivity_credible_lower <- old_credibles$sensitivity_lower
old_data_sensitivity_credible_upper <- old_credibles$sensitivity_upper

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
pdf("Figures/Supplementary/Supp Figure 1 - Old vs New/Supp Figure 1 - Old vs New.pdf", width = 10.32, height = 6.38)
par(mfrow = c(1, 2)) 
xlim <- c(0, 100)
ylim <- c(0, 100)
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
legend("topleft", legend = c("Old Data", "New Data"), col = c("#6BC24E", "#C6D3F3"), pch = 20, pt.cex = 2, cex = 1)

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
       pch = 20, pt.cex = 2, cex = 1, bg = "white")
dev.off()

# Supplementary Figure 2 Plotting - Predicted vs Observed Microscopy Prevalence
pdf("Figures/Supplementary/Supp Figure 2 - Empirical vs Modelled/Supp Figure 2 - Empirical vs Modelled.pdf", width = 7.51, height = 7.33)
observed_PCR_prevalences_full <- full_data$PCR_Prev
matched_PCR_indices_full <-  match.closest(observed_PCR_prevalences_full, PCR_prevalence_full*100)
logit_matched_PCR <- logit(PCR_prevalence_full[matched_PCR_indices_full])
predicted_LM_logit <- (full_delt_mean - (full_beta_mean * mean(logit_matched_PCR))) + ((1 + full_beta_mean)  * logit_matched_PCR)
predicted_LM <- expit(predicted_LM_logit) 
cor(full_data$Micro_Prev, predicted_LM)

observed_PCR_prevalences_old <- old_data$PCR_Prev
matched_PCR_indices_old <-  match.closest(observed_PCR_prevalences_old, PCR_prevalence_full)
logit_matched_PCR_old <- logit(PCR_prevalence_full[matched_PCR_indices_old])
predicted_LM_logit_old <- (full_delt_mean - (full_beta_mean * mean(logit_matched_PCR_old))) + ((1 + full_beta_mean)  * logit_matched_PCR_old)
predicted_LM_old <- expit(predicted_LM_logit_old) 
cor(old_data$Micro_Prev, predicted_LM_old)

observed_PCR_prevalences_new <- new_data$PCR_Prev
matched_PCR_indices_new <-  match.closest(observed_PCR_prevalences_new, PCR_prevalence_full)
logit_matched_PCR_new <- logit(PCR_prevalence_full[matched_PCR_indices_new])
predicted_LM_logit_new <- (full_delt_mean - (full_beta_mean * mean(logit_matched_PCR_new))) + ((1 + full_beta_mean)  * logit_matched_PCR_new)
predicted_LM_new <- expit(predicted_LM_logit_new) 
cor(new_data$Micro_Prev, predicted_LM_new)

par(mfrow = c(1, 1))
plot(old_data$Micro_Prev*100, predicted_LM_old*100, xlim = c(0, 100), ylim = c(0, 100), col = "#C6D3F3", pch = 1, lwd = 2, ylab = "Predicted PCR Prevalence (%)", xlab = "Observed PCR Prevalence (%)", las = 1)
points(new_data$Micro_Prev*100, predicted_LM_new*100, xlim = c(0, 100), ylim = c(0, 100), pch = 1, lwd = 2, col = "#6BC24E")
lines(seq(0,100,0.01), seq(0,100,0.01), lwd = 2, lty = 2)
cor(full_data$Micro_Prev, predicted_LM)
dev.off()

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

