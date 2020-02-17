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
library(tictoc); library(BayesianTools); library(R2jags); library(bayesmix); library(ggplot2); library(cowplot);
library(magrittr); library(dplyr)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/Sub-Patent Malarial Infections/Sub_Patent_Malaria_Analysis/")
source("Functions/Submicroscopic_Analysis_Functions.R")
seed <- 193
fresh_run <- FALSE

# Load in the dataset and subset the data into children aged 0-5 years (young_children), old children, 
# aged 5-15 years (old_children) and data from adults, aged 15+ years (adults)
data_frame <- read.csv("Data/SI_Systematic_Review_Results_R_Import.csv")
data_frame <- data_frame %>%
  mutate(PCR_Prev = 100 * (PCR_N_Positive/PCR_N_Tested)) %>%
  mutate(LM_Prev = 100 * (Microscopy_N_Positive/Microscopy_N_Tested)) %>%
  mutate(Prev_Ratio = LM_Prev/PCR_Prev)
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
if (fresh_run == TRUE) {
  young_child_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, young_children)
  old_child_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, old_children)
  adult_model <- run_rJAGS_model(10000, 4, model_file, params, initial_values_function, adults)
  saveRDS(young_child_model, "Outputs/Young_Children_MCMC_Output.rds")
  saveRDS(old_child_model, "Outputs/Old_Children_MCMC_Output.rds")
  saveRDS(adult_model, "Outputs/Adults_MCMC_Output.rds")
} else {
  young_child_model <- readRDS("Outputs/Young_Children_MCMC_Output.rds")
  old_child_model <- readRDS("Outputs/Old_Children_MCMC_Output.rds")
  adult_model <- readRDS("Outputs/Adults_MCMC_Output.rds")
}

# Supplementary Figure 11 - MCMC Output and Parameter Tables
pdf("Figures/Supplementary/Supp Figure 11 - Age Stratified MCMC Output/Supp Figure 11 - Young Children MCMC Output.pdf", width = 7.33, height = 7.51)
young_children_param_table <- param_table(young_child_model, params)
plot(young_child_model, col = c("yellow"), las = 1)
dev.off()

pdf("Figures/Supplementary/Supp Figure 11 - Age Stratified MCMC Output/Supp Figure 11 - Old Children MCMC Output.pdf", width = 7.33, height = 7.51)
old_children_param_table <- param_table(old_child_model, params)
plot(old_child_model, col = c("orange"), las = 1)
dev.off()

pdf("Figures/Supplementary/Supp Figure 11 - Age Stratified MCMC Output/Supp Figure 11 - Adult MCMC Output.pdf", width = 7.33, height = 7.51)
adult_param_table <- param_table(adult_model, params)
plot(adult_model, col = c("red"), las = 1)
dev.off()

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
Young_plotting <- data.frame(LM_Prev = 100 * young_child_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_young, 
                             lower = 100 * young_children_credible_lower, upper = 100 * young_children_credible_upper,
                             Prev_Ratio = (young_child_fitted_microscopy/PCR_prevalence_young),
                             lower_pr = young_prev_ratio_lower, upper_pr = young_prev_ratio_upper)

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
Old_plotting <- data.frame(LM_Prev = 100 * old_child_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_old, 
                           lower = 100 * old_children_credible_lower, upper = 100 * old_children_credible_upper,
                           Prev_Ratio = (old_child_fitted_microscopy/PCR_prevalence_old),
                           lower_pr = old_prev_ratio_lower, upper_pr = old_prev_ratio_upper)

adult <- rbind(adult_model[[1]], adult_model[[2]], adult_model[[3]], adult_model[[4]])
adult_beta_mean <- mean(as.array(adult[, "beta"]))
adult_delt_mean <- mean(as.array(adult[, "delt"]))
PCR_prevalence_adults <- seq(0.01, 0.8, 0.001)
adult_fitted_microscopy <- mean_output(PCR_prevalence_adults, adult_beta_mean, adult_delt_mean)
adult_credibles <- credible_intervals(PCR_prevalence_adults, adult)
adult_credible_lower <- adult_credibles$credible_lower
adult_credible_upper <- adult_credibles$credible_upper
adult_prev_ratio_lower <- adult_credibles$sensitivity_lower
adult_prev_ratio_upper <- adult_credibles$sensitivity_upper
Adult_plotting <- data.frame(LM_Prev = 100 * adult_fitted_microscopy, PCR_Prev = 100 * PCR_prevalence_adults, 
                             lower = 100 * adult_credible_lower, upper = 100 * adult_credible_upper,
                             Prev_Ratio = (adult_fitted_microscopy/PCR_prevalence_adults),
                             lower_pr = adult_prev_ratio_lower, upper_pr = adult_prev_ratio_upper)

a <- ggplot(young_children, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "#EAE600") +
  geom_ribbon(data = Young_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#EAE600", linetype = 0) + 
  geom_line(data = Young_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "#EAE600", size = 2) +
  labs(y = "Prevalence Ratio", x = "PCR Prevalence (%)") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.15, 1, 0.7),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

b <- ggplot(old_children, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "orange") +
  geom_ribbon(data = Old_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "orange", linetype = 0) + 
  geom_line(data = Old_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "orange", size = 2) +
  labs(y = "", x = "PCR Prevalence (%)") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 0.2, 1, 0),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

c <- ggplot(adults, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "red") +
  geom_ribbon(data = Adult_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "red", linetype = 0) + 
  geom_line(data = Adult_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "red", size = 2) +
  labs(y = "", x = "PCR Prevalence (%)") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(0.5, 1, 1, 0),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

d <- ggplot(adults, aes(x = PCR_Prev, y = LM_Prev)) +
  geom_point(size = 2, col = "red") +
  geom_point(data = young_children, aes(x = PCR_Prev, y = LM_Prev), col = "orange", size = 2) +
  geom_point(data = old_children, aes(x = PCR_Prev, y = LM_Prev), col = "#EAE600", size = 2) +
  geom_ribbon(data = Adult_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "red", linetype = 0) +
  geom_ribbon(data = Old_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "orange", linetype = 0) + 
  geom_ribbon(data = Young_plotting, aes(x = PCR_Prev, ymin = lower, ymax = upper, border = NULL), alpha = 0.2, col = "white", fill = "yellow", linetype = 0) +
  geom_segment(x = 0, y = 0, xend = 100, yend = 100, colour = "black", size = 0.5, linetype = 2) +
  geom_line(data = Adult_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "red", size = 2) +
  geom_line(data = Old_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "orange", size = 2) +
  geom_line(data = Young_plotting, aes(x = PCR_Prev, y = LM_Prev), col = "#EAE600", size = 2) +
  labs(y = "Microscopy\nPrevalence (%)", x = "PCR Prevalence (%)") + 
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(1, 1.5, 0.5, 1),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

combined <- plot_grid(a, b, c, label_size = 30, ncol = 3, rel_widths = c(1.30, 1.05, 1.12))
combined_two <- plot_grid(d, combined, nrow = 2, rel_heights = c(1.2, 0.9))
combined_two +
  draw_plot_label(
    c("A", "B"), 
    c(0, 0),
    c(1, 0.48),
    size = 30) 
ggsave("Figures/Figure 4 - Age/Figure_4.pdf", width = 9.03, height = 7,
       plot = last_plot(), device = NULL, path = NULL, scale = 1, units = c("in", "cm", "mm"), dpi = 300)


age_data <- data_frame[data_frame$Age_Group != "All Ages", ]

age_counts <- age_data %>%
  group_by(Age_Group) %>%
  summarise(count = n(), total_sample_size = sum(PCR_N_Tested), mean_sample_size = mean(PCR_N_Tested), mean = mean(Prev_Ratio), sd = sd(Prev_Ratio), se = sd/sqrt(n()))

age_data <- age_data %>%
  left_join(age_counts, by = "Age_Group") %>%
  mutate(weight = PCR_N_Tested/total_sample_size) %>%
  mutate(weight_for_plotting = PCR_N_Tested/mean_sample_size) %>%
  mutate(weight_for_plotting = ifelse(weight_for_plotting > 3, 3, weight_for_plotting)) %>%
  mutate(weight_for_plotting = ifelse(weight_for_plotting < 1, 1, weight_for_plotting))

summarised <- age_data %>% 
  mutate(weighted_prev_ratio = weight * Prev_Ratio) %>%
  group_by(Age_Group) %>%
  summarise(mean_prev_ratio = sum(weighted_prev_ratio))

age_data <- age_data %>% 
  left_join(summarised, by = "Age_Group")

age_data$Age_Group <- factor(age_data$Age_Group, levels = c("0-5", "5-15years", "15+")) 

a <- ggplot(young_children, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "#CD99EF") +
  geom_ribbon(data = Young_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#CD99EF", linetype = 0) + 
  geom_line(data = Young_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "#CD99EF", size = 2) +
  labs(y = "Prevalence Ratio", x = "") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 15,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(1.5, 0.15, 2.5, 0.7),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

b <- ggplot(old_children, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "#DD7596") +
  geom_ribbon(data = Old_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#DD7596", linetype = 0) + 
  geom_line(data = Old_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "#DD7596", size = 2) +
  labs(y = "", x = "") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(1.5, 0.15, 2.5, 0.05),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

c <- ggplot(adults, aes(x = PCR_Prev, y = Prev_Ratio)) +
  geom_point(size = 2, col = "#6A94B7") +
  geom_ribbon(data = Adult_plotting, aes(x = PCR_Prev, ymin = lower_pr, ymax = upper_pr, border = NULL), alpha = 0.2, col = "white", fill = "#6A94B7", linetype = 0) + 
  geom_line(data = Adult_plotting, aes(x = PCR_Prev, y = Prev_Ratio), col = "#6A94B7", size = 2) +
  labs(y = "", x = "") + 
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 100)) +
  theme_bw() + 
  theme(legend.position = c(0.15, 0.7), axis.text = element_text(size = 12, face = "bold"), axis.title = element_text(size = 12,face="bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.title.y = element_text(size = 15, vjust = +6),
        plot.margin = unit(c(1.5, 0.25, 2.5, 0.05),"cm"), legend.title = element_blank(),
        legend.text = element_text(size = 14), legend.key = element_rect(fill = NA, color = NA),
        legend.spacing.x = unit(0.2, 'cm'), legend.key.size = unit(1, "cm"), legend.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())


d <- ggplot(age_data, aes(x = Age_Group, y = Prev_Ratio, col = Age_Group)) + 
  geom_boxplot(outlier.alpha = 0) +
  labs(x = "", y = "Prevalence Ratio") + 
  geom_jitter(shape = 16, position = position_jitter(0.2), size = age_data$weight_for_plotting, alpha = 1) + 
  scale_x_discrete(labels = c("0-5" = paste0("0-5 Years\nn = (", age_counts$count[1], ")"), 
                              "5-15years" = paste0("5-15 Years\nn = (", age_counts$count[3], ")"),
                              "15+" = paste0("Adults\nn = (", age_counts$count[2], ")")))  +
  scale_color_manual(values = c("#CD99EF", "#DD7596", "#6A94B7")) +
  theme_bw() + 
  theme(legend.position = "none", axis.text = element_text(size = 15, face = "bold"), axis.title = element_text(size = 15,face = "bold"),
        axis.title.x = element_text(size = 15, vjust = -3), axis.text.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 15, vjust = +3),
        plot.margin = unit(c(1,1,1,1),"cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  
  geom_errorbar(aes(ymin = mean_prev_ratio, ymax = mean_prev_ratio, size = 3), linetype = 1, width = 0.75)

combined <- plot_grid(a, b, c, label_size = 30, ncol = 3, rel_widths = c(1.30, 1.05, 1.05))
combined_two <- plot_grid(d, combined, nrow = 1, rel_widths = c(0.7, 1.3))
combined_two +
  draw_plot_label(
    c("A", "B", "PCR Prevalence (%)"), 
    c(0, 0.33, 0.57),
    c(1, 1, 0.22),
    size = c(30, 30, 16)) 
ggsave("Figures/Figure 4 - Age/alternate_Figure_4.pdf", width = 12.58, height = 4.56,
       plot = last_plot(), device = NULL, path = NULL, scale = 1, units = c("in", "cm", "mm"), dpi = 300)


# Statistical Tests Carried Out On The Data
# ANOVA - Testing for Differences in Means
data_frame_ANOVA <- data_frame[(data_frame$Age_Group == "0-5") | (data_frame$Age_Group ==  "5-15years") | (data_frame$Age_Group ==  "15+"), ] # all age-disaggregated data
data_frame_ANOVA$prev_ratio <- (data_frame_ANOVA$Microscopy_N_Positive/data_frame_ANOVA$Microscopy_N_Tested)/
                                (data_frame_ANOVA$PCR_N_Positive/data_frame_ANOVA$PCR_N_Tested)
ANOVA_object <- aov(prev_ratio ~ Age_Group, data = data_frame_ANOVA)
summary(ANOVA_object)
TukeyHSD(ANOVA_object)

