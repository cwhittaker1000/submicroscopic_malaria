# Whittaker et al., 2019 : Variation in the Prevalence of Submicroscopic Malaria Infections: Historical Transmission 
#                          Intensity and Age as Key Determinants

# This paper represents an update on the systematic reviews of submicroscopic malaria infections  published
# by Okell et al in 2009 (https://academic.oup.com/jid/article/200/10/1509/879741) and by Okell et al in 2012
# (https://www.nature.com/articles/ncomms2241). In it, statistical analyses of both data from these previous 
# reviews and new data collected during the updating process are undertaken and results displayed graphically.

# The code below is responsible for the analyses and plotting that produced Figure 1 of the paper. Any questions, 
# queries, comments, or mistakes, please feel free to get in touch at charles.whittaker16@imperial.ac.uk :) 

# Analyses Responsible for Producing Figure 2 In the Paper

# Access and load various required packages for the analyses
library(rjags); library(moments); library(gtools); library(nortest); library(ggplot2); library(LMest); library(ssa); library(binom);
library(tidyverse)

# Load in the dataset and subset the data into:
#        1) Data from young children, aged 0-5 years (young_children)
#        2) Data from old children, aged 5-15 years (old_children)
#        3) Data from adults, aged 15+ years (adults)
data_frame <- Whittaker.et.al.R.Import
young_children <- data_frame[data_frame$Age_Group == "0-5" , ]
old_children <- data_frame[data_frame$Age_Group ==  "5-15years", ]
adults <- data_frame[data_frame$Age_Group == "15+", ]

# Specifying the Data in the Format Required for Input into RJAGS
young_child_microscopy_PCR_comparison <- list(prev_pcr = young_children$PCR_N_Positive, ## number positive by PCR,
                                              prev_microscopy = young_children$Microscopy_N_Positive, ## number positive by microscopy,
                                              total_pcr = young_children$PCR_N_Tested, ## number tested by PCR,
                                              total_microscopy = young_children$Microscopy_N_Tested, ## number tested by microscopy,
                                              N = 40)                                # Total sample data overall

old_child_microscopy_PCR_comparison <- list(prev_pcr = old_children$PCR_N_Positive, ## number positive by PCR,
                                            prev_microscopy = old_children$Microscopy_N_Positive, ## number positive by microscopy,
                                            total_pcr = old_children$PCR_N_Tested, ## number tested by PCR,
                                            total_microscopy = old_children$Microscopy_N_Tested, ## number tested by microscopy,
                                            N = 37 )                                # Total sample data overall

adult_microscopy_PCR_comparison <- list(prev_pcr = adults$PCR_N_Positive, ## number positive by PCR,
                                        prev_microscopy = adults$Microscopy_N_Positive,## number positive by microscopy,
                                        total_pcr = adults$PCR_N_Tested, ## number tested by PCR,
                                        total_microscopy = adults$Microscopy_N_Tested, ## number tested by microscopy,
                                        N = 41)                           


# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
params <- c("delt", "beta", "taud", "gamma")
inits <- list(beta=0.0001,delt=0.0001,taud=0.0001)


# Specifying and initialising the RJAGS model- creates a JAGS model object
young_children_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                                 data = young_child_microscopy_PCR_comparison, n.chains = 1, inits = inits)

old_children_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                               data = old_child_microscopy_PCR_comparison, n.chains = 1, inits = inits)

adult_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # NEW DATA
                                        data = adult_microscopy_PCR_comparison, n.chains = 1, inits = inits)

# Running, updating and iterating the RJAGS model

# YOUNG CHILD DATA
update(young_children_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
young_children_micr_PCR_comp_model <- coda.samples(young_children_micr_PCR_comp_model, params, 
                                                   n.iter = 100000, thin = 70) # Model updating
# Exploring MCMC output
summary(young_children_micr_PCR_comp_model)
plot(young_children_micr_PCR_comp_model)
autocorr.plot(young_children_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
young_child_beta_mean <- mean(as.array(young_children_micr_PCR_comp_model[, 1]))
young_child_delt_mean <- mean(as.array(young_children_micr_PCR_comp_model[, 2]))

# OLD CHILD DATA
update(old_children_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
old_children_micr_PCR_comp_model <- coda.samples(old_children_micr_PCR_comp_model, params, 
                                                 n.iter = 100000, thin = 70) # Model updating
# Exploring MCMC output
summary(old_children_micr_PCR_comp_model)
plot(old_children_micr_PCR_comp_model)
autocorr.plot(old_children_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
old_child_beta_mean <- mean(as.array(old_children_micr_PCR_comp_model[, 1]))
old_child_delt_mean <- mean(as.array(old_children_micr_PCR_comp_model[, 2]))

# ADULT DATA
update(adult_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
adult_micr_PCR_comp_model <- coda.samples(adult_micr_PCR_comp_model, params, 
                                          n.iter = 100000, thin = 70) # Model updating
# Exploring MCMC output
summary(adult_micr_PCR_comp_model)
plot(adult_micr_PCR_comp_model)
autocorr.plot(adult_micr_PCR_comp_model)

# Saving the MCMC estimates of parameter mean values
adult_beta_mean <- mean(as.array(adult_micr_PCR_comp_model[, 1]))
adult_delt_mean <- mean(as.array(adult_micr_PCR_comp_model[, 2]))

# Calculating the Predicted Values of Microscopy Prevalence Given a PCR Prevalence
  # Equation: logit(Microscopy Prevalence) = delta' + (1 + beta) * logit(PCR_prevalence)
PCR_prevalence_young <- seq(0.005,0.9,0.001)
logit_PCR_prevalence_young <- logit(PCR_prevalence_young)
young_child_fitted_logit_microscopy <- young_child_delt_mean + (1 + young_child_beta_mean)  * logit_PCR_prevalence_young
young_child_fitted_microscopy <- expit(young_child_fitted_logit_microscopy)

PCR_prevalence_old <- seq(0.04, 0.97, 0.001)
logit_PCR_prevalence_old <- logit(PCR_prevalence_old)
old_child_fitted_logit_microscopy <- old_child_delt_mean + (1 + old_child_beta_mean)  * logit_PCR_prevalence_old
old_child_fitted_microscopy <- expit(old_child_fitted_logit_microscopy)

PCR_prevalence_adults <- seq(0.02, 0.8, 0.001)
logit_PCR_prevalence_adults <- logit(PCR_prevalence_adults)
adult_fitted_logit_microscopy <- adult_delt_mean + (1 + adult_beta_mean) * logit_PCR_prevalence_adults
adult_fitted_microscopy <- expit(adult_fitted_logit_microscopy)

# 95% Credible Interval Calculation for Prevalence and Sensitivity
  # Young Children
young_children_chains <- young_children_micr_PCR_comp_model[[1]]
young_children_pred_mean_dist <- matrix(NA, nrow = nrow(young_children_chains), ncol = length(PCR_prevalence_young))
for (i in 1:nrow(young_children_pred_mean_dist)){
  young_children_value_logit_scale <- young_children_chains[i, "delt"] + 
    (1 + young_children_chains[i, "beta"]) * logit(PCR_prevalence_young)
  young_children_pred_mean_dist[i, ] <- expit(young_children_value_logit_scale)
}
young_children_credible_lower <- apply(young_children_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
young_children_credible_upper <- apply(young_children_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
young_data_sensitivity_lower <- young_children_credible_lower / PCR_prevalence_young
young_data_sensitivity_upper <- young_children_credible_upper / PCR_prevalence_young

  # Old Children
old_children_chains <- old_children_micr_PCR_comp_model[[1]]
old_children_pred_mean_dist <- matrix(NA, nrow = nrow(old_children_chains), ncol = length(PCR_prevalence_old))
for (i in 1:nrow(old_children_pred_mean_dist)){
  old_children_value_logit_scale <- old_children_chains[i, "delt"] + 
    (1 + old_children_chains[i, "beta"]) * logit(PCR_prevalence_old)
  old_children_pred_mean_dist[i, ] <- expit(old_children_value_logit_scale)
}
old_children_credible_lower <- apply(old_children_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
old_children_credible_upper <- apply(old_children_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
old_data_sensitivity_lower <- old_children_credible_lower / PCR_prevalence_old
old_data_sensitivity_upper <- old_children_credible_upper / PCR_prevalence_old

  # Adults
adults_chains <- adult_micr_PCR_comp_model[[1]]
adults_pred_mean_dist <- matrix(NA, nrow = nrow(adults_chains), ncol = length(PCR_prevalence_adults))
for (i in 1:nrow(adults_pred_mean_dist)){
  adults_value_logit_scale <- adults_chains[i, "delt"] + 
    (1 + adults_chains[i, "beta"]) * logit(PCR_prevalence_adults)
  adults_pred_mean_dist[i, ] <- expit(adults_value_logit_scale)
}
adults_credible_lower <- apply(adults_pred_mean_dist, MARGIN = 2, quantile, prob = 0.025)
adults_credible_upper <- apply(adults_pred_mean_dist, MARGIN = 2, quantile, prob = 0.975)
adult_data_sensitivity_lower <- adults_credible_lower / PCR_prevalence_adults
adult_data_sensitivity_upper <- adults_credible_upper / PCR_prevalence_adults

# Figure 2A Plotting - Microscopy Prevalence Against PCR Prevalence for All 3 Age Groups- Data & Modelled Relationship
plot(young_children$PCR_Prev, young_children$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#EAE600",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(old_children$PCR_Prev, old_children$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "orange")
points(adults$PCR_Prev, adults$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "red")
lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)
lines(PCR_prevalence_young, young_child_fitted_microscopy, col = "#EAE600", lwd = 3)
lines(PCR_prevalence_old, old_child_fitted_microscopy, col = "orange", lwd = 3)
lines(PCR_prevalence_adults, adult_fitted_microscopy, col = "red", lwd = 3)
polygon(x = c(PCR_prevalence_young, rev(PCR_prevalence_young)), 
        y = c(young_children_credible_upper, rev(young_children_credible_lower)), 
        col = adjustcolor("#EAE600", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_old, rev(PCR_prevalence_old)), 
        y = c(old_children_credible_upper, rev(old_children_credible_lower)), 
        col = adjustcolor("orange", alpha.f = 0.5), border = NA)
polygon(x = c(PCR_prevalence_adults, rev(PCR_prevalence_adults)), 
        y = c(adults_credible_upper, rev(adults_credible_lower)), 
        col = adjustcolor("red", alpha.f = 0.5), border = NA)
legend("topleft", legend = c("Young Children", "Older Children", "Adults"), col = c("#EAE600", "orange", "red"),
       lty = 1, lwd = 2, pt.cex = 2, cex = 0.8)

# Figure 2B Plotting - Microscopy Sensitivity Against PCR Prevalence for YOUNG CHILDREN - Raw Data & Modelled Relationship
plot(young_children$PCR_Prev, young_children$Micro_Prev/young_children$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "#EAE600",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(PCR_prevalence_young, young_child_fitted_microscopy/PCR_prevalence_young, col = "#EAE600", lwd = 3)
polygon(x = c(PCR_prevalence_young, rev(PCR_prevalence_young)), 
        y = c(young_data_sensitivity_upper, rev(young_data_sensitivity_lower)), 
        col = adjustcolor("#EAE600", alpha.f = 0.5), border = NA)

# Figure 2C Plotting - Microscopy Sensitivity Against PCR Prevalence for OLD CHILDREN - Raw Data & Modelled Relationship
plot(old_children$PCR_Prev, old_children$Micro_Prev/old_children$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "orange",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(PCR_prevalence_old, old_child_fitted_microscopy/PCR_prevalence_old, col = "orange", lwd = 3)
polygon(x = c(PCR_prevalence_old, rev(PCR_prevalence_old)), 
        y = c(old_data_sensitivity_upper, rev(old_data_sensitivity_lower)), 
        col = adjustcolor("orange", alpha.f = 0.5), border = NA)

# Figure 2D Plotting - Microscopy Sensitivity Against PCR Prevalence for ADULTS - Raw Data & Modelled Relationship
plot(adults$PCR_Prev, adults$Micro_Prev/adults$PCR_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "red",
     xlab = "PCR Prevalence (%)", ylab = "Sensitivity of Microscopy (%)")
lines(PCR_prevalence_adults, adult_fitted_microscopy/PCR_prevalence_adults, col = "red", lwd = 3)
polygon(x = c(PCR_prevalence_adults, rev(PCR_prevalence_adults)), 
        y = c(adult_data_sensitivity_upper, rev(adult_data_sensitivity_lower)), 
        col = adjustcolor("red", alpha.f = 0.5), border = NA)

## Possible Graphs ## 

# Figure 2E Plotting - Empirical Mean Microscopy Sensitivities - Raw Data & Confidence Intervals
young_children$Sensitivity <- (young_children$Microscopy_N_Positive/young_children$Microscopy_N_Tested)/(young_children$PCR_N_Positive/young_children$PCR_N_Tested)
young_mean_sens <- mean(young_children$Sensitivity)
young_mean_se <- sd(young_children$Sensitivity)/sqrt(length(young_children$Sensitivity))

old_children$Sensitivity <- (old_children$Microscopy_N_Positive/old_children$Microscopy_N_Tested)/(old_children$PCR_N_Positive/old_children$PCR_N_Tested)
old_mean_sens <- mean(old_children$Sensitivity)
old_mean_se <- sd(old_children$Sensitivity)/sqrt(length(old_children$Sensitivity))

adults$Sensitivity <- (adults$Microscopy_N_Positive/adults$Microscopy_N_Tested)/(adults$PCR_N_Positive/adults$PCR_N_Tested)
adult_mean_sens <- mean(adults$Sensitivity)
adult_mean_se <- sd(adults$Sensitivity)/sqrt(length(adults$Sensitivity))


plot(0, 0, ylim = c(0, 1), xlim = c(0, 3), cex = 0, xlab = "Age Group", ylab = "Mean Sensitivity", xaxt = "n")
rect(xleft = 0, ybottom = 0, xright = 1, ytop = young_mean_sens, col = "#EAE600")
arrows(x0 = 0.5, y0 = young_mean_sens - 1.96 * young_mean_se, 
       x1 = 0.5, y1 = young_mean_sens + 1.96 * young_mean_se, 
       col=1, angle=90, code=3, length = 0.05)
rect(xleft = 1, ybottom = 0, xright = 2, ytop = old_mean_sens, col = "orange")
arrows(x0 = 1.5, y0 = old_mean_sens - 1.96 * old_mean_se, 
       x1 = 1.5, y1 = old_mean_sens + 1.96 * old_mean_se, 
       col=1, angle=90, code=3, length = 0.05)
rect(xleft = 2, ybottom = 0, xright = 3, ytop = adult_mean_sens, col = "red")
arrows(x0 = 2.5, y0 = adult_mean_sens - 1.96 * adult_mean_se, 
       x1 = 2.5, y1 = adult_mean_sens + 1.96 * adult_mean_se, 
       col=1, angle=90, code=3, length = 0.05)
axis(1, at = c(0.5, 1.5, 2.5), labels = c("Young Children", "Old Children", "Adults"))

# t-tests
t.test(young_children$Sensitivity, old_children$Sensitivity)
t.test(adults$Sensitivity, old_children$Sensitivity)
t.test(young_children$Sensitivity, adults$Sensitivity)


# Figure 2F Plotting - Difference in Mean Sensitivity, Pairwise Comparisons - Modelled Relationships
Fig2E_PCR_prevalence_young <- seq(0.04, 0.8, 0.001)
Fig2E_logit_PCR_prevalence_young <- logit(Fig2E_PCR_prevalence_young)
Fig2E_young_child_fitted_logit_microscopy <- young_child_delt_mean + (1 + young_child_beta_mean)  * Fig2E_logit_PCR_prevalence_young
Fig2E_young_child_fitted_microscopy <- expit(Fig2E_young_child_fitted_logit_microscopy)

Fig2E_PCR_prevalence_old <- seq(0.04, 0.8, 0.001)
Fig2E_logit_PCR_prevalence_old <- logit(Fig2E_PCR_prevalence_old)
Fig2E_old_child_fitted_logit_microscopy <- old_child_delt_mean + (1 + old_child_beta_mean)  * Fig2E_logit_PCR_prevalence_old
Fig2E_old_child_fitted_microscopy <- expit(Fig2E_old_child_fitted_logit_microscopy)

Fig2E_PCR_prevalence_adults <- seq(0.04, 0.8, 0.001)
Fig2E_logit_PCR_prevalence_adults <- logit(Fig2E_PCR_prevalence_adults)
Fig2E_adult_fitted_logit_microscopy <- adult_delt_mean + (1 + adult_beta_mean) * Fig2E_logit_PCR_prevalence_adults
Fig2E_adult_fitted_microscopy <- expit(Fig2E_adult_fitted_logit_microscopy)

Fig2E_adult_young_diff <- (Fig2E_young_child_fitted_microscopy/Fig2E_PCR_prevalence_young) - (Fig2E_adult_fitted_microscopy/Fig2E_PCR_prevalence_adults)
Fig2E_adult_old_diff <- (Fig2E_old_child_fitted_microscopy/Fig2E_PCR_prevalence_old) - (Fig2E_adult_fitted_microscopy/Fig2E_PCR_prevalence_adults)

plot(Fig2E_PCR_prevalence_young, Fig2E_adult_young_diff, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "#17BEBB", lwd = 2)
lines(Fig2E_PCR_prevalence_young, Fig2E_adult_old_diff, type = "l", xlim = c(0, 1), ylim = c(0, 1), col = "#76B041", lwd = 2)

# Figure 2H Plotting - Lollipop Plot of Data From Studies With All 3 Age Categories Only
data_frame_complete_age_studies_only <- read.csv("Q:/Sub-Patent Malarial Infections/Data/Age Data Studies With All Three Categories.csv")

ggplot(data_frame_complete_age_studies_only) +
  scale_colour_manual(values = c("#EAE600", "orange", "red"), labels = c("0-5 Years", "5-15 Years", "15+ Years")) +
  geom_line(aes(x = Sensitivity, y = Author, group = Author)) +
  geom_point(aes(x = Sensitivity, y = Author, color = Age.Group, group = Author), size=3 ) + 
  labs(x = "Sensitivity", y = "Study", color = "Age Group")  

# Figure 2I Plotting - Lollipop Plot of Data From Studies With All 3 Age Categories Only, But Only 2 Categories (Children and Adults) Plotted
data_frame_two_cats <- read.csv("Q:/Sub-Patent Malarial Infections/Data/2 CATS ONLY Age Data Studies With All Three Categories.csv")
data_frame_two_cats <- data_frame_two_cats[1:54, ]

ggplot(data_frame_two_cats) +
  scale_colour_manual(values = c("#EAE600", "red"), labels = c("Children", "Adults")) +
  geom_line(aes(x = Sensitivity, y = Author, group = Author)) +
  geom_point(aes(x = Sensitivity, y = Author, color = Age.Group, group = Author), size=3 ) + 
  labs(x = "Sensitivity", y = "Study", color = "Age Group")  
