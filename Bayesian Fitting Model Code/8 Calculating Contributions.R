# R Script for Assessing the Contribution of the Submicroscopic Population to the Infectious Reservoir

# Ask Hannah 1) about the papers
#       and  2) about the SHINY app
#       and  3) about the figure suggestion she had in our meeting that I've now totally forgotten

# Various Parameters
low_HH <- 0
high_HH <- 1
low_HL <- 0
high_HL <- 1
low_LL <- 0
high_LL <- 1

#just for testing
SOME_TERM_FOR_SENSITIVITY <- seq(0, 1, 0.001)

# High High Transmission Setting 

# Range of High High Survey PCR Prevalences
HH_PCR_Prevalence <- seq(0.03,0.95,0.001) 
# Proportion of individuals with patent infections (equivalent to the sensitivity)
HH_Patent_Percentage <- (High_high_fitted_microscopy/PCR_prevalence_high_high) * 100 
# Proportion of individuals with subpatent infections
HH_Subpatent_Percentage <- (100 - HH_Patent_Percentage) 
# Calculating subpatent contribution to transmission
HH_Subpatent_Contribution <- (HH_Subpatent_Percentage) / ((2.65 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
# Plotting
HH_plot <- plot(HH_PCR_Prevalence, HH_Subpatent_Contribution, xlim = c(0, 1), ylim = c(0, 1))

mean(HH_Subpatent_Contribution)

# High Low Transmission Setting 

# Range of High High Survey PCR Prevalences
HL_PCR_Prevalence <- seq(0.004,0.9,0.001)
# Proportion of individuals with patent infections (equivalent to the sensitivity)
HL_Patent_Percentage <- (High_low_fitted_microscopy/PCR_prevalence_high_low) * 100 
# Proportion of individuals with subpatent infections
HL_Subpatent_Percentage <- (100 - HL_Patent_Percentage) 
# Calculating subpatent contribution to transmission
HL_Subpatent_Contribution <- (HL_Subpatent_Percentage) / ((2.65 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
# Plotting
HL_plot <- plot(HL_PCR_Prevalence, HL_Subpatent_Contribution, xlim = c(0, 1), ylim = c(0, 1))

mean(HL_Subpatent_Contribution)

# Low Low Transmission Setting 

# Range of High High Survey PCR Prevalences
LL_PCR_Prevalence <- seq(0.01,0.55,0.001)
# Proportion of individuals with patent infections (equivalent to the sensitivity)
LL_Patent_Percentage <- (Low_low_fitted_microscopy/PCR_prevalence_low_low) * 100 
# Proportion of individuals with subpatent infections
LL_Subpatent_Percentage <- (100 - LL_Patent_Percentage) 
# Calculating subpatent contribution to transmission
LL_Subpatent_Contribution <- (LL_Subpatent_Percentage) / ((2.65 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
# Plotting
LL_plot <- plot(LL_PCR_Prevalence, LL_Subpatent_Contribution, xlim = c(0, 1), ylim = c(0, 1))

mean(LL_Subpatent_Contribution)





