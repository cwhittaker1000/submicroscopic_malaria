# R Script for Assessing the Contribution of the Submicroscopic Population to the Infectious Reservoir

# Ask Hannah 1) about the papers
#       and  2) about the SHINY app
#       and  3) about the figure suggestion she had in our meeting that I've now totally forgotten

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
plot(0,0, xlim = c(0,1), ylim = c(0,1))
HH_plot <- lines(HH_PCR_Prevalence, HH_Subpatent_Contribution, xlim = c(0, 1), ylim = c(0, 1))

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
plot(0,0, xlim = c(0,1), ylim = c(0,1))
HL_plot <- lines(HL_PCR_Prevalence, HL_Subpatent_Contribution, xlim = c(0, 1), ylim = c(0, 1))

mean(HL_Subpatent_Contribution)

# Low Low Transmission Setting 

# Range of High High Survey PCR Prevalences
LL_PCR_Prevalence <- seq(0.01,0.55,0.001)
# Proportion of individuals with patent infections (equivalent to the sensitivity)
LL_Patent_Percentage <- (Low_low_fitted_microscopy/PCR_prevalence_low_low) * 100 
# Proportion of individuals with subpatent infections
LL_Subpatent_Percentage <- (100 - LL_Patent_Percentage) 
# Calculating subpatent contribution to transmission
LL_Subpatent_Contribution5 <- (LL_Subpatent_Percentage) / ((20 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL_Subpatent_Contribution10 <- (LL_Subpatent_Percentage) / ((5 * LL_Patent_Percentage) + LL_Subpatent_Percentage)
LL_Subpatent_Contribution30 <- (LL_Subpatent_Percentage) / ((2 * LL_Patent_Percentage) + LL_Subpatent_Percentage)

max(LL_Subpatent_Contribution30)

# Plotting
plot(0,0, xlim = c(0,1), ylim = c(0,1))
LL_plot <- lines(LL_PCR_Prevalence, LL_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(LL_PCR_Prevalence, LL_Subpatent_Contribution10, xlim = c(0, 1), ylim = c(0, 1))
lines(LL_PCR_Prevalence, LL_Subpatent_Contribution30, xlim = c(0, 1), ylim = c(0, 1))

zero_line <- rep(0, length(seq(0.01,0.55,0.001)))
plot(0,0, xlim = c(0,1), ylim = c(0,1))
polygon(x = c(seq(0.01,0.55,0.001), rev(seq(0.01,0.55,0.001))), 
        y = c(LL_Subpatent_Contribution5, zero_line))
polygon(x = c(seq(0.01,0.55,0.001), rev(seq(0.01,0.55,0.001))), 
        y = c(LL_Subpatent_Contribution10, zero_line))
polygon(x = c(seq(0.01,0.55,0.001), rev(seq(0.01,0.55,0.001))), 
        y = c(LL_Subpatent_Contribution30, zero_line))


# Low Low Transmission Setting 

# Range of High High Survey PCR Prevalences
HL_PCR_Prevalence <- seq(0.004,0.9,0.001)
# Proportion of individuals with patent infections (equivalent to the sensitivity)
HL_Patent_Percentage <- (High_low_fitted_microscopy/PCR_prevalence_high_low) * 100 
# Proportion of individuals with subpatent infections
HL_Subpatent_Percentage <- (100 - HL_Patent_Percentage) 

# Calculating subpatent contribution to transmission
HL_Subpatent_Contribution5 <- (HL_Subpatent_Percentage) / ((20 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL_Subpatent_Contribution10 <- (HL_Subpatent_Percentage) / ((5 * HL_Patent_Percentage) + HL_Subpatent_Percentage)
HL_Subpatent_Contribution30 <- (HL_Subpatent_Percentage) / ((2 * HL_Patent_Percentage) + HL_Subpatent_Percentage)

# Plotting
plot(0,0, xlim = c(0,1), ylim = c(0,1))
HL_plot <- lines(HL_PCR_Prevalence, HL_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(HL_PCR_Prevalence, HL_Subpatent_Contribution10, xlim = c(0, 1), ylim = c(0, 1))
lines(HL_PCR_Prevalence, HL_Subpatent_Contribution30, xlim = c(0, 1), ylim = c(0, 1))

zero_line <- rep(0, length(seq(0.004,0.9,0.001)))
plot(0,0, xlim = c(0,1), ylim = c(0,1))
HL_plot <- lines(HL_PCR_Prevalence, HL_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
polygon(x = c(seq(0.004,0.9,0.001), rev(seq(0.004,0.9,0.001))), 
        y = c(HL_Subpatent_Contribution5, zero_line))
polygon(x = c(seq(0.004,0.9,0.001), rev(seq(0.004,0.9,0.001))), 
        y = c(HL_Subpatent_Contribution10, zero_line))
polygon(x = c(seq(0.004,0.9,0.001), rev(seq(0.004,0.9,0.001))), 
        y = c(HL_Subpatent_Contribution30, zero_line))


?polygon

polygon(x = c(PCR_prevalence_low_low, rev(PCR_prevalence_low_low)), 
        y = c(low_low_data_sensitivity_upper, rev(low_low_data_sensitivity_lower)), 
        col = adjustcolor("darkgrey", alpha.f = 0.5), border = NA)



# Range of High High Survey PCR Prevalences
HH_PCR_Prevalence <- seq(0.03,0.95,0.001) 
# Proportion of individuals with patent infections (equivalent to the sensitivity)
HH_Patent_Percentage <- (High_high_fitted_microscopy/PCR_prevalence_high_high) * 100 
# Proportion of individuals with subpatent infections
HH_Subpatent_Percentage <- (100 - HH_Patent_Percentage) 

# Calculating subpatent contribution to transmission
HH_Subpatent_Contribution5 <- (HH_Subpatent_Percentage) / ((20 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH_Subpatent_Contribution10 <- (HH_Subpatent_Percentage) / ((5 * HH_Patent_Percentage) + HH_Subpatent_Percentage)
HH_Subpatent_Contribution30 <- (HH_Subpatent_Percentage) / ((2 * HH_Patent_Percentage) + HH_Subpatent_Percentage)

# Plotting
plot(0,0, xlim = c(0,1), ylim = c(0,1))
HL_plot <- lines(HH_PCR_Prevalence, HH_Subpatent_Contribution5, xlim = c(0, 1), ylim = c(0, 1))
lines(HH_PCR_Prevalence, HH_Subpatent_Contribution10, xlim = c(0, 1), ylim = c(0, 1))
lines(HH_PCR_Prevalence, HH_Subpatent_Contribution30, xlim = c(0, 1), ylim = c(0, 1))

zero_line <- rep(0, length(seq(0.03,0.95,0.001)))
plot(0,0, xlim = c(0,1), ylim = c(0,1))
polygon(x = c(seq(0.03,0.95,0.001) , rev(seq(0.03,0.95,0.001))), 
        y = c(HH_Subpatent_Contribution5, zero_line))
polygon(x = c(seq(0.03,0.95,0.001) , rev(seq(0.03,0.95,0.001))), 
        y = c(HH_Subpatent_Contribution10, zero_line))
polygon(x = c(seq(0.03,0.95,0.001) , rev(seq(0.03,0.95,0.001))), 
        y = c(HH_Subpatent_Contribution30, zero_line))


