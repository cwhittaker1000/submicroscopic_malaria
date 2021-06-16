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
##        Supp Figure 4: Comparison of Different Methodology (PCR) Prevalence Ratios             ##                    
##                                                                                               ##
##    Any questions, queries, comments, or mistakes, please feel free to get in touch at:        ##
##        charles.whittaker16@imperial.ac.uk :)                                                  ##
##                                                                                               ##
###################################################################################################
library(rjags); library(ssa); library(binom); library(MALDIquant); library(formattable); 
library(tictoc); library(BayesianTools); library(R2jags); library(bayesmix)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/Sub-Patent Malarial Infections/Sub_Patent_Malaria_Analysis/")
source("Functions/Submicroscopic_Analysis_Functions.R")

# Loading In and Processing the Dataset
data_frame <- read.csv("Data/SI_Systematic_Review_Results_R_Import.csv")
data_frame$PCR_Method <- factor(data_frame$PCR_Method)
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]
full_data$Prev_Ratio <- full_data$Micro_Prev/full_data$PCR_Prev
table(full_data$PCR_Method)

# Figure 6A - Relating Methodology to Prevalence Ratio
variance <- full_data$PCR_N_Tested * full_data$PCR_Prev * (1 - full_data$PCR_Prev)
stdev <- sqrt(variance)
weights <- mean(stdev)/stdev

weighted_LM_fields_model <- lm(Prev_Ratio ~ PCR_Prev + Microscopy_Fields, data = full_data, na.action = na.omit, weights = 1/variance)
summary(weighted_LM_fields_model)
Anova(weighted_LM_fields_model)

weighted_PCR_method_model <- lm(Prev_Ratio ~ PCR_Prev + PCR_Method, data = full_data, na.action = na.omit, weights = 1/variance)
summary(weighted_PCR_method_model)
summary(aov(weighted_PCR_method_model))
TukeyHSD(aov(weighted_PCR_method_model))

unweighted_PCR_method_model <- lm(Prev_Ratio ~ PCR_Method, data = full_data, na.action = na.omit)
summary(unweighted_PCR_method_model)
Anova(unweighted_PCR_method_model)
TukeyHSD(aov(unweighted_PCR_method_model))

# Supplementary Figure 4A - PCR Methodology Prevalence Ratio Boxplots
pdf("Figures/Supplementary/Supp Figure 6 - Diagnostic Quality & Methodology/Supp Figure 6 PCR Diagnostic Quality.pdf", 
    width = 14, height = 7.51, useDingbats = FALSE)
par(mfrow = c(1, 2))
nested <- full_data$Prev_Ratio[full_data$PCR_Method == "Nested"]
LDR <- full_data$Prev_Ratio[full_data$PCR_Method == "PCR-LDR"]
qPCR <- full_data$Prev_Ratio[full_data$PCR_Method == "qPCR"]
semi_nested <- full_data$Prev_Ratio[full_data$PCR_Method == "Semi-Nested"]
RT_PCR <- full_data$Prev_Ratio[full_data$PCR_Method == "RT-PCR"]

colours <- c("#F15025", "#F2328C", "#7ACC70", "#00A7E1", "#EDB21C")
bp <- boxplot(Prev_Ratio ~ PCR_Method, data = full_data, las = 1, col = adjustcolor(colours, alpha.f = 0.4), 
              border = colours, boxlwd = 3, whisklty = 1, staplelwd = 3, whisklwd = 3, ylab = "Prevalence Ratio", xlab = "PCR Method")

nested_weights <- weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "Nested"]
ldr_weights <- weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "PCR-LDR"]
qPCR_weights <- weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "qPCR"]
semi_weights <- weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "Semi-Nested"]

mylevels <- levels(full_data$PCR_Method)
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- full_data$Prev_Ratio[full_data$PCR_Method == thislevel]
  myjitter <- jitter(rep(i, length(thisvalues)), amount = 0.3)
  if (i == 1) {
    points(myjitter, thisvalues, pch = 20, col =  "#F15025", cex = nested_weights)
  } else if (i == 2) {
    points(myjitter, thisvalues, pch = 20, col =  "#F2328C", cex = ldr_weights)
  } else if (i == 3) {
    points(myjitter, thisvalues, pch = 20, col =  "#7ACC70", cex = qPCR_weights)
  } else if (i == 4) {
    points(myjitter, thisvalues, pch = 20, col =  "#00A7E1", cex = semi_weights)
  } else if (i == 5) {
    points(myjitter, thisvalues, pch = 20, col =  "#EDB21C", cex = semi_weights)
  }
}

weighted_mean_nested <- weighted.mean(nested, w = weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "Nested"])
weighted_mean_ldr <- weighted.mean(LDR, w = weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "PCR-LDR"])
weighted_mean_qPCR <- weighted.mean(qPCR, w = weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "qPCR"])
weighted_mean_RT_PCR <- weighted.mean(RT_PCR, w = weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "RT-PCR"])
weighted_mean_semi_nested <- weighted.mean(semi_nested, w = weights[!is.na(full_data$PCR_Method) & as.character(full_data$PCR_Method) == "Semi-Nested"])

lines(seq(0.61, 1.39, 0.001), rep(weighted_mean_nested, length(seq(0.61, 1.39, 0.001))), pch = 20, cex = 5, col = "#F15025", lwd = 5)
lines(seq(1.61, 2.39, 0.001), rep(weighted_mean_ldr, length(seq(1.61, 2.39, 0.001))), pch = 20, cex = 5, col = "#F2328C", lwd = 5)
lines(seq(2.61, 3.39, 0.001), rep(weighted_mean_qPCR, length(seq(1.61, 2.39, 0.001))), pch = 20, cex = 5, col = "#7ACC70", lwd = 5)
lines(seq(3.61, 4.39, 0.001), rep(weighted_mean_RT_PCR, length(seq(1.61, 2.39, 0.001))), pch = 20, cex = 5, col = "#00A7E1", lwd = 5)
lines(seq(4.61, 5.39, 0.001), rep(weighted_mean_semi_nested, length(seq(1.61, 2.39, 0.001))), pch = 20, cex = 5, col = "#EDB21C", lwd = 5)

# Supplementary Figure 4B - Barplot of PCR Methodology Prevalence Ratio Means
nested_mean <- mean(nested)
nested_std_err <- sd(nested)/sqrt(length(nested))
nested_median <- median(nested)
LDR_mean <- mean(LDR)
LDR_std_err <- sd(LDR)/sqrt(length(LDR))
LDR_median <- median(LDR)
qPCR_mean <- mean(qPCR)
qPCR_std_err <- sd(qPCR)/sqrt(length(qPCR))
qPCR_median <- median(qPCR)
semi_nested_mean <- mean(semi_nested)
semi_nested_std_err <- sd(semi_nested)/sqrt(length(semi_nested))
semi_nested_median <- median(semi_nested)
RT_PCR_mean <- mean(RT_PCR)
RT_PCR_std_err <- sd(RT_PCR)/sqrt(length(RT_PCR))
RT_PCR_median <- median(RT_PCR)

medians <- c(nested_median, LDR_median, qPCR_median, RT_PCR_median, semi_nested_median)
means <- c(nested_mean, LDR_mean, qPCR_mean, RT_PCR_mean, semi_nested_mean)
lower_ci <- means - 1.96 * c(nested_std_err, LDR_std_err, qPCR_std_err, RT_PCR_std_err, semi_nested_std_err)
upper_ci <- means + 1.96 * c(nested_std_err, LDR_std_err, qPCR_std_err, RT_PCR_std_err, semi_nested_std_err)

palette(c("#F15025", "#F2328C", "#7ACC70", "#00A7E1", "#EDB21C"))
diagnostic_sensitivity <- barplot(means, names.arg = c("Nested", "LDR", "qPCR", "RT-PCR", "Semi-Nested"), las = 1, ylim = c(0, 1), col = palette(), xlab = "PCR Method", ylab = "Prevalence Ratio")
arrows(x0 = diagnostic_sensitivity, y0 = lower_ci, x1 = diagnostic_sensitivity, y1 = upper_ci, col = "black", angle = 90, code = 3, length = 0.05)
dev.off()
