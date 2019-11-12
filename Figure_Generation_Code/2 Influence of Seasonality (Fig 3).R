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
##        Figure 3: Modelling the Influence of Seasonality on the Prevalence Ratio               ##
##                                                                                               ##
##    Any questions, queries, comments, or mistakes, please feel free to get in touch at:        ##
##        charles.whittaker16@imperial.ac.uk :)                                                  ##
##                                                                                               ##
###################################################################################################
library(rjags); library(ssa); library(binom); library(MALDIquant); library(formattable); 
library(tictoc); library(BayesianTools); library(R2jags); library(bayesmix); library(car)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/Sub-Patent Malarial Infections/Sub_Patent_Malaria_Analysis/")
seed <- 193

# Loading In and Processing the Dataset
data_frame <- read.csv("Data/Submicroscopic_Review_Data_R_Import.csv")
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2 & !is.na(data_frame$Sampling_Season), ]
full_data$prev_ratio <- full_data$Micro_Prev/full_data$PCR_Prev
wet_season <- full_data[!is.na(full_data$Sampling_Season) & as.character(full_data$Sampling_Season) == "Wet", ]
dry_season <- full_data[!is.na(full_data$Sampling_Season) & as.character(full_data$Sampling_Season) == "Dry", ]

###################################################################################################
##                                                                                               ##
##                 Plotting the Prevalence Ratio Stratified By Sampling Season                   ##
##                                                                                               ##
###################################################################################################
colours <- c("#F28D13", "#0FA7C1")
bp <- boxplot(prev_ratio ~ Sampling_Season, data = full_data, las = 1, col = adjustcolor(colours, alpha.f = 0.4), 
              border = colours, boxlwd = 3, whisklty = 1, staplelwd = 3, whisklwd = 3, ylab = "Prevalence Ratio", xlab = "Sampling Season")
mylevels <- levels(full_data$Sampling_Season)
levelProportions <- summary(full_data$Sampling_Season[!is.na(full_data$Sampling_Season)]) / length(full_data$Sampling_Season[!is.na(full_data$Sampling_Season)])
for(i in 1:length(mylevels)){
  thislevel <- mylevels[i]
  thisvalues <- full_data$prev_ratio[full_data$Sampling_Season == thislevel]
  myjitter <- jitter(rep(i, length(thisvalues)), amount = 0.3)
  if (i == 1) {
    points(myjitter, thisvalues, pch = 20, col =  "#F28D13", cex = 1.5)
  } else {
    points(myjitter, thisvalues, pch = 20, col =  "#0FA7C1", cex = 1.5)
  }
}

# Statistical Model Fitting and Results 
variance <- full_data$PCR_N_Tested * full_data$PCR_Prev * (1 - full_data$PCR_Prev)
weighted_seasonal_model <- lm(prev_ratio ~ Sampling_Season, data = full_data, na.action = na.omit, weights = 1/variance)
summary(weighted_seasonal_model)
Anova(weighted_seasonal_model)

weighted_seasonal_model <- lm(prev_ratio ~ PCR_Prev + Sampling_Season, data = full_data, na.action = na.omit, weights = 1/variance)
summary(weighted_seasonal_model)
Anova(weighted_seasonal_model)
