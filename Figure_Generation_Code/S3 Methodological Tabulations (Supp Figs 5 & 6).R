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
##        Supp Figure 5: Tabulation of Different Methodologies by Global Region                  ##                    
##        Supp Figure 6: Tabulation of Different Methodologies by Transmission Archetype         ## 
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
data_frame <- read.csv("Data/Submicroscopic_Review_Data_R_Import.csv")
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2 & !is.na(data_frame$Sampling_Season), ]
full_data$Sensitivity <- full_data$Micro_Prev/full_data$PCR_Prev

# Regional Stratification and Tabulation of Diagnostic Methodologies
data_frame$Microscopy_Fields <- as.factor(data_frame$Microscopy_Fields)
data_frame$Microscopy_Leucocytes <- as.factor(data_frame$Microscopy_Leucocytes)
Asia_Oceania <- data_frame[data_frame$Global_Region == "Asia&Oceania" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
East_Africa <- data_frame[data_frame$Global_Region == "East Africa" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
South_America <- data_frame[data_frame$Global_Region == "South America" & data_frame$Full_Or_Age_Disagg_Data == 2, ]
West_Africa <- data_frame[data_frame$Global_Region == "West Africa" & data_frame$Full_Or_Age_Disagg_Data == 2, ]

Asia_Num_Studies <- sum(table(Asia_Oceania$PCR_Method))
East_Africa_Num_Studies <- sum(table(East_Africa$PCR_Method))
South_America_Num_Studies <- sum(table(South_America$PCR_Method))
West_Africa_Num_Studies <- sum(table(West_Africa$PCR_Method))

par(mfrow = c(2, 2))
PCR_Methodology_By_Continent <- rbind("Asia & Oceania" = table(Asia_Oceania$PCR_Method), "East Africa" = table(East_Africa$PCR_Method), 
                                      "South America" = table(South_America$PCR_Method), "West Africa" = table(West_Africa$PCR_Method))
Normalised_PCR_Methodology_By_Continent <- rbind("Asia & Oceania" = table(Asia_Oceania$PCR_Method)/Asia_Num_Studies, "East Africa" = table(East_Africa$PCR_Method)/East_Africa_Num_Studies,
                                                 "South America" = table(South_America$PCR_Method)/South_America_Num_Studies, "West Africa" = table(West_Africa$PCR_Method)/West_Africa_Num_Studies)
palette(c("#F15025", "#F2328C", "#7ACC70", "#00A7E1", "#8C00FF"))
barplot(t(Normalised_PCR_Methodology_By_Continent), beside = TRUE, ylim = c(0, 1.1), las = 1, ylab = "Proportion of Studies", col = palette())
legend(1, 1.1, c("Nested", "PCR-LDR", "qPCR", "RT-PCR", "Semi-Nested"), fill = palette(), ncol = 5, title = "PCR Methodology")

LM_Fields_By_Continent <- rbind("Asia & Oceania" = table(Asia_Oceania$Microscopy_Fields), "East Africa" = table(East_Africa$Microscopy_Fields), 
                                "South America" = table(South_America$Microscopy_Fields), "West Africa" = table(West_Africa$Microscopy_Fields))
Normalised_LM_Fields_By_Continent <- rbind("Asia & Oceania" = table(Asia_Oceania$Microscopy_Fields)/Asia_Num_Studies, "East Africa" = table(East_Africa$Microscopy_Fields)/East_Africa_Num_Studies,
                                                 "South America" = table(South_America$Microscopy_Fields)/South_America_Num_Studies, "West Africa" = table(West_Africa$Microscopy_Fields)/West_Africa_Num_Studies)
palette(c("#F15025", "#F2328C", "#7ACC70", "#00A7E1", "#8C00FF"))
barplot(t(Normalised_LM_Fields_By_Continent), beside = TRUE, ylim = c(0, 1), las = 1, ylab = "Proportion of Studies", col = palette())
legend(1, 1, c("100", "200", "300", "400", "1000"), fill = palette(), ncol = 5, title = "Number of Fields Checked")

LM_Leucs_By_Continent <- rbind("Asia & Oceania" = table(Asia_Oceania$Microscopy_Leucocytes), "East Africa" = table(East_Africa$Microscopy_Leucocytes), 
                                "South America" = table(South_America$Microscopy_Leucocytes), "West Africa" = table(West_Africa$Microscopy_Leucocytes))
Normalised_LM_Leucs_By_Continent <- rbind("Asia & Oceania" = table(Asia_Oceania$Microscopy_Leucocytes)/Asia_Num_Studies, "East Africa" = table(East_Africa$Microscopy_Leucocytes)/East_Africa_Num_Studies,
                                           "South America" = table(South_America$Microscopy_Leucocytes)/South_America_Num_Studies, "West Africa" = table(West_Africa$Microscopy_Leucocytes)/West_Africa_Num_Studies)
palette(c("#F15025", "#F2328C", "#F9E931", "#7ACC70", "#00A7E1", "#3160F9", "#8C00FF"))
barplot(t(Normalised_LM_Leucs_By_Continent), beside = TRUE, ylim = c(0, 1), las = 1, ylab = "Proportion of Studies", col = palette())
legend(1, 1, c("50", "200", "300", "350", "400", "500", "1000"), fill = palette(), ncol = 4, title = "Number of Leucocytes Counted")

# Archetype Stratification and Tabulation of Diagnostic Methodologies
data_frame$Microscopy_Fields <- as.factor(data_frame$Microscopy_Fields)
data_frame$Microscopy_Leucocytes <- as.factor(data_frame$Microscopy_Leucocytes)
sub_data_frame <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ] #Remove age disaggregated data
sub_data_frame <- sub_data_frame[order(data_frame$Transmission_Setting_10), ] # orders by transmission history status (puts NAs at end)
sub_data_frame <- sub_data_frame[1:165, ] # removes NAs which are surveys not conducted in Africa for which Trans_Hist data was not available 
high_high_subset <- sub_data_frame[sub_data_frame$Transmission_Setting_15 == "High_High", ]
high_low_subset <- sub_data_frame[sub_data_frame$Transmission_Setting_15 == "High_Low", ]
low_low_subset <- sub_data_frame[sub_data_frame$Transmission_Setting_15 == "Low_Low", ]

High_High_Num_Studies <- sum(table(high_high_subset$PCR_Method))
High_Low_Num_Studies <- sum(table(high_low_subset$PCR_Method))
Low_Low_Num_Studies <- sum(table(low_low_subset$PCR_Method))

par(mfrow = c(2, 2))
PCR_Methodology_By_Archetype <- rbind("High_High" = table(high_high_subset$PCR_Method), "High Low" = table(high_low_subset$PCR_Method), "Low Low" = table(low_low_subset$PCR_Method))
Normalised_PCR_Methodology_By_Archetype <- rbind("High_High" = table(high_high_subset$PCR_Method)/High_High_Num_Studies, "High Low" = table(high_low_subset$PCR_Method)/High_Low_Num_Studies,
                                                 "Low Low" = table(low_low_subset$PCR_Method)/Low_Low_Num_Studies)
palette(c("#F15025", "#F2328C", "#7ACC70", "#00A7E1", "#8C00FF"))
barplot(t(Normalised_PCR_Methodology_By_Archetype), beside = TRUE, ylim = c(0, 1.1), las = 1, ylab = "Proportion of Studies", col = palette())
legend(1, 1.1, c("Nested", "PCR-LDR", "qPCR", "RT-PCR", "Semi-Nested"), fill = palette(), ncol = 5, title = "PCR Methodology")

LM_Fields_By_Archetype <- rbind("High_High" = table(high_high_subset$Microscopy_Fields), "High Low" = table(high_low_subset$Microscopy_Fields), "Low Low" = table(low_low_subset$Microscopy_Fields))
Normalised_LM_Fields_By_Archetype <- rbind("High_High" = table(high_high_subset$Microscopy_Fields)/High_High_Num_Studies, "High Low" = table(high_low_subset$Microscopy_Fields)/High_Low_Num_Studies, "Low Low" = table(low_low_subset$Microscopy_Fields)/Low_Low_Num_Studies)
palette(c("#F15025", "#F2328C", "#7ACC70", "#00A7E1"))
barplot(t(Normalised_LM_Fields_By_Archetype), beside = TRUE, ylim = c(0, 1.1), las = 1, ylab = "Proportion of Studies", col = palette())
legend(1, 1, c("100", "200", "300", "400"), fill = palette(), ncol = 4, title = "Number of Fields Checked")

LM_Leucs_By_Archetype <- rbind("High_High" = table(high_high_subset$Microscopy_Leucocytes), "High Low" = table(high_low_subset$Microscopy_Leucocytes), "Low Low" = table(low_low_subset$Microscopy_Leucocytes))
Normalised_LM_Leucs_By_Archetype <- rbind("High_High" = table(high_high_subset$Microscopy_Leucocytes)/High_High_Num_Studies, "High Low" = table(high_low_subset$Microscopy_Leucocytes)/High_Low_Num_Studies, "Low Low" = table(low_low_subset$Microscopy_Leucocytes)/Low_Low_Num_Studies)
palette(c("#F15025", "#F2328C", "#F9E931", "#7ACC70", "#00A7E1", "#8C00FF"))
barplot(t(Normalised_LM_Leucs_By_Archetype), beside = TRUE, ylim = c(0, 1.1), las = 1, ylab = "Proportion of Studies", col = palette())
legend(1, 1, c("50", "200", "300", "350", "500", "1000"), fill = palette(), ncol = 4, title = "Number of Leucocytes Counted")

