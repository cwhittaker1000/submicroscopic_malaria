# Access relevant packages
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

# Assign whichever dataset being used to data_frame
data_frame <- as.data.frame(Data.4th.February)

# Manipulate data_frame in prep for analysis
data_frame$Rainy_Or_Dry <- as.factor(data_frame$Rainy_Or_Dry)
  
# Create individual datasets for each global region
rainy_season_data <- data_frame[data_frame$Rainy_Or_Dry == "rainy", ]
dry_season_data <- data_frame[data_frame$Rainy_Or_Dry == "dry", ]


# Running the linear regression, with PCR regressed onto Microscopy. 
rainy_season_lm <- summary(lm(logit(rainy_season_data$Micro_Prev) ~ logit(rainy_season_data$PCR_Prev)))
dry_season_lm <- summary(lm(logit(dry_season_data$Micro_Prev) ~ logit(dry_season_data$PCR_Prev)))


# Defining logitPCR for plotting purposes
logit_PCR <- logit(seq(0,1,0.001))


# Logit scale plotting
plot(logit(rainy_season_data$PCR_Prev), logit(rainy_season_data$Micro_Prev))
lines(logit_PCR, rainy_season_lm$coefficients[1,1] + rainy_season_lm$coefficients[2,1]*logit_PCR, lwd = 2)

plot(logit(dry_season_data$PCR_Prev), logit(dry_season_data$Micro_Prev))
lines(logit_PCR, dry_season_lm$coefficients[1,1] + dry_season_lm$coefficients[2,1]*logit_PCR, lwd = 2)

# Natural scale plotting
plot(rainy_season_data$PCR_Prev, rainy_season_data$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "black")
lines(expit(logit_PCR), expit(rainy_season_lm$coefficients[1,1] + rainy_season_lm$coefficients[2,1]*logit_PCR), col = "black", lwd = 2)

points(dry_season_data$PCR_Prev, dry_season_data$Micro_Prev, xlim = c(0,1), ylim = c(0,1), pch = 20, col = "orange")
lines(expit(logit_PCR), expit(dry_season_lm$coefficients[1,1] + dry_season_lm$coefficients[2,1]*logit_PCR), col = "orange", lwd = 2)

### LOW TRANSMISSION SETTING EXPLORATION WITH LINEAR RELATIONSHIP

# Create individual datasets for each global region
subset <- 0.2
rainy_season_data <- data_frame[data_frame$Rainy_Or_Dry == "rainy" & data_frame$PCR_Prev < subset, ]
dry_season_data <- data_frame[data_frame$Rainy_Or_Dry == "dry" & data_frame$PCR_Prev < subset, ]

# Running the linear regression, with PCR regressed onto Microscopy. 
rainy_season_lm <- summary(lm(rainy_season_data$Micro_Prev ~ 0 + rainy_season_data$PCR_Prev))
dry_season_lm <- summary(lm(dry_season_data$Micro_Prev ~ 0 + dry_season_data$PCR_Prev))

# Defining logitPCR for plotting purposes
PCR <- seq(0,subset,0.001)

# Natural scale plotting
plot(rainy_season_data$PCR_Prev, rainy_season_data$Micro_Prev, xlim = c(0,subset), pch = 20, col = "black")
lines(PCR, rainy_season_lm$coefficients[1,1]*PCR, col = "black", lwd = 2)

points(dry_season_data$PCR_Prev, dry_season_data$Micro_Prev, xlim = c(0,subset), ylim = c(0,1), pch = 20, col = "orange")
lines(PCR, dry_season_lm$coefficients[1,1]*PCR, col = "orange", lwd = 2)

