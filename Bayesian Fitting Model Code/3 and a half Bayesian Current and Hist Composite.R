# Access relevant packages
library(rjags)
library(moments)
library(gtools)
library(nortest)
library(ggplot2)
library(LMest)
library(ssa)

# Load in the dataset
data_frame <- Data.4th.February
data_frame <- Data.4th.February.Dec.Rounded

# Manipulate data_frame in prep for analysis
data_frame$Old_or_New <- as.factor(data_frame$Old_or_New)
data_frame$Trans_Hist <- data_frame$Transmission_Setting_History
data_frame$Trans_Hist <- as.factor(data_frame$Trans_Hist)
data_frame$Trans_Now <- data_frame$Transmission_Setting_Current
data_frame$Trans_Now <- as.factor(data_frame$Trans_Now)

# OPTIONAL- REMOVE THOSE WHICH WERE INITIALLY ZERO
data_frame <- data_frame[data_frame$Was_Initially_Zero. == "N", ]

# Subsetting data_frame, do if 0s NOT removed
data_frame <- data_frame[order(data_frame$Trans_Hist),]
data_frame <- data_frame[1:234,]

# Subsetting data_frame, do if 0s removed
data_frame <- data_frame[order(data_frame$Trans_Hist),]
data_frame <- data_frame[1:229,]

# Subsettin the data_frame, only new data, checking whether including Lucy's data changes the conclusions
# NEED TO GO BACK AND CHANGE AS THE OLD DATA IS THE WRONG WAY AROUND, SO DEFINITELY PROCEED WITH THIS FOR NOW
data_frame <- data_frame[data_frame$Old_or_New == "New",]

# Add columns for composite components
for (i in 1:length(data_frame[,1])) {
  if (data_frame$Trans_Hist[i] == 1 & data_frame$Trans_Now[i] == 1) {
    data_frame$high_high[i] <- 1 }
  else {
    data_frame$high_high[i] <- 0
  }
  if (data_frame$Trans_Hist[i] == 1 & data_frame$Trans_Now[i] == 0) {
    data_frame$high_low[i] <- 1 }
  else {
    data_frame$high_low[i] <- 0 
  }
  if (data_frame$Trans_Hist[i] == 0 & data_frame$Trans_Now[i] == 1) {
    data_frame$low_high[i] <- 1 }
  else {
    data_frame$low_high[i] <- 0
  }
  if (data_frame$Trans_Hist[i] == 0 & data_frame$Trans_Now[i] == 0) {
    data_frame$low_low[i] <- 1 }
  else {
    data_frame$low_low[i] <- 0 
  }
}

# Create individual datasets for each composite component
high_high_subset <- data_frame[data_frame$high_high == 1,]
high_low_subset <- data_frame[data_frame$high_low == 1,]
low_low_subset <- data_frame[data_frame$low_low == 1,]

## Specifying the data in the format required for input into RJAGS

#### NOTE THAT THERE ARE CURRENTLY DECIMALS IN THE DATA WHERE I USED LUCY'S TRICK- UNSURE HOW TO INCORPORATE
#### THESE CURRENTLY, MUST CHECK WITH HER. 

# High high data
high_high_microscopy_PCR_comparison <- list(prev_pcr = high_high_subset$PCR_N_Positive, ## number positive by PCR,
                                       prev_microscopy = high_high_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                       total_pcr = high_high_subset$PCR_N_Tested, ## number tested by PCR,
                                       total_microscopy = high_high_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                       N = 64) #49) #64) zeroes included  #first is all data, no zeroes, second 2  are new data without and with 0s.

# High low data
high_low_microscopy_PCR_comparison <- list(prev_pcr = high_low_subset$PCR_N_Positive, ## number positive by PCR,
                                              prev_microscopy = high_low_subset$Microscopy_N_Positive,## number positive by microscopy,
                                              total_pcr = high_low_subset$PCR_N_Tested, ## number tested by PCR,
                                              total_microscopy = high_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                              N = 141) #119) #120) zeroes included  #first is all data, no zeroes, second 2  are new data without and with 0s.

# Low low data
low_low_microscopy_PCR_comparison <- list(prev_pcr = low_low_subset$PCR_N_Positive, ## number positive by PCR,
                                                prev_microscopy = low_low_subset$Microscopy_N_Positive, ## number positive by microscopy,
                                                total_pcr = low_low_subset$PCR_N_Tested, ## number tested by PCR,
                                                total_microscopy = low_low_subset$Microscopy_N_Tested, ## number tested by microscopy,
                                                N = 24) #20) #28) zeroes included  #first is all data, no zeroes, second 2  are new data without and with 0s.


# Specifying the parameters of interest that RJAGS will track and output, and the 
# initial values to start the chain with
params <- c("delt", "beta", "taud")
inits <- list(beta=0.0001,delt=0.0001,taud=0.0001)

# Specifying and initialising the RJAGS model- creates a JAGS model object
High_high_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # OLD DATA
                                       data = high_high_microscopy_PCR_comparison, n.chains = 1, inits = inits)

High_low_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # NEW DATA
                                              data = high_low_microscopy_PCR_comparison, n.chains = 1, inits = inits)

Low_low_micr_PCR_comp_model <- jags.model('1_Basic_Model_Old_And_New.txt',   # ALL DATA
                                                data = low_low_microscopy_PCR_comparison, n.chains = 1, inits = inits)

# Running, updating and iterating the RJAGS model

# High high DATA
update(High_high_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
High_high_micr_PCR_comp_model <- coda.samples(High_high_micr_PCR_comp_model, params, 
                                         n.iter = 35000, thin = 70) # Model updating
summary(High_high_micr_PCR_comp_model)
plot(High_high_micr_PCR_comp_model)
autocorr.plot(High_high_micr_PCR_comp_model)

High_high_beta_mean <- mean(as.array(High_high_micr_PCR_comp_model[, 1]))
High_high_delt_mean <- mean(as.array(High_high_micr_PCR_comp_model[, 2]))

# High low DATA
update(High_low_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
High_low_micr_PCR_comp_model <- coda.samples(High_low_micr_PCR_comp_model, params, 
                                                n.iter = 35000, thin = 70) # Model updating
summary(High_low_micr_PCR_comp_model)
plot(High_low_micr_PCR_comp_model)
autocorr.plot(High_low_micr_PCR_comp_model)

High_low_beta_mean <- mean(as.array(High_low_micr_PCR_comp_model[, 1]))
High_low_delt_mean <- mean(as.array(High_low_micr_PCR_comp_model[, 2]))

# Low low DATA
update(Low_low_micr_PCR_comp_model, 5000) # Updates the model 5000 iterations, basically the burn in
Low_low_micr_PCR_comp_model <- coda.samples(Low_low_micr_PCR_comp_model, params, 
                                                  n.iter = 35000, thin = 70) # Model updating
summary(Low_low_micr_PCR_comp_model)
plot(Low_low_micr_PCR_comp_model)
autocorr.plot(Low_low_micr_PCR_comp_model)

Low_low_beta_mean <- mean(as.array(Low_low_micr_PCR_comp_model[, 1]))
Low_low_delt_mean <- mean(as.array(Low_low_micr_PCR_comp_model[, 2]))

# Specifying PCR prevalence to plot against
PCR_prevalence <- seq(0.001,1,0.001)
logit_PCR_prevalence <- logit(PCR_prevalence)

# Calculating the fitted values specifying the best fit line on the logit scale
# DON'T FORGET IT'S delta' + (1 + BETA) * logit_PCR_prevalence
High_high_fitted_logit_microscopy <- High_high_delt_mean + (1 + High_high_beta_mean)  * logit_PCR_prevalence
High_low_fitted_logit_microscopy <- High_low_delt_mean + (1 + High_low_beta_mean) * logit_PCR_prevalence
Low_low_fitted_logit_microscopy <- Low_low_delt_mean + (1 + Low_low_beta_mean)  * logit_PCR_prevalence

# Converting them to the natural scale
High_high_fitted_microscopy <- expit(High_high_fitted_logit_microscopy)
High_low_fitted_microscopy <- expit(High_low_fitted_logit_microscopy)
Low_low_fitted_microscopy <- expit(Low_low_fitted_logit_microscopy)

# Plotting the results
plot(high_high_subset$PCR_Prev, high_high_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "slateblue4",
     xlab = "PCR Prevalence", ylab = "LM Prevalence")
points(high_low_subset$PCR_Prev, high_low_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "tan1")
points(low_low_subset$PCR_Prev, low_low_subset$Micro_Prev, xlim = c(0, 1), ylim = c(0, 1), pch = 20, col = "green3")

lines(seq(0,1,0.01), seq(0,1,0.01))
lines(PCR_prevalence, High_high_fitted_microscopy, col = "slateblue4", lwd = 3)
lines(PCR_prevalence, High_low_fitted_microscopy, col = "tan1", lwd = 3)
lines(PCR_prevalence, Low_low_fitted_microscopy, col = "green3", lwd = 3)

legend("topleft", 
       legend = c("High High", "High Low", "Low low"), 
       col = c("slateblue4", "tan1", "green3"),
       lty = 1,
       lwd = 2,
       pt.cex = 2, 
       cex = 0.8)
