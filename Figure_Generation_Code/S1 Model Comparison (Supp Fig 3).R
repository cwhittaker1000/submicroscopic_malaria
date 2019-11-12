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
##        Supp Figure 3: Comparison of Different Logistic Models                                 ##                    
##                                                                                               ##
##    Any questions, queries, comments, or mistakes, please feel free to get in touch at:        ##
##        charles.whittaker16@imperial.ac.uk :)                                                  ##
##                                                                                               ##
###################################################################################################
library(rjags); library(ssa); library(binom); library(MALDIquant); library(formattable); 
library(tictoc); library(BayesianTools); library(R2jags); library(bayesmix)
setwd("C:/Users/cw1716/Documents/Q_Drive_Copy/Sub-Patent Malarial Infections/Sub_Patent_Malaria_Analysis/")
source("Functions/Submicroscopic_Analysis_Functions.R")

# Load in the dataset and subset the data into:
data_frame <- read.csv("Data/Submicroscopic_Review_Data_R_Import.csv")
full_data <- data_frame[data_frame$Full_Or_Age_Disagg_Data == 2, ]

# Specifying the Data in the Format Required for Input into RJAGS
full_microscopy_PCR_comparison <- list(prev_pcr = full_data$PCR_N_Positive, ## number positive by PCR,
                                       prev_microscopy = full_data$Microscopy_N_Positive, ## number positive by microscopy,
                                       total_pcr = full_data$PCR_N_Tested, ## number tested by PCR,
                                       total_microscopy = full_data$Microscopy_N_Tested, ## number tested by microscopy
                                       N = length(full_data$Name)) # Total sample data overall

PCR_values <- seq(0.001, 0.99, 0.001)
logit_PCR <- logit(PCR_values)
seed <- seq(143, 193, 1) 

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
basic_params <- c("delt", "taud")
basic_jags_inits <- function(){
  list("delt" = 0.001, "taud" = 0.001)
}

# Specifying and initialising the RJAGS model- creates a JAGS model object
Basic_DIC_vector <- c()
for (i in 1:50) {
  LM_basic_model <- jags.parallel(data = full_microscopy_PCR_comparison, inits = basic_jags_inits, parameters.to.save = basic_params, jags.seed = seed[i], model.file = "JAGS_Model/non_centred_LM_Basic_Bayesian_Logit_Linear_Model.txt", n.chains = 4, n.iter = 100000, n.burnin = 20000, n.thin = 20, DIC = TRUE)
  Basic_DIC_vector[i] <- LM_basic_model$BUGSoutput$DIC
  print(i)
}

LM_basic_model_DIC <- LM_basic_model$BUGSoutput$DIC
LM_basic_model <- as.mcmc(LM_basic_model)
autocorr.plot(LM_basic_model)
plot(LM_basic_model)
summary(LM_basic_model)

basic_LM_overall_chain <- rbind(LM_basic_model[[1]], LM_basic_model[[2]], LM_basic_model[[3]], LM_basic_model[[4]])
basic_delt <- mean(as.array(basic_LM_overall_chain[, 1])) # note if using multiple chains, have to calculate this including results from all chains (change subsetting etc)
logit_LM_basic <- logit_PCR + basic_delt 
LM_basic <- expit(logit_LM_basic)

plot(logit_PCR, logit_LM_basic, pch = 20, col = "red", lwd = 3, type = "l", ylim = c(-10, 4))
points(logit(full_data$PCR_Prev), logit(full_data$Micro_Prev + 0.001), pch = 20)
plot(PCR_values, LM_basic, type = "l", col = "red", lwd = 3)
points(full_data$PCR_Prev, full_data$Micro_Prev, pch = 20)

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
standard_params <- c("beta" ,"delt", "taud")
standard_jags_inits <- function(){
  list("beta" = 0.0001, "delt" = 0.001, "taud" = 0.001)
}

Standard_DIC_vector <- c()
for (i in 1:50) {
  LM_standard_model <- jags.parallel(data = full_microscopy_PCR_comparison, inits = standard_jags_inits, parameters.to.save = standard_params, jags.seed = seed[i], model.file = 'JAGS_Model/non_centred_LM_Standard_Bayesian_Logit_Linear_Model.txt', n.chains = 4, n.iter = 100000, n.burnin = 20000, n.thin = 20, DIC = TRUE)
  LM_standard_model_DIC <- LM_standard_model$BUGSoutput$DIC
  Standard_DIC_vector[i] <- LM_standard_model_DIC
  print(i)
}
LM_standard_model <- as.mcmc(LM_standard_model)
autocorr.plot(LM_standard_model)
plot(LM_standard_model)
summary(LM_standard_model)

standard_LM_overall_chain <- rbind(LM_standard_model[[1]], LM_standard_model[[2]], LM_standard_model[[3]], LM_standard_model[[4]])
standard_beta <- mean(as.array(standard_LM_overall_chain[, 1])) 
standard_delt <- mean(as.array(standard_LM_overall_chain[, 2])) 

logit_LM_standard <- logit_PCR + standard_delt + standard_beta * (logit_PCR - mean(logit_PCR))
LM_standard <- expit(logit_LM_standard)

plot(logit_PCR, logit_LM_standard, pch = 20, col = "red", lwd = 3, type = "l", ylim = c(-10, 4))
points(logit(full_data$PCR_Prev), logit(full_data$Micro_Prev + 0.001), pch = 20)
plot(PCR_values, LM_standard, type = "l", col = "red", lwd = 3)
points(full_data$PCR_Prev, full_data$Micro_Prev, pch = 20)

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
quadratic_params <- c("sigma","beta" ,"delt", "taud")
quadratic_jags_inits <- function(){
  list("sigma" = 0.0001, "beta" = 0.0001, "delt" = 0.001, "taud" = 0.001)
}

Quadratic_DIC_vector <- c()
for (i in 1:50) {
  LM_quadratic_model <- jags.parallel(data = full_microscopy_PCR_comparison, inits = quadratic_jags_inits, parameters.to.save = quadratic_params, jags.seed = seed[i], model.file = 'JAGS_Model/non_centred_LM_Quadratic_Bayesian_Logit_Linear_Model.txt', n.chains = 4, n.iter = 100000, n.burnin = 20000, n.thin = 20, DIC = TRUE)
  LM_quadratic_model_DIC <- LM_quadratic_model$BUGSoutput$DIC
  Quadratic_DIC_vector[i] <- LM_quadratic_model_DIC
  print(i)
}
LM_quadratic_model <- as.mcmc(LM_quadratic_model)
autocorr.plot(LM_quadratic_model)
plot(LM_quadratic_model)
summary(LM_quadratic_model)

quadratic_LM_overall_chain <- rbind(LM_quadratic_model[[1]], LM_quadratic_model[[2]], LM_quadratic_model[[3]], LM_quadratic_model[[4]])
quadratic_sigma <- mean(as.array(quadratic_LM_overall_chain[, 4])) 
quadratic_beta <- mean(as.array(quadratic_LM_overall_chain[, 1])) 
quadratic_delt <- mean(as.array(quadratic_LM_overall_chain[, 2])) 

logit_LM_quadratic <- logit_PCR + quadratic_delt + quadratic_beta * (logit_PCR - mean(logit_PCR)) + quadratic_sigma * (logit_PCR - mean(logit_PCR))^2
LM_quadratic <- expit(logit_LM_quadratic)

plot(logit_PCR, logit_LM_quadratic, pch = 20, col = "red", lwd = 3, type = "l", ylim = c(-10, 4))
points(logit(full_data$PCR_Prev), logit(full_data$Micro_Prev + 0.001), pch = 20)
plot(PCR_values, LM_quadratic, type = "l", col = "red", lwd = 3)
points(full_data$PCR_Prev, full_data$Micro_Prev, pch = 20)

# Specifying the parameters of interest that RJAGS will Track and Output, and the initial values to start the chain with
cubic_params <- c("sigma","beta", "gamma", "delt", "taud")
cubic_jags_inits <- function(){
  list("sigma" = 0.0001, "beta" = 0.0001, "gamma" = 0.001, "delt" = 0.001, "taud" = 0.001)
}

Cubic_DIC_vector <- c()
for (i in 1:50) {
  LM_cubic_model <- jags.parallel(data = full_microscopy_PCR_comparison, inits = cubic_jags_inits, parameters.to.save = cubic_params, model.file = 'JAGS_Model/non_centred_LM_Cubic_Bayesian_Logit_Linear_Model.txt', n.chains = 4, jags.seed = seed[i], n.iter = 100000, n.burnin = 20000, n.thin = 20, DIC = TRUE)
  LM_cubic_model_DIC <- LM_cubic_model$BUGSoutput$DIC
  Cubic_DIC_vector[i] <- LM_cubic_model_DIC
  print(i)
}
LM_cubic_model <- as.mcmc(LM_cubic_model)
autocorr.plot(LM_cubic_model)
plot(LM_cubic_model)
summary(LM_cubic_model)

cubic_LM_overall_chain <- rbind(LM_cubic_model[[1]], LM_cubic_model[[2]], LM_cubic_model[[3]], LM_cubic_model[[4]])
cubic_sigma <- mean(as.array(cubic_LM_overall_chain[, 5])) 
cubic_beta <- mean(as.array(cubic_LM_overall_chain[, 1])) 
cubic_delt <- mean(as.array(cubic_LM_overall_chain[, 2])) 
cubic_gamma <- mean(as.array(cubic_LM_overall_chain[, 4])) 

logit_LM_cubic <- logit_PCR + cubic_delt + cubic_beta * (logit_PCR - mean(logit_PCR)) + cubic_sigma * (logit_PCR - mean(logit_PCR))^2 + cubic_gamma * (logit_PCR - mean(logit_PCR))^3
LM_cubic <- expit(logit_LM_cubic)

plot(logit_PCR, logit_LM_cubic, pch = 20, col = "red", lwd = 3, type = "l", ylim = c(-10, 4))
points(logit(full_data$PCR_Prev), logit(full_data$Micro_Prev + 0.0001), pch = 20)
plot(PCR_values, LM_cubic, type = "l", col = "red", lwd = 3)
points(full_data$PCR_Prev, full_data$Micro_Prev, pch = 20)

# Comparing DIC Vectors
Basic_DIC_vector
Standard_DIC_vector
Quadratic_DIC_vector
Cubic_DIC_vector

basic_mean <- mean(Basic_DIC_vector)
standard_mean <- mean(Standard_DIC_vector)
quadratic_mean <- mean(Quadratic_DIC_vector)
cubic_mean <- mean(Cubic_DIC_vector[-23])

basic_sd <- sd(Basic_DIC_vector)
standard_sd <- sd(Standard_DIC_vector)
quadratic_sd <- sd(Quadratic_DIC_vector)
cubic_sd <- sd(Cubic_DIC_vector)

basic_se <- sd(Basic_DIC_vector)/sqrt(length(Basic_DIC_vector))
standard_se <- sd(Standard_DIC_vector)/sqrt(length(Standard_DIC_vector))
quadratic_se <- sd(Quadratic_DIC_vector)/sqrt(length(Quadratic_DIC_vector))
cubic_se <- sd(Cubic_DIC_vector[-23])/sqrt(length(Cubic_DIC_vector))


# Plotting Altogether
plot(logit_PCR, logit_LM_basic, pch = 20, col = "red", lwd = 3, type = "l", ylim = c(-8, 4), xlab = "PCR Prevalence (Logit Scale)", ylab = "LM Prevalence (Logit Scale)", las = 1)
lines(logit_PCR, logit_LM_standard, pch = 20, col = "green", lwd = 3, type = "l", ylim = c(-10, 4))
lines(logit_PCR, logit_LM_quadratic, pch = 20, col = "blue", lwd = 3, type = "l", ylim = c(-10, 4))
lines(logit_PCR, logit_LM_cubic, pch = 20, col = "purple", lwd = 3, type = "l", ylim = c(-10, 4))
points(logit(full_data$PCR_Prev), logit(full_data$Micro_Prev + 0.001), pch = 20)
legend(x = -7, y = 4, legend = c("Basic", "Standard", "Quadratic", "Cubic"), col = c("red", "green", "blue", "purple"), lty = rep(1, 4))

plot(PCR_values, LM_basic, type = "l", col = "red", lwd = 3, xlab = "PCR Prevalence", ylab = "LM Prevalence", las = 1)
lines(PCR_values, LM_standard, type = "l", col = "green", lwd = 3)
lines(PCR_values, LM_quadratic, type = "l", col = "blue", lwd = 3)
lines(PCR_values, LM_cubic, type = "l", col = "purple", lwd = 3)
points(full_data$PCR_Prev, full_data$Micro_Prev, pch = 20)
legend(x = 0, y = 1, legend = c("Basic", "Standard", "Quadratic", "Cubic"), col = c("red", "green", "blue", "purple"), lty = rep(1, 4))
lines(seq(0,1,0.01), seq(0,1,0.01), lwd = 2, lty = 2)

LM_basic_model_DIC
LM_standard_model_DIC
LM_quadratic_model_DIC
LM_cubic_model_DIC

p_D_basic <- var(as.array(LM_basic_model[, 2]))/2 # effective number of parameters
p_D_standard <- var(as.array(LM_standard_model[, 3]))/2 # effective number of parameters
p_D_quadratic <- var(as.array(LM_quadratic_model[, 3]))/2 # effective number of parameters
p_D_cubic <- var(as.array(LM_cubic_model[, 3]))/2 # effective number of parameters

posterior_mean_deviance_basic <- mean(as.array(LM_basic_model[, 2])) 
posterior_mean_deviance_standard <- mean(as.array(LM_standard_model[, 3])) 
posterior_mean_deviance_quadratic <- mean(as.array(LM_quadratic_model[, 3])) 
posterior_mean_deviance_cubic <- mean(as.array(LM_cubic_model[, 3])) 

deviance_basic <- p_D_basic + posterior_mean_deviance_basic
deviance_standard <- p_D_standard + posterior_mean_deviance_standard
deviance_quadratic <- p_D_quadratic + posterior_mean_deviance_quadratic
deviance_cubic <- p_D_cubic + posterior_mean_deviance_cubic

