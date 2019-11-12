# run the rJAGS model 
run_rJAGS_model <- function(number_iterations, number_chains, model_file, params, initial_values_function, data) {
  
  # Specifying the Data in the Format Required for Input into RJAGS
  data_JAGS <- list(prev_pcr = data$PCR_N_Positive, ## number positive by PCR,
                    prev_microscopy = data$Microscopy_N_Positive, ## number positive by microscopy,
                    total_pcr = data$PCR_N_Tested, ## number tested by PCR,
                    total_microscopy = data$Microscopy_N_Tested, ## number tested by microscopy
                    N = length(data$Name)) # Total sample data overall
  
  # Specifying and initialising the RJAGS model- creates a JAGS model object
  JAGS_output <- jags.parallel(data = data_JAGS, inits = initial_values_function, 
                               parameters.to.save = params, jags.seed = seed, 
                               model.file = model_file, n.chains = 4, n.iter = number_iterations, 
                               n.burnin = number_iterations/2, n.thin = 25, DIC = TRUE)
  JAGS_output <- as.mcmc(JAGS_output)
  return(JAGS_output)
  
}

# extract parameter values from JAGS output
param_table <- function(JAGS_output, parameters) {
  param_estimates <- summary(JAGS_output)
  param_estimates <- as.data.frame(cbind(round(param_estimates$statistics, 3)[parameters, 1:2], 
                                         round(param_estimates$quantiles, 3)[parameters, c(1, 3, 5)]))
  colnames(param_estimates) <- c("Mean", "SD", "2.5%", "Median", "97.5%")
  param_table <- formattable(param_estimates, align = c("r", "c", "r", "c", "l"))
  return(param_table)
}

# calculate mean output for loglinear model
mean_output <- function(PCR_range, beta_mean, delt_mean) {
  logit_PCR <- logit(PCR_range)
  logit_microscopy <- (delt_mean - (beta_mean * mean(logit_PCR))) + ((1 + beta_mean)  * logit_PCR)
  fitted_microscopy <- expit(logit_microscopy)
  return(fitted_microscopy)
}

# calculate credible intervals
credible_intervals <- function(PCR_range, JAGS_output) {
  
  logit_PCR <- logit(PCR_range)
  pred_dist <- matrix(NA, nrow = nrow(JAGS_output), ncol = length(PCR_range))
  for (i in 1:nrow(JAGS_output)) {
    delt <- JAGS_output[i, "delt"]
    beta <- JAGS_output[i, "beta"]
    pred_values <- (delt - (beta * mean(logit_PCR))) + ((1 + beta)  * logit_PCR)
    pred_dist[i, ] <- expit(pred_values)
  }
  
  credible_lower <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.025)
  credible_upper <- apply(pred_dist, MARGIN = 2, quantile, prob = 0.975)
  sensitivity_credible_lower <- credible_lower / PCR_range
  sensitivity_credible_upper <- credible_upper / PCR_range
  
  return(list(credible_lower = credible_lower, credible_upper = credible_upper,
              sensitivity_lower = sensitivity_credible_lower, sensitivity_upper = sensitivity_credible_upper))
  
}



