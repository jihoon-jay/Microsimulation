# Jihoon Lim
# Economic Modelling - Microsimulation
# 06 - Model Calibration

# Note: Une partie de la compilation est effectuée à partir de données provenant 
# du © Gouvernement du Québec (année de la publication du Fichier de recherche: 
# 2009-2019). Le © Gouvernement du Québec n’est pas responsable des compilations
# ni de l’interprétation des résultats produits à l’aide du Fichier de recherche.

library(readxl)
library(data.table)
library(flextable)
library(ggplot2)
#library(rriskDistributions)
library(lhs) #latin hypercube sampling
theme_set(theme_bw())
library(mvtnorm)
library(doParallel)
library(parallel)
library(foreach)
library(compiler)
library(future)
library(future.apply)

########## ---- Part A: Data Generation ##########
# | 1. Model Input Parameters ####
source("~/01 - Setup.R")

# | 2. Mean and SD Inputs for Prior ####
t_calib <- read_excel('~/psa_params.xlsx', sheet = "calibration_params")
#Display it nicely
t_calib |>
  flextable() |> #turn into flextable object
  merge_v(j=1) |> #Merge cells in first column with same value (group probabilities, costs, etc.)
  theme_box() |> #Apply a theme for aesthetics
  autofit() #automatically set column widths to reasonable values

calib_param_names <- t_calib$rname
risk_mean_log <- c()
risk_sd_log <- c()
for (param in calib_param_names){
  risk_mult <- suppressMessages(get.lnorm.par(p=c(0.025, 0.5, 0.975),
                                              q=c(t_calib$value_min[t_calib$rname == param],
                                                  t_calib$value[t_calib$rname == param],
                                                  t_calib$value_max[t_calib$rname == param])))
  risk_mean_log[param] <- risk_mult[1]
  risk_sd_log[param] <- risk_mult[2]
}

########## ---- Part B: Functions ##########
# | 1. Sample Prior ####
sample_prior <- cmpfun(function(n) {
  # n: the number of samples desired
  draws0 <- randomLHS(n=n,k=15)
  draws  <- data.frame(rr.nsp.sharing = qlnorm(draws0[,1], risk_mean_log[1], risk_sd_log[1]),
                       rr.nsp.reusing = qlnorm(draws0[,2], risk_mean_log[2], risk_sd_log[2]),
                       rr.sharing.needles = qlnorm(draws0[,3], risk_mean_log[3], risk_sd_log[3]),
                       rr.reusing.needles = qlnorm(draws0[,4], risk_mean_log[4], risk_sd_log[4]),
                       rr.inj.freq.2x = qlnorm(draws0[,5], risk_mean_log[5], risk_sd_log[5]),
                       rr.inj.freq.7x = qlnorm(draws0[,6], risk_mean_log[6], risk_sd_log[6]),
                       rr.num.attempts2 = qlnorm(draws0[,7], risk_mean_log[7], risk_sd_log[7]),
                       rr.num.attempts3 = qlnorm(draws0[,8], risk_mean_log[8], risk_sd_log[8]),
                       rr.num.attempts4 = qlnorm(draws0[,9], risk_mean_log[9], risk_sd_log[9]),
                       rr.years.idu2to5 = qlnorm(draws0[,10], risk_mean_log[10], risk_sd_log[10]),
                       rr.years.idu5to8 = qlnorm(draws0[,11], risk_mean_log[11], risk_sd_log[11]),
                       rr.years.idu8plus = qlnorm(draws0[,12], risk_mean_log[12], risk_sd_log[12]),
                       rr.age25to44 = qlnorm(draws0[,13], risk_mean_log[13], risk_sd_log[13]),
                       rr.age44plus = qlnorm(draws0[,14], risk_mean_log[14], risk_sd_log[14]),
                       rr.male = qlnorm(draws0[,15], risk_mean_log[15], risk_sd_log[15])
  )
  return(as.matrix(draws))
})

# | 2. Log Prior ####
l_prior <- cmpfun(function(par_vector) {
  # par_vector: a vector (or matrix) of model parameters (omits c)
  if(is.null(dim(par_vector))) par_vector <- t(par_vector)
  lprior <- rep(0,nrow(par_vector))
  lprior <- lprior+dlnorm(par_vector[,1], risk_mean_log[1], risk_sd_log[1], log=TRUE)    # NSP: Sharing needles
  lprior <- lprior+dlnorm(par_vector[,2], risk_mean_log[2], risk_sd_log[2], log=TRUE)    # NSP: Reusing needles
  lprior <- lprior+dlnorm(par_vector[,3], risk_mean_log[3], risk_sd_log[3], log=TRUE)    # Sharing needles
  lprior <- lprior+dlnorm(par_vector[,4], risk_mean_log[4], risk_sd_log[4], log=TRUE)    # Reusing needles
  lprior <- lprior+dlnorm(par_vector[,5], risk_mean_log[5], risk_sd_log[5], log=TRUE)    # Injection frequency 2x
  lprior <- lprior+dlnorm(par_vector[,6], risk_mean_log[6], risk_sd_log[6], log=TRUE)    # Injection frequency 7x
  lprior <- lprior+dlnorm(par_vector[,7], risk_mean_log[7], risk_sd_log[7], log=TRUE)    # Number of attempts 2
  lprior <- lprior+dlnorm(par_vector[,8], risk_mean_log[8], risk_sd_log[8], log=TRUE)    # Number of attempts 3
  lprior <- lprior+dlnorm(par_vector[,9], risk_mean_log[9], risk_sd_log[9], log = TRUE)    # Number of attempts 4
  lprior <- lprior+dlnorm(par_vector[,10], risk_mean_log[10], risk_sd_log[10], log = TRUE)    # IDU (2-5 years)
  lprior <- lprior+dlnorm(par_vector[,11], risk_mean_log[11], risk_sd_log[11], log = TRUE)    # IDU (5-8 years)
  lprior <- lprior+dlnorm(par_vector[,12], risk_mean_log[12], risk_sd_log[12], log = TRUE)    # IDU (8+ years)
  lprior <- lprior+dlnorm(par_vector[,13], risk_mean_log[13], risk_sd_log[13], log = TRUE)    # Age 25-44
  lprior <- lprior+dlnorm(par_vector[,14], risk_mean_log[14], risk_sd_log[14], log = TRUE)    # Age > 44
  lprior <- lprior+dlnorm(par_vector[,15], risk_mean_log[15], risk_sd_log[15], log = TRUE)    # Male
  return(lprior)
})

# | 3. Log Likelihood ####
l_likelihood <- cmpfun(function(par_vector) {
  # par_vector: a vector (or matrix) of model parameters
  if(is.null(dim(par_vector))) par_vector <- t(par_vector)
  llik <- rep(0, nrow(par_vector))
  for(j in 1:nrow(par_vector)) {
    jj <- tryCatch( {
      #Run model
      res_j <- mod(as.numeric(par_vector[j,])) 
      
      # Calculate likelihood for OCM
      llik[j] <- llik[j]+sum(dgamma(x = (res_j[["OCM_IR"]]),
                                    shape = 1,
                                    scale = 1 / 32.7646,
                                    log=TRUE))
      
      # Calculate likelihood for SSTVI deaths
      llik[j] <- llik[j]+sum(dgamma(x= (res_j[["SSTVI_IR"]]),
                                    shape = 1,
                                    scale = 1 / 5.5036,
                                    log=TRUE))
      
    }, error = function(e) NA)
    if(is.na(jj)) { llik[j] <- -Inf } 
  }
  return(llik)
})

# | 4. Log Posterior ####
l_post <- cmpfun(function(par_vector) {
  return(l_prior(par_vector) + l_likelihood(par_vector))
})

########## ---- Part C: Model Calibration Setup ##########
## With 1,000 patients per cohort, it takes about 1 hour 10 minutes to run 2000 samples.
## 100 samples: 4 minutes 20 seconds (22 unique samples)
## 500 samples: 16 minutes 40 seconds (131 unique samples)
## 2500 samples: 1 hour 18 minutes (657 unique samples) / 1 hour 14 minutes (659)

# | 1. Draw samples from Prior ####
set.seed(123)   # Set random seed for reproducibility
n.i <- 1000 # Run this first before running calibration (SIR or MAP)!!!
n_samples <- 2500 # Generate 500 draws from prior
m_samp <- sample_prior(n_samples) #matrix of samples from prior distribution

# | 2. Sampling Importance Resampling (SIR) ####
# || a. Set up parallel processing plan with 30 cores ####
sir_llik <- matrix(NA, nrow = n_samples, ncol = 1)
plan(multisession, workers = 40)
# || b. Run model on parallel computing platform ####
time_foreach <- system.time({
  sir_llik <- future_lapply(1:n_samples, function(i) {
    set.seed(i)
    l_likelihood(m_samp[i,])
  })
})
time_foreach[3]
# || c. Clean up and stop the workers ####
future:::plan("sequential")
# || d. Combine the results ####
sir_llik <- do.call(rbind, sir_llik)
# || e. Calculate likelihood of each sample ####
lik_samps <- exp(sir_llik - max(sir_llik)) 
# || f. Resample with replacement ####
m_resamples <- m_samp[sample(length(lik_samps),
                             size = length(lik_samps),
                             prob = lik_samps,
                             replace = TRUE),]
# || g. Calculate the number of unique parameter sets ####
SIR_unique_sets <- nrow(unique(m_resamples))
SIR_unique_sets
# || h. Matrix of unique parameter sets ####
posterior_params <- unique(m_resamples)
risk_multiplier_posterior <- colMeans(posterior_params)
# Calculate column medians for numeric columns
risk_multiplier_posterior <- apply(posterior_params, 2, median); risk_multiplier_posterior
# Calculate the 95% confidence intervals for each parameter
risk_multiplier_ci_lower <- apply(posterior_params, 2, quantile, probs = 0.025); risk_multiplier_ci_lower
risk_multiplier_ci_upper <- apply(posterior_params, 2, quantile, probs = 0.975); risk_multiplier_ci_upper

########## ---- Part D: Model Validation ##########
# | 1. OCM ####
targets <- "OCM"
obs_post <- 33.85; ll_post <- 33.33; ul_post <- 34.38
obs_pre <- 33.85; ll_pre <- 33.33; ul_pre <- 34.38
calib_values <- 32.76; ll_calib <- 31.82; ul_calib <- 33.72
# Create a data frame
df <- data.frame(targets, obs_post, ll_post, ul_post, obs_pre, ll_pre, ul_pre, 
                 calib_values, ll_calib, ul_calib)
# Create a plot
plot(1, type = "n", main = "Model vs. Target: Other-Cause Mortality", ylab = "", 
     xlab = "Mortality Rate (per 1,000 population)", 
     xlim = c(31, 35), ylim = c(0.5, 1.5), axes = FALSE)
# Plot horizontal error bars for each group
plot_error_bars(df$obs_post, 1, df$ll_post, df$ul_post, "black", 19, "Obs Post")
plot_error_bars(df$obs_pre, 1.2, df$ll_pre, df$ul_pre, "blue", 19, "Obs Pre")
plot_error_bars(df$calib_values, 1.4, df$ll_calib, df$ul_calib, "red", 19, "Calib Values")
axis(1, at = NULL)
# Add legend
legend("bottomleft", legend = c("Posterior", "Prior", "Calibration Target"),
       col = c("black", "blue", "red"), pch = 19, cex = 0.5, title = "Legend")

# | 2. SSTVI ####
targets <- c("SSTVI")
obs_post <- 5.57; ll_post <- 5.36; ul_post <- 5.79
obs_pre <- 5.60; ll_pre <- 5.39; ul_pre <- 5.82
calib_values <- 5.50; ll_calib <- 5.12; ul_calib <- 5.91
# Create a data frame
df <- data.frame(targets, obs_post, ll_post, ul_post, obs_pre, ll_pre, ul_pre, 
                 calib_values, ll_calib, ul_calib)
# Create a plot
plot(1, type = "n", main = "Model vs. Target: SSTVI Mortality", ylab = "", 
     xlab = "Mortality Rate (per 1,000 population)", 
     xlim = c(5, 6), ylim = c(0.5, 1.5), axes = FALSE)
# Plot horizontal error bars for each group
plot_error_bars(df$obs_post, 1, df$ll_post, df$ul_post, "black", 19, "Obs Post")
plot_error_bars(df$obs_pre, 1.2, df$ll_pre, df$ul_pre, "blue", 19, "Obs Pre")
plot_error_bars(df$calib_values, 1.4, df$ll_calib, df$ul_calib, "red", 19, "Calib Values")
axis(1, at = NULL)
# Add legend
legend("bottomleft", legend = c("Posterior", "Prior", "Calibration Target"),
       col = c("black", "blue", "red"), pch = 19, cex = 0.5, title = "Legend")

# End Script
