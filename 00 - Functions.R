# Jihoon Lim
# Economic Modelling - Microsimulation
# 00 - Functions

# Note: Une partie de la compilation est effectuée à partir de données provenant 
# du © Gouvernement du Québec (année de la publication du Fichier de recherche: 
# 2009-2019). Le © Gouvernement du Québec n’est pas responsable des compilations
# ni de l’interprétation des résultats produits à l’aide du Fichier de recherche.

library(dampack)
library(profvis)
library("compiler")
library(dplyr)
library("stringr")
library(tidyr)
library(survival)
library(survminer)
library(reda)
library(doParallel)
library(parallel)
library(foreach)
library(ggsurvfit)
library(cmprsk)

########## ---- Part A: Data Generation ##########
# | 1. Cohort Creation ####
data_generation <- cmpfun(function(n.i, nsp_status, par_vector) {
  # Input Parameters
  id <- seq(1:n.i)
  id <- as.numeric(id)
  
  ## Injection frequency in the past month
  inj.freq.dist <- c(0.129, 0.222, 0.156, 0.144, 0.349)
  names(inj.freq.dist) <- c("Never", "Not every week", "1-2 days a week",
                            "3-6 days a week", "Every day")
  v.inj.freq <- sample(x = 0:4, size = n.i, prob = inj.freq.dist, replace = T)
  v.inj.freq <- as.numeric(v.inj.freq)
  table(v.inj.freq)
  
  ## Number of injections in the past month
  inj.num.dist1 <- c(0.543, 0.457); inj.num.dist2 <- c(0.552, 0.448)
  inj.num.dist3 <- c(0.458, 0.542); inj.num.dist4 <- c(0.570, 0.430)
  v.inj.num <- vector("double", length = length(v.inj.freq))
  for (i in 1:length(v.inj.freq)){
    v.inj.num[i] <- ifelse(v.inj.freq[i] == 1, sample(1:2, size = 1, prob = inj.num.dist1, replace = T),
                           ifelse(v.inj.freq[i] == 2, sample(3:4, size = 1, prob = inj.num.dist2, replace = T),
                                  ifelse(v.inj.freq[i] == 3, sample(5:6, size = 1, prob = inj.num.dist3, replace = T),
                                         ifelse(v.inj.freq[i] == 4, sample(7:8, size = 1, prob = inj.num.dist4, replace = T),
                                                0)
                                  )
                           )
    )
  }
  
  ## Injection frequency in the past month
  inj.year.dist <- c(0.4012, 0.3555, 0.1816, 0.0616)
  names(inj.year.dist) <- c("< 2 years", "2-<5 years", "5-<8 years", "8+ years")
  v.yrs.idu <- sample(x = 1:4, size = n.i, prob = inj.year.dist, replace = T)
  v.yrs.idu <- as.numeric(v.yrs.idu)
  
  ## Risk factors for injection behaviour modification
  rr.sharing.nsp <- par_vector[1] # RR = 0.42 if NSP or 2.38 if no NSP
  rr.reusing.nsp <- par_vector[2] # RR = 0.79 if NSP or 1.27 if no NSP
  # Assume that the number of people who share needles/syringes increases by 2.38-fold
  # in a world without NSP. Similarly, 1.27 times the number of individuals who reuse
  # needles/syringes in a world without NSP. Use 95% CI for sensitivity analysis.
  
  if (nsp_status == 1) {
    sharing.prob <- 0.205
    reusing.prob <- 0.297
  } else {
    sharing.prob <- 0.205 * (1/rr.sharing.nsp)
    reusing.prob <- 0.297 * (1/rr.reusing.nsp)
  }
  
  ## Percent of PWID sharing needles/syringe with others
  v.sharing <- ifelse(v.inj.freq != 0, rbinom(n=sum(v.inj.freq >= 1), size=1, p=sharing.prob), 0)
  sum(v.sharing)
  
  ## Percent of PWID reusing their own needles/syringes
  v.reusing <- ifelse(v.inj.freq != 0, rbinom(n=sum(v.inj.freq >= 1), size=1, p=reusing.prob), 0)
  sum(v.reusing)
  
  ## Sex (0 = Female, 1 = Male)
  sex.dist <- c(0.3795, 0.6205)
  v.sex <- sample(x = 0:1, size = n.i, prob = sex.dist, replace = T)
  v.sex <- as.numeric(v.sex)
  
  ## Age
  # Define the desired percentages
  percentages <- c(0.1065, 0.5108, 1)
  # Find the quantiles for the specified percentages
  quantiles <- qnorm(percentages, mean = 42.711, sd = 12.822)
  # Generate ages based on the quantiles
  v.age <- round(rnorm(n.i, 42.711, 12.822), 0)
  v.age <- ifelse(v.age < quantiles[1], runif(sum(v.age < quantiles[1]), 18, 24),
                  ifelse(v.age < quantiles[2], runif(sum(v.age < quantiles[2]), 25, 44),
                         runif(sum(v.age >= quantiles[2]), 45, 65)))
  
  ## Mortality rate
  # Note: The rates for ages 60+ reflect rates combined from the 60-64 and 65-69 
  # groups. Individuals may start at ages between 60 and 64 but die when they are 
  # 65+ years old. This may result in small cell counts in the latter group
  # (i.e., 65-69), and we took this approach to account for this possibility.
  v.mortality <- ifelse(v.age >= 18 & v.age <= 24, 0.0040, 
                        ifelse(v.age >= 25 & v.age <= 29, 0.0068,
                               ifelse(v.age >= 30 & v.age <= 34, 0.0089,
                                      ifelse(v.age >= 35 & v.age <= 39, 0.0130, 
                                             ifelse(v.age >= 40 & v.age <= 44, 0.0181,
                                                    ifelse(v.age >= 45 & v.age <= 49, 0.0336,
                                                           ifelse(v.age >= 50 & v.age <= 54, 0.0526,
                                                                  ifelse(v.age >= 55 & v.age <= 59, 0.0703, 0.0622)
                                                           )
                                                    )
                                             )
                                      )
                               )
                        )
  )
  v.ocm <- round((1 - exp(-1*v.mortality/52)), 6)
  
  ## NSP Status
  if (nsp_status == 1) {
    v.nsp <- rep(1, n.i)
  } else {
    v.nsp <- rep(0, n.i)
  }
  
  # Risk Multipliers
  # Because interaction between these factors is unknown, take the maximum of the risk
  # ratios as the risk multiplier applicable for that particular patient. This will enable
  # us to obtain a more conservative estimate of the costs incurred and the ICER/NMB
  # generated from the data. Use 95% CI of the max RR for sensitivity analysis.
  
  # Dataset
  ## Create a data frame for PWID baseline characteristics
  pwid.baseline <- data.frame(id, v.inj.freq, v.inj.num, v.yrs.idu, v.sharing, v.reusing, v.sex, v.age, v.ocm, v.nsp)
  
  ## Add relative risk variables associated with SSTVI into the dataset
  pwid.baseline <- pwid.baseline |>
    dplyr::mutate(rr_sharing = if_else(v.sharing == 1, par_vector[3], 1),
                  rr_reusing = if_else(v.reusing == 1, par_vector[4], 1),
                  rr_inj_freq = case_when(
                    v.inj.freq %in% c(1, 2) ~ par_vector[5],
                    v.inj.freq %in% c(3, 4) ~ par_vector[6],
                    TRUE ~ 1
                  ),
                  rr_inj_num = case_when(
                    v.inj.num %in% c(3, 4) ~ par_vector[7],
                    v.inj.num %in% c(5, 6) ~ par_vector[8],
                    v.inj.num %in% c(7, 8) ~ par_vector[9],
                    TRUE ~ 1
                  ),
                  rr_inj_year = case_when(
                    v.yrs.idu == 2 ~ par_vector[10],
                    v.yrs.idu == 3 ~ par_vector[11],
                    v.yrs.idu == 4 ~ par_vector[12],
                    TRUE ~ 1
                  ),
                  rr_age = case_when(
                    v.age >= 25 & v.age <= 44 ~ par_vector[13],
                    v.age > 44 ~ par_vector[14],
                    TRUE ~ 1
                  ),
                  rr_sex = if_else(v.sex == 1, par_vector[15], 1)
    )
  
  ## Analysis Dataset
  pwid.baseline <- pwid.baseline |>
    mutate(risk.multiplier(data_name = pwid.baseline, prob = 0.001458, prev = 0.001452))
})

# | 2. Risk Multiplier Function ####
# Function to multiply the largest and three smallest in each row
multiplier_algorithm <- function(variables) {
  sort_desc_order <- sort(variables, decreasing = TRUE)
  result <- sort_desc_order[1] * sort_desc_order[4] * sort_desc_order[5] * sort_desc_order[6] * sort_desc_order[7]
  return(result)
}

risk.multiplier <- cmpfun(function(data_name, prob, prev){
  v_rand <- runif(n=nrow(data_name))
  data_name <- data_name |>
    dplyr::mutate(rr_multiplier = apply(data_name[, c("rr_sharing", "rr_reusing", "rr_inj_freq", "rr_inj_num", "rr_inj_year", "rr_age", "rr_sex")], 1, multiplier_algorithm),
                  prob_infection = prob*rr_multiplier,
                  prevalent = rbinom(nrow(data_name), 1, prev),
                  infection = ifelse(prevalent == 0, as.numeric(v_rand < prob_infection), 1),
                  state = ifelse(infection == 1, "SSTVI", "Healthy"))
}) # Creates different probability of SSTVI for each individual based on baseline characteristics

########## ---- Part B: Microsimulation ##########
# | 1. Sampling Next State from Transition Probabilities ####
samplev <- cmpfun(function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-03))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
})

# | 2. Main Microsimulation Function ####
MicroSim <- cmpfun(function(v.M_1, n.i, n.t, n.years, v.n, v.inf.prob, v.ocm, d.c, d.e, nsp, TS.out = TRUE, seed = 1) {
  
  v.dwc <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))    # Calculate the cost discount weight based on the discount rate d.c 
  v.dwe <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))    # Calculate the QALY discount weight based on the discount rate d.e
  
  # Create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <- matrix(nrow = n.i, ncol = n.t, 
                              dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                              paste("cycle", 0:(n.t-1), sep = " ")))  
  
  m.M[, 1] <- v.M_1                  # Indicate the initial health state   
  set.seed(seed)                  # Set the seed for every individual for the random number generator
  m.C[, 1] <- Costs(m.M[, 1])  # Estimate costs per individual for the initial health state
  m.E[, 1] <- Effs (m.M[, 1])  # Estimate QALYs per individual for the initial health state
  
  for (t in 1:(n.t-1)) { # t <- 3
    m.p <- Probs(m.M[, t], v.inf.prob, v.ocm)           # Calculate the transition probabilities at cycle t
    m.M[, t + 1] <- samplev(prob = m.p, m = 1)  # Sample the next health state and store that state in matrix m.M
    m.C[, t + 1] <- Costs(m.M[, t + 1])   # Estimate costs per individual during cycle t + 1 conditional on treatment
    m.E[, t + 1] <- Effs(m.M[, t + 1])   # Estimate QALYs per individual during cycle t + 1 conditional on treatment
    cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))       # Display the progress of the simulation
  } # close the loop for the time points
  
  v.dwc <- as.matrix(as.data.frame(v.dwc))
  v.dwe <- as.matrix(as.data.frame(v.dwe))
  tc <- m.C %*% v.dwc     # Total (discounted) cost per individual
  te <- m.E %*% v.dwe       # Total (discounted) QALYs per individual 
  
  tc_hat <- mean(tc)        # Average (discounted) cost 
  te_hat <- mean(te)        # Average (discounted) QALYs
  
  if (TS.out == TRUE) {  # Create a matrix of transitions across states
    TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # Transitions from one state to the other
    TS <- matrix(TS, nrow = n.i)
    rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # Name the rows 
    colnames(TS) <- paste("Cycle", 0:(n.t-1), sep = " ")   # Name the columns 
  } else {
    TS <- NULL
  }
  
  results <- list(m.p = m.p, m.M = m.M, m.C = m.C, m.E = m.E, tc = tc, te = te, 
                  tc_hat = tc_hat, te_hat = te_hat, TS = TS) # Store the results from the simulation in a list  
  return(results)  # Return the results
})  # End of the MicroSim function  

# | 3. Transition Probabilities Function ####
# This function updates the transition probabilities at every cycle.
Probs <- cmpfun(function(M_it, v.inf.prob, v.ocm) { 
  # M_it:   health state occupied by individual i at cycle t (character variable)
  m.p.it <- matrix(NA, n.s, n.i)     # Create vector of state transition probabilities
  rownames(m.p.it) <- v.n            # Assign names to the vector
  
  # Use for loop to update the matrix with the appropriate probabilities and account for 
  # different individual probabilities of SSTVI
  for (i in 1:n.i){
    m.p.it[,M_it == "Healthy"] <- t(c((1 - v.inf.prob[i])*(1-v.ocm[i]),
                                      v.inf.prob[i]*(1-v.ocm[i]),
                                      rep(0, 10),
                                      v.ocm[i]))
    m.p.it[,M_it == "SSTVI"] <- c((1 - p_SSTVI_to_SelfTrt_Death - p_SSTVI_to_SSTVI - p_SSTVI_to_OP_Purulent - p_SSTVI_to_ED_Purulent - p_SSTVI_to_OP_NonPurulent - p_SSTVI_to_ED_NonPurulent)*(1-v.ocm[i]),
                                  p_SSTVI_to_SSTVI*(1-v.ocm[i]),
                                  p_SSTVI_to_SelfTrt_Death*(1-v.ocm[i]),
                                  p_SSTVI_to_OP_Purulent*(1-v.ocm[i]),
                                  p_SSTVI_to_ED_Purulent*(1-v.ocm[i]),
                                  rep(0, 2),
                                  p_SSTVI_to_OP_NonPurulent*(1-v.ocm[i]),
                                  p_SSTVI_to_ED_NonPurulent*(1-v.ocm[i]),
                                  rep(0, 3),
                                  v.ocm[i]) # Transition probabilities when SSTVI
    m.p.it[,M_it == "SelfTrt_Death"] <- c(rep(0, 2),
                                          1,
                                          rep(0, 10)) # Transition probabilities when death from SSTVI
    m.p.it[,M_it == "Purulent_OP"] <- c((1 - p_OP_to_IP_Purulent - p_OP_to_IPC_Purulent)*(1-v.ocm[i]),
                                        rep(0, 4),
                                        p_OP_to_IP_Purulent*(1-v.ocm[i]),
                                        p_OP_to_IPC_Purulent*(1-v.ocm[i]),
                                        rep(0, 5),
                                        v.ocm[i]) # Transition probabilities when outpatient setting for purulent SSTVI
    m.p.it[,M_it == "Purulent_ED"] <- c((1 - p_ED_to_SSTVI_Purulent - p_ED_to_IP_Purulent - p_ED_to_IPC_Purulent - p_ED_to_SSTVI_Death_Purulent)*(1-v.ocm[i]),
                                        p_ED_to_SSTVI_Purulent*(1-v.ocm[i]),
                                        rep(0, 3),
                                        p_ED_to_IP_Purulent*(1-v.ocm[i]),
                                        p_ED_to_IPC_Purulent*(1-v.ocm[i]),
                                        rep(0, 4),
                                        p_ED_to_SSTVI_Death_Purulent*(1-v.ocm[i]),
                                        v.ocm[i]) # Transition probabilities when ED setting for purulent SSTVI
    m.p.it[,M_it == "Purulent_IP"] <- c((1 - p_IP_to_SSTVI_Purulent - p_IP_to_IP_Purulent - p_IP_to_SSTVI_Death_Purulent)*(1-v.ocm[i]),
                                        p_IP_to_SSTVI_Purulent*(1-v.ocm[i]),
                                        rep(0, 3),
                                        p_IP_to_IP_Purulent*(1-v.ocm[i]),
                                        rep(0, 5),
                                        p_IP_to_SSTVI_Death_Purulent*(1-v.ocm[i]),
                                        v.ocm[i]) # Transition probabilities when inpatient setting for purulent SSTVI
    m.p.it[,M_it == "Purulent_IPC"] <- c((1 - p_IPC_to_SSTVI_Purulent - p_IPC_to_SSTVI_Death_Purulent)*(1-v.ocm[i]),
                                         p_IPC_to_SSTVI_Purulent*(1-v.ocm[i]),
                                         rep(0, 9),
                                         p_IPC_to_SSTVI_Death_Purulent*(1-v.ocm[i]),
                                         v.ocm[i]) # Transition probabilities when complications in inpatient setting (purulent)
    m.p.it[,M_it == "NonPurulent_OP"] <- c((1 - p_OP_to_IP_NonPurulent - p_OP_to_IPC_NonPurulent)*(1-v.ocm[i]),
                                           rep(0, 8),
                                           p_OP_to_IP_NonPurulent*(1-v.ocm[i]),
                                           p_OP_to_IPC_NonPurulent*(1-v.ocm[i]),
                                           0,
                                           v.ocm[i]) # Transition probabilities when outpatient setting for non-purulent SSTVI
    m.p.it[,M_it == "NonPurulent_ED"] <- c((1 - p_ED_to_SSTVI_NonPurulent - p_ED_to_IP_NonPurulent - p_ED_to_IPC_NonPurulent - p_ED_to_SSTVI_Death_NonPurulent)*(1-v.ocm[i]),
                                           p_ED_to_SSTVI_NonPurulent*(1-v.ocm[i]),
                                           rep(0, 7),
                                           p_ED_to_IP_NonPurulent*(1-v.ocm[i]),
                                           p_ED_to_IPC_NonPurulent*(1-v.ocm[i]),
                                           p_ED_to_SSTVI_Death_NonPurulent*(1-v.ocm[i]),
                                           v.ocm[i]) # Transition probabilities when ED setting for non-purulent SSTVI
    m.p.it[,M_it == "NonPurulent_IP"] <- c((1 - p_IP_to_SSTVI_NonPurulent - p_IP_to_IP_NonPurulent - p_IP_to_SSTVI_Death_NonPurulent)*(1-v.ocm[i]),
                                           p_IP_to_SSTVI_NonPurulent*(1-v.ocm[i]),
                                           rep(0, 7),
                                           p_IP_to_IP_NonPurulent*(1-v.ocm[i]),
                                           0,
                                           p_IP_to_SSTVI_Death_NonPurulent*(1-v.ocm[i]),
                                           v.ocm[i]) # Transition probabilities when inpatient setting for non-purulent SSTVI
    m.p.it[,M_it == "NonPurulent_IPC"] <- c((1 - p_IPC_to_SSTVI_NonPurulent - p_IPC_to_SSTVI_Death_NonPurulent)*(1-v.ocm[i]),
                                            p_IPC_to_SSTVI_NonPurulent*(1-v.ocm[i]),
                                            rep(0, 9),
                                            p_IPC_to_SSTVI_Death_NonPurulent*(1-v.ocm[i]),
                                            v.ocm[i]) # Transition probabilities when complications in inpatient setting (non-purulent)
    m.p.it[,M_it == "SSTVI_Death"] <- c(rep(0, 11), 1, 0) # Transition probabilities when death from OCM
    m.p.it[,M_it == "OCM"] <- c(rep(0, 12), 1) # Transition probabilities when death from OCM
  }
  # Return the transition probabilities or produce an error
  ifelse(colSums(m.p.it) - 1 < 1e-5, return(t(m.p.it)), print("Probabilities do not sum to 1")) 
})  

# | 4. Cost Function ####
# This function estimates the costs at every cycle.
Costs <- cmpfun(function (M_it) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  # Default cost for everyone is zero
  c.it <- 0
  
  # Update the costs as follows
  c.it[M_it == "Healthy"] <- c_Healthy # Healthy
  c.it[M_it == "SSTVI"] <- c_SSTVI # SSTVI
  c.it[M_it == "SelfTrt_Death"] <- c_SelfTrt_Death # Death from SSTVI (including IE)
  c.it[M_it == "Purulent_OP"] <- c_OP_Purulent # Outpatient (purulent)
  c.it[M_it == "Purulent_ED"] <- c_ED_Purulent # ED (purulent)
  c.it[M_it == "Purulent_IP"] <- c_IP_Purulent # Inpatient (purulent)
  c.it[M_it == "Purulent_IPC"] <- c_IPC_Purulent # Inpatient (purulent) involving surgeries
  c.it[M_it == "NonPurulent_OP"] <- c_OP_NonPurulent # Outpatient (non-purulent)
  c.it[M_it == "NonPurulent_ED"] <- c_ED_NonPurulent # ED (non-purulent)
  c.it[M_it == "NonPurulent_IP"] <- c_IP_NonPurulent # Inpatient (non-purulent)
  c.it[M_it == "NonPurulent_IPC"] <- c_IPC_NonPurulent # Inpatient (non-purulent) involving surgeries
  c.it[M_it == "SSTVI_Death"] <- c_SSTVI_Death # Death from SSTVI
  c.it[M_it == "OCM"] <- c_OCM # Death from OCM
  
  # Return the costs
  return(c.it)
})


# | 5. Health Outcomes Function ####
# This function updates the utilities at every cycle.
Effs <- cmpfun(function (M_it) {
  # M_it: health state occupied by individual i at cycle t (character variable)
  
  # Default utility for everyone is zero
  u.it <- 0
  
  # Update the utilities as follows
  u.it[M_it == "Healthy"] <- q_Healthy # Healthy
  u.it[M_it == "SSTVI"] <- q_SSTVI # SSTVI
  u.it[M_it == "SelfTrt_Death"] <- q_SelfTrt_Death # Death from SSTVI including IE
  u.it[M_it == "Purulent_OP"] <- q_OP_Purulent # Outpatient (purulent)
  u.it[M_it == "Purulent_ED"] <- q_ED_Purulent # ED (purulent)
  u.it[M_it == "Purulent_IP"] <- q_IP_Purulent # Inpatient (purulent)
  u.it[M_it == "Purulent_IPC"] <- q_IPC_Purulent # Inpatient (purulent) involving surgeries
  u.it[M_it == "NonPurulent_OP"] <- q_OP_NonPurulent # Outpatient (non-purulent)
  u.it[M_it == "NonPurulent_ED"] <- q_ED_NonPurulent # ED (non-purulent)
  u.it[M_it == "NonPurulent_IP"] <- q_IP_NonPurulent # Inpatient (non-purulent)
  u.it[M_it == "NonPurulent_IPC"] <- q_IPC_NonPurulent # Inpatient (non-purulent) involving surgeries
  u.it[M_it == "SSTVI_Death"] <- q_SSTVI_Death # Death from SSTVI
  u.it[M_it == "OCM"] <- q_OCM # Death from OCM
  
  return(u.it)                  # return the QALYs
})

########## ---- Part C: Data Analysis ##########
# | 1. Calculation of Time of Death ####
death_calculation <- cmpfun(function(simulation_group, nsp_status) {
  # Identify when the person dies
  transition.path <- simulation_group$TS
  transition.path <- as.data.frame(transition.path)
  transition.ocm <- apply(transition.path, 1, function(x) match(TRUE, str_sub(x, -3, -1) == "OCM"))
  transition.sstvi.death <- apply(transition.path, 1, function(x) match(TRUE, str_sub(x, -5, -1) == "Death"))
  
  # Vectorize the death times (other-cause mortality, SSI, and complications)
  transition.ocm <- as.vector(transition.ocm, mode = "numeric")
  transition.sstvi.death <- as.vector(transition.sstvi.death, mode = "numeric")
  
  # Create data frame
  death.dates <- data.frame(OCM = transition.ocm,
                            SSTVI_Death = transition.sstvi.death)
  death.dates$Death <- apply(death.dates[,c(1:2)], 1, min, na.rm = TRUE)
  death.dates$Death <- ifelse(death.dates$Death == Inf, n.t, death.dates$Death)
  death.dates$ocm_status <- ifelse(is.na(death.dates$OCM) == TRUE & is.na(death.dates$SSTVI_Death) == TRUE, 0,
                                   ifelse(is.na(death.dates$OCM) == FALSE & is.na(death.dates$SSTVI_Death) == TRUE, 1, 2))
  death.dates$sstvi_status <- ifelse(is.na(death.dates$OCM) == TRUE & is.na(death.dates$SSTVI_Death) == TRUE, 0,
                                     ifelse(is.na(death.dates$OCM) == FALSE & is.na(death.dates$SSTVI_Death) == TRUE, 2, 1))
  death.dates$stop_t <- ifelse(is.na(death.dates$OCM) == TRUE & is.na(death.dates$SSTVI_Death) == TRUE, 260,
                               ifelse(is.na(death.dates$OCM) == FALSE & is.na(death.dates$SSTVI_Death) == TRUE,
                                      death.dates$OCM, death.dates$SSTVI_Death))
  
  if (nsp_status == 1) {
    death.dates$id <- seq(1, n.i)
    death.dates$nsp <- rep(1, n.i)
  } else {
    death.dates$id <- seq(100001, (100000 + n.i))
    death.dates$nsp <- rep(0, n.i)
  }
  
  # Person-time for each year
  death.dates$yr1 <- pmin(death.dates$Death / 52, 1)
  death.dates$yr2 <- pmin((death.dates$Death - 52) / 52, 1)
  death.dates$yr2[death.dates$Death <= 52] <- 0
  death.dates$yr3 <- pmin((death.dates$Death - 104) / 52, 1)
  death.dates$yr3[death.dates$Death <= 104] <- 0
  death.dates$yr4 <- pmin((death.dates$Death - 156) / 52, 1)
  death.dates$yr4[death.dates$Death <= 156] <- 0
  death.dates$yr5 <- pmin((death.dates$Death - 208) / 52, 1)
  death.dates$yr5[death.dates$Death <= 208] <- 0
  death.dates$yr5[death.dates$Death >= 260] <- 1
  return(death.dates)
})

# | 2a. Counting Process Function - Single Pathway ####
counting.process.single <- cmpfun(function(simulation_group, transition_text, pwid_baseline_data) {
  # Identify when the person dies
  transition.path <- simulation_group$TS
  transition.path <- as.data.frame(transition.path)
  transition.ocm <- apply(transition.path, 1, function(x) match(TRUE, str_sub(x, -3, -1) == "OCM"))
  transition.sstvi.death <- apply(transition.path, 1, function(x) match(TRUE, str_sub(x, -5, -1) == "Death"))
  
  # Vectorize the death times (other-cause mortality, SSI, and complications)
  transition.ocm <- as.vector(transition.ocm, mode = "numeric")
  transition.sstvi.death <- as.vector(transition.sstvi.death, mode = "numeric")
  
  # Create data frame
  death.dates <- data.frame(OCM = transition.ocm,
                            SSTVI_Death = transition.sstvi.death)
  death.dates$Death <- apply(death.dates[,c(1:2)], 1, min, na.rm = TRUE)
  death.dates$Death <- ifelse(death.dates$Death == Inf, n.t, death.dates$Death)
  death.dates$id <- seq(1:n.i)
  
  # Identify occurrences of event (OP, ED, and IP)
  index_transition <- (transition.path == transition_text) + 0
  
  # Put the occurrences of event in a list
  index_transition_time <- lapply(1:n.i, function(i) which(index_transition[i, ] == 1))
  
  # Put event time in a matrix format
  n.obs <- sapply(index_transition_time, length)
  seq.max <- seq_len(max(n.obs))
  recurrent <- t(sapply(index_transition_time, "[", i = seq.max))
  recurrent <- as.data.frame(recurrent)
  colnames(recurrent) <- paste("V", 1:ncol(recurrent), sep = "")
  
  # Transform data into long format
  recurrent_long1 <- recurrent %>% 
    pivot_longer(
      cols = everything(), # This could change depending on the number of events 
      names_to = "event_num",
      values_to = "t1"
    )
  recurrent_long1$id <- rep(1:n.i, each = ncol(recurrent))
  recurrent_long1$event <- ifelse(is.na(recurrent_long1$t1), 0, 1)
  recurrent_long2 <- left_join(recurrent_long1 %>% select(id, t1, event), 
                               death.dates %>% select(id, Death), 
                               by = "id")
  
  # Create t0 and t1 variables
  recurrent_long3 <- recurrent_long2 %>% 
    dplyr::filter(event != 0) %>%  
    group_by(id) %>%
    mutate(t0 = ifelse(row_number() == 1, 0, dplyr::lag(t1))) %>% 
    # mutate(start = lag(t1)) %>% 
    # mutate(start = ifelse(is.na(start), 0, start)) %>%
    ungroup() %>% 
    select(id, t0, t1, event, Death)
  
  # Include last row for each patient before getting censored or died
  recurrent_long4 <- recurrent_long3 %>% 
    bind_rows(recurrent_long3 %>% group_by(id) %>% 
                slice(n()) %>% 
                mutate(t0 = t1,
                       t1 = Death,
                       event = 0) %>% 
                ungroup()) %>% 
    arrange(id, t0, t1, event) %>% 
    select(id, t0, t1, event)
  
  # Include baseline characteristics from data generating step
  recurrent_long5 <- left_join(recurrent_long4, 
                               pwid_baseline_data %>% 
                                 select(id, v.inj.freq, v.inj.num, v.sharing, v.reusing, v.age, v.nsp), 
                               by = "id")
  
  return(recurrent_long5)
})

# | 2b. Counting Process Function - Multiple Pathways ####
counting.process.multiple <- cmpfun(function(simulation_group, transition_text1, transition_text2, pwid_baseline_data) {
  # Identify when the person dies
  transition.path <- simulation_group$TS
  transition.path <- as.data.frame(transition.path)
  transition.ocm <- apply(transition.path, 1, function(x) match(TRUE, str_sub(x, -3, -1) == "OCM"))
  transition.sstvi.death <- apply(transition.path, 1, function(x) match(TRUE, str_sub(x, -5, -1) == "Death"))
  
  # Vectorize the death times (other-cause mortality, SSI, and complications)
  transition.ocm <- as.vector(transition.ocm, mode = "numeric")
  transition.sstvi.death <- as.vector(transition.sstvi.death, mode = "numeric")
  
  # Create data frame
  death.dates <- data.frame(OCM = transition.ocm,
                            SSTVI_Death = transition.sstvi.death)
  death.dates$Death <- apply(death.dates[,c(1:2)], 1, min, na.rm = TRUE)
  death.dates$Death <- ifelse(death.dates$Death == Inf, n.t, death.dates$Death)
  death.dates$id <- seq(1:n.i)
  
  # Identify occurrences of event (OP, ED, and IP)
  index_transition <- (transition.path == transition_text1 | transition.path == transition_text2) + 0
  
  # Put the occurrences of event in a list
  index_transition_time <- lapply(1:n.i, function(i) which(index_transition[i, ] == 1))
  
  # Put event time in a matrix format
  n.obs <- sapply(index_transition_time, length)
  seq.max <- seq_len(max(n.obs))
  recurrent <- t(sapply(index_transition_time, "[", i = seq.max))
  if (length(seq.max) > 1) {
    recurrent <- as.data.frame(recurrent)
  } else {
    recurrent <- as.data.frame(t(recurrent)) 
  }
  colnames(recurrent) <- paste("V", 1:ncol(recurrent), sep = "")
  
  # Transform data into long format
  recurrent_long1 <- recurrent %>% 
    pivot_longer(
      cols = everything(), # This could change depending on the number of events 
      names_to = "event_num",
      values_to = "t1"
    )
  
  recurrent_long1$id <- rep(1:n.i, each = ncol(recurrent))
  recurrent_long1$event <- ifelse(is.na(recurrent_long1$t1), 0, 1)
  recurrent_long2 <- left_join(recurrent_long1 %>% select(id, t1, event), 
                               death.dates %>% select(id, Death), 
                               by = "id")
  
  # Create t0 and t1 variables
  recurrent_long3 <- recurrent_long2 %>% 
    dplyr::filter(event != 0) %>%  
    group_by(id) %>%
    mutate(t0 = ifelse(row_number() == 1, 0, dplyr::lag(t1))) %>% 
    # mutate(start = lag(t1)) %>% 
    # mutate(start = ifelse(is.na(start), 0, start)) %>%
    ungroup() %>% 
    select(id, t0, t1, event, Death)
  
  # Include last row for each patient before getting censored or died
  recurrent_long4 <- recurrent_long3 %>% 
    bind_rows(recurrent_long3 %>% group_by(id) %>% 
                slice(n()) %>% 
                mutate(t0 = t1,
                       t1 = Death,
                       event = 0) %>% 
                ungroup()) %>% 
    arrange(id, t0, t1, event) %>% 
    select(id, t0, t1, event)
  
  # Include baseline characteristics from data generating step
  recurrent_long5 <- left_join(recurrent_long4, 
                               pwid_baseline_data %>% 
                                 select(id, v.inj.freq, v.inj.num, v.sharing, v.reusing, v.age, v.nsp), 
                               by = "id")
  
  return(recurrent_long5)
})

# | 3a. Counting Process for Each Recurrent Outcome - Single Pathway ####
recurrent_event_data1 <- cmpfun(function (trt_output, no_trt_output, transition_text) {
  recurrent_with_NSP <- counting.process.single(trt_output, transition_text, pwid.baseline.with.nsp)
  recurrent_no_NSP <- counting.process.single(no_trt_output, transition_text, pwid.baseline.no.nsp)
  recurrent_no_NSP$id <- recurrent_no_NSP$id + n.i
  recurrent_data <- bind_rows(recurrent_with_NSP, recurrent_no_NSP)
  return(recurrent_data)
})

# | 3b. Counting Process for Each Recurrent Outcome - Multiple Pathways ####
recurrent_event_data2 <- cmpfun(function (trt_output, no_trt_output, transition_text1, transition_text2) {
  recurrent_with_NSP <- counting.process.multiple(trt_output, transition_text1, transition_text2, pwid.baseline.with.nsp)
  recurrent_no_NSP <- counting.process.multiple(no_trt_output, transition_text1, transition_text2, pwid.baseline.no.nsp)
  recurrent_no_NSP$id <- recurrent_no_NSP$id + n.i
  recurrent_data <- bind_rows(recurrent_with_NSP, recurrent_no_NSP)
  return(recurrent_data)
})

# | 4. Mean Cumulative Function Plot ####
mcf_plot <- cmpfun(function (recurrent_data, recurrent_outcome) {
  ## Internal knots are set as 33% and 67% quantiles of time variable
  (splineFit <- rateReg(Recur(t1, id, event) ~ v.nsp, data = recurrent_data,
                        df = 6, degree = 3, spline = "mSplines"))
  newDat <- data.frame(v.nsp = c(1, 0), group = c("Treat", "Contr"))
  estmcf <- mcf(splineFit, newdata = newDat, groupName = "Group",
                groupLevels = c("NSP", "No NSP"))
  mcf_plot <- plot(estmcf) +
    geom_ribbon(data = estmcf@MCF, alpha = 0.2,
                aes(x = time, ymin = lower, ymax = upper, fill = Group)) +
    ggtitle(substitute(paste("Mean Cumulative Count (", recurrent_outcome, ")"))) + xlab("Weeks")
  return(mcf_plot)
})

# | 5. Marginal Means Model ####
marginal_means <- cmpfun(function (recurrent_data) {
  marginal_model <- coxph(Surv(t0, t1, event) ~ as.factor(v.nsp) + cluster(id), 
                          data=recurrent_data[,c(1:4, 10)])
  marginal_model_summary <- summary(marginal_model)
  hr_estimates <- marginal_model_summary$conf.int
  return(hr_estimates)
})

########## ---- Part D: Probabilistic Sensitivity Analysis ##########
# | 1. Beta Distribution Functions ####
# Adapted from prevalence R package
#.  https://github.com/cran/prevalence/blob/master/R/betaExpert.R
betaExpert <- function(best, lower, upper, p = 0.95, method = "mean"){
  ## check presence
  if (missing(best))
    stop("'best' is missing")
  if (missing(lower) & missing(upper))
    stop("at least 'lower' or 'upper' must be specified")
  
  ## check input values: order
  if (!missing(lower))
    if (lower > best) stop("'lower' cannot be greater than 'best'")
  if (!missing(upper))
    if (upper < best) stop("'upper' cannot be smaller than 'best'")
  if (!missing(lower) & !missing(upper)) # useless??
    if (lower > upper) stop("'lower' cannot be greater than 'upper'")
  
  ## functions to optimize ~ mode
  f_mode <-
    function(x, mode, p, target){
      return(
        sum(
          (qbeta(p = p,
                 shape1 = x,
                 shape2 = (x * (1 - mode) + 2 * mode - 1) / mode) -
             target) ^ 2
        ))
    }
  
  f_mode_zero <-
    function(x, p, target){
      return((qbeta(p = p, shape1 = 1, shape2 = x) - target) ^ 2)
    }
  
  f_mode_one <-
    function(x, p, target){
      return((qbeta(p = p, shape1 = x, shape2 = 1) - target) ^ 2)
    }
  
  ## functions to optimize ~ mean
  f_mean <-
    function(x, mean, p, target){
      return(
        sum(
          (qbeta(p = p,
                 shape1 = x,
                 shape2 = (x * (1 - mean)) / mean) -
             target) ^ 2
        ))
    }
  
  ## define 'target' and 'p'
  if (!missing(lower) & missing(upper)){
    target <- lower
    p <- 1 - p
  } else if (!missing(upper) & missing(lower)){
    target <- upper
  } else if (!missing(upper) & !missing(lower)){
    target <- c(lower, upper)
    p <- c(0, p) + (1 - p) / 2
  }
  
  ## derive a and b (=shape1 and shape2)
  if (method == "mode"){
    if (best == 0){
      a <- 1
      b <- optimize(f_mode_zero, c(0, 1000), p = p, target = target)$minimum
    } else if (best == 1) {
      a <- optimize(f_mode_one, c(0, 1000), p = p, target = target)$minimum
      b <- 1
    } else {
      a <- optimize(f_mode, c(0, 1000),
                    mode = best, p = p, target = target)$minimum
      b <- (a * (1 - best) + 2 * best - 1) / best
    }
  } else if (method == "mean"){
    a <- optimize(f_mean, c(0, 1000),
                  mean = best, p = p, target = target)$minimum
    b <- (a * (1 - best)) / best
  }
  
  ## create 'out' dataframe
  out <- list(alpha = a, beta = b)
  class(out) <- "betaExpert"
  
  ## return 'out'
  return(out)
}

rbeta_mean_quantile <- function(mean, lb, ub, n=1){
  params <- betaExpert(mean, lb, ub)
  return(rbeta(n=n, params$alpha, params$beta))
}

# | 2. Gamma Distribution Functions ####
gammaExpert <- function(best, lower, upper, p = 0.95, method = "mean"){
  ## check presence
  if (missing(best))
    stop("'best' is missing")
  if (missing(lower) & missing(upper))
    stop("at least 'lower' or 'upper' must be specified")
  
  ## check input values: order
  if (!missing(lower))
    if (lower > best) stop("'lower' cannot be greater than 'best'")
  if (!missing(upper))
    if (upper < best) stop("'upper' cannot be smaller than 'best'")
  if (!missing(lower) & !missing(upper)) # useless??
    if (lower > upper) stop("'lower' cannot be greater than 'upper'")
  
  
  f_mean <-
    function(x, mean, p, target){
      return(
        sum(
          (qgamma(p = p,
                  shape = x,
                  scale = mean/x) -
             target) ^ 2
        ))
    }
  
  ## define 'target' and 'p'
  if (!missing(lower) & missing(upper)){
    target <- lower
    p <- 1 - p
  } else if (!missing(upper) & missing(lower)){
    target <- upper
  } else if (!missing(upper) & !missing(lower)){
    target <- c(lower, upper)
    p <- c(0, p) + (1 - p) / 2
  }
  
  ## derive a and b (=shape1 and shape2)
  
  if (method == "mean"){
    a <- optimize(f_mean, c(0, 1000),
                  mean = best, p = p, target = target)$minimum #a is shape parameter
    b <- best / a #b is scale parameter
  }
  
  ## create 'out' dataframe
  out <- list(shape = a, scale = b)
  #class(out) <- "betaExpert"
  
  ## return 'out'
  return(out)
}

rgamma_mean_quantile <- function(mean, lb, ub, n=1){
  params <- gammaExpert(mean, lb, ub)
  return(rgamma(n=n, shape=params$shape, scale=params$scale))
}

# | 3. Manual Creation of get.lnorm.par Function ####
is.error <- function(x) inherits(x, "try-error")
get.lnorm.par <- function(p = c(0.025, 0.5, 0.975), q, 
                          show.output = TRUE, plot = TRUE, 
                          tol = 0.001, fit.weights = rep(1, length(p)), 
                          scaleX = c(0.1, 0.9), ...) {
  #-----------------------------------------------------------------------------
  # checking consistency of the input data
  #-----------------------------------------------------------------------------
  if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
  }
  if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
  }
  if (min(p) < 0 | max(p) > 1) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
  }
  if (min(q) <= 0) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, percentiles are out of the domain (0, inf) => Lognormal distribution couldn't be fitted!", call. = FALSE)
  }
  if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
  }
  if (length(q) < 2) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
  }
  if (!is.logical(show.output)) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
  }
  if (!is.logical(plot)) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
  }
  if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
    # on.exit(return(invisible(NA)))
    stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
  }
  
  #-----------------------------------------------------------------------------
  # minimizing procedure
  #-----------------------------------------------------------------------------
  fit.weights.original <- fit.weights
  fit.weights <- fit.weights/sum(fit.weights)
  minimize <- function(theta) {
    summand <- suppressWarnings(stats::plnorm(q = q, meanlog = theta[1], sdlog = theta[2]) - p)
    summand <- summand * fit.weights
    sum(summand^2)
  }
  fit <- c(); fit$value <- tol + 1
  try1 <- try(
    fit <- stats::optim(par = c(1, 0.35), 
                        minimize, method = "L-BFGS-B", 
                        lower = c(-10000, 0.001), 
                        upper = c(10000, 10000)), 
    silent = TRUE
  )
  
  #-----------------------------------------------------------------------------
  # checking results
  #-----------------------------------------------------------------------------
  if (is.error(try1) || fit$value >= tol) {
    warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
    fit <- c(); fit$value <- tol + 1
    try2 <- try(
      fit <- stats::optim(par = c(1, 3),
                          minimize, 
                          method = "Nelder-Mead"), 
      silent = TRUE)
    if (is.error(try2) || fit$value >= tol) { 
      warning("The fitting procedure 'Nelder-Mead' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
      Par <- NA
    } else if (fit$value < tol) {
      message("The fitting procedure 'Nelder-Mead' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
      Par <- fit$par
      names(Par) <- c("meanlog", "sdlog")
      if (show.output) print(fit) 
    }
  } else if (fit$value < tol) {
    message("The fitting procedure 'L-BFGS-B' was successful!") 
    Par <- fit$par
    names(Par) <- c("meanlog", "sdlog")
    if (show.output) print(fit) 
  }
  
  
  
  #-----------------------------------------------------------------------------
  # output
  #-----------------------------------------------------------------------------
  return(Par)
}

# | 4. PSA Cohort Creation ####
psa_cohort <- cmpfun(function(n.i, nsp_status, psa_sharing, psa_reusing) {
  # ID generation
  id <- seq(1:n.i)
  id <- as.numeric(id)
  
  # Injection frequency in the past month
  inj.freq.dist <- c(0.129, 0.222, 0.156, 0.144, 0.349)
  names(inj.freq.dist) <- c("Never", "Not every week", "1-2 days a week",
                            "3-6 days a week", "Every day")
  v.inj.freq <- sample(x = 0:4, size = n.i, prob = inj.freq.dist, replace = TRUE)
  v.inj.freq <- as.numeric(v.inj.freq)
  table(v.inj.freq)
  
  # Number of injections in the past month
  inj.num.dist1 <- c(0.543, 0.457)
  inj.num.dist2 <- c(0.552, 0.448)
  inj.num.dist3 <- c(0.458, 0.542)
  inj.num.dist4 <- c(0.570, 0.430)
  
  v.inj.num <- vector("double", length = length(v.inj.freq))
  for (i in 1:length(v.inj.freq)) {
    v.inj.num[i] <- if (v.inj.freq[i] == 1) {
      sample(1:2, size = 1, prob = inj.num.dist1, replace = TRUE)
    } else if (v.inj.freq[i] == 2) {
      sample(3:4, size = 1, prob = inj.num.dist2, replace = TRUE)
    } else if (v.inj.freq[i] == 3) {
      sample(5:6, size = 1, prob = inj.num.dist3, replace = TRUE)
    } else if (v.inj.freq[i] == 4) {
      sample(7:8, size = 1, prob = inj.num.dist4, replace = TRUE)
    } else {
      0
    }
  }
  
  # Injection frequency in the past year
  inj.year.dist <- c(0.4012, 0.3555, 0.1816, 0.0616)
  names(inj.year.dist) <- c("< 2 years", "2-<5 years", "5-<8 years", "8+ years")
  v.yrs.idu <- sample(x = 1:4, size = n.i, prob = inj.year.dist, replace = TRUE)
  v.yrs.idu <- as.numeric(v.yrs.idu)
  
  # Risk factors for injection behavior modification
  rr.sharing.nsp <- psa_sharing # example value, replace with psa_sharing
  rr.reusing.nsp <- psa_reusing # example value, replace with psa_reusing
  
  set.seed(123)  # Set seed for reproducibility
  sharing_cap <- rnorm(1, mean = 0.5, sd = 0.01)
  reusing_cap <- rnorm(1, mean = 0.5, sd = 0.01)
  
  if (nsp_status == 1) {
    sharing.prob <- 0.205
    reusing.prob <- 0.297
  } else {
    sharing.prob <- min(0.205 * (1 / rr.sharing.nsp), sharing_cap)
    reusing.prob <- min(0.297 * (1 / rr.reusing.nsp), reusing_cap)
  }
  
  # Percent of PWID sharing needles/syringe with others
  v.sharing <- ifelse(v.inj.freq != 0, rbinom(n = length(v.inj.freq), size = 1, p = sharing.prob), 0)
  sum(v.sharing)
  
  # Percent of PWID reusing their own needles/syringes
  v.reusing <- ifelse(v.inj.freq != 0, rbinom(n = length(v.inj.freq), size = 1, p = reusing.prob), 0)
  sum(v.reusing)
  
  ## Sex (0 = Female, 1 = Male)
  sex.dist <- c(0.3795, 0.6205)
  v.sex <- sample(x = 0:1, size = n.i, prob = sex.dist, replace = T)
  v.sex <- as.numeric(v.sex)
  
  ## Age
  # Define the desired percentages
  percentages <- c(0.1065, 0.5108, 1)
  # Find the quantiles for the specified percentages
  quantiles <- qnorm(percentages, mean = 42.711, sd = 12.822)
  # Generate ages based on the quantiles
  v.age <- round(rnorm(n.i, 42.711, 12.822), 0)
  v.age <- ifelse(v.age < quantiles[1], runif(sum(v.age < quantiles[1]), 18, 24),
                  ifelse(v.age < quantiles[2], runif(sum(v.age < quantiles[2]), 25, 44),
                         runif(sum(v.age >= quantiles[2]), 45, 65)))
  
  ## Mortality rate
  v.mortality <- ifelse(v.age >= 18 & v.age <= 24, 0.0040, 
                        ifelse(v.age >= 25 & v.age <= 29, 0.0068,
                               ifelse(v.age >= 30 & v.age <= 34, 0.0089,
                                      ifelse(v.age >= 35 & v.age <= 39, 0.0130, 
                                             ifelse(v.age >= 40 & v.age <= 44, 0.0181,
                                                    ifelse(v.age >= 45 & v.age <= 49, 0.0336,
                                                           ifelse(v.age >= 50 & v.age <= 54, 0.0526,
                                                                  ifelse(v.age >= 55 & v.age <= 59, 0.0703, 0.0898)
                                                           )
                                                    )
                                             )
                                      )
                               )
                        )
  )
  v.ocm <- round((1 - exp(-1*v.mortality/52)), 6)
  
  ## NSP Status
  if (nsp_status == 1) {
    v.nsp <- rep(1, n.i)
  } else {
    v.nsp <- rep(0, n.i)
  }
  
  # Dataset
  ## Create a data frame for PWID baseline characteristics
  pwid.baseline <- data.frame(id, v.inj.freq, v.inj.num, v.yrs.idu, v.sharing, v.reusing, v.sex, v.age, v.ocm, v.nsp)
})

# | 5. PSA Risk Ratios ####
mutate_psa <- function(data, params) {
  data %>%
    mutate(
      rr_sharing = if_else(v.sharing == 1, params['rr.sharing.needles'], 1),
      rr_reusing = if_else(v.reusing == 1, params['rr.reusing.needles'], 1),
      rr_inj_freq = case_when(
        v.inj.freq %in% c(1, 2) ~ params['rr.inj.freq.2x'],
        v.inj.freq %in% c(3, 4) ~ params['rr.inj.freq.7x'],
        TRUE ~ 1
      ),
      rr_inj_num = case_when(
        v.inj.num %in% c(3, 4) ~ params['rr.num.attempts2'],
        v.inj.num %in% c(5, 6) ~ params['rr.num.attempts3'],
        v.inj.num %in% c(7, 8) ~ params['rr.num.attempts4'],
        TRUE ~ 1
      ),
      rr_inj_year = case_when(
        v.yrs.idu == 2 ~ params['rr.years.idu2to5'],
        v.yrs.idu == 3 ~ params['rr.years.idu5to8'],
        v.yrs.idu == 4 ~ params['rr.years.idu8plus'],
        TRUE ~ 1
      ),
      rr_age = case_when(
        v.age >= 25 & v.age <= 44 ~ params['rr.age25to44'],
        v.age > 44 ~ params['rr.age44plus'],
        TRUE ~ 1
      ),
      rr_sex = if_else(v.sex == 1, params['rr.male'], 1)
    )
}

# | 6. PSA Microsimulation ####
split_psa <- cmpfun(function(v.M_1, n.i, n.t, n.years, v.n, v.inf.prob, v.ocm, d.c, d.e, nsp, seed = 1) {
  
  v.dwc <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))    # calculate the cost discount weight based on the discount rate d.c 
  v.dwe <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))    # calculate the QALY discount weight based on the discount rate d.e
  
  # Create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- m.C <- m.E <- matrix(nrow = n.i, ncol = n.t, 
                              dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                              paste("cycle", 0:(n.t-1), sep = " ")))  
  
  m.M[, 1] <- v.M_1                  # indicate the initial health state   
  set.seed(seed)                  # set the seed for every individual for the random number generator
  # create the dur variable that stores the number of consecutive cycles the individual occupies either when sick or sicker
  m.C[, 1] <- Costs(m.M[, 1])  # estimate costs per individual for the initial health state
  m.E[, 1] <- Effs (m.M[, 1])  # estimate QALYs per individual for the initial health state
  
  for (t in 1:(n.t-1)) { # t <- 3
    m.p <- Probs(m.M[, t], v.inf.prob, v.ocm)           # calculate the transition probabilities at cycle t 
    m.M[, t + 1] <- samplev(prob = m.p, m = 1)  # sample the next health state and store that state in matrix m.M 
    m.C[, t + 1] <- Costs(m.M[, t + 1])   # estimate costs per individual during cycle t + 1 conditional on treatment
    m.E[, t + 1] <- Effs(m.M[, t + 1])   # estimate QALYs per individual during cycle t + 1 conditional on treatment
    cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))       # display the progress of the simulation
  } # close the loop for the time points 
  
  v.dwc <- as.matrix(as.data.frame(v.dwc))
  v.dwe <- as.matrix(as.data.frame(v.dwe))
  tc <- m.C %*% v.dwc     # total (discounted) cost per individual
  te <- m.E %*% v.dwe       # total (discounted) QALYs per individual 
  # tc_hat <- mean(tc) # Adding this part here will reduce the ICER by around 10-fold
  # te_hat <- mean(te)
  
  results <- list(tc = tc, te = te) # store the results from the simulation in a list  
  return(results)  # return the results
})

# | 7. Store PSA Results ####
psa_results <- function(num_sub_data, psa_split_list) {
  # num_sub_data: This is the number of sub-datasets per PSA iteration.
  # psa_split_list: This is the list that contains tc_hat and te_hat once the PSA
  # is completed using the psa_run() function.
  num_groups <- nrow(psa_split_list) %/% num_sub_data
  # Loop through the list by groups
  tc_hat <- vector(mode = "numeric", length = n_psa_iter)
  te_hat <- vector(mode = "numeric", length = n_psa_iter)
  for (i in 1:num_groups) {
    # Calculate the start and end indices for each group
    start_index <- (i - 1) * num_sub_data + 1
    end_index <- i * num_sub_data
    # Extract the rows for each PSA iteration
    tc_psa <- lapply(start_index:end_index, function(index) psa_split_list[[index]])
    te_psa <- lapply((start_index + nrow(psa_split_list)):(end_index + nrow(psa_split_list)), 
                     function(index) psa_split_list[[index]])
    tc_hat[i] <- mean(unlist(tc_psa))
    te_hat[i] <- mean(unlist(te_psa))
  }
  psa_output <- cbind(tc_hat, te_hat)
  return(psa_output)
}

########## ---- Part E: Model Calibration ##########
# | 1. Model Function ####
mod <- cmpfun(function(par_vector) {
  # | 1. Input Parameters ####
  id <- seq(1:n.i)
  id <- as.numeric(id)
  
  ## | a. Injection frequency in the past month ####
  inj.freq.dist <- c(0.129, 0.222, 0.156, 0.144, 0.349)
  names(inj.freq.dist) <- c("Never", "Not every week", "1-2 days a week",
                            "3-6 days a week", "Every day")
  v.inj.freq <- sample(x = 0:4, size = n.i, prob = inj.freq.dist, replace = T)
  v.inj.freq <- as.numeric(v.inj.freq)
  table(v.inj.freq)
  
  ## Number of injections in the past month
  inj.num.dist1 <- c(0.543, 0.457); inj.num.dist2 <- c(0.552, 0.448)
  inj.num.dist3 <- c(0.458, 0.542); inj.num.dist4 <- c(0.570, 0.430)
  v.inj.num <- vector("double", length = length(v.inj.freq))
  for (i in 1:length(v.inj.freq)){
    v.inj.num[i] <- ifelse(v.inj.freq[i] == 1, sample(1:2, size = 1, prob = inj.num.dist1, replace = T),
                           ifelse(v.inj.freq[i] == 2, sample(3:4, size = 1, prob = inj.num.dist2, replace = T),
                                  ifelse(v.inj.freq[i] == 3, sample(5:6, size = 1, prob = inj.num.dist3, replace = T),
                                         ifelse(v.inj.freq[i] == 4, sample(7:8, size = 1, prob = inj.num.dist4, replace = T),
                                                0)
                                  )
                           )
    )
  }
  
  ## | c. Injection frequency in the past month ####
  inj.year.dist <- c(0.4012, 0.3555, 0.1816, 0.0616)
  names(inj.year.dist) <- c("< 2 years", "2-<5 years", "5-<8 years", "8+ years")
  v.yrs.idu <- sample(x = 1:4, size = n.i, prob = inj.year.dist, replace = T)
  v.yrs.idu <- as.numeric(v.yrs.idu)
  
  ## | d. Risk factors for injection behaviour modification ####
  rr.sharing.nsp <- par_vector[1] # RR = 0.42 if NSP or 2.38 if no NSP
  rr.reusing.nsp <- par_vector[2] # RR = 0.79 if NSP or 1.27 if no NSP
  # Assume that the number of people who share needles/syringes increases by 2.38-fold
  # in a world without NSP. Similarly, 1.27 times the number of individuals who reuse
  # needles/syringes in a world without NSP. Use 95% CI for sensitivity analysis.
  
  ## | e. Percent of PWID sharing needles/syringe with others ####
  sharing.prob <- 0.205; reusing.prob <- 0.297
  v.sharing <- ifelse(v.inj.freq != 0, rbinom(n=sum(v.inj.freq >= 1), size=1, p=sharing.prob), 0)
  sum(v.sharing)
  
  ## | f. Percent of PWID reusing their own needles/syringes ####
  v.reusing <- ifelse(v.inj.freq != 0, rbinom(n=sum(v.inj.freq >= 1), size=1, p=reusing.prob), 0)
  sum(v.reusing)
  
  ## | g. Sex (0 = Female, 1 = Male) ####
  sex.dist <- c(0.3795, 0.6205)
  v.sex <- sample(x = 0:1, size = n.i, prob = sex.dist, replace = T)
  v.sex <- as.numeric(v.sex)
  
  ## | h. Age ####
  # Define the desired percentages
  percentages <- c(0.1065, 0.5108, 1)
  # Find the quantiles for the specified percentages
  quantiles <- qnorm(percentages, mean = 42.711, sd = 12.822)
  # Generate ages based on the quantiles
  v.age <- round(rnorm(n.i, 42.711, 12.822), 0)
  v.age <- ifelse(v.age < quantiles[1], runif(sum(v.age < quantiles[1]), 18, 24),
                  ifelse(v.age < quantiles[2], runif(sum(v.age < quantiles[2]), 25, 44),
                         runif(sum(v.age >= quantiles[2]), 45, 65)))
  
  ## Mortality rate
  v.mortality <- ifelse(v.age >= 18 & v.age <= 24, 0.0040, 
                        ifelse(v.age >= 25 & v.age <= 29, 0.0068,
                               ifelse(v.age >= 30 & v.age <= 34, 0.0089,
                                      ifelse(v.age >= 35 & v.age <= 39, 0.0130, 
                                             ifelse(v.age >= 40 & v.age <= 44, 0.0181,
                                                    ifelse(v.age >= 45 & v.age <= 49, 0.0336,
                                                           ifelse(v.age >= 50 & v.age <= 54, 0.0526,
                                                                  ifelse(v.age >= 55 & v.age <= 59, 0.0703, 0.0898)
                                                           )
                                                    )
                                             )
                                      )
                               )
                        )
  )
  v.ocm <- round((1 - exp(-1*v.mortality/52)), 6)
  
  ## | j. NSP Status ####
  v.nsp <- rep(1, n.i)
  
  # | 2. Dataset ####
  ## | a. Create a data frame for PWID baseline characteristics ####
  pwid.baseline <- data.frame(id, v.inj.freq, v.inj.num, v.yrs.idu, v.sharing, v.reusing, v.sex, v.age, v.ocm, v.nsp)
  
  ## | b. Add relative risk variables associated with SSTVI into the dataset
  pwid.baseline <- pwid.baseline |>
    dplyr::mutate(rr_sharing = if_else(v.sharing == 1, par_vector[3], 1),
                  rr_reusing = if_else(v.reusing == 1, par_vector[4], 1),
                  rr_inj_freq = case_when(
                    v.inj.freq %in% c(1, 2) ~ par_vector[5],
                    v.inj.freq %in% c(3, 4) ~ par_vector[6],
                    TRUE ~ 1
                  ),
                  rr_inj_num = case_when(
                    v.inj.num %in% c(3, 4) ~ par_vector[7],
                    v.inj.num %in% c(5, 6) ~ par_vector[8],
                    v.inj.num %in% c(7, 8) ~ par_vector[9],
                    TRUE ~ 1
                  ),
                  rr_inj_year = case_when(
                    v.yrs.idu == 2 ~ par_vector[10],
                    v.yrs.idu == 3 ~ par_vector[11],
                    v.yrs.idu == 4 ~ par_vector[12],
                    TRUE ~ 1
                  ),
                  rr_age = case_when(
                    v.age >= 25 & v.age <= 44 ~ par_vector[13],
                    v.age > 44 ~ par_vector[14],
                    TRUE ~ 1
                  ),
                  rr_sex = if_else(v.sex == 1, par_vector[15], 1)
    )
  
  pwid.baseline <- pwid.baseline |>
    mutate(risk.multiplier(data_name = pwid.baseline, prob = 0.001458, prev = 0.001452))
  
  v.M_1 <- pwid.baseline$state         # Everyone starts as "healthy" or "SSTVI" 
  v.inf.prob <- pwid.baseline$prob_infection # Probability of SSTVI for each individual
  v.ocm <- pwid.baseline$v.ocm # Probability of other-cause mortality for each individual
  
  # | 3. Simulation ####
  # Create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
  m.M <- matrix(nrow = n.i, ncol = n.t, 
                dimnames = list(paste("ind", 1:n.i, sep = " "), 
                                paste("cycle", 0:(n.t-1), sep = " ")))  
  
  m.M[, 1] <- v.M_1                  # indicate the initial health state   
  
  for (t in 1:(n.t-1)) { # t <- 3
    m.p <- Probs(m.M[, t], v.inf.prob, v.ocm)           # calculate the transition probabilities at cycle t 
    m.M[, t + 1] <- samplev(prob = m.p, m = 1)  # sample the next health state and store that state in matrix m.M 
    cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))       # display the progress of the simulation
  } # close the loop for the time points 
  
  # Create a matrix of transitions across states
  TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other
  TS <- matrix(TS, nrow = n.i)
  rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
  colnames(TS) <- paste("Cycle", 0:(n.t-1), sep = " ")   # name the columns 
  
  # | 4. Mortality - Data Setup ####
  transition.path <- TS
  transition.path <- as.data.frame(transition.path)
  transition.ocm <- apply(transition.path, 1, function(x) match(TRUE, stringr::str_sub(x, -3, -1) == "OCM"))
  transition.sstvi.death <- apply(transition.path, 1, function(x) match(TRUE, stringr::str_sub(x, -5, -1) == "Death"))
  
  # Vectorize the death times (other-cause mortality, SSI, and complications)
  transition.ocm <- as.vector(transition.ocm, mode = "numeric")
  transition.sstvi.death <- as.vector(transition.sstvi.death, mode = "numeric")
  
  # Create data frame
  death.dates <- data.frame(OCM = transition.ocm, SSTVI_Death = transition.sstvi.death)
  death.dates$Death <- apply(death.dates[,c(1:2)], 1, min, na.rm = TRUE)
  death.dates$Death <- ifelse(death.dates$Death == Inf, n.t, death.dates$Death)
  death.dates$id <- seq(1:n.i)
  
  # Person-time for each year
  death.dates$yr1 <- pmin(death.dates$Death / 52, 1)
  death.dates$yr2 <- pmin((death.dates$Death - 52) / 52, 1)
  death.dates$yr2[death.dates$Death <= 52] <- 0
  death.dates$yr3 <- pmin((death.dates$Death - 104) / 52, 1)
  death.dates$yr3[death.dates$Death <= 104] <- 0
  death.dates$yr4 <- pmin((death.dates$Death - 156) / 52, 1)
  death.dates$yr4[death.dates$Death <= 156] <- 0
  death.dates$yr5 <- pmin((death.dates$Death - 208) / 52, 1)
  death.dates$yr5[death.dates$Death <= 208] <- 0
  death.dates$yr5[death.dates$Death >= 260] <- 1
  
  # Population (denominator)
  pop_yr1 <- sum(death_dates_with_nsp$yr1 != 0)
  pop_yr2 <- sum(death_dates_with_nsp$yr2 != 0)
  pop_yr3 <- sum(death_dates_with_nsp$yr3 != 0)
  pop_yr4 <- sum(death_dates_with_nsp$yr4 != 0)
  pop_yr5 <- sum(death_dates_with_nsp$yr5 != 0)
  pop_yrs <- c(pop_yr1, pop_yr2, pop_yr3, pop_yr4, pop_yr5)
  
  ## Counts per year: Other-cause mortality
  ocm.1yr <- length(which(death.dates$OCM > 0 & death.dates$OCM <= 52))
  ocm.2yr <- length(which(death.dates$OCM > 52 & death.dates$OCM <= 104))
  ocm.3yr <- length(which(death.dates$OCM > 104 & death.dates$OCM <= 156))
  ocm.4yr <- length(which(death.dates$OCM > 156 & death.dates$OCM <= 208))
  ocm.5yr <- length(which(death.dates$OCM > 208 & death.dates$OCM <= 260))
  ocm.yrs <- c(ocm.1yr, ocm.2yr, ocm.3yr, ocm.4yr, ocm.5yr)
  ocm.ir <- (sum(ocm.yrs) / sum(pop_yrs)) * 1000
  
  # Counts per year: SSTVI mortality
  sstvi.1yr <- length(which(death.dates$SSTVI_Death > 0 & death.dates$SSTVI_Death <= 52))
  sstvi.2yr <- length(which(death.dates$SSTVI_Death > 52 & death.dates$SSTVI_Death <= 104))
  sstvi.3yr <- length(which(death.dates$SSTVI_Death > 104 & death.dates$SSTVI_Death <= 156))
  sstvi.4yr <- length(which(death.dates$SSTVI_Death > 156 & death.dates$SSTVI_Death <= 208))
  sstvi.5yr <- length(which(death.dates$SSTVI_Death > 208 & death.dates$SSTVI_Death <= 260))
  sstvi.yrs <- c(sstvi.1yr, sstvi.2yr, sstvi.3yr, sstvi.4yr, sstvi.5yr)
  sstvi.ir <- (sum(sstvi.yrs) / sum(pop_yrs)) * 1000
  
  results <- list(OCM_IR = ocm.ir, SSTVI_IR = sstvi.ir) # store the results from the simulation in a list  
  return(results)
})

# Horizontal error bars
plot_error_bars <- function(x, y, lower, upper, col, pch, legend_text) {
  arrows(lower, y, upper, y, angle = 90, code = 3, length = 0.02, col = col)
  points(x, y, pch = pch, col = col)
}
