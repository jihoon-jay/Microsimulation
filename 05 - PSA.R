# Jihoon Lim
# Economic Modelling - Microsimulation
# 05 - Probabilistic Sensitivity Analysis

# Note: Une partie de la compilation est effectuée à partir de données provenant 
# du © Gouvernement du Québec (année de la publication du Fichier de recherche: 
# 2009-2019). Le © Gouvernement du Québec n’est pas responsable des compilations
# ni de l’interprétation des résultats produits à l’aide du Fichier de recherche.

########## ---- Part A: PSA Parameter Setup ##########
library(readxl)
library(flextable)
library(ggplot2)
library(dampack)
library(tibble)
library(dplyr)
#library(rriskDistributions)
library(doParallel)
library(parallel)
library(foreach)
theme_set(theme_bw()) #Makes ggplots look better

# | 1. Read in Cost and Utility Parameters from Excel ####
t_rewards <- read_excel('~/psa_params.xlsx', sheet = "psa")
#Display it nicely
t_rewards |>
  flextable() |> #turn into flextable object
  merge_v(j=1) |> #Merge cells in first column with same value (group probabilities, costs, etc.)
  theme_box() |> #Apply a theme for aesthetics
  autofit() #automatically set column widths to reasonable values

# | 2. Read in Risk Ratio Parameters from Excel ####
t_risk <- read_excel('~/psa_params.xlsx', sheet = "risk_multiplier")
#Display it nicely
t_risk |>
  flextable() |> #turn into flextable object
  merge_v(j=1) |> #Merge cells in first column with same value (group probabilities, costs, etc.)
  theme_box() |> #Apply a theme for aesthetics
  autofit() #automatically set column widths to reasonable values

########## ---- Part B: Sampling from Parameters ##########
n_psa_iter <- 1000 # Number of PSA iterations
# Get names of parameters
reward_param_names <- t_rewards$rname
risk_param_names <- t_risk$rname

# | 1. Sample risk ratio multipliers from log normal distribution ####
# Pre-allocate matrix
m_risk_params_psa <- matrix(NA, nrow=n_psa_iter, 
                            ncol=length(risk_param_names),
                            dimnames=list(seq(1:n_psa_iter), risk_param_names))
for (param in risk_param_names){
  risk_mult <- get.lnorm.par(p=c(0.025, 0.5, 0.975),
                             q=c(t_risk$value_min[t_risk$rname == param],
                                 t_risk$value[t_risk$rname == param],
                                 t_risk$value_max[t_risk$rname == param]))
  m_risk_params_psa[, param] <- rlnorm(n_psa_iter, meanlog = risk_mult[1], sdlog = risk_mult[2])
}

# | 2. Sample costs (gamma) and QALYs (beta) ####
# Pre-allocate matrix
m_reward_params_psa <- matrix(NA, nrow=n_psa_iter, ncol=length(reward_param_names),
                              dimnames=list(seq(1:n_psa_iter), reward_param_names))
# Sample cost and QALY parameters
for (param in reward_param_names){
  if (isTRUE(t_rewards$distribution[t_rewards$rname == param] == "beta")) {
    m_reward_params_psa[, param] <- rbeta_mean_quantile(t_rewards$value[t_rewards$rname == param], 
                                                        t_rewards$value_min[t_rewards$rname == param], 
                                                        t_rewards$value_max[t_rewards$rname == param], 
                                                        n=n_psa_iter)
  } else if(isTRUE(t_rewards$distribution[t_rewards$rname == param] == "gamma")) {
    m_reward_params_psa[, param] <- rgamma(n_psa_iter, 
                                           shape = t_rewards$shape[t_rewards$rname == param], 
                                           scale = t_rewards$scale[t_rewards$rname == param])
  } else if (isTRUE(is.na(t_rewards$distribution[t_rewards$rname == param]))) {
    m_reward_params_psa[, param] <- rep(t_rewards$value[t_rewards$rname == param], n_psa_iter)
  }
}

# | 3. Create PSA datasets ####
# || a. With NSP ####
set.seed(1234)
psa.cohort.with.nsp <- psa_cohort(n.i = 10000, nsp_status = 1)
# psa.cohort.with.nsp <- pwid.baseline.with.nsp[1:10000,] # Create a cohort with 10,000 patients from the original cohort
psa_list1 <- list()
for (i in 1:n_psa_iter) {
  psa_list1[[i]] <- psa.cohort.with.nsp |>
    dplyr::mutate(rr_sharing = if_else(v.sharing == 1, m_risk_params_psa[i, 'rr.sharing.needles'], 1),
                  rr_reusing = if_else(v.reusing == 1, m_risk_params_psa[i, 'rr.reusing.needles'], 1),
                  rr_inj_freq = case_when(
                    v.inj.freq %in% c(1, 2) ~ m_risk_params_psa[i, 'rr.inj.freq.2x'],
                    v.inj.freq %in% c(3, 4) ~ m_risk_params_psa[i, 'rr.inj.freq.7x'],
                    TRUE ~ 1
                  ),
                  rr_inj_num = case_when(
                    v.inj.num %in% c(3, 4) ~ m_risk_params_psa[i, 'rr.num.attempts2'],
                    v.inj.num %in% c(5, 6) ~ m_risk_params_psa[i, 'rr.num.attempts3'],
                    v.inj.num %in% c(7, 8) ~ m_risk_params_psa[i, 'rr.num.attempts4'],
                    TRUE ~ 1
                  ),
                  rr_inj_year = case_when(
                    v.yrs.idu == 2 ~ m_risk_params_psa[i, 'rr.years.idu5to7'],
                    v.yrs.idu == 3 ~ m_risk_params_psa[i, 'rr.years.idu8to10'],
                    TRUE ~ 1
                  ),
                  rr_age = case_when(
                    v.age >= 25 & v.age <= 44 ~ m_risk_params_psa[i, 'rr.age25to44'],
                    v.age > 44 ~ m_risk_params_psa[i, 'rr.age44plus'],
                    TRUE ~ 1
                  ),
                  rr_sex = if_else(v.sex == 1, m_risk_params_psa[i, 'rr.male'], 1)
    )
}

for (i in 1:n_psa_iter) {
  psa_list1[[i]] <- psa_list1[[i]] |>
    mutate(risk.multiplier(data_name = psa_list1[[i]], prob = 0.001458, prev = 0.0014523))
}

# || b. No NSP ####
set.seed(1234)
psa.cohort.no.nsp <- psa_cohort(n.i = 10000, nsp_status = 0)
# psa.cohort.no.nsp <- pwid.baseline.no.nsp[1:10000,] # Create a cohort with 10,000 patients from the original cohort
psa_list0 <- list()
for (i in 1:n_psa_iter) {
  psa_list0[[i]] <- psa.cohort.no.nsp |>
    dplyr::mutate(rr_sharing = if_else(v.sharing == 1, m_risk_params_psa[i, 'rr.sharing.needles'], 1),
                  rr_reusing = if_else(v.reusing == 1, m_risk_params_psa[i, 'rr.reusing.needles'], 1),
                  rr_inj_freq = case_when(
                    v.inj.freq %in% c(1, 2) ~ m_risk_params_psa[i, 'rr.inj.freq.2x'],
                    v.inj.freq %in% c(3, 4) ~ m_risk_params_psa[i, 'rr.inj.freq.7x'],
                    TRUE ~ 1
                  ),
                  rr_inj_num = case_when(
                    v.inj.num %in% c(3, 4) ~ m_risk_params_psa[i, 'rr.num.attempts2'],
                    v.inj.num %in% c(5, 6) ~ m_risk_params_psa[i, 'rr.num.attempts3'],
                    v.inj.num %in% c(7, 8) ~ m_risk_params_psa[i, 'rr.num.attempts4'],
                    TRUE ~ 1
                  ),
                  rr_inj_year = case_when(
                    v.yrs.idu == 2 ~ m_risk_params_psa[i, 'rr.years.idu5to7'],
                    v.yrs.idu == 3 ~ m_risk_params_psa[i, 'rr.years.idu8to10'],
                    TRUE ~ 1
                  ),
                  rr_age = case_when(
                    v.age >= 25 & v.age <= 44 ~ m_risk_params_psa[i, 'rr.age25to44'],
                    v.age > 44 ~ m_risk_params_psa[i, 'rr.age44plus'],
                    TRUE ~ 1
                  ),
                  rr_sex = if_else(v.sex == 1, m_risk_params_psa[i, 'rr.male'], 1)
    )
}

for (i in 1:n_psa_iter) {
  psa_list0[[i]] <- psa_list0[[i]] |>
    mutate(risk.multiplier(data_name = psa_list0[[i]], prob = 0.001458, prev = 0.0014523))
}

########## ---- Part C: Microsimulation Run: With NSP ##########
# With 10,000 patients, it should take around 5 hours to run 1,000 iterations.
# | 1. Initialize number of sub-datasets and number of patients in sub-datasets ####
# Fix n.i at 100!!!
n.i <- 100 # Number of pat in each split dataset (has to match with # pat in each split data)
nr <- nrow(psa.cohort.with.nsp) # Number of rows in each dataset (same as original PSA cohort size)
n_sub <- nr/n.i # Number of sub-datasets you want in each PSA run

# | 2. Parallel computing ####
psa1_split_list <- list() # Initialize list
cl <- parallel::makeCluster(detectCores() - 16) # Detect the number of available cores and create cluster
doParallel::registerDoParallel(cl) # Activate cluster for foreach library
# Run model
set.seed(1234)
time_foreach <- system.time({
  psa1_split_list <- foreach(i = 1:n_psa_iter, .combine = rbind) %:% 
    foreach::foreach(j = 1:n_sub, .combine = rbind) %dopar% {
      # Costs
      c_Healthy <- m_reward_params_psa[i, "c_Healthy"] + c_NSP
      c_SSTVI <- m_reward_params_psa[i, "c_SSTVI"]  + c_NSP # SSTVI
      c_SelfTrt_Death <- m_reward_params_psa[i, "c_SelfTrt_Death"] # Death from IE
      c_OP_Purulent <- m_reward_params_psa[i, "c_OP_Purulent"] + c_NSP # Outpatient (purulent)
      c_ED_Purulent <- m_reward_params_psa[i, "c_ED_Purulent"] + c_NSP # ED (purulent)
      c_IP_Purulent <- m_reward_params_psa[i, "c_IP_Purulent"] + c_NSP # Inpatient (purulent)
      c_IPC_Purulent <- m_reward_params_psa[i, "c_IPC_Purulent"] + c_NSP # Inpatient (purulent) involving surgeries
      c_OP_NonPurulent <- m_reward_params_psa[i, "c_OP_NonPurulent"] + c_NSP # Outpatient (non-purulent)
      c_ED_NonPurulent <- m_reward_params_psa[i, "c_ED_NonPurulent"] + c_NSP # ED (non-purulent)
      c_IP_NonPurulent <- m_reward_params_psa[i, "c_IP_NonPurulent"] + c_NSP # Inpatient (non-purulent)
      c_IPC_NonPurulent <- m_reward_params_psa[i, "c_IPC_NonPurulent"] + c_NSP # Inpatient (non-purulent) involving surgeries
      c_SSTVI_Death <- m_reward_params_psa[i, "c_SSTVI_Death"] # Death from SSTVI
      c_OCM <- m_reward_params_psa[i, "c_OCM"]
      
      # QALYs
      q_Healthy <- m_reward_params_psa[i, "q_Healthy"] / 52
      q_SSTVI <- m_reward_params_psa[i, "q_SSTVI"] / 52 # SSTVI
      q_SelfTrt_Death <- m_reward_params_psa[i, "q_SelfTrt_Death"] # Death from IE
      q_OP_Purulent <- m_reward_params_psa[i, "q_OP_Purulent"] / 52 # Outpatient (purulent)
      q_ED_Purulent <- m_reward_params_psa[i, "q_ED_Purulent"] / 52 # ED (purulent)
      q_IP_Purulent <- m_reward_params_psa[i, "q_IP_Purulent"] / 52 # Inpatient (purulent)
      q_IPC_Purulent <- m_reward_params_psa[i, "q_IPC_Purulent"] / 52 # Inpatient (purulent) involving surgeries
      q_OP_NonPurulent <- m_reward_params_psa[i, "q_OP_NonPurulent"] / 52 # Outpatient (non-purulent)
      q_ED_NonPurulent <- m_reward_params_psa[i, "q_ED_NonPurulent"] / 52 # ED (non-purulent)
      q_IP_NonPurulent <- m_reward_params_psa[i, "q_IP_NonPurulent"] / 52 # Inpatient (non-purulent)
      q_IPC_NonPurulent <- m_reward_params_psa[i, "q_IPC_NonPurulent"] / 52 # Inpatient (non-purulent) involving surgeries
      q_SSTVI_Death <- m_reward_params_psa[i, "q_SSTVI_Death"] # Death from SSI
      q_OCM <- m_reward_params_psa[i, "q_OCM"]
      
      # each = : represents how many people should be in each sub-dataset
      # n.i x n_sub should equal the number of patients in each PSA dataset
      # Ex) n.i = 100, n_sub = 10, then the total number of patients in each PSA
      # run before splitting is 1000.
      tmp_list <- list()
      tmp_list <- split(psa_list1[[i]], rep(1:n_sub, each=n.i, length.out=nr))
      
      # Create input variables to run split_psa() below.
      v.M_1 <- tmp_list[[j]]$state
      v.inf.prob <- tmp_list[[j]]$prob_infection
      v.ocm <- tmp_list[[j]]$v.ocm
      
      # Run microsimulation
      split_psa(v.M_1, n.i = n.i, n.t = n.t, n.years = n.years, 
                v.n, v.inf.prob = v.inf.prob, v.ocm = v.ocm, 
                d.c, d.e, nsp = 1, seed = (i*j))
    }
})
time_foreach[3]
parallel::stopCluster(cl) # Stop cluster to free up resources

########## ---- Part D: Microsimulation Run: No NSP ##########
# | 1. Initialize number of sub-datasets and number of patients in sub-datasets ####
# Fix n.i at 100!!!
n.i <- 100 # Number of pat in each split dataset (has to match with # pat in each split data)
nr <- nrow(psa.cohort.no.nsp) # Number of rows in each dataset (same as original PSA cohort size)
n_sub <- nr/n.i # Number of sub-datasets you want in each PSA run

# | 2. Parallel computing ####
psa0_split_list <- list()
cl <- parallel::makeCluster(detectCores() - 16) # Detect the number of available cores and create cluster
doParallel::registerDoParallel(cl) # Activate cluster for foreach library
# Run model
set.seed(1234)
time_foreach <- system.time({
  psa0_split_list <- foreach(i = 1:n_psa_iter, .combine = rbind) %:% 
    foreach::foreach(j = 1:n_sub, .combine = rbind) %dopar% {
      # Costs
      c_Healthy <- m_reward_params_psa[i, "c_Healthy"]
      c_SSTVI <- m_reward_params_psa[i, "c_SSTVI"]
      c_SelfTrt_Death <- m_reward_params_psa[i, "c_SelfTrt_Death"]
      c_OP_Purulent <- m_reward_params_psa[i, "c_OP_Purulent"]
      c_ED_Purulent <- m_reward_params_psa[i, "c_ED_Purulent"]
      c_IP_Purulent <- m_reward_params_psa[i, "c_IP_Purulent"]
      c_IPC_Purulent <- m_reward_params_psa[i, "c_IPC_Purulent"]
      c_OP_NonPurulent <- m_reward_params_psa[i, "c_OP_NonPurulent"]
      c_ED_NonPurulent <- m_reward_params_psa[i, "c_ED_NonPurulent"]
      c_IP_NonPurulent <- m_reward_params_psa[i, "c_IP_NonPurulent"]
      c_IPC_NonPurulent <- m_reward_params_psa[i, "c_IPC_NonPurulent"]
      c_SSTVI_Death <- m_reward_params_psa[i, "c_SSTVI_Death"]
      c_OCM <- m_reward_params_psa[i, "c_OCM"]
      
      # QALYs
      q_Healthy <- m_reward_params_psa[i, "q_Healthy"] / 52
      q_SSTVI <- m_reward_params_psa[i, "q_SSTVI"] / 52 # SSTVI
      q_SelfTrt_Death <- m_reward_params_psa[i, "q_SelfTrt_Death"] # Death from IE
      q_OP_Purulent <- m_reward_params_psa[i, "q_OP_Purulent"] / 52 # Outpatient (purulent)
      q_ED_Purulent <- m_reward_params_psa[i, "q_ED_Purulent"] / 52 # ED (purulent)
      q_IP_Purulent <- m_reward_params_psa[i, "q_IP_Purulent"] / 52 # Inpatient (purulent)
      q_IPC_Purulent <- m_reward_params_psa[i, "q_IPC_Purulent"] / 52 # Inpatient (purulent) involving surgeries
      q_OP_NonPurulent <- m_reward_params_psa[i, "q_OP_NonPurulent"] / 52 # Outpatient (non-purulent)
      q_ED_NonPurulent <- m_reward_params_psa[i, "q_ED_NonPurulent"] / 52 # ED (non-purulent)
      q_IP_NonPurulent <- m_reward_params_psa[i, "q_IP_NonPurulent"] / 52 # Inpatient (non-purulent)
      q_IPC_NonPurulent <- m_reward_params_psa[i, "q_IPC_NonPurulent"] / 52 # Inpatient (non-purulent) involving surgeries
      q_SSTVI_Death <- m_reward_params_psa[i, "q_SSTVI_Death"] # Death from SSI
      q_OCM <- m_reward_params_psa[i, "q_OCM"]
      
      # each = : represents how many people should be in each sub-dataset
      # n.i x n_sub should equal the number of patients in each PSA dataset
      # Ex) n.i = 100, n_sub = 10, then the total number of patients in each PSA
      # run before splitting is 1000.
      tmp_list <- list()
      tmp_list <- split(psa_list0[[i]], rep(1:n_sub, each=n.i, length.out=nr))
      
      # Create input variables to run split_psa() below.
      v.M_1 <- tmp_list[[j]]$state
      v.inf.prob <- tmp_list[[j]]$prob_infection
      v.ocm <- tmp_list[[j]]$v.ocm
      
      # Run microsimulation
      split_psa(v.M_1, n.i = n.i, n.t = n.t, n.years = n.years, 
                v.n, v.inf.prob = v.inf.prob, v.ocm = v.ocm, 
                d.c, d.e, nsp = 0, seed = (i*j))
    }
})
time_foreach[3]
parallel::stopCluster(cl) # Stop cluster to free up resources

########## ---- Part E: Cost-Effectiveness Plane ##########
# | 1. Generate PSA outputs ####
psa1_output <- psa_results(num_sub_data = n_sub, psa_split_list = psa1_split_list)
psa0_output <- psa_results(num_sub_data = n_sub, psa_split_list = psa0_split_list)

# | 2. Incremental Costs and QALYs ####
inc_c <- inc_e <- icer <- vector(mode = "numeric", length = n_psa_iter)
for (i in 1:n_psa_iter) {
  inc_c[i] <- psa1_output[i, 1] - psa0_output[i, 1]
  inc_e[i] <- psa1_output[i ,2] - psa0_output[i ,2]
  icer[i] <- inc_c[i] / inc_e[i]
}

# | 3. Cost-Effectiveness Plane ####
icer_psa <- data.frame(inc_c, inc_e, icer)
ggplot(icer_psa, aes(x=inc_e, y=inc_c)) + 
  geom_point() + labs(title = "Cost-Effectiveness Plane", 
                      x = "Incremental Effects", 
                      y = "Incremental Costs") + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0)
sum(icer_psa$inc_c > 0 & icer_psa$inc_e > 0) # NE quandrant
sum(icer_psa$inc_c < 0 & icer_psa$inc_e > 0) # SE quandrant
sum(icer_psa$inc_c > 0 & icer_psa$inc_e < 0) # NW quandrant
sum(icer_psa$inc_c < 0 & icer_psa$inc_e < 0) # SW quandrant

########## ---- Part F: Cost-Effectiveness Acceptability Curve ##########
# | 1. Create Threshold Vector ####
wtp_threshold <- seq(0, 300000, 5000) # Modify values here

# | 2. Create data frame ####
id <- seq(1, n_psa_iter)
df <- matrix(NA, nrow = n_psa_iter, ncol = length(wtp_threshold))

# | 3. Calculate incremental net benefits ####
for(i in 1:length(wtp_threshold)){
  df[,i] <- wtp_threshold[i]*inc_e-inc_c  
}

# | 4. Probability of cost-effectiveness ####
ProbCE <- c()
for (i in 1:length(wtp_threshold)) {
  ProbCE[i] <- sum(df[,i] >= 0) / n_psa_iter
}

# | 5. CEAC data frame ####
CEAC_Data <- data.frame(wtp_threshold, ProbCE); CEAC_Data

# | 6. Plot CEAC ####
ggplot(CEAC_Data, aes(y = ProbCE)) +
  geom_line(aes(x = wtp_threshold), color = "blue", line=1) +
  geom_point(aes(x = wtp_threshold), color = "black", size=2) +
  ylim(0, 1) +
  xlab("WTP per QALY") +
  ylab("Probability of Cost Effectiveness") +
  ggtitle("Cost Effectiveness Acceptability Curve")
