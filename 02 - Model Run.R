# Jihoon Lim
# Economic Modelling - Microsimulation
# 02 - Model Run

# Note: Une partie de la compilation est effectuée à partir de données provenant 
# du © Gouvernement du Québec (année de la publication du Fichier de recherche: 
# 2009-2019). Le © Gouvernement du Québec n’est pas responsable des compilations
# ni de l’interprétation des résultats produits à l’aide du Fichier de recherche.

########## ---- Part A: Data Generation ##########
source("~/00 - Functions.R")
source("~/01 - Setup.R")
p1 = Sys.time()
risk_multiplier_prior <- c(3.31, 2.1, 2.1, 3.1, 2.6, 3.7, 3.8, 3.95, 4.84, 1.19, 1.07, 1.14)
set.seed(10000)
pwid.baseline.with.nsp <- data_generation(n.i = 100000, nsp_status = 1, par_vector = risk_multiplier_posterior)
set.seed(10000)
pwid.baseline.no.nsp <- data_generation(n.i = 100000, nsp_status = 0, par_vector = risk_multiplier_posterior)
comp.time1 <- Sys.time() - p1; comp.time1 # Approximately 15 seconds

# Before running this part, make sure that the transition probability, cost, and QALY
# parameters are run and that they show up in the global environment!!!
########## ---- Part B: Microsimulation: With NSP ##########
# | 0. Update cost variables ####
# Note: The '01 - Setup' contains cost variables with NSP. In a scenario without NSP,
# costs associated with each health state has to be updated!
c_Healthy <- 0 + c_NSP; c_SSTVI <- 0 + c_NSP; c_SelfTrt_Death <- 0
c_OP_Purulent <- 98.58 + c_NSP; c_ED_Purulent <- 1134.59 + c_NSP
c_IP_Purulent <- 6045.12 + c_NSP; c_IPC_Purulent <- 42526.20 + c_NSP
c_OP_NonPurulent <- 90.67 + c_NSP; c_ED_NonPurulent <- 1143.70 + c_NSP
c_IP_NonPurulent <- 6818.38 + c_NSP; c_IPC_NonPurulent <- 41016.17 + c_NSP
c_SSTVI_Death <- 0; c_OCM <- 0

# | 1. Split data into equal parts ####
n_split <- 100 # Number of patients in each split dataset
nr <- nrow(pwid.baseline.with.nsp)
n_sub <- nr/n_split # Number of sub-datasets you want
# Each dataset should have NR/N_split rows
# Create a sequence from 1 to n_sub with each number repeating 100 times (100 rows each)
sim_trt_split <- split(pwid.baseline.with.nsp, rep(1:ceiling(n_sub), each=n_split, length.out=nr))

# | 2. Parallel computing ####
# Initialize list
sim_trt_split_list <- list()
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores() - 16)
# Activate cluster for foreach library
doParallel::registerDoParallel(cl)
# Run model
n.i <- 100 # This number has to match with the number of individuals in each of the split datasets
#set.seed(100)
time_foreach <- system.time({
  sim_trt_split_list <- foreach::foreach(i = 1:n_sub, .combine = rbind) %dopar% {
    v.M_1.with.nsp <- sim_trt_split[[i]]$state
    v.inf.prob.with.nsp <- sim_trt_split[[i]]$prob_infection
    v.ocm.with.nsp <- sim_trt_split[[i]]$v.ocm
    
    MicroSim(v.M_1.with.nsp, n.i, n.t = n.t, n.years = n.years, 
             v.n, v.inf.prob = v.inf.prob.with.nsp, v.ocm = v.ocm.with.nsp, 
             d.c, d.e, nsp = 1, seed = i)
  }
})
time_foreach[3] # 73 seconds to run (about 1 minute 13 seconds) with 24 cores
# Stop cluster to free up resources
parallel::stopCluster(cl)

# | 3. Use lapply to extract the specified matrices ####
selected_matrices_m.p <- lapply(1:1000, function(index) sim_trt_split_list[[index]])
selected_matrices_m.M <- lapply(1001:2000, function(index) sim_trt_split_list[[index]])
selected_matrices_m.C <- lapply(2001:3000, function(index) sim_trt_split_list[[index]])
selected_matrices_m.E <- lapply(3001:4000, function(index) sim_trt_split_list[[index]])
selected_matrices_tc <- lapply(4001:5000, function(index) sim_trt_split_list[[index]])
selected_matrices_te <- lapply(5001:6000, function(index) sim_trt_split_list[[index]])
selected_matrices_tc_hat <- lapply(6001:7000, function(index) sim_trt_split_list[[index]])
selected_matrices_te_hat <- lapply(7001:8000, function(index) sim_trt_split_list[[index]])
selected_matrices_TS <- lapply(8001:9000, function(index) sim_trt_split_list[[index]])

# | 4. Use do.call and rbind to append the selected matrices into a single matrix ####
m.p <- do.call(rbind, selected_matrices_m.p)
m.M <- do.call(rbind, selected_matrices_m.M)
m.C <- do.call(rbind, selected_matrices_m.C)
m.E <- do.call(rbind, selected_matrices_m.E)
tc <- do.call(rbind, selected_matrices_tc)
te <- do.call(rbind, selected_matrices_te)
tc_hat <- mean(tc)
te_hat <- mean(te)
TS <- do.call(rbind, selected_matrices_TS)

# | 5. Finalize list ####
appended_matrix_list <- list(m.p, m.M, m.C, m.E, tc, te, tc_hat, te_hat, TS)
sim_trt <- lapply(appended_matrix_list, as.matrix)
sim_trt <- setNames(sim_trt, c("m.p", "m.M", "m.C", "m.E", "tc", "te", 
                               "tc_hat", "te_hat", "TS"))
n.i <- 100000
rownames(sim_trt$TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
colnames(sim_trt$TS) <- paste("Cycle", 0:(n.t-1), sep = " ")   # name the columns 
rm(list = ls()[grep("^selected_matrices_", ls())])
rm(m.p, m.M, m.C, m.E, tc, te, tc_hat, te_hat, TS)
rm(sim_trt_split_list, appended_matrix_list)

# Before, running this part, re-run the code to update the cost parameters. Because
# NSP costs are added in the With NSP cohort, we have to re-define the cost variables
# in the global environment (this time, without the NSP costs).
########## ---- Part C: Microsimulation: No NSP ##########
# | 0. Update cost variables ####
# Note: The '01 - Setup' contains cost variables with NSP. In a scenario without NSP,
# costs associated with each health state has to be updated!
c_Healthy <- 0; c_SSTVI <- 0; c_SelfTrt_Death <- 0
c_OP_Purulent <- 98.58; c_ED_Purulent <- 1134.59
c_IP_Purulent <- 6045.12; c_IPC_Purulent <- 42526.20
c_OP_NonPurulent <- 90.67; c_ED_NonPurulent <- 1143.70
c_IP_NonPurulent <- 6818.38; c_IPC_NonPurulent <- 41016.17
c_SSTVI_Death <- 0; c_OCM <- 0

# | 1. Split data into equal parts ####
n_split <- 100 # Split into N smaller datasets
nr <- nrow(pwid.baseline.no.nsp) 
# Each dataset should have NR/N_split rows
# Create a sequence from 1 to 100 with each number repeating 100 times (100 rows each)
sim_no_trt_split <- split(pwid.baseline.no.nsp, rep(1:ceiling(n_sub), each=n_split, length.out=nr))

# | 2. Parallel computing ####
# Initialize list
sim_no_trt_split_list <- list()
# Detect the number of available cores and create cluster
cl <- parallel::makeCluster(detectCores() - 16)
# Activate cluster for foreach library
doParallel::registerDoParallel(cl)
# Run model
n.i <- 100 # This number has to match with the number of individuals in each of the split datasets
#set.seed(100)
time_foreach <- system.time({
  sim_no_trt_split_list <- foreach::foreach(i = 1:n_sub, .combine = rbind) %dopar% {
    v.M_1.no.nsp <- sim_no_trt_split[[i]]$state
    v.inf.prob.no.nsp <- sim_no_trt_split[[i]]$prob_infection
    v.ocm.no.nsp <- sim_no_trt_split[[i]]$v.ocm
    
    MicroSim(v.M_1.no.nsp, n.i, n.t = n.t, n.years = n.years, 
             v.n, v.inf.prob = v.inf.prob.no.nsp, v.ocm = v.ocm.no.nsp, 
             d.c, d.e, nsp = 0, seed = i)
  }
})
time_foreach[3] # 72 seconds to run (about 1 minute 12 seconds)
# Stop cluster to free up resources
parallel::stopCluster(cl)

# | 3. Use lapply to extract the specified matrices ####
selected_matrices_m.p <- lapply(1:1000, function(index) sim_no_trt_split_list[[index]])
selected_matrices_m.M <- lapply(1001:2000, function(index) sim_no_trt_split_list[[index]])
selected_matrices_m.C <- lapply(2001:3000, function(index) sim_no_trt_split_list[[index]])
selected_matrices_m.E <- lapply(3001:4000, function(index) sim_no_trt_split_list[[index]])
selected_matrices_tc <- lapply(4001:5000, function(index) sim_no_trt_split_list[[index]])
selected_matrices_te <- lapply(5001:6000, function(index) sim_no_trt_split_list[[index]])
selected_matrices_tc_hat <- lapply(6001:7000, function(index) sim_no_trt_split_list[[index]])
selected_matrices_te_hat <- lapply(7001:8000, function(index) sim_no_trt_split_list[[index]])
selected_matrices_TS <- lapply(8001:9000, function(index) sim_no_trt_split_list[[index]])

# | 4. Use do.call and rbind to append the selected matrices into a single matrix ####
m.p <- do.call(rbind, selected_matrices_m.p)
m.M <- do.call(rbind, selected_matrices_m.M)
m.C <- do.call(rbind, selected_matrices_m.C)
m.E <- do.call(rbind, selected_matrices_m.E)
tc <- do.call(rbind, selected_matrices_tc)
te <- do.call(rbind, selected_matrices_te)
tc_hat <- mean(tc)
te_hat <- mean(te)
TS <- do.call(rbind, selected_matrices_TS)

# | 5. Finalize list ####
appended_matrix_list <- list(m.p, m.M, m.C, m.E, tc, te, tc_hat, te_hat, TS)
sim_no_trt <- lapply(appended_matrix_list, as.matrix)
sim_no_trt <- setNames(sim_no_trt, c("m.p", "m.M", "m.C", "m.E", "tc", "te", 
                                     "tc_hat", "te_hat", "TS"))
n.i <- 100000
rownames(sim_no_trt$TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
colnames(sim_no_trt$TS) <- paste("Cycle", 0:(n.t-1), sep = " ")   # name the columns 
rm(list = ls()[grep("^selected_matrices_", ls())])
rm(m.p, m.M, m.C, m.E, tc, te, tc_hat, te_hat, TS)
rm(sim_no_trt_split_list, appended_matrix_list)

########## ---- Part D: Epidemiology: With NSP ##########
# | 1. Cumulative Incidence ####
# Number of people in each state (Use this to obtain cumulative incidence)
## Faster function but categories in alphabetical order
p3 = Sys.time()
prev.with.nsp <- c()
for (i in v.n) {
  prev.with.nsp[i] <- length(which(sim_trt$m.M[,-1] == i))
}
prev.with.nsp
comp.time3 <- Sys.time() - p3; comp.time3 # Approximately 8 seconds

# Group of texts to check
num_group <- c("Healthy->SSTVI", "SSTVI->SSTVI", 
               "SSTVI->Purulent_OP", "SSTVI->Purulent_ED",
               "SSTVI->NonPurulent_OP", "SSTVI->NonPurulent_ED", 
               "Purulent_OP->Purulent_IP", "Purulent_OP->Purulent_IPC",
               "Purulent_ED->Purulent_IP", "Purulent_ED->Purulent_IPC", "Purulent_ED->SSTVI",
               "Purulent_IP->Purulent_IP", "Purulent_IP->SSTVI",
               "Purulent_IPC->SSTVI", # Purulent SSTVI
               "NonPurulent_OP->NonPurulent_IP", "NonPurulent_OP->NonPurulent_IPC",
               "NonPurulent_ED->NonPurulent_IP", "NonPurulent_ED->NonPurulent_IPC", "NonPurulent_ED->SSTVI",
               "NonPurulent_IP->NonPurulent_IP", "NonPurulent_IP->SSTVI",
               "NonPurulent_IPC->SSTVI", # Purulent SSTVI
               "SSTVI->Healthy", "Purulent_OP->Healthy", "Purulent_ED->Healthy",
               "Purulent_IP->Healthy", "Purulent_IPC->Healthy",
               "NonPurulent_OP->Healthy", "NonPurulent_ED->Healthy",
               "NonPurulent_IP->Healthy", "NonPurulent_IPC->Healthy")

# Group of texts to exclude
den_group <- c("Healthy->OCM", "SSTVI->OCM", "SSTVI->SelfTrt_Death",
               "Purulent_OP->OCM", "Purulent_ED->OCM", "Purulent_IP->OCM", "Purulent_IPC->OCM", 
               "NonPurulent_OP->OCM", "NonPurulent_ED->OCM", "NonPurulent_IP->OCM", "NonPurulent_IPC->OCM", 
               "Purulent_ED->SSTVI_Death", "Purulent_IP->SSTVI_Death", "Purulent_IPC->SSTVI_Death",
               "NonPurulent_ED->SSTVI_Death", "NonPurulent_IP->SSTVI_Death", "NonPurulent_IPC->SSTVI_Death",
               "SSTVI_Death->SSTVI_Death", "SelfTrt_Death->SelfTrt_Death", "OCM->OCM")

# Compute prevalence at each week
col_num <- c(); col_den <- c(); col_prop <- c()
for (i in 1:ncol(sim_trt$TS)) {
  col_num[i] <- length(which(sim_trt$TS[,i] %in% num_group))
  col_den[i] <- length(which(!sim_trt$TS[,i] %in% den_group))
  col_prop[i] <- col_num[i] / col_den[i]
}
prev_data <- mean(col_prop); prev_data

# | 2. Person-Time Calculation ####
## Sum up the total number of weeks contributed by each patient
death_dates_with_nsp <- death_calculation(sim_trt, nsp_status = 1)
pw_total_with_nsp <- sum(death_dates_with_nsp$Death) # Person-weeks
py_total_with_nsp <- pw_total_with_nsp / 52 # Person-years
py_with_nsp <- c(sum(death_dates_with_nsp$yr1), sum(death_dates_with_nsp$yr2),
                 sum(death_dates_with_nsp$yr3), sum(death_dates_with_nsp$yr4),
                 sum(death_dates_with_nsp$yr5))

# | 3. Incidence Rate - Recurrent Events ####
options(scipen=999)
## || a. Individual components ####
person_years <- rep(py_total_with_nsp, length(prev.with.nsp))
tmp <- as.matrix(cbind(prev.with.nsp, person_years)); N. <- 1 - ((1 - 0.95)/2)
ir_with_nsp <- round((prev.with.nsp / py_total_with_nsp) * 1000, 3); ir_with_nsp
ir_with_nsp_lb <- ((qchisq(p = 1 - N., df = 2 * tmp[,1])/2)/tmp[,2])*1000; ir_with_nsp_lb
ir_with_nsp_ub <- ((qchisq(p = N., df = 2 * (tmp[,1] + 1))/2)/tmp[,2])*1000; ir_with_nsp_ub

## || b. OP, ED, IP combined ####
op_ir <- ((prev.with.nsp["Purulent_OP"] + prev.with.nsp["NonPurulent_OP"]) / py_total_with_nsp)*1000; op_ir
op_ir_ci <- prop.test(x = (prev.with.nsp["Purulent_OP"] + prev.with.nsp["NonPurulent_OP"]), 
                      n = py_total_with_nsp, conf.level = 0.95)$conf.int
op_ir_ci*1000
ed_ir <- ((prev.with.nsp["Purulent_ED"] + prev.with.nsp["NonPurulent_ED"]) / py_total_with_nsp)*1000; ed_ir
ed_ir_ci <- prop.test(x = (prev.with.nsp["Purulent_ED"] + prev.with.nsp["NonPurulent_ED"]), 
                      n = py_total_with_nsp, conf.level = 0.95)$conf.int
ed_ir_ci*1000
ip_ir <- ((prev.with.nsp["Purulent_IP"] + prev.with.nsp["NonPurulent_IP"] + 
             prev.with.nsp["Purulent_IPC"] + prev.with.nsp["NonPurulent_IPC"]) / py_total_with_nsp)*1000; ip_ir
ip_ir_ci <- prop.test(x = (prev.with.nsp["Purulent_IP"] + prev.with.nsp["NonPurulent_IP"]
                           + prev.with.nsp["Purulent_IPC"] + prev.with.nsp["NonPurulent_IPC"]), 
                      n = py_total_with_nsp, conf.level = 0.95)$conf.int
ip_ir_ci*1000

# | 4. Incidence Rate - Mortality ####
pop_yr1 <- sum(death_dates_with_nsp$yr1 != 0)
pop_yr2 <- sum(death_dates_with_nsp$yr2 != 0)
pop_yr3 <- sum(death_dates_with_nsp$yr3 != 0)
pop_yr4 <- sum(death_dates_with_nsp$yr4 != 0)
pop_yr5 <- sum(death_dates_with_nsp$yr5 != 0)

## Counts per year: other-cause mortality
ocm.1yr.with.nsp <- length(which(death_dates_with_nsp$OCM > 0 & death_dates_with_nsp$OCM <= 52))
ocm.2yr.with.nsp <- length(which(death_dates_with_nsp$OCM > 52 & death_dates_with_nsp$OCM <= 104))
ocm.3yr.with.nsp <- length(which(death_dates_with_nsp$OCM > 104 & death_dates_with_nsp$OCM <= 156))
ocm.4yr.with.nsp <- length(which(death_dates_with_nsp$OCM > 156 & death_dates_with_nsp$OCM <= 208))
ocm.5yr.with.nsp <- length(which(death_dates_with_nsp$OCM > 208 & death_dates_with_nsp$OCM <= 260))

ocm.with.nsp <- c(ocm.1yr.with.nsp, ocm.2yr.with.nsp, ocm.3yr.with.nsp, ocm.4yr.with.nsp, ocm.5yr.with.nsp)
ocm.pop <- c(pop_yr1, pop_yr2, pop_yr3, pop_yr4, pop_yr5)
ocm.rate.with.nsp <- (sum(ocm.with.nsp) / sum(ocm.pop)) * 1000; ocm.rate.with.nsp
ocm.rate.with.nsp.ci <- prop.test(sum(ocm.with.nsp), sum(ocm.pop), conf.level = 0.95)$conf.int
ocm.rate.with.nsp.ci*1000

## Counts per year: SSTVI mortality
sstvi.mort.1yr.with.nsp <- length(which(death_dates_with_nsp$SSTVI_Death > 0 & death_dates_with_nsp$SSTVI_Death <= 52))
sstvi.mort.2yr.with.nsp <- length(which(death_dates_with_nsp$SSTVI_Death > 52 & death_dates_with_nsp$SSTVI_Death <= 104))
sstvi.mort.3yr.with.nsp <- length(which(death_dates_with_nsp$SSTVI_Death > 104 & death_dates_with_nsp$SSTVI_Death <= 156))
sstvi.mort.4yr.with.nsp <- length(which(death_dates_with_nsp$SSTVI_Death > 156 & death_dates_with_nsp$SSTVI_Death <= 208))
sstvi.mort.5yr.with.nsp <- length(which(death_dates_with_nsp$SSTVI_Death > 208 & death_dates_with_nsp$SSTVI_Death <= 260))

sstvi.mort.with.nsp <- c(sstvi.mort.1yr.with.nsp, sstvi.mort.2yr.with.nsp, 
                         sstvi.mort.3yr.with.nsp, sstvi.mort.4yr.with.nsp, 
                         sstvi.mort.5yr.with.nsp)
sstvi.pop <- c(pop_yr1, pop_yr2, pop_yr3, pop_yr4, pop_yr5)
sstvi.rate.with.nsp <- (sum(sstvi.mort.with.nsp) / sum(sstvi.pop)) * 1000; sstvi.rate.with.nsp
sstvi.rate.with.nsp.ci <- prop.test(sum(sstvi.mort.with.nsp), sum(sstvi.pop), conf.level = 0.95)$conf.int
sstvi.rate.with.nsp.ci*1000

########## ---- Part E: Epidemiology: No NSP ##########
# | 1. Cumulative Incidence ####
# Number of people in each state (Use this to obtain prevalence)
## Faster function but categories in alphabetical order
p3 = Sys.time()
prev.no.nsp <- c()
for (i in v.n) {
  prev.no.nsp[i] <- length(which(sim_no_trt$m.M[,-1] == i))
}
prev.no.nsp
comp.time3 <- Sys.time() - p3; comp.time3

# | 2. Person-Time Calculation ####
## Sum up the total number of weeks contributed by each patient
death_dates_no_nsp <- death_calculation(sim_no_trt, nsp_status = 0)
pw_total_no_nsp <- sum(death_dates_no_nsp$Death) # Person-weeks
py_total_no_nsp <- pw_total_no_nsp / 52 # Person-years
py_no_nsp <- c(sum(death_dates_no_nsp$yr1), sum(death_dates_no_nsp$yr2),
               sum(death_dates_no_nsp$yr3), sum(death_dates_no_nsp$yr4),
               sum(death_dates_no_nsp$yr5))

# | 3. Incidence Rate - Recurrent Events ####
## || a. Individual components ####
options(scipen=999)
person_years <- rep(py_total_no_nsp, length(prev.no.nsp))
tmp <- as.matrix(cbind(prev.no.nsp, person_years)); N. <- 1 - ((1 - 0.95)/2)
ir_no_nsp <- round((prev.no.nsp / py_total_no_nsp) * 1000, 3); ir_no_nsp
ir_no_nsp_lb <- ((qchisq(p = 1 - N., df = 2 * tmp[,1])/2)/tmp[,2])*1000; ir_no_nsp_lb
ir_no_nsp_ub <- ((qchisq(p = N., df = 2 * (tmp[,1] + 1))/2)/tmp[,2])*1000; ir_no_nsp_ub

## || b. OP, ED, IP combined ####
op_ir <- ((prev.no.nsp["Purulent_OP"] + prev.no.nsp["NonPurulent_OP"]) / py_total_no_nsp)*1000; op_ir
op_ir_ci <- prop.test(x = (prev.no.nsp["Purulent_OP"] + prev.no.nsp["NonPurulent_OP"]), 
                      n = py_total_with_nsp, conf.level = 0.95)$conf.int
op_ir_ci*1000
ed_ir <- ((prev.no.nsp["Purulent_ED"] + prev.no.nsp["NonPurulent_ED"]) / py_total_no_nsp)*1000; ed_ir
ed_ir_ci <- prop.test(x = (prev.no.nsp["Purulent_ED"] + prev.no.nsp["NonPurulent_ED"]), 
                      n = py_total_with_nsp, conf.level = 0.95)$conf.int
ed_ir_ci*1000
ip_ir <- ((prev.no.nsp["Purulent_IP"] + prev.no.nsp["NonPurulent_IP"] + 
             prev.no.nsp["Purulent_IPC"] + prev.no.nsp["NonPurulent_IPC"]) / py_total_no_nsp)*1000; ip_ir
ip_ir_ci <- prop.test(x = (prev.no.nsp["Purulent_IP"] + prev.no.nsp["NonPurulent_IP"]
                           + prev.no.nsp["Purulent_IPC"] + prev.no.nsp["NonPurulent_IPC"]), 
                      n = py_total_with_nsp, conf.level = 0.95)$conf.int
ip_ir_ci*1000

# | 4. Incidence Rate - Mortality ####
pop_yr1 <- sum(death_dates_no_nsp$yr1 != 0)
pop_yr2 <- sum(death_dates_no_nsp$yr2 != 0)
pop_yr3 <- sum(death_dates_no_nsp$yr3 != 0)
pop_yr4 <- sum(death_dates_no_nsp$yr4 != 0)
pop_yr5 <- sum(death_dates_no_nsp$yr5 != 0)

## Counts per year: other-cause mortality
ocm.1yr.no.nsp <- length(which(death_dates_no_nsp$OCM > 0 & death_dates_no_nsp$OCM <= 52))
ocm.2yr.no.nsp <- length(which(death_dates_no_nsp$OCM > 52 & death_dates_no_nsp$OCM <= 104))
ocm.3yr.no.nsp <- length(which(death_dates_no_nsp$OCM > 104 & death_dates_no_nsp$OCM <= 156))
ocm.4yr.no.nsp <- length(which(death_dates_no_nsp$OCM > 156 & death_dates_no_nsp$OCM <= 208))
ocm.5yr.no.nsp <- length(which(death_dates_no_nsp$OCM > 208 & death_dates_no_nsp$OCM <= 260))

ocm.no.nsp <- c(ocm.1yr.no.nsp, ocm.2yr.no.nsp, ocm.3yr.no.nsp, ocm.4yr.no.nsp, ocm.5yr.no.nsp)
ocm.pop <- c(pop_yr1, pop_yr2, pop_yr3, pop_yr4, pop_yr5)
ocm.rate.no.nsp <- (sum(ocm.no.nsp) / sum(ocm.pop)) * 1000; ocm.rate.no.nsp
ocm.rate.no.nsp.ci <- prop.test(sum(ocm.no.nsp), sum(ocm.pop), conf.level = 0.95)$conf.int
ocm.rate.no.nsp.ci*1000

## Counts per year: SSTVI mortality
sstvi.mort.1yr.no.nsp <- length(which(death_dates_no_nsp$SSTVI_Death > 0 & death_dates_no_nsp$SSTVI_Death <= 52))
sstvi.mort.2yr.no.nsp <- length(which(death_dates_no_nsp$SSTVI_Death > 52 & death_dates_no_nsp$SSTVI_Death <= 104))
sstvi.mort.3yr.no.nsp <- length(which(death_dates_no_nsp$SSTVI_Death > 104 & death_dates_no_nsp$SSTVI_Death <= 156))
sstvi.mort.4yr.no.nsp <- length(which(death_dates_no_nsp$SSTVI_Death > 156 & death_dates_no_nsp$SSTVI_Death <= 208))
sstvi.mort.5yr.no.nsp <- length(which(death_dates_no_nsp$SSTVI_Death > 208 & death_dates_no_nsp$SSTVI_Death <= 260))

sstvi.mort.no.nsp <- c(sstvi.mort.1yr.no.nsp, sstvi.mort.2yr.no.nsp, 
                       sstvi.mort.3yr.no.nsp, sstvi.mort.4yr.no.nsp, 
                       sstvi.mort.5yr.no.nsp)
sstvi.pop <- c(pop_yr1, pop_yr2, pop_yr3, pop_yr4, pop_yr5)
sstvi.rate.no.nsp <- (sum(sstvi.mort.no.nsp) / sum(sstvi.pop)) * 1000; sstvi.rate.no.nsp
sstvi.rate.no.nsp.ci <- prop.test(sum(sstvi.mort.no.nsp), sum(sstvi.pop), conf.level = 0.95)$conf.int
sstvi.rate.no.nsp.ci*1000

# End of Script