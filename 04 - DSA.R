# Jihoon Lim
# Economic Modelling - Microsimulation
# 04 - Deterministic Sensitivity Analysis

# Note: Une partie de la compilation est effectuée à partir de données provenant 
# du © Gouvernement du Québec (année de la publication du Fichier de recherche: 
# 2009-2019). Le © Gouvernement du Québec n’est pas responsable des compilations
# ni de l’interprétation des résultats produits à l’aide du Fichier de recherche.

########## ---- Part A: Discounting Rates ##########
# | 1. At 3% ####
d.c <- d.e <- 0.03          # Equal discounting of costs and QALYs by 3%
v.dwc <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))
v.dwe <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))
v.dwc <- as.matrix(as.data.frame(v.dwc))
v.dwe <- as.matrix(as.data.frame(v.dwe))
## | a. With NSP ####
tc003_trt <- sim_trt$m.C %*% v.dwc     # total (discounted) cost per individual
te003_trt <- sim_trt$m.E %*% v.dwe       # total (discounted) QALYs per individual 
tc003_hat_trt <- mean(tc003_trt)        # average (discounted) cost 
te003_hat_trt <- mean(te003_trt)        # average (discounted) QALYs
## | b. No NSP ####
tc003_no_trt <- sim_no_trt$m.C %*% v.dwc     # total (discounted) cost per individual
te003_no_trt <- sim_no_trt$m.E %*% v.dwe       # total (discounted) QALYs per individual 
tc003_hat_no_trt <- mean(tc003_no_trt)        # average (discounted) cost 
te003_hat_no_trt <- mean(te003_no_trt)        # average (discounted) QALYs

# | 2. At 0% ####
d.c <- d.e <- 0          # Equal discounting of costs and QALYs by 0%
v.dwc <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))
v.dwe <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))
v.dwc <- as.matrix(as.data.frame(v.dwc))
v.dwe <- as.matrix(as.data.frame(v.dwe))
## | a. With NSP ####
tc000_trt <- sim_trt$m.C %*% v.dwc     # total (discounted) cost per individual
te000_trt <- sim_trt$m.E %*% v.dwe       # total (discounted) QALYs per individual 
tc000_hat_trt <- mean(tc000_trt)        # average (discounted) cost 
te000_hat_trt <- mean(te000_trt)        # average (discounted) QALYs
## | b. No NSP ####
tc000_no_trt <- sim_no_trt$m.C %*% v.dwc     # total (discounted) cost per individual
te000_no_trt <- sim_no_trt$m.E %*% v.dwe       # total (discounted) QALYs per individual 
tc000_hat_no_trt <- mean(tc000_no_trt)        # average (discounted) cost 
te000_hat_no_trt <- mean(te000_no_trt)        # average (discounted) QALYs

########## ---- Part B: DSA CEA for Discount Rates ##########
# | 1. Discounting rate at 3% ####
# || a. Calculate costs and QALYs ####
# Mean costs (and the MCSE) of each strategy
v.C  <- c(tc003_hat_no_trt, tc003_hat_trt) 
sd.C <- c(sd(tc003_no_trt), sd(tc003_trt)) / sqrt(n.i)
# Mean QALYs (and the MCSE) of each strategy
v.E  <- c(te003_hat_no_trt, te003_hat_trt)
sd.E <- c(sd(te003_no_trt), sd(te003_trt)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # Incremental costs
delta.E <- v.E[2] - v.E[1]                   # Incremental QALYs
sd.delta.C <- sd(tc003_trt - tc003_no_trt) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental costs
sd.delta.E <- sd(te003_trt - te003_no_trt) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # Calculate the ICER
results <- c(delta.C, delta.E, ICER)         # Store the values in a new variable

# || b. Incremental cost-effectiveness analysis table ####
table_micro003 <- data.frame(
  c(round(v.C, 0),  ""),           # Costs per arm
  c(round(sd.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # Health outcomes per arm
  c(round(sd.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # Incremental costs
  c("", round(sd.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),  # Incremental QALYs 
  c("", round(sd.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
v.Trt <- c("No Treatment", "Treatment") # Store the strategy names
rownames(table_micro003) <- c(v.Trt, "* are MCSE values")  # Name the rows
colnames(table_micro003) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro003  # Print the table

# | 2. Discounting rate at 0% ####
# || a. Calculate costs and QALYs ####
# Mean costs (and the MCSE) of each strategy
v.C  <- c(tc000_hat_no_trt, tc000_hat_trt) 
sd.C <- c(sd(tc000_no_trt), sd(tc000_trt)) / sqrt(n.i)
# Mean QALYs (and the MCSE) of each strategy
v.E  <- c(te000_hat_no_trt, te000_hat_trt)
sd.E <- c(sd(te000_no_trt), sd(te000_trt)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # Incremental costs
delta.E <- v.E[2] - v.E[1]                   # Incremental QALYs
sd.delta.C <- sd(tc000_trt - tc000_no_trt) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental costs
sd.delta.E <- sd(te000_trt - te000_no_trt) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # Calculate the ICER
results <- c(delta.C, delta.E, ICER)         # Store the values in a new variable

# || b. Incremental cost-effectiveness analysis table ####
table_micro000 <- data.frame(
  c(round(v.C, 0),  ""),           # Costs per arm
  c(round(sd.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # Health outcomes per arm
  c(round(sd.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # Incremental costs
  c("", round(sd.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),  # Incremental QALYs 
  c("", round(sd.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
v.Trt <- c("No Treatment", "Treatment") # Store the strategy names
rownames(table_micro000) <- c(v.Trt, "* are MCSE values")  # Name the rows
colnames(table_micro000) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro000  # Print the table

########## ---- Part C: NSP Costs ##########
# Average cost of NSP per week per person (undiscounted)
c_NSP300 <- round(300 / cycle, 2); c_NSP350 <- round(350 / cycle, 2)
c_NSP375 <- round(375 / cycle, 2); c_NSP400 <- round(400 / cycle, 2)
m.C.new <- sim_trt$m.C - c_NSP
# Discounting rate at 1.5%
d.c <- d.e <- 0.015
v.dwc <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))
v.dwe <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))
v.dwc <- as.matrix(as.data.frame(v.dwc))
v.dwe <- as.matrix(as.data.frame(v.dwe))

# | 1. Assume annual cost is $300 (in 2022 CAD) ####
m.C.300 <- m.C.new + c_NSP300
tc.nsp.300 <- m.C.300 %*% v.dwc; te.nsp.300 <- sim_trt$m.E %*% v.dwc
tc_hat.nsp.300 <- mean(tc.nsp.300); te_hat.nsp.300 <- mean(te.nsp.300)
# || a. Calculate costs and QALYs ####
# Mean costs (and the MCSE) of each strategy
v.C  <- c(sim_no_trt$tc_hat, tc_hat.nsp.300) 
sd.C <- c(sd(sim_no_trt$tc), sd(tc.nsp.300)) / sqrt(n.i)
# Mean QALYs (and the MCSE) of each strategy
v.E  <- c(sim_no_trt$te_hat, te_hat.nsp.300)
sd.E <- c(sd(sim_no_trt$te), sd(te.nsp.300)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # Incremental costs
delta.E <- v.E[2] - v.E[1]                   # Incremental QALYs
sd.delta.C <- sd(tc.nsp.300 - sim_no_trt$tc) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental costs
sd.delta.E <- sd(te.nsp.300 - sim_no_trt$te) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # Calculate the ICER
results <- c(delta.C, delta.E, ICER)         # Store the values in a new variable

# || b. Incremental cost-effectiveness analysis table ####
table_micro_nsp_300 <- data.frame(
  c(round(v.C, 0),  ""),           # Costs per arm
  c(round(sd.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # Health outcomes per arm
  c(round(sd.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # Incremental costs
  c("", round(sd.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),  # Incremental QALYs 
  c("", round(sd.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
v.Trt <- c("No Treatment", "Treatment") # Store the strategy names
rownames(table_micro_nsp_300) <- c(v.Trt, "* are MCSE values")  # Name the rows
colnames(table_micro_nsp_300) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro_nsp_300  # Print the table

# | 2. Assume annual cost is $350 (in 2022 CAD) ####
m.C.350 <- m.C.new + c_NSP350
tc.nsp.350 <- m.C.350 %*% v.dwc; te.nsp.350 <- sim_trt$m.E %*% v.dwc
tc_hat.nsp.350 <- mean(tc.nsp.350); te_hat.nsp.350 <- mean(te.nsp.350)
# || a. Calculate costs and QALYs ####
# Mean costs (and the MCSE) of each strategy
v.C  <- c(sim_no_trt$tc_hat, tc_hat.nsp.350) 
sd.C <- c(sd(sim_no_trt$tc), sd(tc.nsp.350)) / sqrt(n.i)
# Mean QALYs (and the MCSE) of each strategy
v.E  <- c(sim_no_trt$te_hat, te_hat.nsp.350)
sd.E <- c(sd(sim_no_trt$te), sd(te.nsp.350)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # Incremental costs
delta.E <- v.E[2] - v.E[1]                   # Incremental QALYs
sd.delta.C <- sd(tc.nsp.350 - sim_no_trt$tc) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental costs
sd.delta.E <- sd(te.nsp.350 - sim_no_trt$te) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # Calculate the ICER
results <- c(delta.C, delta.E, ICER)         # Store the values in a new variable

# || b. Incremental cost-effectiveness analysis table ####
table_micro_nsp_350 <- data.frame(
  c(round(v.C, 0),  ""),           # Costs per arm
  c(round(sd.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # Health outcomes per arm
  c(round(sd.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # Incremental costs
  c("", round(sd.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),  # Incremental QALYs 
  c("", round(sd.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
v.Trt <- c("No Treatment", "Treatment") # Store the strategy names
rownames(table_micro_nsp_350) <- c(v.Trt, "* are MCSE values")  # Name the rows
colnames(table_micro_nsp_350) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro_nsp_350  # Print the table

# | 3. Assume annual cost is $375 (in 2022 CAD) ####
m.C.375 <- m.C.new + c_NSP375
tc.nsp.375 <- m.C.375 %*% v.dwc; te.nsp.375 <- sim_trt$m.E %*% v.dwc
tc_hat.nsp.375 <- mean(tc.nsp.375); te_hat.nsp.375 <- mean(te.nsp.375)
# || a. Calculate costs and QALYs ####
# Mean costs (and the MCSE) of each strategy
v.C  <- c(sim_no_trt$tc_hat, tc_hat.nsp.375) 
sd.C <- c(sd(sim_no_trt$tc), sd(tc.nsp.375)) / sqrt(n.i)
# Mean QALYs (and the MCSE) of each strategy
v.E  <- c(sim_no_trt$te_hat, te_hat.nsp.375)
sd.E <- c(sd(sim_no_trt$te), sd(te.nsp.375)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # Incremental costs
delta.E <- v.E[2] - v.E[1]                   # Incremental QALYs
sd.delta.C <- sd(tc.nsp.375 - sim_no_trt$tc) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental costs
sd.delta.E <- sd(te.nsp.375 - sim_no_trt$te) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # Calculate the ICER
results <- c(delta.C, delta.E, ICER)         # Store the values in a new variable

# || b. Incremental cost-effectiveness analysis table ####
table_micro_nsp_375 <- data.frame(
  c(round(v.C, 0),  ""),           # Costs per arm
  c(round(sd.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # Health outcomes per arm
  c(round(sd.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # Incremental costs
  c("", round(sd.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),  # Incremental QALYs 
  c("", round(sd.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
v.Trt <- c("No Treatment", "Treatment") # Store the strategy names
rownames(table_micro_nsp_375) <- c(v.Trt, "* are MCSE values")  # Name the rows
colnames(table_micro_nsp_375) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro_nsp_375  # Print the table

# | 4. Assume annual cost is $400 (in 2022 CAD) ####
m.C.400 <- m.C.new + c_NSP400
tc.nsp.400 <- m.C.400 %*% v.dwc; te.nsp.400 <- sim_trt$m.E %*% v.dwc
tc_hat.nsp.400 <- mean(tc.nsp.400); te_hat.nsp.400 <- mean(te.nsp.400)
# || a. Calculate costs and QALYs ####
# Mean costs (and the MCSE) of each strategy
v.C  <- c(sim_no_trt$tc_hat, tc_hat.nsp.400) 
sd.C <- c(sd(sim_no_trt$tc), sd(tc.nsp.400)) / sqrt(n.i)
# Mean QALYs (and the MCSE) of each strategy
v.E  <- c(sim_no_trt$te_hat, te_hat.nsp.400)
sd.E <- c(sd(sim_no_trt$te), sd(te.nsp.400)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # Incremental costs
delta.E <- v.E[2] - v.E[1]                   # Incremental QALYs
sd.delta.C <- sd(tc.nsp.400 - sim_no_trt$tc) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental costs
sd.delta.E <- sd(te.nsp.400 - sim_no_trt$te) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # Calculate the ICER
results <- c(delta.C, delta.E, ICER)         # Store the values in a new variable

# || b. Incremental cost-effectiveness analysis table ####
table_micro_nsp_400 <- data.frame(
  c(round(v.C, 0),  ""),           # Costs per arm
  c(round(sd.C, 0), ""),           # MCSE for costs
  c(round(v.E, 3),  ""),           # Health outcomes per arm
  c(round(sd.E, 3), ""),           # MCSE for health outcomes
  c("", round(delta.C, 0),   ""),  # Incremental costs
  c("", round(sd.delta.C, 0),""),  # MCSE for incremental costs
  c("", round(delta.E, 3),   ""),  # Incremental QALYs 
  c("", round(sd.delta.E, 3),""),  # MCSE for health outcomes (QALYs) gained
  c("", round(ICER, 0),      "")   # ICER
)
v.Trt <- c("No Treatment", "Treatment") # Store the strategy names
rownames(table_micro_nsp_400) <- c(v.Trt, "* are MCSE values")  # Name the rows
colnames(table_micro_nsp_400) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro_nsp_400  # Print the table
