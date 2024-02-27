# Jihoon Lim
# Economic Modelling - Microsimulation
# 03 - Analysis

# Note: Une partie de la compilation est effectuée à partir de données provenant 
# du © Gouvernement du Québec (année de la publication du Fichier de recherche: 
# 2009-2019). Le © Gouvernement du Québec n’est pas responsable des compilations
# ni de l’interprétation des résultats produits à l’aide du Fichier de recherche.

########## ---- Preface ##########
# Run the microsimulation and create the sim_trt and sim_no_trt variables in the
# 'Analysis' scripts before conducting economic evaluation and survival analysis here.

########## ---- Part A: Cost-Effectiveness Analysis ##########
# | 1. Calculate costs and QALYs ####
# Mean costs (and the MCSE) of each strategy
v.C  <- c(sim_no_trt$tc_hat, sim_trt$tc_hat) 
sd.C <- c(sd(sim_no_trt$tc), sd(sim_trt$tc)) / sqrt(n.i)
# Mean QALYs (and the MCSE) of each strategy
v.E  <- c(sim_no_trt$te_hat, sim_trt$te_hat)
sd.E <- c(sd(sim_no_trt$te), sd(sim_trt$te)) / sqrt(n.i)

delta.C <- v.C[2] - v.C[1]                   # Incremental costs
delta.E <- v.E[2] - v.E[1]                   # Incremental QALYs
sd.delta.C <- sd(sim_trt$tc - sim_no_trt$tc) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental costs
sd.delta.E <- sd(sim_trt$te - sim_no_trt$te) / sqrt(n.i) # Monte Carlo Squared Error (MCSE) of incremental QALYs
ICER    <- delta.C / delta.E                 # Calculate the ICER
results <- c(delta.C, delta.E, ICER)         # Store the values in a new variable

# | 2. Incremental cost-effectiveness analysis table ####
table_micro <- data.frame(
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
rownames(table_micro) <- c(v.Trt, "* are MCSE values")  # Name the rows
colnames(table_micro) <- c("Costs", "*",  "QALYs", "*", "Incremental Costs", "*", "QALYs Gained", "*", "ICER") # name the columns
table_micro  # Print the table

########## ---- Part B: Survival Analysis ##########
# | 1. Create counting process data ####
# Self-Treatment
p1 = Sys.time()
# Purulent SSTVI
recurrent_purulent_OP <- recurrent_event_data1(sim_trt, sim_no_trt, 'SSTVI->Purulent_OP')
recurrent_purulent_ED <- recurrent_event_data1(sim_trt, sim_no_trt, 'SSTVI->Purulent_ED')
recurrent_purulent_IP <- recurrent_event_data2(sim_trt, sim_no_trt, 'Purulent_OP->Purulent_IP', 'Purulent_ED->Purulent_IP')
recurrent_purulent_IPC <- recurrent_event_data2(sim_trt, sim_no_trt, 'Purulent_OP->Purulent_IPC', 'Purulent_ED->Purulent_IPC')

# Non-purulent SSTVI
recurrent_nonpurulent_OP <- recurrent_event_data1(sim_trt, sim_no_trt, 'SSTVI->NonPurulent_OP')
recurrent_nonpurulent_ED <- recurrent_event_data1(sim_trt, sim_no_trt, 'SSTVI->NonPurulent_ED')
recurrent_nonpurulent_IP <- recurrent_event_data2(sim_trt, sim_no_trt, 'NonPurulent_OP->NonPurulent_IP', 'NonPurulent_ED->NonPurulent_IP')
recurrent_nonpurulent_IPC <- recurrent_event_data2(sim_trt, sim_no_trt, 'NonPurulent_OP->NonPurulent_IPC', 'NonPurulent_ED->NonPurulent_IPC')

# Departure against Medical Advice
recurrent_purulent_DAA <- recurrent_event_data2(sim_trt, sim_no_trt, 'Purulent_IP->SSTVI', 'Purulent_IPC->SSTVI')
recurrent_nonpurulent_DAA <- recurrent_event_data2(sim_trt, sim_no_trt, 'NonPurulent_IP->SSTVI', 'NonPurulent_IPC->SSTVI')
recurrent_IP_DAA <- recurrent_event_data2(sim_trt, sim_no_trt, 'Purulent_IP->SSTVI', 'NonPurulent_IP->SSTVI')
recurrent_IPC_DAA <- recurrent_event_data2(sim_trt, sim_no_trt, 'Purulent_IPC->SSTVI', 'NonPurulent_IPC->SSTVI')

comp.time1 <- Sys.time() - p1; comp.time1 # Approximately 7 minutes

# | 2. Generate mean cumulative function plot ####
p2 = Sys.time()
mcf_plot(recurrent_purulent_OP, "Purulent Outpatient") # Purulent OP
mcf_plot(recurrent_purulent_ED, "Purulent Emergency") # Purulent ED
mcf_plot(recurrent_purulent_IP, "Purulent Inpatient") # Purulent IP
mcf_plot(recurrent_purulent_IPC, "Purulent Inpatient/Complications") # Purulent IPC
mcf_plot(recurrent_nonpurulent_OP, "Non-Purulent Outpatient") # Non-purulent OP
mcf_plot(recurrent_nonpurulent_ED, "Non-Purulent Emergency") # Non-purulent ED
mcf_plot(recurrent_nonpurulent_IP, "Non-Purulent Inpatient") # Non-purulent IP
mcf_plot(recurrent_nonpurulent_IPC, "Non-Purulent Inpatient/Complications") # Non-purulent IPC
mcf_plot(recurrent_purulent_DAA, "Departure Against Advice (Purulent)") # Purulent DAA
mcf_plot(recurrent_nonpurulent_DAA, "Departure Against Advice (Non-Purulent)") # Non-Purulent DAA
mcf_plot(recurrent_IP_DAA, "Departure Against Advice (Inpatient)") # Non-Purulent DAA
mcf_plot(recurrent_IPC_DAA, "Departure Against Advice (Inpatient/Complications)") # Non-Purulent DAA
comp.time2 <- Sys.time() - p2; comp.time2 # Approximately 2.5 minutes

# | 3. Marginal means model ####
p3 = Sys.time()
marginal_means(recurrent_purulent_OP) # Purulent OP
marginal_means(recurrent_purulent_ED) # Purulent ED
marginal_means(recurrent_purulent_IP) # Purulent IP
marginal_means(recurrent_purulent_IPC) # Purulent IPC
marginal_means(recurrent_nonpurulent_OP) # Non-Purulent OP
marginal_means(recurrent_nonpurulent_ED) # Non-Purulent ED
marginal_means(recurrent_nonpurulent_IP) # Non-Purulent IP
marginal_means(recurrent_nonpurulent_IPC) # Non-Purulent IPC
marginal_means(recurrent_purulent_DAA) # Purulent DAA
marginal_means(recurrent_nonpurulent_DAA) # Non-Purulent DAA
marginal_means(recurrent_IP_DAA) # IP DAA
marginal_means(recurrent_IPC_DAA) # IPC DAA
comp.time3 <- Sys.time() - p3; comp.time3 # Approximately 3.5 seconds

# | 4. Competing risks analysis ####
cr_data <- bind_rows(death_dates_with_nsp, death_dates_no_nsp)
cr_data_new <- cr_data[,(5:8)]

# Cumulative incidence function
cif_result <- cuminc(cr_data_new$stop_t, cr_data_new$sstvi_status, group = cr_data_new$nsp)
plot(cif_result$`0 1`$time, cif_result$`0 1`$est, ylim = c(0, 0.04),
     type = "s", main = "Cumulative Incidence Function (SSTVI Mortality)", 
     xlab = "Weeks", ylab = "Cumulative Incidence", col = "blue", lty = 1, lwd = 1)
lines(cif_result$`1 1`$time, cif_result$`1 1`$est, col="red", lwd = 1)
legend('topleft', legend=c("NSP", "No NSP"), col=c("red", "blue"), lty=c(1, 1), cex = 0.5)

# Competing risks regression (Fine-Gray Model)
p4 = Sys.time()
fine_gray_model <- crr(ftime = cr_data$stop_t, fstatus = cr_data$sstvi_status, cov1 = cr_data$nsp)
summary(fine_gray_model)
comp.time4 <- Sys.time() - p4; comp.time4 # Approximately 2 minutes
# Allows assessment of impact of covariates on the cumulative incidence function (CIF) 
# for a specific event while accounting for the presence of competing events.

