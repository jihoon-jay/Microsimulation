# Jihoon Lim
# Economic Modelling - Microsimulation
# 01 - Setup

# Note: Une partie de la compilation est effectuée à partir de données provenant 
# du © Gouvernement du Québec (année de la publication du Fichier de recherche: 
# 2009-2019). Le © Gouvernement du Québec n’est pas responsable des compilations
# ni de l’interprétation des résultats produits à l’aide du Fichier de recherche.

########## ---- Part A: Time Horizon and States for the Model ##########
cycle <- 52 # Cycle: 1 week
n.years   <- 5 # 5 years
n.t <- n.years*cycle # Total number of time points: 260
v.n   <- c("Healthy","SSTVI", "SelfTrt_Death", "Purulent_OP", "Purulent_ED", 
           "Purulent_IP", "Purulent_IPC", "NonPurulent_OP", "NonPurulent_ED", 
           "NonPurulent_IP", "NonPurulent_IPC", "SSTVI_Death", "OCM") # Model states: 20 in total
n.s   <- length(v.n) # Number of health states
d.c <- d.e <- 0.015           # Equal discounting of costs and QALYs by 1.5%
v.dwc <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))
v.dwe <- 1 / (1 + d.c) ^ (rep(0:(n.years-1), each=(n.t/n.years)))

########## ---- Part B: Transition Probabilities (per cycle) ##########
# Transition probabilities (per cycle)
# | 1. Self-Treatment of SSTVI ####
p_SelfTreat <- 0.75 #0.474 #0.46 #0.857  #0.72 #0.3272#0.63 # Proportion of self-treatment
p_SSTVI_to_SelfTrt_Death <- p_SelfTreat*0.011012 #0.011012, 0.018091, 0.001049
p_SSTVI_to_SSTVI <- p_SelfTreat*0.65 # p_SelfTreat*0.71 # Accounts for PWID delaying treatment/not improving on their own (not seeking trt 0.3208)

# | 2. Transition from SSTVI to Healthcare settings ####
p_HealthUtilization <- 1-p_SelfTreat
p_Purulent <- 0.171; p_NonPurulent <- 0.829
p_SSTVI_to_OP_Purulent <- p_Purulent*p_HealthUtilization*0.8565
p_SSTVI_to_ED_Purulent <- p_Purulent*p_HealthUtilization*0.1435
p_SSTVI_to_OP_NonPurulent <- p_NonPurulent*p_HealthUtilization*0.8720
p_SSTVI_to_ED_NonPurulent <- p_NonPurulent*p_HealthUtilization*0.1280

# | 3. Purulent SSTVI ####
p_OP_to_IP_Purulent <- 0.0362; p_OP_to_IPC_Purulent <- 0.0223
p_ED_to_IP_Purulent <- 0.0333; p_ED_to_IPC_Purulent <- 0.0530
p_ED_to_SSTVI_Purulent <- 0.002468 # ED departure against advice (Assume same as non-purulent!)
p_ED_to_SSTVI_Death_Purulent <- 0.000048 # SSTVI death (Assume same as non-purulent!)

p_IP_to_IP_Purulent <- 0.0124 # Re-admission (small cell count but non-trivial percentage! Assume same IP->IP as non-purulent!)
p_IP_to_SSTVI_Purulent <- 0.0909 # Hospital departure against advice
p_IP_to_SSTVI_Death_Purulent <- 0.000212 # SSTVI death

p_IPC_to_SSTVI_Purulent <- 0.037 # Hospital departure against advice (Assume same as non-purulent!)
p_IPC_to_SSTVI_Death_Purulent <- 0.000359 # SSTVI death (Assume same as non-purulent!)

# | 4. Non-Purulent SSTVI ####
p_OP_to_IP_NonPurulent <- 0.0411; p_OP_to_IPC_NonPurulent <- 0.0042
p_ED_to_IP_NonPurulent <- 0.0887; p_ED_to_IPC_NonPurulent <- 0.0198
p_ED_to_SSTVI_NonPurulent <- 0.002468 # ED departure against advice
p_ED_to_SSTVI_Death_NonPurulent <- 0.000048 # SSTVI death

p_IP_to_IP_NonPurulent <- 0.0124 # Re-admission
p_IP_to_SSTVI_NonPurulent <- 0.0579 # Hospital departure against advice
p_IP_to_SSTVI_Death_NonPurulent <- 0.000387 # SSTVI death

p_IPC_to_SSTVI_NonPurulent <- 0.037 # Hospital departure against advice
p_IPC_to_SSTVI_Death_NonPurulent <- 0.000359 # SSTVI death

########## ---- Part D: Cost Parameters ##########
# | 1. NSP Costs ####
# Cost of NSP (Create vectors and account for CPI)
cpi.2016 <- 122.4;  cpi.2022 <- 139.3

# Subject to change
c_NSP.2016 <- 283 # Ontario total in 2015 but assume this is true for Quebec
# CPI adjusted costs (2022 CAD)
c_NSP.2022 <- c_NSP.2016 * cpi.2022 / cpi.2016
# Average cost of NSP per week per person (undiscounted)
c_NSP <- round(c_NSP.2022 / cycle, 2)

# | 2. Healthcare Costs ####
c_Healthy <- 0 + c_NSP; c_SSTVI <- 0 + c_NSP; c_SelfTrt_Death <- 0
c_OP_Purulent <- 100.78 + c_NSP; c_ED_Purulent <- 1135.70 + c_NSP
c_IP_Purulent <- 6022.35 + c_NSP; c_IPC_Purulent <- 42526.40 + c_NSP
c_OP_NonPurulent <- 93.73 + c_NSP; c_ED_NonPurulent <- 1144.64 + c_NSP
c_IP_NonPurulent <- 6821.06 + c_NSP; c_IPC_NonPurulent <- 41016.23 + c_NSP
c_SSTVI_Death <- 0; c_OCM <- 0

########## ---- Part D: Utility Parameters ##########
pwid_multiplier <- 0.9 # Otherwise, go with 0.9
# "Evaluating the cost-effectiveness of needle and syringe programmes in 
# preventing Hepatitis C transmission in people who inject drugs" (this paper
# has the PWID multiplier.)
q_Healthy <- pwid_multiplier/52; q_SSTVI <- pwid_multiplier*0.97607/52
q_OP_Purulent <- q_OP_NonPurulent <- pwid_multiplier*0.97607/52
q_ED_Purulent <- q_ED_NonPurulent <- pwid_multiplier*0.97607/52
q_IP_Purulent <- q_IP_NonPurulent <- pwid_multiplier*0.642/52
q_IPC_Purulent <- q_IPC_NonPurulent <- pwid_multiplier*0.642/52
q_SelfTrt_Death <- q_SSTVI_Death <- q_OCM <- 0
