################################################################
# Dynamic Transmission Model of MRSA in an ICU      	       #
# This model keeps nurses randomly mixing but separates the MD #
# Queue-based Steady State Populations             	           #
# Author: Matthew Mietchen (matthew.mietchen@wsu.edu)          #
################################################################

# Descriptive Information for PML File
Modelname: MRSA MD Cohort 
Description: PML Implementation of MRSA transmission model 

# Set model to run with numbers of individuals
Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####

# Reactions Governing Movement of Nurses (N) #
R1:
	N_u > N_c
	rho_N * sigma * N_u * (P_c / (P_c + P_u))

R2:
	N_c > N_u
	N_c * iota_N

R3:
	N_c > N_u
	N_c * tau_N * (P_c / (P_c + P_u))
	

########################

# Reactions Governing Movement of the Doctor (D) 
R4:
	D_u > D_c
	rho_D * sigma * D_u * (P_c / (P_c + P_u))

R5:
	D_c > D_u
	D_c * iota_D

R6:
	D_c > D_u
	D_c * tau_D * (P_c / (P_c + P_u))
	

########################

# Reactions Involving Uncontaminated Patients (P_u) #
R7:
	P_u > P_c + Acquisition
	rho_N * psi * P_u * (N_c / (N_c + N_u))
	
R8:
    P_u > P_c + Acquisition
    rho_D * psi * P_u * (D_c / (D_c + D_u))
	
R9:
	P_u > P_u
	theta * P_u * (1-nu)

R10:	
	P_u > P_c
	theta * P_u * nu

########################

# Reactions Involving Contaminated Patients (P_c) #	
R11:
    P_c > P_u
    mu * P_c
    
R12:
	P_c > P_c
	theta * P_c * nu

R13:
	P_c > P_u
	theta * P_c * (1-nu)
	

########################

### Parameter Values ###

## Time Values are in HOURS ##
# Compartments #

N_u = 6
N_c = 0

D_u = 1
D_c = 0

P_u = 18
P_c = 0

Acquisition = 0

### Contact Rates and Contamination Probabilities #
rho_N = 3.973   # nurse direct care tasks per patient per hour
rho_D = 0.181   # doctor direct care tasks per patient per hour 
sigma = 0.054   # hand contamination probability
psi = 0.0236    # new fitted probability of a successful colonization of an uncolonized patient

### Exit (death/discharge) rates
theta = 0.00949   # probability of death/discharge

### Admission Proportions
nu = 0.0779   # proportion of admissions of colonized with MRSA

### Handwashing and Gown/Glove Change Rates
iota_N = 6.404   # 11.92 nurse direct care tasks per hour with 56.55% compliance and 95% efficacy
iota_D = 1.748   # 3.25 doctor direct care tasks per hour with 56.55% compliance and 95% efficacy
tau_N = 2.728    # 3.30 nurse gown/glove changes per hour with 82.66% compliance
tau_D = 0.744    # 0.90 doctor gown/glove changes per hour with 82.66% compliance

mu = 0.002083   # natural decolonization rate median 20 days per Star*ICU trial
