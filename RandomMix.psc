#######################################################
# Dynamic Transmission Model of MRSA in an ICU		  #
# Queue-based Steady State Population                 #
# Author: Matthew Mietchen (matthew.mietchen@wsu.edu) #
#######################################################

# Descriptive Information for PML File
Modelname: MRSA ICU Random Mixing
Description: PML Implementation of MRSA transmission model 

# Set model to run with numbers of individuals
Species_In_Conc: False
Output_In_Conc: False

### Model Reactions ###

# Reactions Governing Movement of Contaminated Staff (S_c) and Uncontaminated Staff (S_u)#
R1:
	S_c > S_u
	S_c * iota

R2:
	S_c > S_u
	S_c * tau * (P_c/(P_c+P_u))
	
R3:
	S_u > S_c
	rho * sigma * S_u * (P_c/(P_c + P_u))

# Reactions Involving Uncontaminated, Low Risk Patients (P_u) #
R4:
	P_u > P_c + Acquisition
	rho * psi * P_u * (S_c/(S_c + S_u))
	
R5:
	P_u > P_u
	theta * P_u * (1-nu)

R6:	
	P_u > P_c
	theta * P_u * nu

	
# Reactions Involving Contaminated, Low Risk Patients (P_c) #
	
R7:
	P_c > P_c
	theta * P_c * nu

R8:
	P_c > P_u
	theta * P_c * (1-nu)
	
R9:
    P_c > P_u
    mu * P_c
    

### Parameter Values ###
## Time Values are in HOURS ##
# Compartments #

S_u = 7
S_c = 0
P_u = 18
P_c = 0
Acquisition = 0

# Contact Rates and Contamination Probabilities
rho = 4.154  	# direct care tasks per patient per hour
sigma = 0.054
psi = 0.0236    # new fitted probability of a successful colonization of an uncolonized patient

# Exit (death/discharge) rates
theta = 0.00949

# Admission Proportions
nu = 0.0779

# Handwashing and Gown/Glove Change Rates
iota = 5.74  	# 10.682 direct care tasks per hour with 56.55% compliance and 95% efficacy
tau = 2.445

mu = 0.002083  	# natural decolonization rate median 20 days per Star*ICU trial
