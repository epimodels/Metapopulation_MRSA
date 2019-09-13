############################################################
# Dynamic MRSA Transmission Model in ICU            	   #
# Metapopulation Model                                     #
# Cohorts depict a group of patients assigned to one nurse #
############################################################

# Descriptive Information for PML File
Modelname: Metapopulation Model
Description: PML Implementation of MRSA transmission model 

# Set model to run with numbers of individuals
Species_In_Conc: False
Output_In_Conc: False

#### Model Reactions ####

# Reactions Governing Movement of Nurses (N) Cohort 1 #
R1:
	N_u1 > N_c1
	rho_N * sigma * N_u1 * (P_c1 / (P_c1 + P_u1))

R2:
	N_c1 > N_u1
	N_c1 * iota_N

R3:
	N_c1 > N_u1
	N_c1 * tau_N * (P_c1 / (P_c1 + P_u1))
	

# Reactions Governing Movement of Nurses (N) Cohort 2 #
R4:
	N_u2 > N_c2
	rho_N * sigma * N_u2 * (P_c2 / (P_c2 + P_u2))

R5:
	N_c2 > N_u2
	N_c2 * iota_N

R6:
	N_c2 > N_u2
	N_c2 * tau_N * (P_c2 / (P_c2 + P_u2))

# Reactions Governing Movement of Nurses (N) Cohort 3 #
R7:
	N_u3 > N_c3
	rho_N * sigma * N_u3 * (P_c3 / (P_c3 + P_u3))

R8:
	N_c3 > N_u3
	N_c3 * iota_N

R9:
	N_c3 > N_u3
	N_c3 * tau_N * (P_c3 / (P_c3 + P_u3))

# Reactions Governing Movement of Nurses (N) Cohort 4 #
R10:
	N_u4 > N_c4
	rho_N * sigma * N_u4 * (P_c4 / (P_c4 + P_u4))

R11:
	N_c4 > N_u4
	N_c4 * iota_N

R12:
	N_c4 > N_u4
	N_c4 * tau_N * (P_c4 / (P_c4 + P_u4))

# Reactions Governing Movement of Nurses (N) Cohort 5 #
R13:
	N_u5 > N_c5
	rho_N * sigma * N_u5 * (P_c5 / (P_c5 + P_u5))

R14:
	N_c5 > N_u5
	N_c5 * iota_N

R15:
	N_c5 > N_u5
	N_c5 * tau_N * (P_c5 / (P_c5 + P_u5))

# Reactions Governing Movement of Nurses (N) Cohort 6 #
R16:
	N_u6 > N_c6
	rho_N * sigma * N_u6 * (P_c6 / (P_c6 + P_u6))

R17:
	N_c6 > N_u6
	N_c6 * iota_N

R18:
	N_c6 > N_u6
	N_c6 * tau_N * (P_c6 / (P_c6 + P_u6))

########################

# Reactions Governing Movement of the Doctor 
R19:
	D_u > D_c
	rho_D * sigma * D_u * (P_c1+P_c2+P_c3+P_c4+P_c5+P_c6 / (P_c1+P_c2+P_c3+P_c4+P_c5+P_c6 + P_u1+P_u2+P_u3+P_u4+P_u5+P_u6))

R20:
	D_c > D_u
	D_c * iota_D

R21:
	D_c > D_u
	D_c * tau_D * (P_c1+P_c2+P_c3+P_c4+P_c5+P_c6 / (P_c1+P_c2+P_c3+P_c4+P_c5+P_c6 + P_u1+P_u2+P_u3+P_u4+P_u5+P_u6))
	

########################

# Reactions Involving Uncontaminated Patients (P_u) Cohort 1 #
R22:
	P_u1 > P_c1 + Acquisition
	rho_N * psi * P_u1 * (N_c1 / (N_c1 + N_u1))

R23:
	P_u1 > P_c1 + Acquisition
	rho_D * psi * P_u1 * (D_c / (D_c + D_u))
	
R24:
	P_u1 > P_u1
	theta * P_u1 * (1-nu)

R25:	
	P_u1 > P_c1
	theta * P_u1 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 2 #
R26:
	P_u2 > P_c2 + Acquisition
	rho_N * psi * P_u2 * (N_c2 / (N_c2 + N_u2))

R27:
	P_u2 > P_c2 + Acquisition
	rho_D * psi * P_u2 * (D_c / (D_c + D_u))
	
R28:
	P_u2 > P_u2
	theta * P_u2 * (1-nu)

R29:	
	P_u2 > P_c2
	theta * P_u2 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 3 #
R30:
	P_u3 > P_c3 + Acquisition
	rho_N * psi * P_u3 * (N_c3 / (N_c3 + N_u3))

R31:
	P_u3 > P_c3 + Acquisition
	rho_D * psi * P_u3 * (D_c / (D_c + D_u))
	
R32:
	P_u3 > P_u3
	theta * P_u3 * (1-nu)

R33:	
	P_u3 > P_c3
	theta * P_u3 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 4 #
R34:
	P_u4 > P_c4 + Acquisition
	rho_N * psi * P_u4 * (N_c4 / (N_c4 + N_u4))

R35:
	P_u4 > P_c4 + Acquisition
	rho_D * psi * P_u4 * (D_c / (D_c + D_u))
	
R36:
	P_u4 > P_u4
	theta * P_u4 * (1-nu)

R37:	
	P_u4 > P_c4
	theta * P_u4 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 5 #
R38:
	P_u5 > P_c5 + Acquisition
	rho_N * psi * P_u5 * (N_c5 / (N_c5 + N_u5))

R39:
	P_u5 > P_c5 + Acquisition
	rho_D * psi * P_u5 * (D_c / (D_c + D_u))
	
R40:
	P_u5 > P_u5
	theta * P_u5 * (1-nu)

R41:	
	P_u5 > P_c5
	theta * P_u5 * nu

# Reactions Involving Uncontaminated Patients (P_u) Cohort 6 #
R42:
	P_u6 > P_c6 + Acquisition
	rho_N * psi * P_u6 * (N_c6 / (N_c6 + N_u6))

R43:
	P_u6 > P_c6 + Acquisition
	rho_D * psi * P_u6 * (D_c / (D_c + D_u))
	
R44:
	P_u6 > P_u6
	theta * P_u6 * (1-nu)

R45:	
	P_u6 > P_c6
	theta * P_u6 * nu

########################

# Reactions Involving Contaminated Patients (P_c) Cohort 1 #	
R46:
    P_c1 > P_u1
    mu * P_c1
    
R47:
	P_c1 > P_c1
	theta * P_c1 * nu

R48:
	P_c1 > P_u1
	theta * P_c1 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 2 #
R49:
    P_c2 > P_u2
    mu * P_c2

R50:
	P_c2 > P_c2
	theta * P_c2 * nu

R51:
	P_c2 > P_u2
	theta * P_c2 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 3 #
R52:
    P_c3 > P_u3
    mu * P_c3
    
R53:
	P_c3 > P_c3
	theta * P_c3 * nu

R54:
	P_c3 > P_u3
	theta * P_c3 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 4 #
R55:
    P_c4 > P_u4
    mu * P_c4
    
R56:
	P_c4 > P_c4
	theta * P_c4 * nu

R57:
	P_c4 > P_u4
	theta * P_c4 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 5 #
R58:
    P_c5 > P_u5
    mu * P_c5
    
R59:
	P_c5 > P_c5
	theta * P_c5 * nu

R60:
	P_c5 > P_u5
	theta * P_c5 * (1-nu)

# Reactions Involving Contaminated Patients (P_c) Cohort 6 #
R61:
    P_c6 > P_u6
    mu * P_c6
    
R62:
	P_c6 > P_c6
	theta * P_c6 * nu

R63:
	P_c6 > P_u6
	theta * P_c6 * (1-nu)

########################

### Parameter Values ###

## Time Values are in HOURS ##
# Compartments #
N_u1 = 1
N_u2 = 1
N_u3 = 1
N_u4 = 1
N_u5 = 1
N_u6 = 1

N_c1 = 0
N_c2 = 0
N_c3 = 0
N_c4 = 0
N_c5 = 0
N_c6 = 0

D_u = 1
D_c = 0

P_u1 = 3
P_u2 = 3
P_u3 = 3
P_u4 = 3
P_u5 = 3
P_u6 = 3

P_c1 = 0
P_c2 = 0
P_c3 = 0
P_c4 = 0
P_c5 = 0
P_c6 = 0

Acquisition = 0

# Contact Rates and Contamination Probabilities #
rho_N = 3.973 	# nurse direct care tasks per patient per hour
rho_D = 0.181 	# doctor direct care tasks per patient per hour 
sigma = 0.054 	# hand contamination probability
psi = 0.0236    # new fitted probability of a successful colonization of an uncolonized patient

# Exit (death/discharge) rates
theta = 0.00949 # probability of death/discharge

# Admission Proportions
nu = 0.0779 	# proportion of admissions of colonized with MRSA

# Handwashing and Gown/Glove Change Rates
iota_N = 6.404 	# 11.92 nurse direct care tasks per hour with 56.55% compliance and 95% efficacy
iota_D = 1.748 	# 3.25 doctor direct care tasks per hour with 56.55% compliance and 95% efficacy
tau_N = 2.728 	# 3.30 nurse gown/glove changes per hour with 82.66% compliance
tau_D = 0.744 	# 0.90 doctor gown/glove changes per hour with 82.66% compliance

mu = 0.002083 	# natural decolonization rate median 20 days per Star*ICU trial
