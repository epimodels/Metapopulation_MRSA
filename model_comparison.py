######################################################################
# Model Outcome Comparisons                                          #
# MRSA Colonization Acquisitions among various Population Structures #
######################################################################

# Module Imports
import os
import stochpy 
import numpy as numpy
import random

workingdir = os.getcwd()

# General simulation parameters
random.seed(80209)
start_time = 0.0
end_time = 8760
n_runs = 1000

# Model output storage arrays
acquisitions = numpy.empty([n_runs,3])

# Single Staff Type (SST) Model Run
def SST(iteration):
    model = stochpy.SSA()
    model.Model(model_file='SST.psc', dir=workingdir)
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[4][0][-1] 
    acquisitions[iteration,0] = Incident

# Nurse_MD Model Run
def Nurse_MD(iteration):
    model = stochpy.SSA()
    model.Model(model_file='Nurse_MD.psc', dir=workingdir)
    model.Endtime(end_time)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    Incident = outcomes[6][0][-1]
    acquisitions[iteration,1] = Incident

# Metapopulation Model Run
def MetaPop(iteration):
	model = stochpy.SSA()
	model.Model(model_file='MetaPop.psc', dir=workingdir)
	model.Endtime(end_time)
	model.DoStochSim()
	model.GetRegularGrid(n_samples=end_time)
	outcomes = model.data_stochsim_grid.species
	Incident = outcomes[16][0][-1]
	acquisitions[iteration,2] = Incident
		

for i in range(0,n_runs):
    print("*** Iteration %i of %i ***" % (i+1,n_runs))
    Random(i)
    Nurse_MD(i)
    MetaPop(i)
    
numpy.savetxt('model_data.csv',acquisitions,delimiter=',',header="SST,Nurse_MD,MetaPop",comments='')

print("*************************")
print("***** Runs Complete *****")
print("*************************")
