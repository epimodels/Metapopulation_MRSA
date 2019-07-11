import os
import stochpy
import random
import numpy as numpy
from scipy import stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--suffix", help="suffix for files")
args = parser.parse_args()
suffix = args.suffix

workingdir = os.getcwd()

# Simulation parameters
start_time = 0.0
end_time = 8760
n_runs = 250

# Run is a single run of the MetaPop model that returns the number of incident cases
    
def Leakyrun(pdict):
    model = stochpy.SSA()
    model.Model(model_file='Gamma.psc', dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('gamma',pdict['gamma'])
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[16][0][-1]
    return cases

def Leakysens(iterations):
    parameters = numpy.zeros([iterations,2])
    for k in range(0,iterations):
        pdict = {'gamma':random.uniform(0.1667,1.000)}
         
        leakycases = Leakyrun(pdict=pdict)
        parameters[k,0] = leakycases
        parameters[k,1] = pdict['gamma']
    return parameters

print("Starting Gamma Sweep")
gamma_sweep = Leakysens(n_runs)

print("Saving Files")
numpy.savetxt(''.join(['gamma_sweep_',suffix,'.csv']),gamma_sweep,delimiter=',',comments=',')

