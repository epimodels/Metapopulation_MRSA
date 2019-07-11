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
n_runs = 2500

# Run is a single run of the model that returns the number of incident cases
def Baserun(pdict):
    model = stochpy.SSA()
    model.Model(model_file='SST.psc', dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('rho',pdict['rho']*4.154)
    model.ChangeParameter('sigma',pdict['sigma']*0.054)
    model.ChangeParameter('psi',pdict['psi']*0.0236)
    model.ChangeParameter('theta',pdict['theta']*0.00949)
    model.ChangeParameter('nu',pdict['nu']*0.0779)
    model.ChangeParameter('iota',pdict['iota']*5.74)
    model.ChangeParameter('tau',pdict['tau']*2.445)
    model.ChangeParameter('mu',pdict['mu']*0.002083)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[4][0][-1]
    return cases
    
def Seperatedrun(filename,pdict,type):
    model = stochpy.SSA()
    model.Model(model_file=filename, dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('rho_N',pdict['rho_N']*3.973)
    model.ChangeParameter('rho_D',pdict['rho_D']*0.181)
    model.ChangeParameter('sigma',pdict['sigma']*0.054)
    if type==1:
        model.ChangeParameter('psi',pdict['psi']*0.0285) # Nurse_MD Model Psi
    elif type==2:
        model.ChangeParameter('psi',pdict['psi']*0.0464) # MetaPop Model Psi
    model.ChangeParameter('theta',pdict['theta']*0.00949)
    model.ChangeParameter('nu',pdict['nu']*0.0779)
    model.ChangeParameter('iota_N',pdict['iota_N']*6.404)
    model.ChangeParameter('iota_D',pdict['iota_D']*1.748)
    model.ChangeParameter('tau_N',pdict['tau_N']*2.728)
    model.ChangeParameter('tau_D',pdict['tau_D']*0.744)
    model.ChangeParameter('mu',pdict['mu']*0.002083)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    if type==1:
        cases = outcomes[6][0][-1]
    elif type==2:
        cases = outcomes[16][0][-1]
    return cases
    
def Basesens(iterations):
    parameters = numpy.zeros([iterations,9])
    for k in range(0,iterations):
        pdict = {'rho':random.uniform(0.5,1.5),
             'sigma':random.uniform(0.5,1.5),
             'psi':random.uniform(0.5,1.5),
             'theta':random.uniform(0.5,1.5),
             'nu':random.uniform(0.5,1.5),
             'iota':random.uniform(0.5,1.5),
             'tau':random.uniform(0.5,1.5),
             'mu':random.uniform(0.5,1.5)
            } 
        baseline = Baserun(pdict=pdict)
        parameters[k,0] = baseline
        parameters[k,1] = pdict['rho']
        parameters[k,2] = pdict['sigma']
        parameters[k,3] = pdict['psi']
        parameters[k,4] = pdict['theta']
        parameters[k,5] = pdict['nu']
        parameters[k,6] = pdict['iota']
        parameters[k,7] = pdict['tau']
        parameters[k,8] = pdict['mu']
    return parameters

def Seperatesens(iterations):
    parameters = numpy.zeros([iterations,13])
    for k in range(0,iterations):
        pdict = {
        'rho_N':random.uniform(0.5,1.5),
        'rho_D':random.uniform(0.5,1.5),
        'sigma':random.uniform(0.5,1.5),
        'psi':random.uniform(0.5,1.5),
        'theta':random.uniform(0.5,1.5),
        'nu':random.uniform(0.5,1.5),
        'iota_N':random.uniform(0.5,1.5),
        'iota_D':random.uniform(0.5,1.5), 
        'tau_N':random.uniform(0.5,1.5),
        'tau_D':random.uniform(0.5,1.5),
        'mu':random.uniform(0.5,1.5)
            } 
        md = Seperatedrun('Nurse_MD.psc',pdict=pdict,type=1)
        meta = Seperatedrun('MetaPop.psc',pdict=pdict,type=2)
        parameters[k,0] = md
        parameters[k,1] = meta
        parameters[k,2] = pdict['rho_N']
        parameters[k,3] = pdict['rho_D']
        parameters[k,4] = pdict['sigma']
        parameters[k,5] = pdict['psi']
        parameters[k,6] = pdict['theta']
        parameters[k,7] = pdict['nu']
        parameters[k,8] = pdict['iota_N']
        parameters[k,9] = pdict['iota_D']
        parameters[k,10] = pdict['tau_N']
        parameters[k,11] = pdict['tau_D']
        parameters[k,12] = pdict['mu']
    return parameters

base_sweep = Basesens(n_runs)
print("Baseline Sweep Complete, Starting Nurse_MD and MetaPop Sweeps")
sep_sweep = Seperatesens(n_runs)

print("Saving Files")
numpy.savetxt(''.join(['base_sweep_',suffix,'.csv']),base_sweep,delimiter=',',comments=',')
numpy.savetxt(''.join(['sep_sweep_',suffix,'.csv']),sep_sweep,delimiter=',',comments=',')

