import numpy as np
import matplotlib.pyplot as plt
import sys
import configparser
from utils import calc_rw
from dictionaries import func_dict
from testing_utils import khrgian, analytical_solution_1
Av = 6.022140857e23 #Avogadro's #

# General Settings
analytical_solution = False
sim_name = sys.argv[1]
config = configparser.ConfigParser()
config.read('experiments/' + sim_name + '.ini')

# Time settings
run_time = float(config.get('OPTIONS', 'RUN_TIME'))
dt = float(config.get('OPTIONS', 'DT'))
steps = int(run_time/dt)

# Monomer Settings
inf_monomer = config.getboolean('OPTIONS', 'INF_MONOMER') # If true, monomer field remains unchanged
f1 = 1.1e30 # Initial Monomer Count
R = 0 # Monomer Source
S = 0. # Monomer Sink
v1 = 18.4/Av # Molecular Volume

# Construct Moment Object
u = np.zeros(6)
u[0] = 1.0 # Number
u[1] = 5.0 # First moment (micron)
u[2] = 33.3333 # Second moment (micron^2)
u[3] = 277.7778 # Third moment (micron^3)
u[4] = 2777.7778 # Fourth moment (micron^4)
u[5] = 32407.4000 # Fifth moment (micron^5)

# Assign moment moment functions
growth_func = func_dict[config.get('OPTIONS', 'GRW_FUNC')]
sedimentation_func = func_dict[config.get('OPTIONS', 'SED_FUNC')]
nucleation_func = func_dict[config.get('OPTIONS', 'NUC_FUNC')]

# Construct history objects
f1hist = np.zeros(steps)
reffhist = np.zeros(steps)
uhist = np.zeros((steps, 6))
rhist = np.zeros((steps, 3))
whist = np.zeros((steps, 3))
nuchist = np.zeros((steps, 6))
sedhist = np.zeros((steps, 6))
grwhist = np.zeros((steps, 6))
if analytical_solution is True:
    ansohist = np.zeros((steps, 6))

# Initiate time loop
k = 0
t = 0.
while t < run_time:
    # Check for zeroes
    if f1 < 0.:
        f1 = 0.
    for moment in range(len(u)):
        if u[moment] < 0.:
            u[moment] = 0.

    # Record in history
    f1hist[k] = f1
    uhist[k, :] = u[:]
    if u[2] != 0.0:
        reffhist[k] = u[3]/u[2]

    # Calculate abscissas, weights
    r, w = calc_rw(u/u[0])
    rhist[k, :] = r[:]
    whist[k, :] = w[:]

    # Assign moment equation functions
    growth = np.zeros(6)
    nucleation = np.zeros(6)
    sedimentation = np.zeros(6)

    # Calculate moment equation changes with abscissas, weights
    for i in range(6):
        growth[i] =  i * np.sum(r**(i-1) * growth_func(r, f1) * w)
        #growth[i] =  np.sum(r**(i-1) * growth_func(r, f1) * w)
        sedimentation[i] = np.sum(r**(i+1) * sedimentation_func(r) * w)
    nucleation = nucleation_func(u, f1)

    # Record process histories
    grwhist[k, :] = growth[:]
    sedhist[k, :] = sedimentation[:]
    nuchist[k, :] = nucleation[:]

    # Change array values
    if not inf_monomer:
        f1 += R - S - (4*np.pi/v1) * (nucleation[3] + growth[3])
    for i in range(6):
        u[i] += nucleation[i] + growth[i] - sedimentation[i]

    # Check if analytical solution is requested, calculate for timestep
    if analytical_solution is True:
        radii = np.linspace(0,100,100000)
        anso = np.zeros(6)
        ndf = analytical_solution_1(radii, t, khrgian, 0.78)
        for i in range(6):
            anso[i] = np.nansum(radii**i * ndf)/1000
            anso[i] = anso[i] / anso[0]
        ansohist[k, :] = anso[:]

    # Move time forward
    t = t + dt
    k += 1

# Save histories
save_prefix = 'data/' + sim_name + '_'
np.save(save_prefix + 'uhist.npy', uhist)
np.save(save_prefix + 'f1hist.npy', f1hist)
np.save(save_prefix + 'rhist.npy', rhist)
np.save(save_prefix + 'whist.npy', whist)
np.save(save_prefix + 'reffhist.npy', reffhist)
np.save(save_prefix + 'grwhist.npy', grwhist)
np.save(save_prefix + 'sedhist.npy', sedhist)
np.save(save_prefix + 'nuchist.npy', nuchist)
if analytical_solution is True:
    np.save(save_prefix + 'ansohist.npy', ansohist)
