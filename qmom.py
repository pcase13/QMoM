import numpy as np
import matplotlib.pyplot as plt
from utils import calc_rw
from moment_functions import test_growth_1, test_growth_2
from moment_functions import sulfuric_acid_water_growth
from moment_functions import sulfuric_acid_water_nucleation
from moment_functions import sed_func
from testing_utils import khrgian, analytical_solution_1
Av = 6.022140857e23 #Avogadro's #

# Time settings
run_time = 101#2000. * 60.
dt = 1.

# Monomer Settings
R = 1.0e25 # Source
S = 0. # Sink
v1 = 18.4/Av # Molecular Volume

# Construct Moment Object
f1 = 1.1e30
u = np.zeros(6)
u[0] = 1.0 # Number
u[1] = 5.0 # nm
u[2] = 33.3333 # nm^2
u[3] = 277.7778
u[4] = 2777.7778
u[5] = 32407.4000
#u = u * 1.0e4

f1hist = np.zeros(run_time/dt)
reffhist = np.zeros(run_time/dt)
uhist = np.zeros((run_time/dt, 6))
ansohist = np.zeros((run_time/dt, 6))
rhist = np.zeros((run_time/dt, 3))
whist = np.zeros((run_time/dt, 3))
nuchist = np.zeros((run_time/dt, 6))
sedhist = np.zeros((run_time/dt, 6))
grwhist = np.zeros((run_time/dt, 6))
k = 0
t = 0.
while t < run_time:
    #print 'run time ' + str(t) + 's'
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

    r, w = calc_rw(u/u[0])
    rhist[k, :] = r[:]
    whist[k, :] = w[:]

    #growth_func = sulfuric_acid_water_growth
    growth_func = test_growth_2
    growth = np.zeros(6)
    #nucleation_func = sulfuric_acid_water_nucleation
    nucleation = np.zeros(6)
    #sed_func = sed_func
    sedimentation = np.zeros(6)
    for i in range(6):
        #growth[i] =  i * np.sum(r**(i+1) * growth_func(r, f1) * w)
        growth[i] =  i * np.sum(r**(i-1) * growth_func(r) * w)
        #sedimentation[i] = np.sum(r**(i+1) * sed_func(r) * w)
        sedimentation[i] = 0.
    #nucleation = nucleation_func(u, f1)
    nucleation[:] = 0.
    grwhist[k, :] = growth[:]
    sedhist[k, :] = sedimentation[:]
    nuchist[k, :] = nucleation[:]
    #f1 += R - S - (4*np.pi/v1) * (nucleation[3] + growth[3])
    for i in range(6):
        u[i] += nucleation[i] + growth[i] - sedimentation[i]
    radii = np.linspace(0,100,100000)
    anso = np.zeros(6)
    ndf = analytical_solution_1(radii, t, khrgian, 0.78)
    if t % 25 == 0:
        plt.plot(radii, ndf, label='t = ' + str(t) + 's')
    for i in range(6):
        anso[i] = np.nansum(radii**i * ndf)/1000
        anso[i] = anso[i] / anso[0]
    ansohist[k, :] = anso[:]
    t = t + dt
    k += 1

plt.legend(loc='upper-left')
plt.xlim((0,20))
plt.ylabel('Normalized Probability Density')
plt.xlabel('r ($\mu m$)')
plt.savefig('distributions.png')
print ansohist




fig, axs = plt.subplots(5, figsize=(10,10))
axs[0].plot(uhist[:,0], color='red')
axs[0].plot(ansohist[:,0], linewidth = 3.0, linestyle='--', color='#6C71C3')
axs[0].set_ylabel('Radial Moment 0')
axs[0].set_yticklabels(axs[0].get_yticks())
axs[1].plot(uhist[:,1], color='red')
axs[1].plot(ansohist[:,1], linewidth = 3.0, linestyle='--', color='#6C71C3')
axs[1].set_ylabel('Radial Moment 1')
axs[1].set_yticklabels(axs[1].get_yticks())
axs[2].plot(uhist[:,2], color='red')
axs[2].plot(ansohist[:,2], linewidth = 3.0, linestyle='--', color='#6C71C3')
axs[2].set_ylabel('Radial Moment 2')
axs[2].set_yticklabels(axs[2].get_yticks())
axs[3].plot(uhist[:,3], color='red')
axs[3].plot(ansohist[:,3], linewidth = 3.0, linestyle='--', color='#6C71C3')
axs[3].set_ylabel('Radial Moment 3')
axs[3].set_yticklabels(axs[3].get_yticks())
axs[4].plot(reffhist, color='red')
axs[4].set_ylabel('Effective Radius (micron)')
axs[4].set_yticklabels(axs[4].get_yticks())
plt.ticklabel_format()
plt.tight_layout()
plt.savefig('moments.png')

fig, axs = plt.subplots(3, figsize=(10,10))
axs[0].plot(uhist[:,0], label='k=0')
axs[0].set_ylabel('Radial Moments')
axs[0].plot(uhist[:,1], label='k=1')
axs[0].plot(uhist[:,2], label='k=2')
axs[0].plot(uhist[:,3], label='k=3')
axs[0].legend(loc='center')
axs[1].plot(whist)
axs[1].set_ylabel('Weights')
axs[2].plot(rhist)
axs[2].set_ylabel('Abscissas')
plt.tight_layout()
plt.savefig('u.png')

fig, axs = plt.subplots(4, figsize=(10,10))
axs[0].plot(sedhist)
axs[0].set_ylabel('Sedimentation')
axs[1].plot(nuchist)
axs[1].set_ylabel('Nucleation')
axs[2].plot(grwhist)
axs[2].set_ylabel('Growth')
axs[3].plot(reffhist)
axs[3].set_ylabel('Effective Radius')
plt.tight_layout()
plt.savefig('process_rates.png')
