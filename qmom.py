import numpy as np
import matplotlib.pyplot as plt
Av = 6.022140857e23 #Avogadro's #


def calc_rw(u):
    # Construct P
    P = np.zeros((7,7))
    P[0, 0] = 1.
    for i in range(1,7):
        P[i-1, 1] = ((-1.)**(i-1)) * u[i-1]
    for j in range(3,8):
        for i in range(1,7):
            P[i-1, j-1] = P[0, j-2]*P[i, j-3] - P[0, j-3]*P[i,j-2]

    # Construct alpha
    alpha = np.zeros(6)
    alpha[0] = 0
    alpha[1] = u[1]
    for n in range(3, 7):
        alpha[n-1] = P[0, n] / (P[0, n-1] * P[0, n-2])

    # Construct a, b
    a = np.zeros(3)
    b = np.zeros(2)
    for n in range(1, 4):
        a[n-1] = alpha[2*n-1] + alpha[2*n-2]
        if n < 3:
            b[n-1] = np.sqrt(np.abs(alpha[2*n]*alpha[2*n-1]))

    # Construct Jacobian
    J = np.zeros((3, 3))
    J[0, 0] = a[0]
    J[1, 1] = a[1]
    J[2, 2] = a[2]
    J[1, 0] = b[0]
    J[0, 1] = b[0]
    J[1, 2] = b[1]
    J[2, 1] = b[1]

    # Find Absicassas and Wieghts
    r, evectors = np.linalg.eig(J)
    r = r[::-1]
    w = np.zeros(3)
    for i in range(3):
        w[i] = u[0] * evectors[0,i]**2
    w = w[::-1]
    return r, w

def sulfuric_acid_water_growth(r, f1):
    k = 0.0000000001
    return f1 * k * r

def sulfuric_acid_water_nucleation(u, f1):
    k = 1.0
    r = 1.
    nuc = np.zeros(len(u))
    for moment in range(len(u)):
        nuc[moment] = f1/100000. * k * r ** moment
    return nuc


# Time settings
run_time = 100. #10. * 60.
dt = 1.

# Monomer Settings
R = 500. # Source
S = 0. # Sink
nucl = 1000. # Monomer nucleation
v1 = 100. # Molecular Volume

# Construct Moment Object
f1 = 1000000.
u = np.zeros(6)
u[0] = 1. # Number
u[1] = 5. # nm
u[2] = 33.3333 # nm^2
u[3] = 277.7778
u[4] = 2777.7778
u[5] = 32407.4000

f1hist = np.zeros(run_time/dt)
reffhist = np.zeros(run_time/dt)
uhist = np.zeros((run_time/dt, 6))
rhist = np.zeros((run_time/dt, 3))
whist = np.zeros((run_time/dt, 3))
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
    reffhist[k] = u[3]/u[2]

    r, w = calc_rw(u)
    rhist[k, :] = r[:]
    whist[k, :] = w[:]

    growth_func = sulfuric_acid_water_growth
    growth = np.zeros(6)
    nucleation_func = sulfuric_acid_water_nucleation
    nucleation = np.zeros(6)
    for i in range(6):
        growth[i] = i * np.sum(r**(i-1) * growth_func(r, f1) * w)
    nucleation = nucleation_func(u, f1)
    print 'nucleation'
    print nucleation
    print 'growth'
    print growth
    f1 += R - S - nucl - (4*np.pi/v1) * (nucleation[3] + growth[3])
    for i in range(6):
        u[i] += nucleation[i] + growth[i]
    t = t + dt
    k += 1
fig, axs = plt.subplots(6, figsize=(10,10))
axs[0].plot(f1hist)
axs[0].set_ylabel('Monomers')
axs[1].semilogy(uhist[:,0])
axs[1].set_ylabel('Radial Moment 0')
axs[2].semilogy(uhist[:,1])
axs[2].set_ylabel('Radial Moment 1')
axs[3].semilogy(uhist[:,2])
axs[3].set_ylabel('Radial Moment 2')
axs[4].semilogy(uhist[:,3])
axs[4].set_ylabel('Radial Moment 3')
axs[5].semilogy(reffhist)
axs[5].set_ylabel('Effective Radius')
plt.tight_layout()
plt.show()
fig, axs = plt.subplots(4, figsize=(10,10))
axs[0].plot(f1hist)
axs[0].set_ylabel('Monomers')
axs[1].semilogy(uhist[:,0], label='k=0')
axs[1].set_ylabel('Radial Moments')
axs[1].semilogy(uhist[:,1], label='k=1')
axs[1].semilogy(uhist[:,2], label='k=2')
axs[1].semilogy(uhist[:,3], label='k=3')
axs[1].legend()
axs[2].plot(whist)
axs[2].set_ylabel('Weights')
axs[3].plot(rhist)
axs[3].plot(whist)
axs[3].set_ylabel('Abscissas')
plt.tight_layout()
plt.show()
