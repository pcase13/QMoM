import numpy as np
import matplotlib.pyplot as plt
Av = 6.022140857e23 #Avogadro's #

# Settings
infmonomer = True

# Values
R = 1.0e12 # Source per dt
S = 0. # Sink per dt
D = 1.0e-10 # Eddy diffusion coefficient
ny = 10 # Number of vertical levels
dt = 1. # Time per timestep
runtime = 10.*60. # Total time
v1 = 1.8e-5/Av # Molecular Volume
v = 1.0e-5 # Windspeed
nucl = 0.0
a = 1.0e-11
b = -1.0e-6

f1 = np.zeros(ny) # monomer count
f1 += 1.0e12
f1start = np.mean(f1)
u0 = np.zeros(ny) # Number
u0 += 0.0
u1 = np.zeros(ny) # Radius
u1 += 0.0
u2 = np.zeros(ny) # Surface area
u2 += 0.0
u3 = np.zeros(ny) # Mass
u3 += 0.0
f1hist = np.zeros((ny, int(runtime/dt)+1)) # monomer count
u0hist = np.zeros((ny, int(runtime/dt)+1)) # Number
u1hist = np.zeros((ny, int(runtime/dt)+1)) #
u2hist = np.zeros((ny, int(runtime/dt)+1)) # Surface area
u3hist = np.zeros((ny, int(runtime/dt)+1)) # Mass
reffhist = np.zeros((ny, int(runtime/dt)+1)) # Effective Radius
Rhist = np.zeros(int(runtime/dt)+1) # Effective Radius
J = np.zeros(4)
Jmodifier = 1.0e-5
J += 0

running = True
t = 0.
j = 0.
while running:
    # Check for 0s
    f1[f1 < 0.] = 0.0
    u0[u0 < 0.] = 0.0
    u1[u1 < 0.] = 0.0
    u2[u2 < 0.] = 0.0
    u3[u3 < 0.] = 0.0
    # Record in History
    f1hist[:,j] = f1[:]
    u0hist[:,j] = u0[:]
    u1hist[:,j] = u1[:]
    u2hist[:,j] = u2[:]
    u3hist[:,j] = u3[:]
    Rhist[j] = R
    # Create temporary versions
    f1temp = np.copy(f1)
    u0temp = np.copy(u0)
    u1temp = np.copy(u1)
    u2temp = np.copy(u2)
    u3temp = np.copy(u3)
    # Iterate over vertical levels
    for i in range(ny):
        # Calculate Effective Radius
        if u3temp[i] > 0. and u2temp[i] > 0.:
            reff = u3temp[i]/u2temp[i]
            reffhist[i,j] = reff
        else:
            reff = 1.e-9
        # Set Nucleation Rate
        for k in range(len(J)):
            J[k] = Jmodifier * f1temp[i] * (1.e-9) ** k
        nucl = -1*J[3]/v1
        # Update Modes
        if not infmonomer:
            f1[i] += R - S - np.gradient(D*np.gradient(f1))[i] - \
                np.gradient(f1*v)[i] + nucl - (4*np.pi/v1)*(3*a*u2temp[i] + \
                3*b*u3temp[i])
        u0[i] += -np.gradient(D*np.gradient(u0))[i] - np.gradient(u0*v)[i] + \
            J[0]
        u1[i] += -np.gradient(D*np.gradient(u1))[i] - np.gradient(u1*v)[i] + \
            J[1] + a*u0temp[i] + b*u1temp[i]
        u2[i] += -np.gradient(D*np.gradient(u2))[i] - np.gradient(u2*v)[i] + \
            J[2] + 2*a*u1temp[i] + 2*b*u2temp[i]
        u3[i] += -np.gradient(D*np.gradient(u3))[i] - np.gradient(u3*v)[i] + \
            J[3] + 3*a*u2temp[i] + 3*b*u3temp[i]
    t += dt
    j += 1
    if t > runtime:
        running = False
fig, axs = plt.subplots(6, figsize=(10,10))
axs[0].plot(f1hist[0,:])
axs[0].plot(f1hist[1,:])
axs[0].set_ylabel('Monomers')
axs[1].semilogy(u0hist[0,:])
axs[1].semilogy(u0hist[1,:])
axs[1].set_ylabel('Radial Moment 0')
axs[2].semilogy(u1hist[0,:])
axs[2].semilogy(u1hist[1,:])
axs[2].set_ylabel('Radial Moment 1')
axs[3].semilogy(u2hist[0,:])
axs[3].semilogy(u2hist[1,:])
axs[3].set_ylabel('Radial Moment 2')
axs[4].semilogy(u3hist[0,:])
axs[4].semilogy(u3hist[1,:])
axs[4].set_ylabel('Radial Moment 3')
axs[5].semilogy(reffhist[0,:])
axs[5].semilogy(reffhist[1,:])
axs[5].set_ylabel('Effective Radius (m)')
plt.tight_layout()
plt.show()
plt.contourf(reffhist)
plt.colorbar()
plt.xlabel('timestep')
plt.ylabel('y')
plt.tight_layout()
plt.show()
