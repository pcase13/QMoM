import numpy as np

def test_growth_1(r):
    # Simple version of sulfuric acid/water system
    # phi(r) = k * r
    k = 1.0e-5
    return k * r

# NULL functions: return 0. Use if you want to disable any of the processes.
def null_growth(r, f1):
    return np.zeros(u.shape)
def null_nucleation(u, f1):
    return np.zeros(u.shape)
def null_sedimentation(r):
    return 0

# Growth functions: different particle growth functions
def sulfuric_acid_water_growth(r, f1):
    # r * dr/dt = numerator / (Fk + Fd)
    # numerator = (S-1-(a/r) + (b/r^3))
    # Fk = (L/(RvT)-1)(L*density)/(KT)
    # Fd = (density * Rv * T)/(D*es(T))
    f1norm = 1.0e30
    a = 1.0
    b = 1.0
    S = f1/f1norm
    numerator = (S - 1 - (a/r) + (b/r**3))
    Fk = 1.0
    Fd = 1.0
    return numerator * (Fk + Fd)**-1 #* r**-1

def diffusion_controlled_growth(r, f1):
    # Diffusion controlled growth example
    # phi(r) = k / r
    k = 0.78
    return k / r

# Nucleation functions: different particle nucleation functions
def sulfuric_acid_water_nucleation(u, f1):
    k = 1.0e-9 # molecule / monomer
    r = 0.0001 #micron
    f1norm = 1.0e30
    nuc = np.zeros(len(u))
    for moment in range(len(u)):
        nuc[moment] = (f1/f1norm) * k * r ** moment
    return nuc

# Sedimentation functions: different particle loss functions.
def sedimentation(r):
    k = 1.0e3
    g = 9.8e9 #nm/s2
    density = 1.84e-24 #kg/nm3
    return k * (density * g * r**3)
