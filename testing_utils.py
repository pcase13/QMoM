import numpy as np

def khrgian(r):
    a = 0.108
    b = 0.6
    return a * r**2 * np.exp(-1 * b * r)

def analytical_solution_1(r, t, f0, k2):
    term1 = np.sqrt(r**2 - 2 * k2 * t)
    return (r / term1) * f0(term1)
