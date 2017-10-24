import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm
from scipy.optimize import curve_fit

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

# Construct distributions
norm = False
N, mu, sigma = 100., 3.4, 1.8
r = np.linspace(0.1, 100000, 1000000)
ndf = (N / (2 * np.pi)**0.5 * r * np.log(sigma)) * np.exp(- (np.log(r) - np.log(mu))**2 / (2 * np.log(sigma)**2))
sadf = ndf * 2 * np.pi * r**2
vdf = ndf * 4./3. * np.pi * r**3
reff = np.sum(vdf)/np.sum(sadf)

# Try to invert
mean = reff
variance = 100.0
mu_inverted = np.log(mean / (np.sqrt(1 + variance / mean**2)))
print 1 + variance/ mean **2
sigma_inverted = np.log(1 + variance / mean**2)
print mu_inverted
print sigma_inverted
sadf_inverted = (N / (2 * np.pi)**0.5 * r * np.log(sigma_inverted)) * np.exp(- (np.log(r) - np.log(mu_inverted))**2 / (2 * np.log(sigma_inverted)**2))

# Construct Moments
mus = np.zeros(6)
i = 0
while i < len(mus):
    mus[i] = np.sum(ndf * r**i)
    i += 1

# Product Difference
mus = mus / mus[0]
abscissas, weights = calc_rw(mus)
i = 0
calcmus = np.zeros(6)
while i < len(mus):
    calcmus[i] = np.sum(weights * abscissas**i)
    i += 1

# Normalize
if norm:
    ndf = ndf/np.max(ndf)
    sadf = sadf/np.max(sadf)
    sadf_inverted = sadf_inverted/np.max(sadf_inverted)
    vdf = vdf/np.max(vdf)
    weights = weights/np.max(weights)


# Plot
plt.semilogx(r, ndf, linewidth=2, color='r', label='Number')
plt.semilogx(r, sadf_inverted, linewidth=2, color='k', label='Inverted Surface Area')
#plt.semilogx(r, sadf, linewidth=2, color='y', label='Surface Area')
#plt.semilogx(r, vdf, linewidth=2, color='g', label='Volume')
plt.scatter(abscissas, weights)
plt.axvline(x=mu, color='b', label='$\mu$')
plt.axvline(x=reff, color='k', label='$r_{eff}$')
plt.legend(loc=9, ncol=2, bbox_to_anchor=(0.5,-0.2))
plt.title('Lognormal Distribution')
ymove = 0.15
plt.figtext(0.7, ymove + 0.675, '$N = ' + str(N) + '$')
plt.figtext(0.7, ymove + 0.625, '$\mu = ' + str(mu) + '$')
plt.figtext(0.7, ymove + 0.575, '$\sigma = ' + str(sigma) + '$')
plt.figtext(0.7, ymove + 0.525, '$Normalized = ' + str(norm) + '$')
plt.figtext(0.7, ymove + 0.475, '$r_{eff} = ' + str(reff) + '$')
plt.xlabel('Radius ($\mu m$)')
plt.ylabel('$dN/dr$')
plt.subplots_adjust(bottom=0.35)
plt.show()
