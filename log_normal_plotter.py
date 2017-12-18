import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm
from scipy.optimize import curve_fit

def calc_rw(u):
    # Set u0 to 1
    utemp = np.copy(u)
    utemp = utemp/u[0]

    # Construct P
    P = np.zeros((7,7))
    P[0, 0] = 1.
    for i in range(1,7):
        P[i-1, 1] = ((-1.)**(i-1)) * utemp[i-1]
    for j in range(3,8):
        for i in range(1,7):
            P[i-1, j-1] = P[0, j-2]*P[i, j-3] - P[0, j-3]*P[i,j-2]

    # Construct alpha
    alpha = np.zeros(6)
    alpha[0] = 0
    alpha[1] = utemp[1]
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
norm = True
N, mu, sigma = 100., 5.0, 1.8
#mu = np.log(mu) - sigma**2 / 2
r = np.logspace(-2, 3, 1000000)
ndf = (N / (2 * np.pi)**0.5 * r * np.log(sigma)) * np.exp(- (np.log(r) - np.log(mu))**2 / (2 * np.log(sigma)**2))
sadf = ndf * 2 * np.pi * r**2
vdf = ndf * 4./3. * np.pi * r**3
reff = np.sum(vdf)/np.sum(sadf)

# CARMA distributions
carma_x1 = [0.0424942, 0.05574907, 0.07313842, 0.0959519,  0.1258814,
            0.16514656, 0.21665939, 0.2842402, 0.37290094, 0.4892169,
            0.64181437, 0.84201032, 1.10465177, 1.44921684, 1.9012593,
            2.4943037, 3.27233163, 4.2930435, 5.63213775, 7.38892482,
            9.69369223, 12.71736705]
carma_w1 = [1.32548708e-06, 1.73893552e-06, 2.28134758e-06, 2.99294983e-06,
            3.92651640e-06, 5.15128283e-06, 6.75808072e-06, 8.86607404e-06,
            1.16315966e-05, 1.52597463e-05, 2.00195952e-05, 2.62641452e-05,
            3.44565071e-05, 4.52042459e-05, 5.93044398e-05, 7.78027929e-05,
            1.02071187e-04, 1.33909425e-04, 1.75678706e-04, 2.30476741e-04,
            3.02367482e-04, 3.96682519e-04]
carma_w1 = np.asarray(carma_w1)
carma_w1 = carma_w1 * 10**4
carma_y1 = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.0001,
            0.0008, 0.0031, 0.0103, 0.0278, 0.0609, 0.1082,
            0.1557, 0.1817, 0.1719, 0.1318, 0.0820]

carma_x2 = [0.42494199,   0.55749069,   0.73138425,   0.959519,     1.25881399,
            1.65146563,   2.16659391,   2.84240198,   3.72900939,   4.89216905,
            6.41814367,   8.42010319,  11.04651771,  14.49216842,  19.01259301,
            24.94303699,  32.72331629,  42.93043502,  56.32137753,  73.88924818,
            96.93692227, 127.17367047]
carma_w2 = [1.32548708e-05,  1.73893552e-05,  2.28134758e-05,  2.99294983e-05,
            3.92651640e-05,  5.15128283e-05,  6.75808072e-05,  8.86607404e-05,
            1.16315966e-04,  1.52597463e-04,  2.00195952e-04,  2.62641452e-04,
            3.44565071e-04,  4.52042459e-04,  5.93044398e-04,  7.78027929e-04,
            1.02071187e-03,  1.33909425e-03,  1.75678706e-03,  2.30476741e-03,
            3.02367482e-03,  3.96682519e-03]
carma_w2 = np.asarray(carma_w2)
carma_w2 = carma_w2 * 10**4
carma_y2 = np.zeros(22)
carma_y2 = carma_y2 + 1

"""
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
"""
mus = np.zeros(6)
prenorm = np.zeros(6)
i = 0
while i < len(mus):
    mus[i] = np.sum(ndf * r**i)
    i += 1
abscissas, weights = calc_rw(mus)
print 'abscissas'
print abscissas
print 'weights'
print weights
i=0
while i < len(mus):
    prenorm[i] = np.sum(weights * abscissas**i)
    i += 1
print 'prenorm mus'
print mus
print 'prenorm calcmus'
print prenorm


# Product Difference
mus = mus / mus[0]
abscissas, weights = calc_rw(mus)
i = 0
calcmus = np.zeros(6)
while i < len(mus):
    calcmus[i] = np.sum(weights * abscissas**i)
    i += 1
print 'abscissas'
print abscissas
print 'weights'
print weights
print 'calcmus'
print calcmus
print 'mus'
print mus

"""

# Normalize
if norm:
    ndf = ndf/np.max(ndf)
    sadf = sadf/np.max(sadf)
    #sadf_inverted = sadf_inverted/np.max(sadf_inverted)
    vdf = vdf/np.max(vdf)
    #weights = weights/np.max(weights)
    carma_y1 = carma_y1/np.max(carma_y1)


# Plot
plt.semilogx(r, ndf, linewidth=2, color='k', label='Niemeier 2009')
plt.bar(carma_x1, carma_y2, carma_w1, color='r', alpha=0.5, label='Current Bins')
plt.bar(carma_x2, carma_y2, carma_w2, color='b', alpha=0.5, label='rmin=0.5um Bins')
plt.bar(carma_x1, carma_y1, carma_w1, color='y', label='Current Distribution')
#plt.semilogx(r, sadf_inverted, linewidth=2, color='k', label='Inverted Surface Area')
#plt.semilogx(r, sadf, linewidth=2, color='y', label='Surface Area')
#plt.semilogx(r, vdf, linewidth=2, color='g', label='Volume')
#plt.scatter(abscissas, weights)
#plt.axvline(x=mu, color='b', label='$\mu$')
#plt.axvline(x=reff, color='k', label='$r_{eff}$')
plt.legend(loc=9, ncol=2, bbox_to_anchor=(0.5,-0.2))
plt.title('dMash')
ymove = 0.2
#plt.figtext(0.7, ymove + 0.675, '$N = ' + str(N) + '$')
#plt.figtext(0.7, ymove + 0.625, '$\mu = ' + str(mu) + '$')
#plt.figtext(0.7, ymove + 0.575, '$\sigma = ' + str(sigma) + '$')
#plt.figtext(0.7, ymove + 0.525, '$Normalized = ' + str(norm) + '$')
#plt.figtext(0.7, ymove + 0.475, '$r_{eff} = ' + str(reff) + '$')
plt.xlabel('Radius ($\mu m$)')
plt.ylabel('$dM/dr$')
plt.subplots_adjust(bottom=0.35)
plt.show()
"""
