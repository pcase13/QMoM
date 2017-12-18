import numpy as np
from scipy import optimize as opt
import matplotlib.pyplot as plt


r = np.logspace(-2, 3, 10000)
mu, sigma = 5.0, 1.8
ndf = (1 / (2 * np.pi)**0.5 * r * np.log(sigma)) * np.exp(- (np.log(r) - \
       np.log(mu))**2 / (2 * np.log(sigma)**2))
sadf = ndf * r**2

# Create moments object
u = np.zeros(6)
for i in range(6):
    u[i] = np.sum(ndf * r**i)

# Calculate effective radius
reff = u[3]/u[2]

# Instantiate lognormal parameters
def sadf_fit(v, m, r, u):
    mu = np.log(m / np.sqrt(1 + v/m**2))
    sigma = np.sqrt(np.log(1 + v/m**2))
    sadf = (r**2 / (2 * np.pi)**0.5 * r * np.log(sigma)) * np.exp(- (np.log(r) - \
            np.log(mu))**2 / (2 * np.log(sigma)**2))
#    sadf = sadf * r**2
    testu = np.zeros(6)
    for i in range(-1,5):
        testu[i] = np.sum(sadf * r**i)
    residuals = (testu - u)**2
    return np.sum(residuals/u)

x0 = 0.
xopt = opt.fmin(sadf_fit, x0, args=(reff, r, u))
print xopt

v = xopt
m = reff
mu = np.log(m / np.sqrt(1 + v/m**2))
sigma = np.sqrt(np.log(1 + v/m**2))
print 'mu = ' + str(mu)
print 'sigma = ' + str(sigma)
sadf_inverted = (1 / (2 * np.pi)**0.5 * r * np.log(sigma)) * np.exp(- (np.log(r) - \
        np.log(mu))**2 / (2 * np.log(sigma)**2))
sadf_inverted = sadf_inverted * r**2
#plt.semilogx(r, sadf, color='r', label='original')
plt.semilogx(r, ndf, color='r', label='ndf')
plt.semilogx(r, sadf_inverted, color='b', label='inverted')
plt.axvline(x=reff, color='k', label='$r_{eff}$')
plt.legend()
plt.show()
