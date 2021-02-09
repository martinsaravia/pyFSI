import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt
# Equation 2.29 of Tosi

DP = 100

L = 0.2
b = 5E-3 - 1E-4
h0 = b
hL = 4E-3 - 1E-4
a = (hL - h0) / L
x0 = 0
xL = L
rho = 1.204
zetaIn = 1
zetaOut = 0
xix = 6 / 5.0
f0 = 0.0167

# Viscous term
vt = 0
vt =  -0.25 * f0 * ( (1 / (2 * a * (a * xL + b)**2)) - 1 / (2 * a * (a * x0 + b)**2) )
# Convective term
ct =  rho * xix * ( 1/ (2 *hL **2) - (1 / (2 *h0 **2)) )
# Inlet and outlet term
iot = rho * (zetaIn/(2*h0**2) + zetaOut/(2*hL**2))

Q0 = np.sqrt( DP / (iot + vt + ct))

print('c: ', ct, 'v: ', vt, 'io: ', iot )
print(Q0)
#
x = np.linspace(0, 0.2, 100)
h = ((hL - h0) / L) * x + h0
# anal = -1 / (2*h**2) - (-1 / (2*h0**2) )
nume = si.cumtrapz((1 / h**3) * np.gradient( h,
                                            x,
                                            edge_order=2),
                   x,
                   initial=0)
# plt.plot(x, anal)
# plt.plot(x, nume)



