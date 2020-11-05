import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si


plt.figure()
plt.plot(fsi[0].flow().px[0])
plt.plot(fsi[0].flow().px[1])
dp = fsi[0].flow().px[0] - fsi[0].flow().px[1]
plt.figure()
plt.plot(dp)
v = fsi[0].solid().eigen.vectors
x = fsi[0].solid().mesh().x
 # Modal pressure loads
lm1 = dp * v[0]
lm2 = dp * v[1]
plt.plot(lm1)
plt.plot(lm2)
plt.figure()
print('i1: ', si.simps(lm1, x))
print('i2: ', si.simps(lm2, x))

