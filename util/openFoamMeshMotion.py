import matplotlib.pyplot as plt
import numpy as np

# Diffusivity calculations

# Distance from boundary to cell
x = np.arange(1E-3, 2.50E-3, 1.0E-5)

# Linear and quadratic diffusivities
gammaLinear = (1 / x)
gammaLinear /= gammaLinear[0]
plt.plot(x, gammaLinear)

gammaQuad = gammaLinear**0.5
gammaQuad /= gammaQuad[0]
plt.plot(x, gammaQuad)

# Exponential
alphaList = [1E3, 2E3]
for alpha in alphaList:
    gammaExp = np.exp(-alpha / (1 / x))
    gammaExp /= gammaExp[0]
    plt.plot(x, gammaExp)



plt.show()



# plt.figure()
#