import matplotlib.pyplot as plt
import numpy as np

rho = 8780
b = 0.1
h = 4.0E-4
L = 0.2
E = 1.1E11
Iz = b * h**3 /12
A = b * h
g = 9.8E-2
qw = rho * g * A

x = L
# Cantilever displacement
v = (-qw * x**2) * (6 * L**2 - 4 * L * x + x**2) / (24 * E * Iz)
vMax = qw * L**4 / (8 * E * Iz)
print(v, vMax)

