import numpy as np

w0 = 36.0  # First natural frequency (rad/seg)
w1 = 1225.0  # Last Natural frequency

seda = 0.02   # Desired damping coefficient

S = np.array([[1/(2.0*w0), 0.5*w0],
             [1/(2.0*w1), 0.5*w1]])

f = np.array([seda, seda])

coeff = np.linalg.solve(S,f)

print(coeff)