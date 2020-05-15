import numpy as np
k = 0*10000
c = 10
m = 2
w = np.sqrt(k/m)
J = np.array([[0, 1], [-k/m, -c/m]])

val, vec = np.linalg.eig(J)

print("Omega: " + str(w))
print("lambda 1: " + str(val[0]) + "--> Eigenvector: " + str(vec[:,0]))
print("lambda 2: " + str(val[1]) + "--> Eigenvector: " + str(vec[:,1]))