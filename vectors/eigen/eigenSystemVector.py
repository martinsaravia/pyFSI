# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: BE Beam class
#    I    #     return: solidModel object
# --------------------------------------------------------------------------- #
# Notes:
#   This class creates an eigenSystem as a vector of eigenvectors
#
#
# --------------------------------------------------------------------------- #
import numpy as np
import scipy.integrate as si

from pyFSI.vectors.field import field

class eigenSystemVector(object):

    def __repr__(self):
        return 'eigenSystemVector'

    def __init__(self, values, vectors, name=None):

        # ----- Public attributes ----- #
        self.name = name
        self.size = len(values)
        self.values = values  # Eigenvalues
        self.vectors = np.empty(self.size, dtype=object) # Eigenvectors
        self.normalization = vectors[0].normalization
        self.d1 = np.empty(self.size, dtype=object)  # First derivative
        self.d2 = np.empty(self.size, dtype=object)
        self.d3 = np.empty(self.size, dtype=object)
        self.d4 = np.empty(self.size, dtype=object)
        self.ix = np.empty(self.size, dtype=object)  # indefinite integrals
        self.iL = np.empty(self.size)  # definite integrals in L
        self.N = np.empty((self.size, self.size), dtype=object)  # Norm

        # ----- Private attributes ----- #
        self._mesh = vectors[0].mesh()

        # ----- Procedures ----- #
        # Assemble the eigensystem and its derivatives
        for i in range(self.size):
            self.vectors[i] = vectors[i]
            self.d1[i] = vectors[i].d1
            self.d2[i] = vectors[i].d2
            self.d3[i] = vectors[i].d3
            self.d4[i] = vectors[i].d4
            self.ix[i] = vectors[i].ix
            self.iL[i] = vectors[i].iL

    def correct(self):  # Correct the integral after setting the BCs
        for i in range(self.size):
            self.vectors[i].correct()
            self.N[i, i] = si.simps(self.vectors[i]**2, self._mesh.x)  # OrthoNorm

    # Reconstruct displacements from the state vector and the modes
    def reconstruct(self, state):
        size = self.vectors[0].size  # Number of elements of eigenvectors
        fx = field(np.zeros(size))  # Create a field of zeros
        for i in range(self.size):
            fx += state[i] * self.vectors[i]
        return fx

    # Setters for boundary conditions
    def setV(self, rg, mode, val):
        self.vectors[mode][rg] = val

    def setD1(self, rg, mode, val):
        self.d1[mode][rg] = val

    def setD2(self, rg, mode, val):
        self.d2[mode][rg] = val

    def setD3(self, rg, mode, val):
        self.d3[mode][rg] = val

    def setD4(self, rg, mode, val):
        self.d4[mode][rg] = val

    # Getters
    # def __call__(self):
    #     return self.vectors

