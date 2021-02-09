# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: eigen vector class
#    I    #     return: eigen vector object
# --------------------------------------------------------------------------- #
# Notes:
#   This class creates an eigenvector object inheriting from a numpy array
#   Note the this class inherits from field, which inherits from np.ndarray
#   Not strictly neccessary, but helps to understand inheritance (sorry)
# --------------------------------------------------------------------------- #
from vectors.field import field
import numpy as np
import scipy.integrate as si


class eigenVector(field):
    def __new__(cls, nparray, mesh, info=None, normalize='mass', mass=None):

        # ----- Casting  ----- #
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        # Choose Normalization
        if normalize == 'length':
            factor = abs((mesh.L))**0.5
            obj = np.asarray(nparray).view(cls) / factor
        if normalize == 'mass':
            factor = (1 / si.simps(mass*nparray**2, mesh.x))**0.5
            obj = np.asarray(nparray).view(cls) * factor
        else:
            obj = np.asarray(nparray).view(cls)

        # ----- Public attributes ----- #
        obj.info = info
        obj.normFactor = factor
        edgeOrder = 2
        obj.d1 = np.gradient(obj, mesh.x, edge_order=edgeOrder)
        obj.d2 = np.gradient(obj.d1, mesh.x, edge_order=edgeOrder)
        obj.d3 = np.gradient(obj.d2, mesh.x, edge_order=edgeOrder)
        obj.d4 = np.gradient(obj.d3, mesh.x, edge_order=edgeOrder)
        obj.ix = si.cumtrapz(obj, mesh.x, initial=0.0)
        obj.iL = si.simps(obj, mesh.x)

        # ----- Private attributes ----- #
        obj._mesh = mesh

        return obj

    # Re-execute the integral after application of boundary conditions
    def correct(self):
        self.ix = si.cumtrapz(self, self._mesh.x, initial=0.0)
        self.iL = si.simps(self, self._mesh.x)

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self.info = getattr(obj, 'info', None)

    # Getters
    def mesh(self):
        return self._mesh

