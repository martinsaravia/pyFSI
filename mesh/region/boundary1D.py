# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: Boundary 1D Base class
#    I    #     return:
# --------------------------------------------------------------------------- #
# Notes:
#   Interface for 1D boundaries
#
#
# --------------------------------------------------------------------------- #
from abc import ABCMeta, abstractmethod
import numpy as np
import scipy.integrate as si
from scipy.interpolate import interp1d


class boundary1D(metaclass=ABCMeta):

    def __init__(self, mesh, name, coupled=False):
        # ----- Public attributes ----- #
        self.name = name
        self.size = len(mesh.x)
        self.y = None  # Location of the boundary
        self.d1 = None  # First derivative
        self.ix = None  # Indefinite integral
        self.iL = None  # Definite integral
        self.dy = None  # First time derivative
        self.ddy = None  # Second time derivative
        self.dyi = None  # Indefinite spatial integral of the first time derivative
        self.ddyi = None  # Indefinite spatial integral of the second time deriv
        self.isCoupled = coupled  # if the field is coupled with precice or other

        # ----- Private attributes ----- #
        self._mesh = mesh

    def calculate(self):
        # Ver que en viga estan recalculadas, eliminar? [xxx]
        self.d1 = np.gradient(self.y, self._mesh.x, edge_order=2)
        self.ix = si.cumtrapz(self.y, self._mesh.x, initial=0.0)
        self.iL = si.simps(self.y, self._mesh.x)

    # ----- Abstract methods ----- #
    @abstractmethod
    def update(self):
        pass

    @abstractmethod
    def isFlexible(self):
        # Es necesaria? Esta agregada como atributo tambien [xxx]
        pass

    @classmethod
    def fromLineByTwoPoints(cls, mesh, dict, name=None):
        hi = dict['hi']
        hf = dict['hf']
        xvalues = mesh.x
        yvalues = np.zeros(len(xvalues))
        slope = (hf - hi) / (xvalues[-1] - xvalues[0])
        for i, x in enumerate(xvalues):
            yvalues[i] = slope * x + hi
        return cls(mesh, name, yvalues)

    @classmethod
    def fromPoints(cls, mesh, dict, name=None):
        points = np.array(dict['points'])
        if 'interpolation' in dict:
            interpKind = dict['interpolation']
        else:
            interKind = 'linear'
        interpolator = interp1d(points[:, 0], points[:, 1], kind=interpKind)

        return cls(mesh, name, interpolator(mesh.x))


    # ----- Getters ----- #
    def __call__(self):
        return self.y

    def mesh(self, name):
        return self._mesh

    def __repr__(self):
        return 'boundary1D'

    def write(self):
        pass