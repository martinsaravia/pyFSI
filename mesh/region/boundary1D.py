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

class boundary1D(metaclass=ABCMeta):

    def __init__(self, mesh, name):
        # ----- Public attributes ----- #
        self.name = name
        self.y = None  # Location of the boundary
        self.d1 = None  # First derivative
        self.ix = None  # Indefinite integral
        self.iL = None  # Definite integral
        self.dy = None  # First time derivative
        self.ddy = None  # Second time derivative
        self.dyi = None  # Indefinite spatial integral of the first time derivative
        self.ddyi = None  # Indefinite spatial integral of the second time deriv

        # ----- Private attributes ----- #
        self._mesh = mesh

    def calculate(self):
        self.d1 = np.gradient(self.y, self._mesh.x, edge_order=2)
        self.ix = si.cumtrapz(self.y, self._mesh.x, initial=0.0)
        self.iL = si.simps(self.y, self._mesh.x)

    # ----- Abstract methods ----- #
    @abstractmethod
    def update(self):
        pass

    @abstractmethod
    def isFlexible(self):
        pass

    # ----- Getters ----- #
    def __call__(self):
        return self.y

    def mesh(self, name):
        return self._mesh

    def __repr__(self):
        return 'boundary1D'
