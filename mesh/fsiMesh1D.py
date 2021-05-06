# --------------------------------------------------------------------------- #
#    p    #     version: 0.2.0
#    y    #     date: 02/07/2020
#    F    #     author: Martin Saravia
#    S    #     description: 1D Mesh
#    I    #     return: 1D mesh object
# --------------------------------------------------------------------------- #
# Notes:
#   This class creates a 1 dimensional mesh
#
#
# --------------------------------------------------------------------------- #
import numpy as np
from abc import ABCMeta, abstractmethod
# 1D Mesh creation.  We can use different creation methods (@classmethod)

class fsiMesh1D(object):
    def __init__(self, iDict, x=None):

        # ----- Public attributes ----- #
        if 'name' in iDict:
            self.name = iDict['name']
        else:
            self.name = None
        self.x = x  # Mesh coordinates
        self.size = len(self.x)  # number of sampling points
        self.L = self.x[-1] - self.x[0]  # Length of the 1D domain

        # ----- Private attributes ----- #
        self._debug = False

    # Factory Methods:
    # Here the @classmethods act as different constructors, Note that we also
    # use default arguments in __init__ to flexibilize the constructors
    # Construct from 3 parameters

    @classmethod
    def from3Parameters(cls, iDict):
        xi = iDict['xi']
        xf = iDict['xf']
        dx = iDict['dx']
        xarray = np.arange(xi, xf+dx, dx)
        return cls(iDict, x=xarray)

    # Construct from Calculix mesh
    @classmethod
    def fromCalculix(cls, iDict):
        return 0


    @abstractmethod
    def write(self):
        pass

    def __repr__(self):
        return 'fsiMesh1D '

    def debug(self):
        return self._debug


