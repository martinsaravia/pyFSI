import numpy as np
# 1D Mesh creation.  We can use different creation methods (@classmethod)
class fsiMesh1D(object):

    def __repr__(self):
        return 'fsiMesh1D '

    def __init__(self, iDict, name='None', x=None):
        #print('--> Creating ' + self.__repr__())
        self._name = name
        self._x = x
        self._dim = len(self._x)
        self.L = self._x[-1] - self._x[0]  # Length of the 1D domain

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

    # Getters
    def x(self):
        return self._x

    def dim(self):
        return self._dim

    def name(self):
        return self._name
