import numpy as np

class fsiMesh1D(object):

    def __repr__(self):
        return 'fsiMesh1D '

    def __init__(self, input, name='None', x=None):
        print('--> Creating ' + self.__repr__())
        self.__name = name
        self.__x = x
        self.__dim = len(self.__x)

    # Factory Methods:
    # Here the @classmethods act as different constructors, Note that we also
    # use default arguments in __init__ to flexibilize the constructors
    # Construct from parameters
    @classmethod
    def from3Parameters(cls, input):
        xi = input['xi']
        xf = input['xf']
        dx = input['dx']
        xarray = np.arange(xi, xf+dx, dx)
        return cls(input, x=xarray)

    # Construct from Calculix mesh
    @classmethod
    def fromCalculix(cls, input):
        return 0

    # Getters
    def x(self):
        return self.__x

    def dim(self):
        return self.__dim

    def name(self):
        return self.__name
