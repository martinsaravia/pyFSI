import numpy as np
import scipy.integrate as si


class fsiRegion1D(object):

    def __repr__(self):
        return 'fsiRegion1D'

    def __init__(self, mesh, bbot, btop, name=None):
        self.__name = name
        # Channel gap
        self.__bot = bbot.y() # values of the bottom
        self.__top = btop.y() # values of the top
        self.__he = self.__top - self.__bot

        # Eigen system
        if bbot.isFlexible():
            flexBoundary = bbot
        else:
            flexBoundary = btop

        self.__G = flexBoundary.eigen()


        self.__xi = mesh.x()[0] # x initial limit
        self.__xf = mesh.x()[-1] # x final limit
        self.__bot = bbot.y() # values of the bottom
        self.__top = btop.y() # values of the top
        self.__he = self.__top - self.__bot

        # Calculate common integrals and store in dict
        he = self.__he
        ix = {}
        ix['i1'] = si.cumtrapz(1.0/he, mesh.x(), initial=0.0)
        ix['i2'] = si.simps(1.0/he, mesh.x())
        self.__integrals = ix

    # Getters
    def name(self):
        return self.__name

    def he(self):
        return self.__he

    def bot(self):
        return self.__bot

    def top(self):
        return self.__top

    def mesh(self):
        return self.__mesh # Return reference to the mesh

    # Return integral by key
    def it(self, key):
        return self.__integrals[key]
