from vectors.field import field
import numpy as np
from scipy import integrate

# Note the this class inherits from field, which inherits from np.ndarray
# Not strictly neccessary, but helps to understand inheritance (sorry)

class eigenVector(field):
    def __new__(cls, nparray, mesh, info=None, normalize=True):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type



        # Normalization
        if normalize:
            #print("---> Normalizing eigenvector...")
            factor = abs((mesh.x()[-1]-mesh.x()[0]))**0.5
        else:
            factor = 1

        obj = np.asarray(nparray).view(cls) / factor


        obj.__mesh = mesh

            # obj = obj / max([obj.max(), obj.min()], key=abs)

        # add the new attribute to the created instance
        obj.info = info

        # Attemtp using mirrrored ghost boundary (does not work)
        # temp = np.zeros(obj.size + 3)
        # temp[3:] = obj
        # temp[2] = obj[1]
        # temp[1] = obj[2]
        # temp[0] = obj[3]
        # tempd1 = np.gradient(temp, edge_order=edgeOrder)
        # tempd2 = np.gradient(tempd1 , edge_order=edgeOrder)
        # tempd3 = np.gradient(tempd2 , edge_order=edgeOrder)
        # tempd4 = np.gradient(tempd3,  edge_order=edgeOrder)
        # obj.__d1 = tempd1[3:]
        # obj.__d2 = tempd2[3:]
        # obj.__d3 = tempd3[3:]
        # obj.__d4 = tempd4[3:]

        edgeOrder = 2
        obj.__d1 = np.gradient(obj, mesh.x(), edge_order=edgeOrder)
        obj.__d2 = np.gradient(obj.__d1, obj.__mesh.x(), edge_order=edgeOrder)
        obj.__d3 = np.gradient(obj.__d2, obj.__mesh.x(), edge_order=edgeOrder)
        obj.__d4 = np.gradient(obj.__d3, obj.__mesh.x(), edge_order=edgeOrder)
        # obj.__i1i = integrate.cumtrapz(obj, mesh.x(), initial=0.0)
        # obj.__i1d = integrate.simps(obj, mesh.x())

        # Frequency
        obj.freq = 0

        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.info = getattr(obj, 'info', None)

    # Getters
    def mesh(self):
        return self.__mesh

    def d1(self):
        return self.__d1

    def d2(self):
        return self.__d2

    def d3(self):
        return self.__d3

    def d4(self):
        return self.__d4

    # def i1i(self):
    #     return self.__i1i
    #
    # def i1d(self):
    #     return self.__i1d

    # Operators
    # def __mul__(self, other):
    #     return eigenVector(self[:] * other[:])


# class eigenVector(field):
#     def __new__(cls, field):
#         # shape = (dim,)
#         obj = field
#         return obj
#
#     def __init__(self):
#         print('verga')
#         self.d1 = np.gradient(self)

    # @property
    # def d1(self):
    #     return self._d1
    #
    # @d1.setter
    # def d1(self):
    #     self._d1 = np.gradient(self.view())
