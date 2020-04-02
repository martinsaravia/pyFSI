import numpy as np

# ------------------------------------------------------------------------------
# This class creates an eigenSystem as a matrix system; i.e. the eigenvectors
# are columns of the system matrix. This is in contrast to the eigenSystemV
# class in which the system is a vector of eigenvectors objects, which is suitable
# for non FEM solutions.
# ------------------------------------------------------------------------------


class eigenSystemMatrix(object):

    def __repr__(self):
        return 'eigenSystemMatrix'

    def __init__(self, values, vectors, name=None):
        # Create the parameters dictionary
        pdict = {}
        pdict['name'] = name
        pdict['values'] = values
        pdict['vectors'] = vectors
        pdict['dim'] = (vectors[0].size,len(values))
        self.__pdict = pdict
        # Init and fill the eigen vector matrices
        self.__vMatrix = np.zeros(pdict['dim'])
        self.__vMatrixD1 = np.zeros(pdict['dim'])
        self.__vMatrixD2 = np.zeros(pdict['dim'])
        self.__vMatrixD3 = np.zeros(pdict['dim'])
        self.__vMatrixD4 = np.zeros(pdict['dim'])
        self.__vMatrixI1 = np.zeros(pdict['dim'])
        for i in range(0,pdict['dim'][1]):
            self.__vMatrix[i] = vectors[i]
            self.__vMatrixD1[:,i] = vectors[i].d1()
            self.__vMatrixD2[:,i] = vectors[i].d2()
            self.__vMatrixD3[:,i] = vectors[i].d3()
            self.__vMatrixD4[:,i] = vectors[i].d4()
            self.__vMatrixI1[:,i] = vectors[i].i1()

    # Setters for boundary conditions
    def setV(self, rg, mode, val):
        self.__vMatrix[rg,mode] = val

    def setD1(self, rg, mode, val):
        self.__vMatrixD1[rg,mode] = val

    def setD2(self, rg, mode, val):
        self.__vMatrixD2[rg,mode] = val

    def setD3(self, rg, mode, val):
        self.__vMatrixD3[rg,mode] = val

    def setD4(self, rg, mode, val):
        self.__vMatrixD4[rg,mode] = val

    # Getters
    def parDict(self):
        return self.__pdict

    def v(self):
        return self.__vMatrix

    def d1(self):
        return self.__vMatrixD1

    def d2(self):
        return self.__vMatrixD2

    def d3(self):
        return self.__vMatrixD3

    def d4(self):
        return self.__vMatrixD4

    def i1(self):
        return self.__vMatrixI1
