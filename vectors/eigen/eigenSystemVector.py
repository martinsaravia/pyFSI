import numpy as np
import scipy.integrate as si
# ------------------------------------------------------------------------------
# This class creates an eigenSystem as a vector system; i.e. the eigenvectors
# are objects of the system vector.
# ------------------------------------------------------------------------------

class eigenSystemVector(object):

    def __repr__(self):
        return 'eigenSystemVector'

    def __init__(self, values, vectors, name=None):

        # Flags
        self.calculated = False # Calculation of mode integrals and other stuff

        # Create the parameters dictionary
        pdict = {}
        pdict['name'] = name
        pdict['values'] = values
        pdict['vectors'] = vectors
        pdict['size'] = len(values)
        self.pdict = pdict

        # Init the vector of eigenvectors and its derivatives
        self.vMatrix    = np.empty(pdict['size'], dtype=object)
        self.vMatrixD1  = np.empty(pdict['size'], dtype=object)
        self.vMatrixD2  = np.empty(pdict['size'], dtype=object)
        self.vMatrixD3  = np.empty(pdict['size'], dtype=object)
        self.vMatrixD4  = np.empty(pdict['size'], dtype=object)
        for i in range(0, pdict['size']):
            self.vMatrix[i]   = vectors[i]
            self.vMatrixD1[i] = vectors[i].d1()
            self.vMatrixD2[i] = vectors[i].d2()
            self.vMatrixD3[i] = vectors[i].d3()
            self.vMatrixD4[i] = vectors[i].d4()


    def calculate(self):
        dict = self.dict()
        mesh = self.vMatrix[0].mesh()  # Reference to the mesh fron any eigenvector
        # Calculate integrals
        if not self.calculated:
            self.vix = np.empty(dict['size'], dtype=object)  # vector of indefinite integrals
            self.viL = np.empty(dict['size'])
            self.N   = np.empty((dict['size'], dict['size']), dtype=object)
            for i in range(0, dict['size']):
                self.vix[i] = si.cumtrapz(self.vMatrix[i], mesh.x(), initial=0.0)
                self.viL[i] = si.simps(self.vMatrix[i], mesh.x())
                self.N[i, i] = si.simps(self.vMatrix[i]**2, mesh.x()) # Norm assuming orthogonality

        self.calculated = True

    # Setters for boundary conditions
    def setV(self, rg, mode, val):
        self.vMatrix[mode][rg] = val

    def setD1(self, rg, mode, val):
        self.vMatrixD1[mode][rg] = val

    def setD2(self, rg, mode, val):
        self.vMatrixD2[mode][rg] = val

    def setD3(self, rg, mode, val):
        self.vMatrixD3[mode][rg] = val

    def setD4(self, rg, mode, val):
        self.vMatrixD4[mode][rg] = val

    # Getters
    def __call__(self):
        return self.pdict

    def dict(self):
        return self.pdict

    def values(self):
        return self.pdict['values']

    def v(self):
        return self.vMatrix

    def d1(self):
        return self.vMatrixD1

    def d2(self):
        return self.vMatrixD2

    def d3(self):
        return self.vMatrixD3

    def d4(self):
        return self.vMatrixD4

    def ix(self):
        return self.vix

    def iL(self):
        return self.viL

    def norm(self):
        return self.N
