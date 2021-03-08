#---------------------------------------------------------------------------#
#    p    #     version: 0.0
#    y    #     date: 20/03/2020
#    F    #     author: Martin Saravia
#    S    #     description: Class holding eigen system.
#    I    #     return: eigen system list of dictionaries
#---------------------------------------------------------------------------#
# Notes:
#   The class returns a list with dictionaries, each subdictionary contains
#   a eigenvalue (evalue key) and an eigenvector (evector key)
#---------------------------------------------------------------------------#

import numpy as np
from pyFSI.vectors.eigen import eigenValue, eigenVector

class eigenSystem(object):

    def __init__(self, evalues=[], evectors=[], sort=False):
        self._eigen = []

        # Sort in ascending order of the eigenvalues imaginary part
        if sort is True:
            idx = np.imag(evalues).argsort()[::1]
            for i, val in enumerate(evalues):
                self.add(evalues[idx[i]], evectors[idx[i]])
        # Not sorted
        else:
            for i, val in enumerate(evalues):
                self.add(evalues[i], evectors[i])

    # main function for adding values a vectors
    def add(self, evalue, evector):
        # Create a dictionary for the value and the vector
        es = {}
        es['eval'] = evalue
        es['evec'] = evector

        # Add to eigensystem list
        self._eigen.append(es)

    # Getter for the eigensystem
    def eigen(self):
        return self._eigen

    # Getter for the eigenvelues (gives them as a list)
    def evalues(self):
        return [es['eval'] for es in self._eigen]

    # Getter for the eigenvectors (gives them as a list)
    def evectors(self):
        return [es['evec'] for es in self._eigen]
    
    def size(self):
        return len(self._eigen)
