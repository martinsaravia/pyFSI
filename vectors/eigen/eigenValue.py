#---------------------------------------------------------------------------#
#    p    #     version: 0.0
#    y    #     date: 20/03/2020
#    F    #     author: Martin Saravia
#    S    #     description: Class holding eigenvalues.
#    I    #     return: eigenValue Object
#---------------------------------------------------------------------------#
# Notes:
#   The class inherits from fields, I think it could be a better idea to inherit
#   from float.
#---------------------------------------------------------------------------#

from vectors.field import field
import numpy as np

class eigenValue(field):
    def __new__(cls, nparray, info=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(nparray).view(cls)
        # add the new attribute to the created instance
        obj.info = info
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.info = getattr(obj, 'info', None)

    def hz(self):
        return self /( 2 * np.pi)
