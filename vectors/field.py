import numpy as np

class field(np.ndarray):

    def __new__(cls, nparray, info=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(nparray).view(cls)
        # add the new attribute to the created instance
        obj.info = info
        obj.size
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.info = getattr(obj, 'info', None)

#use the @property decorator to set an atribute called everytime is invoked
