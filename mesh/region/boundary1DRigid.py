from mesh.region.boundary1D import *

class boundary1DRigid(boundary1D):
    def __repr__(self):
        return 'boundary1D: boundary1DRigid'

    def __init__(self, mesh, name, y):
        self._y = y
        super(boundary1DRigid, self).__init__(mesh, name)


    # constructors
    # Construct from stragith line by two points
    @classmethod
    def fromLineByTwoPoints(cls, mesh, dict, name=None):
        hi = dict['hi']
        hf = dict['hf']
        xvalues = mesh.x()
        yvalues = np.zeros(len(xvalues))
        slope = (hf - hi) / (xvalues[-1] - xvalues[0])
        for i, x in enumerate(xvalues):
            yvalues[i] = slope * x + hi
        return cls(mesh, name, yvalues)

    # Abstract methods
    def update(self):
        pass #rigid boundary are not updated

    def isFlexible(self): # The boundary is rigid
        return False

    def y(self):
        return self._y

    def ybot(self):
        return self.y()

    def ytop(self):
        return self.y()
