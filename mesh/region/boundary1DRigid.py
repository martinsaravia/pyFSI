from mesh.region.boundary1D import *

class boundary1DRigid(boundary1D):
    def __repr__(self):
        return 'boundary1D: boundary1DRigid'

    def __init__(self, mesh, name, y):
        super().__init__(mesh, name)

        # ----- Public attributes ----- #
        self.y = y

        # ----- Procedures ----- #
        self.calculate()  # Calculate derivatives
        self.dy = np.zeros(mesh.size)  # Rigid boundaries have zero time derivative
        self.ddy = np.zeros(mesh.size)
        self.dyi = np.zeros(mesh.size)  # Indefinite integral of first time derivative
        self.ddyi = np.zeros(mesh.size)  # Indefinite integral of second time deriv

    # Construct from stragith line by two points
    @classmethod
    def fromLineByTwoPoints(cls, mesh, dict, name=None):
        hi = dict['hi']
        hf = dict['hf']
        xvalues = mesh.x
        yvalues = np.zeros(len(xvalues))
        slope = (hf - hi) / (xvalues[-1] - xvalues[0])
        for i, x in enumerate(xvalues):
            yvalues[i] = slope * x + hi
        return cls(mesh, name, yvalues)

    # Abstract methods
    def update(self):
        pass  # rigid boundaries are not updated

    def isFlexible(self):  # The boundary is rigid
        return False
