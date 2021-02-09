from mesh.region.boundary1D import *
from abc import ABCMeta, abstractmethod

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

    # Abstract methods
    def update(self):
        pass  # rigid boundaries are not updated

    def isFlexible(self):  # The boundary is rigid
        return False