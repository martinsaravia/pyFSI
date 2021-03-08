from pyFSI.mesh.region.boundary1D import *

class boundary1DBeam(boundary1D):
    def __repr__(self):
        return 'boundary1D: boundary1DBeam'

    def __init__(self, beam, control):
        super().__init__(beam.mesh(), control['name'])

        # ----- Public attributes ----- #
        self.name = control['name']

        # ----- Private attributes ----- #
        self._beam = beam
        self._control = control  # Reference to the control dictionary

    # ----- Abstract methods ----- #
    def update(self):
        # Update the position of the boundary
        self.y = self._beam.y[self._control['surface']]
        # Update the time derivatives
        self.dy = self._beam.dy['mid']
        self.ddy = self._beam.ddy['mid']
        self.dyi = si.cumtrapz(self.dy, self._mesh.x, initial=0.0)
        self.ddyi = si.cumtrapz(self.ddy, self._mesh.x, initial=0.0)
        # Re calculate the spatial derivatives and integrals
        self.calculate()

    def isFlexible(self):
        return True

    # ----- Getters ----- #
    def beam(self):
        return self._beam

    def control(self):
        return self._control

    def eigen(self):
        return self._beam.eigen


