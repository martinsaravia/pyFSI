from mesh.region.boundary1D import *

class boundary1DBeam(boundary1D):
    def __repr__(self):
        return 'boundary1D: boundary1DBeam'

    def __init__(self, beam, dict):
        self._beam = beam
        super(boundary1DBeam, self).__init__(beam.mesh(), beam.name())

    # Abstract Methods
    def update(self):
        pass

    def isFlexible(self):
        return True

    def y(self):
        return self._beam.y()

    def ybot(self):
        return self._beam.ybot()

    def ytop(self):
        return self._beam.ytop()

    # Getters
    def eigen(self):
        return self._beam.eigen()

    def N(self):
        return self._beam.N()
