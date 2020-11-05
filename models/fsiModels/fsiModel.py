from properties.dimensionlessNumbers import makeDimensionlessNumbers
from abc import ABCMeta, abstractmethod

class fsiModel(metaclass=ABCMeta):
    def __init__(self, execution, control, beam, flow):
        # ----- Public attributes ----- #
        self.dof = None  # Number of degrees of freedom
        self.name = control['name']
        # Store objects
        self._beam = beam
        self._flow = flow
        self._execution = execution
        self._control = control

        # Dimensionless numbers
        self.dimNumbers = None
        # self.calcNumbers()

    def finish(self):
        self._beam.finish()
        self._flow.finish()

    def write(self):
        self._beam.write()
        self._flow.write()

    def calcNumbers(self):
        self.dimNumbers = makeDimensionlessNumbers(fsi=self)

    # Getters
    def control(self):
        return self._control

    def execution(self):
        return self._execution

    def solid(self):
        return self._beam

    def flow(self):
        return self._flow

