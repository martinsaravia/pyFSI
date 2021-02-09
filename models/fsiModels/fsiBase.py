from models.properties.dimensionlessNumbers import makeDimensionlessNumbers
from abc import ABC, abstractmethod

class fsiBase(ABC):
    def __init__(self, execution, control, solid, flow):
        # ----- Public attributes ----- #
        self.dof = None  # Number of degrees of freedom
        self.name = control['name']
        self.output = []
        # Dimensionless numbers
        self.dimNumbers = None

        # ----- Private attributes ----- #
        self._solid = solid
        self._flow = flow
        self._execution = execution
        self._control = control

    @abstractmethod
    def finish(self):
        self._solid.finish()
        self._flow.finish()
        for i in self.output:
            i.close()  # Close all files

    @abstractmethod
    def write(self):
        pass

    def calcNumbers(self):
        self.dimNumbers = makeDimensionlessNumbers(fsi=self)

    # Getters
    def control(self):
        return self._control

    def execution(self):
        return self._execution

    def solid(self):
        return self._solid

    def flow(self):
        return self._flow

